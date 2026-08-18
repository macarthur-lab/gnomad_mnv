"""Single-pass scan-based MNV discovery for gnomAD v4.

Operates on the unsplit sparse VDS variant_data in local allele space
(LGT/LPGT/LA/PID), extracting biallelic alleles via LA indexing so no separate
``sparse_split_multi`` step is needed. The scan tracks a window of recent
non-ref entries per sample using ``hl.scan.fold`` with ``hl.case()`` branching.

Before the scan, genotypes are adjusted in gnomAD QC's canonical order: flag
het_non_ref, apply the hom-alt depletion hotfix, adjust sex ploidy, then
annotate adj. Each non-ref entry is then reduced to one record per distinct
carried alt allele, so het_non_ref (``1/2``) carriers contribute both alts.

For each sample, finds pairs of SNVs within ``max_distance`` bp that co-occur
on the same haplotype. Alt pairs are classified as:

    - het-het: both het, both phased, same PID, same haplotype (in cis)
    - hom-hom: both homozygous
    - het-hom: exactly one homozygous (always in cis)

Results are aggregated per variant pair across all samples, with both raw counts
(all carriers) and ``_adj`` counts (both constituent genotypes pass adj), and
written as a Hail Table.

Usage::

    python -m v4.discover_mnv --discover --overwrite
    python -m v4.discover_mnv --annotate
    python -m v4.discover_mnv --discover --annotate --test
    python -m v4.discover_mnv --discover --annotate --intervals \\
        PCSK9=chr1:55039447-55064852 PCNT=chr21:46324141-46445769
    python -m v4.discover_mnv --discover --chr 21 --overwrite
    python -m v4.discover_mnv --discover --ukb-only --overwrite
"""

import argparse
import logging
import time
from typing import List, Optional, Tuple

import hail as hl
from gnomad.resources.grch38.gnomad import public_release
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import get_adj_expr
from gnomad_qc.v4.resources.annotations import get_info, get_vep
from gnomad_qc.v4.resources.meta import meta

from v4.resources import (
    MNV_ENTRIES_TO_KEEP,
    get_annotated_mnvs,
    get_discovered_mnvs,
    get_gnomad_v4_ukb_vds,
    get_gnomad_v4_vds,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("discover_mnv")
logger.setLevel(logging.INFO)

MAX_MNV_DISTANCE = 2
"""Default maximum distance in bp between the two SNVs of an MNV pair.

2 keeps both SNVs within one codon frame, matching the v2 pipeline and the
source paper; override per run with ``--max-distance``.
"""

HOM_ALT_FIX_AF_CUTOFF = 0.01
"""AF above which a high-AB het-ref call is reclassified to hom-alt.

Part of the GATK hom-alt depletion hotfix (see :func:`discover_mnv`); mirrors the
gnomAD default (``gnomad_qc.v3.utils.hom_alt_depletion_fix``).
"""

HOM_ALT_FIX_AB_CUTOFF = 0.9
"""Allele balance above which a common-variant het-ref call is treated as hom-alt."""

CLASSIFY_N_PARTITIONS = 500
"""Suggested value for ``--classify-n-partitions`` when needed (the CLI arg
defaults to None/off): a reasonable starting point if a classification-step
off-heap OOM needs relieving."""

DEFAULT_TEST_GENES = [
    # PCNT on GRCh38.
    ("PCNT", "chr21:46324141-46445769"),
]
"""Default gene(s)/interval(s) used when ``--intervals`` is given with no values."""


# ---------------------------------------------------------------------------
# Expression helpers
# ---------------------------------------------------------------------------


def _is_nonref(entry: hl.expr.StructExpression) -> hl.expr.BooleanExpression:
    """
    Check if an entry has a defined non-ref local genotype.

    :param entry: Entry expression with an ``LGT`` field.
    :return: Boolean expression.
    """
    return hl.is_defined(entry.LGT) & entry.LGT.is_non_ref()


def _carries_snv(entry: hl.expr.StructExpression) -> hl.expr.BooleanExpression:
    """
    Check if an entry carries at least one SNV alt, i.e. can form an MNV pair.

    Decides whether an entry enters the scan window. :func:`_is_nonref` isn't
    enough because a site can carry both SNV and indel alts and a sample may carry only the
    indel — non-ref, yet its ``_alts`` is empty because :func:`_get_carried_alts`
    keeps SNV alts only. Such an entry can never yield an output pair
    (classification keeps SNP-SNP pairs only), so excluding it here prevents
    bloating the window state, the ``prev`` arrays, and the checkpoint.

    :param entry: Entry expression with the ``_alts`` field built by
        :func:`_get_carried_alts` (missing for ref entries).
    :return: Boolean expression.
    """
    return hl.is_defined(entry._alts) & (hl.len(entry._alts) > 0)


def _get_carried_alts(entry: hl.expr.StructExpression) -> hl.expr.ArrayExpression:
    """
    Get one record per distinct carried alt of a non-ref genotype.

    An array (not one struct) so het_non_ref (``1/2``) contributes both alts;
    het-ref/hom-var yield one. Alleles are min-repped so output keys join the
    split/min-repped release HT — a raw multiallelic alt (e.g.
    ``["TCCGGG","GCCGGG"]``) would miss the join. Built at localization, not per
    pair, so the row ``alleles`` array isn't carried through the scan window on
    every sample's entry.

    :param entry: Entry with ``LGT``, ``LPGT``, ``LA``, ``_alleles``, ``_locus``.
        ``LGT`` may be haploid post sex-ploidy, so extraction is ploidy-safe.
    :return: Array of ``{locus, alleles, is_hom, hap}`` (locus/alleles
        min-repped; ``hap`` is the int32 haplotype index or missing).
    """
    carried = hl.array(
        hl.set(hl.range(entry.LGT.ploidy).map(lambda i: entry.LGT[i])).filter(
            lambda a: a > 0
        )
    )
    return carried.map(
        lambda li: hl.bind(
            lambda mr: hl.struct(
                locus=mr.locus,
                alleles=mr.alleles,
                # is_hom describes the whole call: every alt record from one
                # call gets the same value (1/2 -> False on both). Haploid calls
                # are hom-var in Hail, and the predicate only tests whether the
                # call's allele indices are equal and non-zero, so working in
                # local allele space does not change the answer is_hom_var()
                # gives. Hemizygous alts get is_hom=True and count in n_homhom
                # (correctly cis: one haplotype). het-hom's unconditional-cis
                # rule leans on this flag; see KNOWN_ISSUES.md #2 for the one
                # gap.
                is_hom=entry.LGT.is_hom_var(),
                # ``find`` (not a fixed index) so it's ploidy-safe and yields
                # missing — never a wrong index — if LPGT doesn't carry li.
                hap=hl.or_missing(
                    entry.LPGT.phased,
                    hl.range(entry.LPGT.ploidy).find(lambda h: entry.LPGT[h] == li),
                ),
            ),
            hl.min_rep(
                entry._locus,
                hl.array([entry._alleles[0], entry._alleles[entry.LA[li]]]),
            ),
        )
        # SNP-only: classification keeps SNP-SNP pairs exclusively, so an indel
        # alt can never reach the output. Dropping it here (rather than in the
        # pair filter) keeps it out of the scan window and the checkpoint.
    ).filter(lambda r: hl.is_snp(r.alleles[0], r.alleles[1]))


def _classify_alt_pair(
    current_alt: hl.expr.StructExpression,
    prev_alt: hl.expr.StructExpression,
    same_phase_set: hl.expr.BooleanExpression,
) -> hl.expr.StructExpression:
    """
    Classify a carried-alt pair as het-het, hom-hom, or het-hom (all cis).

    Valid in local allele space: het/hom status and the haplotype index are the
    same as post-split. het-het needs both phased into the same phase set and
    haplotype to confirm cis; het-hom is unconditionally cis because a hom
    occupies both haplotypes. A het pair whose cis status is unconfirmed matches
    none and is dropped downstream — that covers genuinely trans pairs, but also
    cis pairs that are unphased, in different phase sets, or PID-missing, so
    het-het is a conservative undercount rather than an exact measure.

    :param current_alt: Current carried-alt record (with ``is_hom``, ``hap``).
    :param prev_alt: Previous carried-alt record (with ``is_hom``, ``hap``).
    :param same_phase_set: Whether the two entries share a PID (phasing block).
    :return: Struct with ``is_hethet``, ``is_homhom``, ``is_hethom`` booleans.
    """
    return hl.struct(
        is_hethet=(
            ~current_alt.is_hom
            & ~prev_alt.is_hom
            & hl.is_defined(current_alt.hap)
            & hl.is_defined(prev_alt.hap)
            & same_phase_set
            & (current_alt.hap == prev_alt.hap)
        ),
        is_homhom=current_alt.is_hom & prev_alt.is_hom,
        is_hethom=current_alt.is_hom != prev_alt.is_hom,
    )


# ---------------------------------------------------------------------------
# Pipeline steps
# ---------------------------------------------------------------------------


def _scan_for_candidates(ht: hl.Table, max_distance: int) -> hl.Table:
    """
    Annotate each non-ref entry with its window of nearby prior non-ref
    entries (``prev``), then keep rows where some sample has a candidate pair.

    :param ht: Localized Table; entries carry ``LGT`` plus the precomputed
        ``_alts``/``PID``/``adj`` used downstream.
    :param max_distance: Maximum bp distance between two SNVs.
    :return: Table with a ``_cands`` array of candidate carrier entries (each
        with its ``prev`` window), rows without candidates dropped, checkpointed
        to GCS. Narrowed before the checkpoint so the widened all-sample
        ``__entries`` array is never materialized.
    """
    # missing (not empty array) so the fold can tell "no non-ref seen yet" from
    # "one seen but now out of range". Window = array of (locus, entry) tuples.
    missing_window = hl.missing(
        hl.tarray(hl.ttuple(ht.locus.dtype, ht.__entries.dtype.element_type))
    )

    def _in_window(x):
        # NOTE: windows on the SITE locus, while the pipeline emits min-repped
        # loci — these can differ when a SNP is left-padded at a multiallelic
        # site, so a valid pair can be missed / an out-of-range pair emitted.
        # Known, un-fixed edge case; see v4/KNOWN_ISSUES.md.
        return (ht.locus.contig == x[0].contig) & (
            (ht.locus.position - x[0].position) <= max_distance
        )

    # Per-sample fold: maintain a sliding window of recent non-ref entries.
    scan_result = hl.scan.array_agg(
        lambda entry: hl.scan.fold(
            missing_window,
            lambda window: (
                hl.case()
                # Nothing seen yet, current is ref: still nothing.
                .when(hl.is_missing(window) & ~_carries_snv(entry), missing_window)
                # Nothing seen yet, current is non-ref: start the window.
                .when(hl.is_missing(window) & _carries_snv(entry), [(ht.locus, entry)])
                # Window has an in-range entry: prune to in-range, then append the
                # current entry if it is non-ref.
                .when(
                    hl.any(window.map(_in_window)),
                    hl.if_else(
                        _carries_snv(entry),
                        window.filter(_in_window).append((ht.locus, entry)),
                        window.filter(_in_window),
                    ),
                )
                # Window exists but all entries are now out of range: restart from
                # the current entry, or drop to missing if it is ref.
                .default(
                    hl.if_else(
                        _carries_snv(entry),
                        [(ht.locus, entry)],
                        missing_window,
                    )
                )
            ),
            # combop: merge per-partition windows (coalesce handles a missing side).
            lambda a, b: hl.coalesce(hl.flatten([a, b]), a, b),
        ),
        ht.__entries,
    )

    # First row has no scan_result (missing) -> fall back to per-sample missing.
    missing_entry = hl.range(ht.__cols.length()).map(lambda _: missing_window)
    ht = ht.select(
        __entries=hl.zip(ht.__entries, hl.coalesce(scan_result, missing_entry)).map(
            lambda pair: pair[0].annotate(
                prev=hl.or_missing(
                    _carries_snv(pair[0]),
                    hl.bind(
                        # Re-filter the window by the CURRENT row's locus,
                        # since the scan result reflects the window state after
                        # the previous row (which may be further than
                        # max_distance away).
                        lambda p: hl.or_missing(p.length() > 0, p),
                        pair[1].filter(_in_window),
                    ),
                )
            )
        ),
    )

    # Narrow to candidate carriers *before* materializing. The widened table has
    # one entry slot per sample, so checkpointing it whole would write mostly
    # non-candidates; ``_cands`` keeps only entries that have a window. The
    # length test is the same predicate as "some entry has a defined prev", so
    # this filter subsumes the row filter that used to precede the checkpoint.
    ht = ht.select(_cands=ht.__entries.filter(lambda entry: hl.is_defined(entry.prev)))
    ht = ht.filter(hl.len(ht._cands) > 0)
    return ht.checkpoint(hl.utils.new_temp_file("mnv_scan_candidates", "ht"))


def _classify_mnv_pairs(ht: hl.Table) -> hl.Table:
    """
    Explode to one row per (candidate carrier, classified alt pair).

    Carriers are exploded to one row each *before* their
    ``current_alt × prev_alt`` pairs are built (the cross product of the current
    entry's ``_alts`` and each windowed prev entry's ``_alts``), bounding per-task
    memory to one carrier's pairs rather than a per-row array across all samples
    (the off-heap hotspot) — so no classification-step repartition is needed at
    scale. Keeps SNP-SNP pairs matching a classification; a pair is adj iff both
    genotypes are.

    :param ht: Output of :func:`_scan_for_candidates`: a ``_cands`` array of
        candidate carrier entries, each carrying ``_alts``, ``PID``, ``adj`` and
        ``prev``.
    :return: Table with one ``_mnv`` struct per row.
    """

    def _build_pair_record(entry, prev_tuple, current_alt, prev_alt):
        # Loci/alleles come from the (min-repped) carried-alt records, not the
        # scan's site locus, so output keys join the release HT.
        prev_e = prev_tuple[1]
        return hl.bind(
            lambda gt_class: hl.struct(
                prev_locus=prev_alt.locus,
                prev_alleles=prev_alt.alleles,
                cur_locus=current_alt.locus,
                cur_alleles=current_alt.alleles,
                is_hethet=gt_class.is_hethet,
                is_homhom=gt_class.is_homhom,
                is_hethom=gt_class.is_hethom,
                adj=entry.adj & prev_e.adj,
            ),
            _classify_alt_pair(current_alt, prev_alt, entry.PID == prev_e.PID),
        )

    def _classify_entry(entry):
        """Classify all MNV alt pairs for one entry with a defined prev window."""
        # Build every (prev-window entry) x (current alt) x (prev alt) combo,
        # then the two hl.flatten calls collapse the nesting to one flat array of
        # pair records.
        return hl.flatten(
            # For each prior nearby non-ref entry in this sample's window
            # (prev_tuple = (site_locus, prev_entry)):
            entry.prev.map(
                lambda prev_tuple: hl.flatten(
                    # For each carried alt of the current entry:
                    entry._alts.map(
                        # cross it with each carried alt of the prev entry,
                        # building one pair record per combination.
                        lambda current_alt: prev_tuple[1]._alts.map(
                            lambda prev_alt: _build_pair_record(
                                entry, prev_tuple, current_alt, prev_alt
                            )
                        )
                    )
                )
            )
        ).filter(
            # Keep only pairs that matched a classification and are SNP-SNP.
            lambda p: (
                (p.is_hethet | p.is_homhom | p.is_hethom)
                & hl.is_snp(p.cur_alleles[0], p.cur_alleles[1])
                & hl.is_snp(p.prev_alleles[0], p.prev_alleles[1])
            )
        )

    # ``_cands`` was narrowed and checkpointed by :func:`_scan_for_candidates`;
    # explode straight from it rather than re-materializing.
    logger.info("Exploding candidate carriers...")
    ht = ht.explode("_cands")
    ht = ht.select(_mnv=_classify_entry(ht._cands))
    ht = ht.filter(hl.len(ht._mnv) > 0)
    return ht.explode("_mnv")


def _aggregate_mnv_pairs(ht: hl.Table) -> hl.Table:
    """
    Count occurrences per variant pair: raw (all carriers) and ``_adj`` (both
    genotypes pass adj), split into het-het / hom-hom / het-hom.

    :param ht: Table output of :func:`_classify_mnv_pairs` (one ``_mnv`` per row).
    :return: Table keyed by the later SNV of the pair (``locus``, ``alleles``)
        then the earlier (``prev_locus``, ``prev_alleles``), with the eight
        ``n_*`` counts and ``dist``, the min-repped-locus distance
        (``locus.position - prev_locus.position``). ``dist`` can fall outside
        ``[1, max_distance]`` when min-rep shifted a locus (see
        KNOWN_ISSUES.md #1); the post-write guard in ``main`` flags — but does
        not drop — such rows.
    """
    per_pair = ht.group_by(
        locus=ht._mnv.cur_locus,
        alleles=ht._mnv.cur_alleles,
        prev_locus=ht._mnv.prev_locus,
        prev_alleles=ht._mnv.prev_alleles,
    ).aggregate(
        n_hethet=hl.agg.count_where(ht._mnv.is_hethet),
        n_homhom=hl.agg.count_where(ht._mnv.is_homhom),
        n_hethom=hl.agg.count_where(ht._mnv.is_hethom),
        n_total=hl.agg.count(),
        n_hethet_adj=hl.agg.count_where(ht._mnv.is_hethet & ht._mnv.adj),
        n_homhom_adj=hl.agg.count_where(ht._mnv.is_homhom & ht._mnv.adj),
        n_hethom_adj=hl.agg.count_where(ht._mnv.is_hethom & ht._mnv.adj),
        n_total_adj=hl.agg.count_where(ht._mnv.adj),
    )

    return per_pair.annotate(
        dist=per_pair.locus.position - per_pair.prev_locus.position,
    )


def _parse_intervals(values: List[str]) -> List[Tuple[str, str]]:
    """
    Parse ``--intervals`` values into (gene_name, interval) pairs.

    :param values: List of strings in the form ``GENE_NAME=INTERVAL``, e.g.
        ``"PCSK9=chr1:55039447-55064852"``.
    :return: List of (gene_name, interval) tuples.
    """
    parsed = []
    for value in values:
        if "=" not in value:
            raise ValueError(
                f"Invalid --intervals value {value!r}; expected GENE_NAME=INTERVAL,"
                " e.g. 'PCSK9=chr1:55039447-55064852'."
            )
        name, interval = value.split("=", 1)
        # _interval_suffix drops an empty name silently (and a name that is
        # empty alongside others yields a stray "_"), so the run's scope would
        # not show up in the output path. Reject it rather than mislabel.
        if not name:
            raise ValueError(
                f"Invalid --intervals value {value!r}; GENE_NAME must not be empty."
            )
        parsed.append((name, interval))
    return parsed


def _normalize_chrom(value: str) -> str:
    """
    Normalize a chromosome argument to a ``chrN`` contig name.

    :param value: Chromosome as a number/letter (e.g. ``"21"``, ``"X"``) or an
        already-prefixed contig name (e.g. ``"chr21"``).
    :return: GRCh38 contig name, e.g. ``"chr21"`` or ``"chrX"``.
    """
    return value if value.lower().startswith("chr") else f"chr{value}"


def _per_alt_af(mt: hl.MatrixTable, release_ht: hl.Table) -> hl.Table:
    """
    Build the per-alt release-AF row array for the hom-alt hotfix: ``_alt_af[g-1]``
    is the release AF of global alt ``g``.

    :param mt: MatrixTable keyed by ``(locus, alleles)``.
    :param release_ht: Split release HT keyed by ``(locus, alleles)``; ``freq[0].AF``
        is the adj AF.
    :return: Table keyed by ``(locus, alleles)`` with an ``_alt_af`` array.
    """
    r = mt.rows().select()
    r = r.annotate(_ai=hl.range(1, hl.len(r.alleles)))
    r = r.explode("_ai")
    # min-rep so the key matches the min-repped release HT (else AF missing).
    r = r.annotate(_mr=hl.min_rep(r.locus, hl.array([r.alleles[0], r.alleles[r._ai]])))
    r = r.annotate(_af=release_ht[r._mr.locus, r._mr.alleles].freq[0].AF)
    return r.group_by(r.locus, r.alleles).aggregate(
        _alt_af=hl.sorted(hl.agg.collect((r._ai, r._af)), key=lambda t: t[0]).map(
            lambda t: t[1]
        )
    )


def _drop_all_lowqual_sites(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Drop sites where every alt is AS_lowqual — release-excluded, can't be a
    real MNV, and dropping them shrinks the table fed to the scan.

    The unsplit info HT (``AS_lowqual`` is ``array<bool>`` over alts) is re-read
    with the variant_data's partition intervals so the join needs no shuffle.

    :param mt: MatrixTable of variant data keyed by ``(locus, alleles)``.
    :return: ``mt`` with all-AS_lowqual sites removed.
    """
    intervals = mt._calculate_new_partitions(mt.n_partitions())
    info_ht = hl.read_table(get_info(split=False).path, _intervals=intervals)
    info_ht = info_ht.select(_all_lowqual=hl.all(info_ht.AS_lowqual))
    mt = mt.annotate_rows(
        _all_lowqual=hl.coalesce(info_ht[mt.row_key]._all_lowqual, False)
    )
    return mt.filter_rows(~mt._all_lowqual).drop("_all_lowqual")


# ---------------------------------------------------------------------------
# Top-level discovery + annotation
# ---------------------------------------------------------------------------


def discover_mnv(
    vds: hl.vds.VariantDataset,
    release_ht: hl.Table,
    max_distance: int = MAX_MNV_DISTANCE,
    entries_to_keep: List[str] = MNV_ENTRIES_TO_KEEP,
    classify_n_partitions: Optional[int] = None,
) -> hl.Table:
    """
    Run single-pass MNV discovery on an unsplit VDS.

    Stays in local allele space (LGT/LPGT/LA/PID). Pipeline: drop all-lowqual sites →
    adjust genotypes (het_non_ref flag → hom-alt hotfix → sex ploidy → adj,
    gnomAD QC's order) → localize → build ``_alts`` → :func:`_scan_for_candidates` →
    :func:`_classify_mnv_pairs` → :func:`_aggregate_mnv_pairs`.

    .. note::

        Scans a localized Table via ``hl.scan.array_agg``/``fold``, not an MT
        entry scan — the latter hits ``KeyError: 'agg_capability'`` in Hail's
        ``CSEAnalysisPass``. Sex-ploidy adjustment is applied
        (``adjusted_sex_ploidy_expr``) and the code is ploidy-safe, so
        ``--chr X``/``Y`` runs execute, but their output is incomplete for XY
        samples: non-PAR XY het calls are dropped before the scan, leaving only
        hom-hom pairs. See DOCUMENTATION.md §3 ("Sex chromosomes") for the
        measured rates and validation scope, and KNOWN_ISSUES.md #2 for the
        PAR-boundary case.

    :param vds: Unsplit gnomAD v4 VariantDataset.
    :param release_ht: Release sites HT supplying the per-alt AF the hom-alt
        depletion hotfix compares against.
    :param max_distance: Max bp between the two SNVs. Default 2 (codon frame).
    :param entries_to_keep: Entry fields to keep (global/post-split names).
    :param classify_n_partitions: If set, repartition before classification.
        Only pass at ``--intervals`` scale — see the CLI help for why.
    :return: Hail Table of MNV pairs with per-pair counts.
    """
    mt = vds.variant_data
    mt = _drop_all_lowqual_sites(mt)
    mt = mt.filter_rows(hl.any(lambda a: hl.is_snp(mt.alleles[0], a), mt.alleles[1:]))

    # entries_to_keep uses global names; L-prefix the ones the split remaps.
    pre_split_entries = entries_to_keep + ["LA"]
    pre_split_entries = [
        "L" + e if e in {"GT", "AD", "PL", "PGT"} else e for e in pre_split_entries
    ]
    mt = mt.select_entries(*pre_split_entries)

    # Row/col inputs the hotfix + sex-ploidy need (per-alt release AF, sex
    # karyotype, and the model flag that exempts samples from the hotfix).
    af_ht = _per_alt_af(mt, release_ht)
    mt = mt.select_rows(_alt_af=af_ht[mt.locus, mt.alleles]._alt_af)
    meta_ht = meta().ht()
    mt = mt.select_cols(
        sex_karyotype=meta_ht[mt.s].sex_imputation.sex_karyotype,
        fixed_homalt_model=meta_ht[mt.s].project_meta.fixed_homalt_model,
    )

    # Genotype adjustments in gnomAD QC's canonical order (see
    # gnomad_qc/v3/create_release/create_hgdp_tgp_subset.py):
    # het_non_ref flag -> hom-alt hotfix -> sex ploidy -> adj.
    mt = mt.annotate_entries(_het_non_ref=mt.LGT.is_het_non_ref())

    # hom-alt hotfix (local-space analog of hom_alt_depletion_fix): only het-ref
    # is affected, so hl.call(LGT[1], LGT[1]) is right where the stock
    # biallelic hl.call(1, 1) would not be at multiallelic sites.
    high_ab_het_ref = (
        mt.LGT.is_het_ref()
        & ~mt._het_non_ref
        & ~mt.fixed_homalt_model
        & (mt._alt_af[mt.LA[mt.LGT[1]] - 1] > HOM_ALT_FIX_AF_CUTOFF)
        & (hl.float64(mt.LAD[mt.LGT[1]]) / mt.DP > HOM_ALT_FIX_AB_CUTOFF)
    )
    mt = mt.annotate_entries(
        LGT=hl.if_else(
            high_ab_het_ref, hl.call(mt.LGT[1], mt.LGT[1]), mt.LGT, missing_false=True
        )
    )
    mt = mt.annotate_entries(
        LGT=adjusted_sex_ploidy_expr(mt.locus, mt.LGT, mt.sex_karyotype)
    )
    mt = mt.annotate_entries(adj=get_adj_expr(mt.LGT, mt.GQ, mt.DP, mt.LAD))

    # No densify: the scan treats a missing entry like 0/0 (both "not non-ref"),
    # so hom-ref entries stay sparse instead of materializing a 0/0 per sample
    # per site (verified equivalent).
    logger.info("Localizing and scanning for MNV candidates (single-pass)...")
    ht = mt._localize_entries("__entries", "__cols")
    ht = ht.annotate(
        __entries=ht.__entries.map(
            lambda entry: entry.annotate(
                _alts=hl.or_missing(
                    _is_nonref(entry),
                    _get_carried_alts(
                        entry.annotate(_alleles=ht.alleles, _locus=ht.locus)
                    ),
                )
            ).select("LGT", "PID", "adj", "_alts")
        )
    )

    ht = _scan_for_candidates(ht, max_distance)

    if classify_n_partitions is not None:
        logger.info(
            "Repartitioning to %d partitions before classification...",
            classify_n_partitions,
        )
        ht = ht.repartition(classify_n_partitions, shuffle=True)

    logger.info("Classifying MNV pairs...")
    ht = _classify_mnv_pairs(ht)

    return _aggregate_mnv_pairs(ht)


def annotate_mnv(mnv_ht: hl.Table, release_ht: hl.Table, vep_ht: hl.Table) -> hl.Table:
    """
    Add AC/AF/filters (exomes release) and VEP to both SNVs of each pair.

    Exomes release, not joint: the joint table lacks ``filters`` and
    per-data-type frequency arrays (only ``joint_freq``/``joint_faf``).

    .. note::

        Discovery runs on the VDS, so a discovered SNV can be absent from the
        release entirely; its annotation fields come back missing. Missing
        ``filters`` is not PASS — PASS is an *empty set* — so distinguish
        not-in-release (missing) from passing (empty) when consuming this table.

    :param mnv_ht: MNV Hail Table output from :func:`discover_mnv`.
    :param release_ht: Release sites HT supplying AC/AF/filters for both SNVs. Within
        one invocation this is the same table :func:`discover_mnv` is given, so the
        hotfix AF and the reported frequencies come from one release.
    :param vep_ht: VEP HT supplying transcript consequences for both SNVs.
    :return: MNV Table with AC, AF, filters, and VEP for both SNVs.
    """
    snv2_data = release_ht[mnv_ht.locus, mnv_ht.alleles]
    snv1_data = release_ht[mnv_ht.prev_locus, mnv_ht.prev_alleles]
    mnv_ht = mnv_ht.annotate(
        filters=snv2_data.filters,
        AC=snv2_data.freq[0].AC,
        AF=snv2_data.freq[0].AF,
        prev_filters=snv1_data.filters,
        prev_AC=snv1_data.freq[0].AC,
        prev_AF=snv1_data.freq[0].AF,
    )

    return mnv_ht.annotate(
        vep=vep_ht[mnv_ht.locus, mnv_ht.alleles].vep,
        prev_vep=vep_ht[mnv_ht.prev_locus, mnv_ht.prev_alleles].vep,
    )


def main(args: argparse.Namespace) -> None:
    """
    Execute MNV discovery and/or annotation based on parsed arguments.

    :param args: Parsed argparse namespace.
    """
    hl.init(
        log="discover_mnv.log",
        tmp_dir="gs://gnomad-tmp-4day/discover_mnv",
    )

    gene_names = None
    intervals = None
    chrom = None

    if args.chr is not None:
        chrom = _normalize_chrom(args.chr)
        intervals = [chrom]
    elif args.intervals is not None:
        parsed = (
            _parse_intervals(args.intervals) if args.intervals else DEFAULT_TEST_GENES
        )
        gene_names = [name for name, _ in parsed]
        intervals = [interval for _, interval in parsed]

    # A gene-interval subset is always a test run.
    test_enabled = args.test or args.intervals is not None

    release_ht = public_release("exomes").ht()

    if args.discover:
        if chrom is not None:
            logger.info("Starting MNV discovery on chromosome %s...", chrom)
        elif gene_names is not None:
            logger.info(
                "Subsetting to intervals: %s...",
                ", ".join(
                    f"{name} ({interval})"
                    for name, interval in zip(gene_names, intervals)
                ),
            )
        else:
            logger.info("Starting MNV discovery.")
        if test_enabled:
            logger.info("Writing output to the test bucket.")
        if args.ukb_only:
            logger.info("Using the UKB-only VDS subset.")
        start = time.time()

        vds_loader = get_gnomad_v4_ukb_vds if args.ukb_only else get_gnomad_v4_vds
        vds = vds_loader(
            filter_intervals=intervals,
            high_quality_only=args.high_quality_only,
            release_only=args.release_only,
        )

        # classify_n_partitions is --intervals-only; the CLI rejects it otherwise.
        mnv_ht = discover_mnv(
            vds,
            release_ht,
            max_distance=args.max_distance,
            entries_to_keep=MNV_ENTRIES_TO_KEEP,
            classify_n_partitions=args.classify_n_partitions,
        )

        discovery_resource = get_discovered_mnvs(
            test=test_enabled,
            interval_names=gene_names,
            ukb_only=args.ukb_only,
            chrom=chrom,
            suffix=args.output_suffix,
        )
        logger.info("Writing MNV discovery results to %s.", discovery_resource.path)
        mnv_ht.write(discovery_resource.path, overwrite=args.overwrite)

        # Guard: ``hl.min_rep`` can shift a SNP's locus when it sits inside a
        # padded multi-allelic ref, moving it outside the site-locus pairing
        # window. Read back the written (small, aggregated) HT and flag any
        # emitted pair whose distance is outside [1, max_distance]. Write those
        # rows to a sibling HT so they can be triaged later (see
        # v4/KNOWN_ISSUES.md); they are NOT dropped from the discovery output.
        # Written even when empty, so a clean rerun replaces any stale outliers
        # table from a previous run at the same path.
        written = discovery_resource.ht()
        outliers = written.filter(
            (written.dist < 1) | (written.dist > args.max_distance)
        )
        n_bad = outliers.count()
        outlier_path = discovery_resource.path.removesuffix(".ht") + ".dist_outliers.ht"
        outliers.write(outlier_path, overwrite=args.overwrite)
        if n_bad:
            logger.warning(
                "%d MNV pair(s) have dist outside [1, %d] — a SNP's locus likely"
                " shifted inside a padded multi-allelic ref during min-rep. Wrote"
                " them to %s for triage (see v4/KNOWN_ISSUES.md; note this guard"
                " only flags emitted pairs, not silently-missed ones).",
                n_bad,
                args.max_distance,
                outlier_path,
            )
        else:
            logger.info("No dist outliers (wrote empty %s).", outlier_path)

        elapsed = time.time() - start
        logger.info("Finished MNV discovery in %.1f seconds.", elapsed)

    if args.annotate:
        logger.info("Annotating MNV pairs with frequency and VEP data...")
        mnv_ht = get_discovered_mnvs(
            test=test_enabled,
            interval_names=gene_names,
            ukb_only=args.ukb_only,
            chrom=chrom,
            suffix=args.output_suffix,
        ).ht()
        vep_ht = get_vep(data_type="exomes").ht()
        mnv_ht = annotate_mnv(mnv_ht, release_ht, vep_ht)

        annotated_resource = get_annotated_mnvs(
            test=test_enabled,
            interval_names=gene_names,
            ukb_only=args.ukb_only,
            chrom=chrom,
            suffix=args.output_suffix,
        )
        logger.info("Writing annotated MNV results to %s.", annotated_resource.path)
        mnv_ht.write(annotated_resource.path, overwrite=args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Scan-based MNV discovery for gnomAD v4. Identifies pairs of biallelic"
            " SNVs within a configurable distance that co-occur on the same haplotype."
        ),
    )

    step_args = parser.add_argument_group(
        "Pipeline steps",
        "Specify which pipeline steps to run. At least one is required.",
    )
    step_args.add_argument(
        "--discover",
        help="Run MNV discovery and write the raw pair counts HT.",
        action="store_true",
    )
    step_args.add_argument(
        "--annotate",
        help=(
            "Annotate the discovery HT with AC, AF, filters (exomes release)"
            " and VEP consequences (v4 exomes), then write the annotated HT."
        ),
        action="store_true",
    )

    parser.add_argument(
        "--overwrite",
        help="Whether to overwrite existing output files.",
        action="store_true",
    )
    parser.add_argument(
        "--output-suffix",
        help=(
            "Free-form suffix appended last to the discovery/annotated output"
            " filenames (e.g. your username), so a run doesn't overwrite anyone"
            " else's output. E.g. '--output-suffix mwilson' writes"
            " '...mnv_discovery.mwilson.ht'."
        ),
        type=str,
        default=None,
    )
    parser.add_argument(
        "--test",
        help="Write output to the test bucket instead of the production bucket.",
        action="store_true",
    )
    default_test_genes_str = " ".join(
        f"{name}={interval}" for name, interval in DEFAULT_TEST_GENES
    )
    parser.add_argument(
        "--intervals",
        help=(
            "Subset discovery to one or more 'GENE_NAME=INTERVAL' pairs, e.g."
            " 'PCSK9=chr1:55039447-55064852 PCNT=chr21:46324141-46445769'. Implies"
            " --test; gene names are appended to output paths so runs don't"
            f" overwrite each other. Defaults to {default_test_genes_str} if given"
            " with no values. Cannot be combined with --chr."
        ),
        nargs="*",
        metavar="GENE_NAME=INTERVAL",
        default=None,
    )
    parser.add_argument(
        "--chr",
        help=(
            "Run discovery on a single whole chromosome. Pass a chromosome"
            " number/letter, e.g. '--chr 21' or '--chr X'; it is normalized to a"
            " 'chrN' contig name. Output gets a '.<chrom>' suffix so each"
            " chromosome lands at a distinct path (submit one job per chromosome)."
            " Writes to the production bucket unless --test is also given. Cannot"
            " be combined with --intervals."
        ),
        type=str,
        metavar="CHROM",
        default=None,
    )

    sample_filter_args = parser.add_argument_group(
        "Sample filtering",
        "Arguments controlling which samples are included (--discover only).",
    )
    sample_filter_args.add_argument(
        "--high-quality-only",
        help="Filter to only high-quality samples.",
        action="store_true",
    )
    sample_filter_args.add_argument(
        "--release-only",
        help="Filter to only release samples.",
        action="store_true",
    )
    sample_filter_args.add_argument(
        "--ukb-only",
        help=(
            "Run discovery on the UKB-only VDS subset instead of the full v4 VDS."
            " Output filenames get a '.ukb_only' component so they don't collide"
            " with full-cohort output."
        ),
        action="store_true",
    )

    detection_args = parser.add_argument_group(
        "Detection parameters",
        "Arguments controlling MNV detection behavior (--discover only).",
    )
    detection_args.add_argument(
        "--max-distance",
        help=(
            "Maximum distance in bp between two SNVs to consider as an MNV."
            f" Default is {MAX_MNV_DISTANCE}."
        ),
        type=int,
        default=MAX_MNV_DISTANCE,
    )
    detection_args.add_argument(
        "--classify-n-partitions",
        help=(
            "If set, repartition to this many partitions before classification."
            " Requires --intervals: it shuffles the fully-widened table,"
            " too costly for a full-cohort run. Pass it (e.g."
            f" {CLASSIFY_N_PARTITIONS}) if an --intervals run hits a Hail off-heap"
            " memory error during classification."
        ),
        type=int,
        default=None,
    )

    args = parser.parse_args()
    if not args.discover and not args.annotate:
        parser.error("At least one of --discover or --annotate is required.")
    if args.chr is not None and args.intervals is not None:
        parser.error(
            "--chr and --intervals both scope which intervals to load; pass only"
            " one."
        )
    if args.classify_n_partitions is not None and args.intervals is None:
        parser.error(
            "--classify-n-partitions requires --intervals: the repartition"
            " shuffles the fully-widened table, too costly for --chr or"
            " full-cohort runs."
        )
    main(args)
