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
    get_gnomad_v4_ukb_vds,
    get_gnomad_v4_vds,
    mnv_annotated,
    mnv_discovery,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("discover_mnv")
logger.setLevel(logging.INFO)

MAX_MNV_DISTANCE = 2

HOM_ALT_FIX_AF_CUTOFF = 0.01
"""AF above which a high-AB het-ref call is reclassified to hom-alt.

Part of the GATK hom-alt depletion hotfix (see :func:`discover_mnv`); mirrors the
gnomAD default (``gnomad_qc.v3.utils.hom_alt_depletion_fix``).
"""

HOM_ALT_FIX_AB_CUTOFF = 0.9
"""Allele balance above which a common-variant het-ref call is treated as hom-alt."""

CLASSIFY_N_PARTITIONS = 500
"""Suggested number of partitions to repartition to before classification.

Repartitioning is off by default; pass ``--classify-n-partitions`` to enable it
(this is a reasonable starting value). Each row carries an entry per sample in
the cohort, so the classification step (which builds nested per-entry MNV-pair
arrays) is memory-heavy per row. Spreading the (already scan-filtered, typically
small) row count across more partitions keeps any single task from exceeding
Hail's off-heap memory limit.
"""

DEFAULT_TEST_GENES = [
    # PCNT on GRCh38.
    ("PCNT", "chr21:46324141-46445769"),
]
"""Default gene(s)/interval(s) used when ``--intervals`` is given with no values."""


# ---------------------------------------------------------------------------
# Expression helpers
# ---------------------------------------------------------------------------


def _is_nonref(e: hl.expr.StructExpression) -> hl.expr.BooleanExpression:
    """Check if an entry has a defined non-ref local genotype.

    :param e: Entry expression with an ``LGT`` field.
    :return: Boolean expression.
    """
    return hl.is_defined(e.LGT) & e.LGT.is_non_ref()


def _get_carried_alts(entry: hl.expr.StructExpression) -> hl.expr.ArrayExpression:
    """Build one haplotype-resolved record per distinct carried alt allele.

    Generalizes a single biallelic representation to an array so het_non_ref
    (``1/2``) carriers — which hold two different alt alleles on opposite
    haplotypes — contribute both alts. het-ref (``0/1``) and hom-var (``1/1``)
    each yield a single record, identical to the previous single-struct logic;
    ``1/2`` yields two.

    Each record maps the local alt index to the global allele index via LA and
    carries the haplotype the alt sits on (for cis het-het classification).

    .. note::

        Called once per non-ref entry at localization time (see
        :func:`discover_mnv`), not per candidate pair, so the full (possibly
        multiallelic) row ``alleles`` array never needs to be carried on every
        sample's entry through the scan window.

    :param entry: Entry expression with ``LGT``, ``LPGT``, ``LA``, and
        ``_alleles`` fields. ``LGT`` may be haploid (post sex-ploidy adjustment
        on sex chromosomes), so allele extraction is ploidy-safe.
    :return: Array of structs, one per carried alt, each with ``alleles``
        (array<str>), ``is_snp`` (bool), ``is_hom`` (bool), and ``hap``
        (int32 or missing).
    """
    # Distinct carried alt local indices, ploidy-safe (LGT may be haploid on
    # non-PAR X/Y after sex-ploidy adjustment).
    carried = hl.array(
        hl.set(hl.range(entry.LGT.ploidy).map(lambda i: entry.LGT[i])).filter(
            lambda a: a > 0
        )
    )
    return carried.map(
        lambda li: hl.bind(
            lambda alt: hl.struct(
                alleles=hl.array([entry._alleles[0], alt]),
                is_snp=hl.is_snp(entry._alleles[0], alt),
                is_hom=entry.LGT.is_hom_var(),
                # Haplotype index carrying this alt: the position h where
                # LPGT[h] == li. ``find`` over the ploidy range is ploidy-safe
                # and yields missing when LPGT doesn't carry li (e.g. a GT/PGT
                # allele mismatch), never a wrong index.
                hap=hl.or_missing(
                    entry.LPGT.phased,
                    hl.range(entry.LPGT.ploidy).find(lambda h: entry.LPGT[h] == li),
                ),
            ),
            entry._alleles[entry.LA[li]],
        )
    )


def _classify_alt_pair(
    ca: hl.expr.StructExpression,
    pa: hl.expr.StructExpression,
    same_phase_set: hl.expr.BooleanExpression,
) -> hl.expr.StructExpression:
    """Classify a pair of carried alts as het-het, hom-hom, or het-hom.

    Operates on the per-alt records built by :func:`_get_carried_alts` (each
    with ``is_hom`` and ``hap``), which is valid in local allele space:
    het/hom status is the same in local and global representations, and the
    haplotype index is always relative to local ref/alt.

    Classification rules:

    - **het-het**: Both alts are het (not hom), both have a defined haplotype
      (phased), share the same phase set (``same_phase_set``), and sit on the
      same haplotype (``ca.hap == pa.hap``), confirming they are in cis. This
      also covers het_non_ref (``1/2``), where each alt is its own het record on
      its own haplotype.
    - **hom-hom**: Both alts are homozygous.
    - **het-hom**: Exactly one alt is homozygous. A hom occupies both
      haplotypes, so the pair is always in cis.

    :param ca: Current carried-alt record (with ``is_hom``, ``hap``).
    :param pa: Previous carried-alt record (with ``is_hom``, ``hap``).
    :param same_phase_set: Whether the two entries share a PID (phasing block).
    :return: Struct with ``is_hethet``, ``is_homhom``, ``is_hethom`` booleans.
    """
    return hl.struct(
        is_hethet=(
            ~ca.is_hom
            & ~pa.is_hom
            & hl.is_defined(ca.hap)
            & hl.is_defined(pa.hap)
            & same_phase_set
            & (ca.hap == pa.hap)
        ),
        is_homhom=ca.is_hom & pa.is_hom,
        is_hethom=ca.is_hom != pa.is_hom,
    )


# ---------------------------------------------------------------------------
# Pipeline steps
# ---------------------------------------------------------------------------


def _scan_for_candidates(ht: hl.Table, max_distance: int) -> hl.Table:
    """Scan a localized Table for MNV candidates using a sliding window.

    For each sample, tracks a window of recent non-ref entries using
    ``hl.scan.fold`` with ``hl.case()`` branching. The window stores
    ``(locus, entry)`` tuples and uses missing state to represent "no recent
    non-ref entries seen". After the scan, each entry is annotated with a
    ``prev`` field containing the filtered window of nearby previous non-ref
    variants. The table is filtered to rows where at least one sample has a
    candidate MNV pair.

    :param ht: Localized Table (from ``mt._localize_entries``). Entries must
        contain ``LGT`` (and the precomputed ``_alts`` array plus ``PID``/``adj``,
        carried opaquely through the window for use in classification).
    :param max_distance: Maximum bp distance between two SNVs.
    :return: Filtered localized Table with ``prev`` annotated on each entry,
        checkpointed to a temporary path (spilled to GCS rather than held in
        memory, since the widened one-entry-per-sample table can be large).
    """
    # Typed missing value representing "no window state". The window is an
    # array of (locus, entry) tuples. Using missing (rather than empty array)
    # lets the fold distinguish "never seen a non-ref entry" from "saw one
    # but it fell out of range".
    missing_prev = hl.missing(
        hl.tarray(hl.ttuple(ht.locus.dtype, ht.__entries.dtype.element_type))
    )

    def _in_window(x):
        return (ht.locus.contig == x[0].contig) & (
            (ht.locus.position - x[0].position) <= max_distance
        )

    # Per-sample fold: maintain a sliding window of recent non-ref entries.
    scan_result = hl.scan.array_agg(
        lambda e: hl.scan.fold(
            missing_prev,
            lambda acc: (
                hl.case()
                .when(hl.is_missing(acc) & ~_is_nonref(e), missing_prev)
                .when(hl.is_missing(acc) & _is_nonref(e), [(ht.locus, e)])
                .when(
                    hl.any(acc.map(_in_window)),
                    hl.if_else(
                        _is_nonref(e),
                        acc.filter(_in_window).append((ht.locus, e)),
                        acc.filter(_in_window),
                    ),
                )
                .default(
                    hl.if_else(
                        _is_nonref(e),
                        [(ht.locus, e)],
                        missing_prev,
                    )
                )
            ),
            # combop: merge windows from different partitions. hl.flatten
            # concatenates two defined arrays; hl.coalesce handles the case
            # where one or both sides are missing (no state in that partition).
            lambda a, b: hl.coalesce(hl.flatten([a, b]), a, b),
        ),
        ht.__entries,
    )

    # Zip entries with scan results, annotate each entry with its prev window.
    # On the very first row, scan_result is missing (no previous rows), so
    # coalesce provides a per-sample array of missing windows as the fallback.
    missing_entry = hl.range(ht.__cols.length()).map(lambda _: missing_prev)
    ht = ht.select(
        __entries=hl.zip(ht.__entries, hl.coalesce(scan_result, missing_entry)).map(
            lambda pair: pair[0].annotate(
                prev=hl.or_missing(
                    _is_nonref(pair[0]),
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

    # Checkpoint (rather than .cache()) the widened, scan-filtered table: it
    # carries one entry per cohort sample, so holding it in executor memory
    # risks OOM at chromosome/full-cohort scale. Spilling to GCS trades a
    # write/read round-trip for bounded memory.
    ht = ht.filter(hl.any(ht.__entries.map(lambda x: hl.is_defined(x.prev))))
    return ht.checkpoint(hl.utils.new_temp_file("mnv_scan_candidates", "ht"))


def _classify_mnv_pairs(ht: hl.Table) -> hl.Table:
    """Classify MNV pairs from scan results and explode into one row per pair.

    For each entry with a defined ``prev`` window, classifies every
    ``cur._alts × prev._alts`` alt pair using :func:`_classify_alt_pair` and
    reads the biallelic alleles from each carried-alt record (built once per
    non-ref entry at localization time in :func:`discover_mnv`). het_non_ref
    (``1/2``) entries carry two alts, so a single entry can contribute multiple
    distinct pairs. A pair is adj-passing when both constituent entry genotypes
    pass adj. Filters to pairs where both carried alleles are SNPs and at least
    one classification matches.

    :param ht: Localized Table with ``prev`` annotated on entries (output of
        :func:`_scan_for_candidates`). Entries must carry ``_alts`` (array of
        carried-alt records), ``PID``, and ``adj``.
    :return: Table with one ``_mnv`` struct per row (after checkpoint + explode).
    """
    _mnv_pair_type = hl.tstruct(
        prev_locus=hl.tlocus("GRCh38"),
        prev_alleles=hl.tarray(hl.tstr),
        cur_alleles=hl.tarray(hl.tstr),
        is_hethet=hl.tbool,
        is_homhom=hl.tbool,
        is_hethom=hl.tbool,
        adj=hl.tbool,
    )

    def _build_pair_record(entry, prev_tuple, ca, pa):
        """Build an MNV pair struct for one (current alt, previous alt) pair."""
        prev_e = prev_tuple[1]
        return hl.bind(
            lambda gt_class: hl.struct(
                prev_locus=prev_tuple[0],
                prev_alleles=pa.alleles,
                cur_alleles=ca.alleles,
                is_hethet=gt_class.is_hethet,
                is_homhom=gt_class.is_homhom,
                is_hethom=gt_class.is_hethom,
                adj=entry.adj & prev_e.adj,
            ),
            _classify_alt_pair(ca, pa, entry.PID == prev_e.PID),
        )

    def _classify_entry(entry):
        """Classify all MNV alt pairs for one entry with a defined prev window."""
        return hl.flatten(
            entry.prev.map(
                lambda prev_tuple: hl.flatten(
                    entry._alts.map(
                        lambda ca: prev_tuple[1]._alts.map(
                            lambda pa: _build_pair_record(entry, prev_tuple, ca, pa)
                        )
                    )
                )
            )
        ).filter(
            lambda p: (
                (p.is_hethet | p.is_homhom | p.is_hethom)
                & hl.is_snp(p.cur_alleles[0], p.cur_alleles[1])
                & hl.is_snp(p.prev_alleles[0], p.prev_alleles[1])
            )
        )

    # Flatten classified pairs across all samples into one array per row.
    _mnv = hl.flatten(
        ht.__entries.map(
            lambda entry: hl.if_else(
                hl.is_defined(entry.prev),
                _classify_entry(entry),
                hl.empty_array(_mnv_pair_type),
            )
        )
    )
    ht = ht.select(_mnv=_mnv)

    # Filter, checkpoint, and explode to one row per MNV pair.
    ht = ht.filter(hl.len(ht._mnv) > 0)
    ht = ht.checkpoint(hl.utils.new_temp_file("mnv_scan_results", "ht"))
    logger.info("Scan checkpoint complete, exploding MNV pairs...")
    return ht.explode("_mnv")


def _aggregate_mnv_pairs(ht: hl.Table) -> hl.Table:
    """Aggregate MNV pairs per variant pair across all samples.

    Groups by biallelic (locus, alleles) for both the current and previous SNV,
    then counts het-het, hom-hom, het-hom, and total occurrences. Each count has
    a raw variant (all carriers) and an ``_adj`` variant restricted to pairs
    where both constituent genotypes pass adj.

    :param ht: Table with one ``_mnv`` struct per row (output of
        :func:`_classify_mnv_pairs`).
    :return: Aggregated Table keyed by (locus, alleles, prev_locus, prev_alleles).
    """
    per_pair = ht.group_by(
        locus=ht.locus,
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
    """Parse ``--intervals`` values into (gene_name, interval) pairs.

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
        parsed.append((name, interval))
    return parsed


def _normalize_chrom(value: str) -> str:
    """Normalize a chromosome argument to a ``chrN`` contig name.

    :param value: Chromosome as a number/letter (e.g. ``"21"``, ``"X"``) or an
        already-prefixed contig name (e.g. ``"chr21"``).
    :return: GRCh38 contig name, e.g. ``"chr21"`` or ``"chrX"``.
    """
    return value if value.lower().startswith("chr") else f"chr{value}"


def _per_alt_af(mt: hl.MatrixTable, freq_ht: hl.Table) -> hl.Table:
    """Build a per-alt released-AF row array for the hom-alt hotfix.

    Returns a Table keyed by ``(locus, alleles)`` with an ``_alt_af`` array
    holding the release AF of each biallelic ``[ref, alt_i]`` (ordered by alt
    index, so ``_alt_af[g - 1]`` is the AF of global alt allele ``g``).

    .. note::

        Built by explode → join → group rather than a join inside an
        ``hl.range(...).map(...)`` lambda. A table join whose key depends on the
        map loop variable produces IR that Hail's renderer cannot lower
        (``KeyError`` in ``CSEAnalysisPass`` once the row feeds entry
        expressions carried through the scan). SNP alts are min-repped and hit
        the split release HT; a non-min-repped indel alt may miss (AF → missing
        → hotfix skipped), which is safe and irrelevant to SNP-only MNV output.

    :param mt: MatrixTable keyed by ``(locus, alleles)``.
    :param freq_ht: Split release HT keyed by ``(locus, alleles)`` with a
        ``freq`` array (``freq[0].AF`` is the adj AF).
    :return: Table keyed by ``(locus, alleles)`` with an ``_alt_af`` array.
    """
    r = mt.rows().select()
    r = r.annotate(_ai=hl.range(1, hl.len(r.alleles)))
    r = r.explode("_ai")
    r = r.annotate(
        _af=freq_ht[r.locus, hl.array([r.alleles[0], r.alleles[r._ai]])].freq[0].AF
    )
    return r.group_by(r.locus, r.alleles).aggregate(
        _alt_af=hl.sorted(hl.agg.collect((r._ai, r._af)), key=lambda t: t[0]).map(
            lambda t: t[1]
        )
    )


def _drop_all_lowqual_sites(mt: hl.MatrixTable) -> hl.MatrixTable:
    """Drop unsplit sites where every alt allele is AS_lowqual.

    Such sites are excluded from the gnomAD release and cannot form a real MNV,
    so removing them before the scan shrinks the widened one-entry-per-sample
    table. This is a row-only filter — no genotypes are touched.

    The unsplit info HT (``AS_lowqual`` is an ``array<bool>`` over alts) is
    re-read with the variant_data's own partition intervals so the row join is
    co-partitioned and needs no shuffle (the idiom the v4 VDS loader itself
    uses; see ``gnomad_qc.v4.resources.basics``).

    :param mt: variant_data MatrixTable keyed by ``(locus, alleles)``.
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
    max_distance: int = MAX_MNV_DISTANCE,
    entries_to_keep: List[str] = MNV_ENTRIES_TO_KEEP,
    classify_n_partitions: Optional[int] = None,
) -> hl.Table:
    """Run single-pass MNV discovery on an unsplit VDS.

    Works in local allele space (LGT/LPGT/LA/PID), extracting biallelic alleles
    via LA indexing, so no separate unlocalize + split step is needed. Scans
    across rows using ``hl.scan.fold`` with ``hl.case()`` branching to track a
    window of recent non-ref entries per sample.

    Pipeline:

    0. **Drop junk sites + adjust genotypes** (in-line, before localizing):
       drop sites where every alt is AS_lowqual (:func:`_drop_all_lowqual_sites`),
       then adjust genotypes in gnomAD QC's canonical order: flag het_non_ref →
       hom-alt depletion hotfix (high-AB het-ref at common variants → hom-alt,
       skipping het_non_ref and fixed-hom-alt-model samples) → sex-ploidy
       adjustment → adj annotation.
    1. **Scan** (:func:`_scan_for_candidates`): Localize (no densify — the scan
       treats missing and 0/0 entries alike), build each non-ref entry's
       per-carried-alt records (``_alts``), scan
       using ``hl.scan.fold`` tracking a sliding window of non-ref entries per
       sample. Filter to rows with candidate MNV pairs.
    2. **Classify** (:func:`_classify_mnv_pairs`): For each candidate, classify
       every ``cur._alts × prev._alts`` alt pair (het-het / hom-hom / het-hom)
       and read the precomputed biallelic alleles for the aggregation key.
    3. **Aggregate** (:func:`_aggregate_mnv_pairs`): Group by variant pair
       across all samples, emitting raw and ``_adj`` counts.

    .. note::

        Sex-ploidy adjustment is a no-op on autosomes; on ``--chr X/Y`` it makes
        adj ploidy-correct, but hemizygous cis-MNV classification is not
        specially modeled (the source paper restricted to autosomes).

    .. note::

        Uses ``_localize_entries`` + ``hl.scan.array_agg`` + ``hl.scan.fold``
        instead of MT entry scans. Direct scans in MT entry context cause
        ``KeyError: 'agg_capability'`` in Hail's ``CSEAnalysisPass``. The
        localized pattern (from gnomAD's ``compute_last_ref_block_end``) works.

    :param vds: Unsplit gnomAD v4 VariantDataset.
    :param max_distance: Maximum bp distance between SNVs to consider as an MNV.
        Default is 2 (codon reading frame).
    :param entries_to_keep: Entry fields needed for MNV discovery.
    :param classify_n_partitions: If set, repartition (with a shuffle) to this
        many partitions before the classification step. Default is None (no
        repartition).

        .. warning::

            This triggers a shuffle of the fully-widened, one-entry-per-sample
            table. Safe for small (e.g. ``--test``) row counts, but likely too
            expensive/memory-heavy for a full-cohort, genome/exome-wide run —
            leave unset there.
    :return: Hail Table of MNV pairs with per-pair counts.
    """
    mt = vds.variant_data

    # Drop sites where every alt is AS_lowqual (release-excluded, can't form a
    # real MNV) to shrink the widened table fed to the scan.
    mt = _drop_all_lowqual_sites(mt)

    # Filter to rows with at least one SNP alt allele.
    mt = mt.filter_rows(hl.any(lambda a: hl.is_snp(mt.alleles[0], a), mt.alleles[1:]))

    # Select entries needed for classification (pre-split names) + LA.
    pre_split_entries = entries_to_keep + ["LA"]
    pre_split_entries = [
        "L" + e if e in {"GT", "AD", "PL", "PGT"} else e for e in pre_split_entries
    ]
    mt = mt.select_entries(*pre_split_entries)

    # Per-alt released AF (indexed by alt, i.e. global allele index - 1) for the
    # hom-alt depletion hotfix, via :func:`_per_alt_af`.
    freq_ht = public_release("exomes").ht()
    af_ht = _per_alt_af(mt, freq_ht)
    mt = mt.select_rows(_alt_af=af_ht[mt.locus, mt.alleles]._alt_af)

    # Sample sex karyotype (sex-ploidy adjustment) and the hom-alt model flag
    # (samples with the fixed model don't need the depletion hotfix).
    meta_ht = meta().ht()
    mt = mt.select_cols(
        sex_karyotype=meta_ht[mt.s].sex_imputation.sex_karyotype,
        fixed_homalt_model=meta_ht[mt.s].project_meta.fixed_homalt_model,
    )

    # --- Genotype adjustments, in gnomAD QC's canonical order (see
    # gnomad_qc/v3/create_release/create_hgdp_tgp_subset.py): mark het_non_ref ->
    # hom-alt depletion hotfix -> adjust sex ploidy -> annotate adj. ---

    # (1) Flag het_non_ref (1/2) before any GT edit; the hotfix must skip these.
    mt = mt.annotate_entries(_het_non_ref=mt.LGT.is_het_non_ref())

    # (2) Hom-alt depletion hotfix: reclassify a high-AB het-ref call at a common
    # variant to hom-alt. Local-space analog of
    # gnomad_qc.v3.utils.hom_alt_depletion_fix, whose biallelic hl.call(1, 1)
    # would be wrong here — but only het-ref calls are ever affected (het_non_ref
    # and homs are skipped), so hl.call of the single carried local alt (LGT[1])
    # is correct.
    high_ab_het_ref = (
        mt.LGT.is_het_ref()
        & ~mt._het_non_ref
        & ~mt.fixed_homalt_model
        & (mt._alt_af[mt.LA[mt.LGT[1]] - 1] > HOM_ALT_FIX_AF_CUTOFF)
        & (mt.LAD[mt.LGT[1]] / mt.DP > HOM_ALT_FIX_AB_CUTOFF)
    )
    mt = mt.annotate_entries(
        LGT=hl.if_else(
            high_ab_het_ref,
            hl.call(mt.LGT[1], mt.LGT[1]),
            mt.LGT,
            missing_false=True,
        )
    )

    # (3) Adjust sex ploidy on the hotfixed GT (no-op on autosomes; haploidizes
    # non-PAR X/Y in XY samples, sets XX-on-Y to missing).
    mt = mt.annotate_entries(
        LGT=adjusted_sex_ploidy_expr(mt.locus, mt.LGT, mt.sex_karyotype)
    )

    # (4) Annotate adj on the fully adjusted GT (ploidy-aware DP thresholds).
    mt = mt.annotate_entries(adj=get_adj_expr(mt.LGT, mt.GQ, mt.DP, mt.LAD))

    # Localize to a Table for scan-based processing. No unfilter_entries /
    # densify step is needed: the scan treats a missing entry identically to
    # 0/0 (both "not non-ref"), so leaving hom-ref entries sparse yields the
    # same MNV results while avoiding materializing a 0/0 for every sample at
    # every site (verified equivalent on synthetic data).
    logger.info("Localizing and scanning for MNV candidates (single-pass)...")
    ht = mt._localize_entries("__entries", "__cols")
    # Build per-carried-alt records once per non-ref entry (guarded to non-ref,
    # missing otherwise), reading the row alleles here so the multiallelic
    # alleles array isn't duplicated onto every sample's entry and carried
    # through the scan window. Keep only LGT (scan non-ref test), PID + adj
    # (classification), and _alts; drop the rest.
    ht = ht.annotate(
        __entries=ht.__entries.map(
            lambda e: e.annotate(
                _alts=hl.or_missing(
                    _is_nonref(e),
                    _get_carried_alts(e.annotate(_alleles=ht.alleles)),
                )
            ).select("LGT", "PID", "adj", "_alts")
        )
    )

    # Scan → repartition (test mode only) → classify → aggregate.
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


def annotate_mnv(mnv_ht: hl.Table) -> hl.Table:
    """Annotate MNV pairs with frequency and VEP data for both SNVs.

    Joins AC/AF/filters from the gnomAD v4 exomes release
    (``gnomad.resources.grch38.gnomad.public_release("exomes")``) and VEP consequences
    from the v4 exomes VEP HT (``gnomad_qc.v4.resources.annotations.get_vep``) onto
    each variant pair.

    .. note::

        Uses the exomes release (not joint) because the joint release table only
        contains ``joint_freq``/``joint_faf`` fields — it lacks ``filters`` and
        per-data-type frequency arrays.

    :param mnv_ht: MNV Hail Table output from :func:`discover_mnv`.
    :return: MNV Table annotated with AC, AF, filters, and VEP for both SNVs.
    """
    # --- Frequency and filter annotations from exomes release. ---
    freq_ht = public_release("exomes").ht()

    # Annotate SNV2 (current variant).
    snv2_data = freq_ht[mnv_ht.locus, mnv_ht.alleles]
    mnv_ht = mnv_ht.annotate(
        filters=snv2_data.filters,
        AC=snv2_data.freq[0].AC,
        AF=snv2_data.freq[0].AF,
    )

    # Annotate SNV1 (previous variant).
    snv1_data = freq_ht[mnv_ht.prev_locus, mnv_ht.prev_alleles]
    mnv_ht = mnv_ht.annotate(
        prev_filters=snv1_data.filters,
        prev_AC=snv1_data.freq[0].AC,
        prev_AF=snv1_data.freq[0].AF,
    )

    # --- VEP annotations from v4 exomes. ---
    vep_ht = get_vep(data_type="exomes").ht()

    # Annotate SNV2 (current variant).
    mnv_ht = mnv_ht.annotate(vep=vep_ht[mnv_ht.locus, mnv_ht.alleles].vep)

    # Annotate SNV1 (previous variant).
    mnv_ht = mnv_ht.annotate(
        prev_vep=vep_ht[mnv_ht.prev_locus, mnv_ht.prev_alleles].vep
    )

    return mnv_ht


def main(args: argparse.Namespace) -> None:
    """Execute MNV discovery and/or annotation based on parsed arguments.

    :param args: Parsed argparse namespace.
    """
    hl.init(
        log="discover_mnv.log",
        tmp_dir="gs://gnomad-tmp-4day/discover_mnv",
    )

    # Interval scope is set by --chr (a whole chromosome) or --intervals (gene
    # intervals); ``args.intervals`` is None when omitted, [] when given with no
    # values (use the default gene(s)), or a list of "GENE_NAME=INTERVAL" strings.
    gene_names: Optional[List[str]] = None
    intervals: Optional[List[str]] = None
    chrom: Optional[str] = None
    if args.chr is not None:
        chrom = _normalize_chrom(args.chr)
        intervals = [chrom]
    elif args.intervals is not None:
        parsed: List[Tuple[str, str]] = (
            _parse_intervals(args.intervals) if args.intervals else DEFAULT_TEST_GENES
        )
        gene_names = [name for name, _ in parsed]
        intervals = [interval for _, interval in parsed]

    # A gene-interval subset is always a test run, so --intervals implies the
    # test bucket even without --test.
    test_enabled = args.test or args.intervals is not None

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

        # --- Load unsplit VDS, optionally filtering to test intervals. ---
        vds_loader = get_gnomad_v4_ukb_vds if args.ukb_only else get_gnomad_v4_vds
        vds = vds_loader(
            filter_intervals=intervals,
            high_quality_only=args.high_quality_only,
            release_only=args.release_only,
        )

        # Repartitioning before classification shuffles the fully-widened
        # (one-entry-per-sample) table, so it's only safe at --intervals (gene
        # interval) scale. A --chr run is a spatial subset of the full run and,
        # like it, relies on native partitioning instead of a shuffle.
        mnv_ht = discover_mnv(
            vds,
            max_distance=args.max_distance,
            entries_to_keep=MNV_ENTRIES_TO_KEEP,
            classify_n_partitions=(
                args.classify_n_partitions if args.intervals is not None else None
            ),
        )

        discovery_resource = mnv_discovery(
            test=test_enabled,
            interval_names=gene_names,
            ukb_only=args.ukb_only,
            chrom=chrom,
            suffix=args.output_suffix,
        )
        logger.info("Writing MNV discovery results to %s.", discovery_resource.path)
        mnv_ht.write(discovery_resource.path, overwrite=args.overwrite)

        elapsed = time.time() - start
        logger.info("Finished MNV discovery in %.1f seconds.", elapsed)

    if args.annotate:
        logger.info("Annotating MNV pairs with frequency and VEP data...")
        mnv_ht = mnv_discovery(
            test=test_enabled,
            interval_names=gene_names,
            ukb_only=args.ukb_only,
            chrom=chrom,
            suffix=args.output_suffix,
        ).ht()
        mnv_ht = annotate_mnv(mnv_ht)

        annotated_resource = mnv_annotated(
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
            "Annotate the discovery HT with AC, AF, filters (joint release)"
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
            "Subset discovery to one or more gene intervals (implies --test). Pass"
            " one or more 'GENE_NAME=INTERVAL' pairs, e.g. '--intervals"
            " PCSK9=chr1:55039447-55064852 PCNT=chr21:46324141-46445769', to run"
            " different genes in separate job submissions or together in one. Gene"
            " names are appended to output filenames so runs don't overwrite each"
            f" other. If given with no values, defaults to {default_test_genes_str}."
            " Cannot be combined with --chr."
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
            "If set, repartition to this many partitions before the classification"
            " step, which builds nested per-entry MNV-pair arrays and is"
            " memory-heavy per row (each row carries an entry per cohort sample)."
            " Only applied when --intervals is set, since it shuffles the"
            " fully-widened table and is too expensive for a full-cohort run. Left"
            " unset there is no repartition; pass it (e.g."
            f" {CLASSIFY_N_PARTITIONS}) if you hit a Hail off-heap memory error"
            " during an --intervals run's classification step."
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
    main(args)
