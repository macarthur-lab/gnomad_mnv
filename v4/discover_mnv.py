"""Single-pass scan-based MNV discovery for gnomAD v4.

Uses ``hl.vds.lgt_to_gt()`` to convert local genotypes to global allele indices
inline, eliminating the need for a separate unlocalize + split step. The scan
operates on the unsplit sparse VDS variant_data, tracking a window of recent
non-ref entries per sample using ``hl.scan.fold`` with ``hl.case()`` branching.

For each sample, finds pairs of SNVs within ``max_distance`` bp that co-occur
on the same haplotype. Variant pairs are classified as:

    - het-het: both het-ref, both phased, same PID, same phase orientation
    - hom-hom: both hom-var
    - het-hom: one het-ref and one hom-var (either direction)

Results are aggregated per variant pair across all samples and written as a
Hail Table.

Usage::

    python -m v4.discover_mnv --discover --overwrite
    python -m v4.discover_mnv --annotate
    python -m v4.discover_mnv --discover --annotate --test
"""

import argparse
import logging
import time
from typing import List

import hail as hl

from v4.resources import (
    MNV_ENTRIES_TO_KEEP,
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
# PCSK9 on GRCh38 (chr1:55,039,447-55,064,852).
TEST_INTERVAL = "chr1:55039447-55064852"


# ---------------------------------------------------------------------------
# Expression helpers
# ---------------------------------------------------------------------------


def _is_nonref(e: hl.expr.StructExpression) -> hl.expr.BooleanExpression:
    """Check if an entry has a defined non-ref local genotype.

    :param e: Entry expression with an ``LGT`` field.
    :return: Boolean expression.
    """
    return hl.is_defined(e.LGT) & e.LGT.is_non_ref()


def _get_biallelic_snv(entry: hl.expr.StructExpression) -> hl.expr.StructExpression:
    """Extract biallelic alleles and SNP status from a local-genotype entry.

    Uses the entry's LA (local alleles) array to map from the local alt allele
    index to the global allele index, then constructs a biallelic ``[ref, alt]``
    representation.

    :param entry: Entry expression with ``LGT``, ``LA``, and ``_alleles`` fields.
    :return: Struct with ``alleles`` (array<str>) and ``is_snp`` (bool).
    """
    # For het-ref (LGT=0/1), max gives the alt index (1). For hom-var
    # (LGT=1/1), max also gives the alt index (1). LA then maps from
    # local → global allele index.
    alt_global = entry.LA[hl.max(entry.LGT[0], entry.LGT[1])]
    return hl.struct(
        alleles=hl.array([entry._alleles[0], entry._alleles[alt_global]]),
        is_snp=hl.is_snp(entry._alleles[0], entry._alleles[alt_global]),
    )


def _classify_genotype_pair(
    cur: hl.expr.StructExpression, prev: hl.expr.StructExpression
) -> hl.expr.StructExpression:
    """Classify a pair of entries as het-het, hom-hom, or het-hom.

    Uses LGT/LPGT in local allele space, which is valid without splitting
    multiallelics because het-ref/hom-var status is the same in local and
    global representations, and phase orientation (LPGT 0|1 vs 1|0) is always
    relative to local ref/alt.

    Classification rules:

    - **het-het**: Both entries are het-ref, both have phased LPGT, share the
      same PID (phasing block), and have the same phase orientation
      (``LPGT[0] == prev.LPGT[0]``), confirming the alt alleles are on the
      same haplotype.
    - **hom-hom**: Both entries are homozygous alt.
    - **het-hom**: One entry is het-ref and the other is hom-var (either
      direction). The alt allele is guaranteed to be on at least one shared
      haplotype.

    :param cur: Current entry expression (with LGT, LPGT, PID).
    :param prev: Previous entry expression (with LGT, LPGT, PID).
    :return: Struct with ``is_hethet``, ``is_homhom``, ``is_hethom`` booleans.
    """
    return hl.struct(
        is_hethet=(
            cur.LGT.is_het_ref()
            & prev.LGT.is_het_ref()
            & cur.LPGT.phased
            & prev.LPGT.phased
            & (cur.PID == prev.PID)
            & (cur.LPGT[0] == prev.LPGT[0])
            & (cur.LPGT[1] == prev.LPGT[1])
        ),
        is_homhom=cur.LGT.is_hom_var() & prev.LGT.is_hom_var(),
        is_hethom=(
            (cur.LGT.is_het_ref() & prev.LGT.is_hom_var())
            | (cur.LGT.is_hom_var() & prev.LGT.is_het_ref())
        ),
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
        contain ``LGT`` and ``_alleles``.
    :param max_distance: Maximum bp distance between two SNVs.
    :return: Filtered localized Table with ``prev`` annotated on each entry.
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

    return ht.filter(hl.any(ht.__entries.map(lambda x: hl.is_defined(x.prev)))).cache()


def _classify_mnv_pairs(ht: hl.Table) -> hl.Table:
    """Classify MNV pairs from scan results and explode into one row per pair.

    For each entry with a defined ``prev`` window, classifies all pairs using
    :func:`_classify_genotype_pair` and extracts biallelic alleles using
    :func:`_get_biallelic_snv`. Filters to pairs where both carried alleles are
    SNPs and at least one classification matches.

    :param ht: Localized Table with ``prev`` annotated on entries (output of
        :func:`_scan_for_candidates`).
    :return: Table with one ``_mnv`` struct per row (after checkpoint + explode).
    """
    _mnv_pair_type = hl.tstruct(
        prev_locus=hl.tlocus("GRCh38"),
        prev_alleles=hl.tarray(hl.tstr),
        cur_alleles=hl.tarray(hl.tstr),
        is_hethet=hl.tbool,
        is_homhom=hl.tbool,
        is_hethom=hl.tbool,
    )

    def _build_pair_record(entry, prev_tuple):
        """Build an MNV pair struct for one (current, previous) entry pair."""
        cur_snv = _get_biallelic_snv(entry)
        # Two levels of hl.bind: outer binds prev_tuple[1] (the entry struct
        # from the window tuple) to prev_e; inner binds the biallelic SNV and
        # genotype classification so they're computed once and shared across
        # the struct fields.
        return hl.bind(
            lambda prev_e: hl.bind(
                lambda prev_snv, gt_class: hl.struct(
                    prev_locus=prev_tuple[0],
                    prev_alleles=prev_snv.alleles,
                    cur_alleles=cur_snv.alleles,
                    is_hethet=gt_class.is_hethet,
                    is_homhom=gt_class.is_homhom,
                    is_hethom=gt_class.is_hethom,
                ),
                _get_biallelic_snv(prev_e),
                _classify_genotype_pair(entry, prev_e),
            ),
            prev_tuple[1],
        )

    def _classify_entry(entry):
        """Classify all MNV pairs for one entry with a defined prev window."""
        cur_is_snp = _get_biallelic_snv(entry).is_snp
        return entry.prev.map(
            lambda prev_tuple: _build_pair_record(entry, prev_tuple)
        ).filter(
            lambda p: (
                (p.is_hethet | p.is_homhom | p.is_hethom)
                & cur_is_snp
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
    then counts het-het, hom-hom, het-hom, and total occurrences.

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
    )

    return per_pair.annotate(
        dist=per_pair.locus.position - per_pair.prev_locus.position,
    )


# ---------------------------------------------------------------------------
# Top-level discovery + annotation
# ---------------------------------------------------------------------------


def discover_mnv(
    vds: hl.vds.VariantDataset,
    max_distance: int = MAX_MNV_DISTANCE,
    entries_to_keep: List[str] = MNV_ENTRIES_TO_KEEP,
) -> hl.Table:
    """Run single-pass MNV discovery on an unsplit VDS.

    Uses ``hl.vds.lgt_to_gt`` to convert local genotypes to global allele
    indices inline, avoiding a separate unlocalize + split step. Scans across
    rows using ``hl.scan.fold`` with ``hl.case()`` branching to track a window
    of recent non-ref entries per sample.

    Pipeline:

    1. **Scan** (:func:`_scan_for_candidates`): Unfilter entries, localize,
       annotate with row alleles, scan using ``hl.scan.fold`` tracking a
       sliding window of non-ref entries per sample. Filter to rows with
       candidate MNV pairs.
    2. **Classify** (:func:`_classify_mnv_pairs`): For each candidate, use
       LGT/LPGT directly for classification (het-het / hom-hom / het-hom).
       Extract biallelic alleles via LA indexing for the aggregation key.
    3. **Aggregate** (:func:`_aggregate_mnv_pairs`): Group by variant pair
       across all samples.

    .. note::

        Uses ``_localize_entries`` + ``hl.scan.array_agg`` + ``hl.scan.fold``
        instead of MT entry scans. Direct scans in MT entry context cause
        ``KeyError: 'agg_capability'`` in Hail's ``CSEAnalysisPass``. The
        localized pattern (from gnomAD's ``compute_last_ref_block_end``) works.

    :param vds: Unsplit gnomAD v4 VariantDataset.
    :param max_distance: Maximum bp distance between SNVs to consider as an MNV.
        Default is 2 (codon reading frame).
    :param entries_to_keep: Entry fields needed for MNV discovery.
    :return: Hail Table of MNV pairs with per-pair counts.
    """
    mt = vds.variant_data

    # Filter to rows with at least one SNP alt allele.
    mt = mt.filter_rows(hl.any(lambda a: hl.is_snp(mt.alleles[0], a), mt.alleles[1:]))

    # Select entries needed for classification (pre-split names) + LA.
    pre_split_entries = entries_to_keep + ["LA"]
    pre_split_entries = [
        "L" + e if e in {"GT", "AD", "PL", "PGT"} else e for e in pre_split_entries
    ]
    mt = mt.select_entries(*pre_split_entries)
    mt = mt.select_rows()
    mt = mt.select_cols()

    # Unfilter sparse entries (fills in ref-block-covered positions with
    # LGT=0/0), then localize to a Table for scan-based processing.
    logger.info("Localizing and scanning for MNV candidates (single-pass)...")
    mt = mt.unfilter_entries()
    ht = mt._localize_entries("__entries", "__cols")
    # Carry the row's alleles into each entry struct so that previous entries
    # stored in the scan window retain their original row's alleles (needed
    # for biallelic allele extraction in the classification step).
    ht = ht.annotate(
        __entries=ht.__entries.map(lambda e: e.annotate(_alleles=ht.alleles))
    )

    # Scan → classify → aggregate.
    ht = _scan_for_candidates(ht, max_distance)

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
    from gnomad.resources.grch38.gnomad import public_release
    from gnomad_qc.v4.resources.annotations import get_vep

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

    if args.discover:
        if args.test:
            logger.info(
                "Running in test mode, subsetting to PCSK9 (%s)...",
                TEST_INTERVAL,
            )
        else:
            logger.info("Starting MNV discovery.")
        start = time.time()

        # --- Load unsplit VDS, optionally filtering to test interval. ---
        filter_intervals = [TEST_INTERVAL] if args.test else None
        vds = get_gnomad_v4_vds(
            filter_intervals=filter_intervals,
            high_quality_only=args.high_quality_only,
            release_only=args.release_only,
        )

        # --- Run discovery (single-pass scan + classify + aggregate). ---
        mnv_ht = discover_mnv(
            vds,
            max_distance=args.max_distance,
            entries_to_keep=MNV_ENTRIES_TO_KEEP,
        )

        # --- Write discovery output. ---
        discovery_resource = mnv_discovery(test=args.test)
        logger.info("Writing MNV discovery results to %s.", discovery_resource.path)
        mnv_ht.write(discovery_resource.path, overwrite=args.overwrite)

        elapsed = time.time() - start
        logger.info("Finished MNV discovery in %.1f seconds.", elapsed)

    if args.annotate:
        logger.info("Annotating MNV pairs with frequency and VEP data...")
        mnv_ht = mnv_discovery(test=args.test).ht()
        mnv_ht = annotate_mnv(mnv_ht)

        annotated_resource = mnv_annotated(test=args.test)
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
        "--test",
        help="Run on the PCSK9 locus only (chr1:55039447-55064852) for quick testing.",
        action="store_true",
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

    args = parser.parse_args()
    if not args.discover and not args.annotate:
        parser.error("At least one of --discover or --annotate is required.")
    main(args)
