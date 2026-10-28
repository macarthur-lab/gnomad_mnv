"""Resource paths and data loading utilities for gnomAD v4 MNV discovery."""

from typing import List, Optional, Tuple, Union

import hail as hl
from gnomad.resources.resource_utils import TableResource, VariantDatasetResource

CURRENT_VERSION = "4.1"

MNV_ENTRIES_TO_KEEP = ["GT", "PGT", "PID", "GQ", "DP", "AD"]
"""Entry fields needed for MNV discovery (global names, post-split).

GQ/DP/AD are used for the hom-alt depletion hotfix (allele balance) and adj
genotype-quality filtering. The L-prefix loop in :func:`discover_mnv` maps
GT->LGT, AD->LAD, PGT->LPGT and leaves GQ/DP/PID unchanged.
"""

UKB_VDS = VariantDatasetResource("gs://gnomad/v4.0/raw/exomes/gnomad_v4.0.ukb.vds")
"""UKB-only subset of the gnomAD v4 exomes VDS.

Written by ``gnomad_qc.v4.split_ukb_vds`` from
``get_gnomad_v4_vds(remove_hard_filtered_samples=False)`` filtered to the columns
whose sample ID (``s``) starts with ``"UKB"``.
"""

AUTOSOME_TEST_GENES = {
    "AHNAK2": ("chr14", 104937244, 104978374),
    "ANO5": ("chr11", 21782659, 22283567),
    "CAPN3": ("chr15", 42359498, 42412949),
    "DYSF": ("chr2", 71453561, 71686763),
    "SGCA": ("chr17", 50164214, 50175928),
}
"""Autosomal test genes as GRCh38 ``(contig, start, end)`` gene spans.

The five genes the variant co-occurrence exome pipeline is tested on; the spans
match ``TEST_INTERVALS`` in gnomad_chets ``v4/resources.py`` (GENCODE gene spans
from the gnomAD browser).
"""

SEX_CHR_TEST_GENES = {
    # PAR1: diploid in XY samples.
    "SHOX": ("chrX", 624344, 659411),
    # Straddles the PAR1 boundary (PAR1 ends at chrX:2781479), so XY samples
    # switch from diploid to hemizygous within the gene.
    "XG": ("chrX", 2752040, 2816500),
    # Non-PAR X: hemizygous in XY samples.
    "EMD": ("chrX", 154379273, 154381574),
    # Y: haploid in XY samples, missing in XX samples.
    "DDX3Y": ("chrY", 12904108, 12920478),
}
"""Sex-chromosome test genes as GRCh38 ``(contig, start, end)`` gene spans.

One gene per ploidy case, plus one across the PAR1 boundary (GENCODE gene spans
from the gnomAD browser).
"""

TEST_GENE_SETS = {
    "default": {"SGCA": AUTOSOME_TEST_GENES["SGCA"]},
    "autosomes": AUTOSOME_TEST_GENES,
    "sex_chr": SEX_CHR_TEST_GENES,
    "all": {**AUTOSOME_TEST_GENES, **SEX_CHR_TEST_GENES},
}
"""Named test gene sets selectable with ``--test-genes``.

``"default"`` is SGCA, at 11.7 kb the smallest of the autosomal test genes, so a
default run is quick. It's also what ``--intervals`` uses when given no values.
"""

TEST_GENE_PADDING = 100
"""bp added to each side of a test gene's span when building its interval.

Wider than the largest pair distance, so pairs at a gene's edges aren't cut off,
but close enough to the gene that the region still matches the variant
co-occurrence test intervals.
"""


def get_test_gene_intervals(
    gene_set: str, padding: int = TEST_GENE_PADDING
) -> List[Tuple[str, str]]:
    """
    Get the padded (gene_name, interval) pairs for a named test gene set.

    :param gene_set: Key of :data:`TEST_GENE_SETS`: ``"default"``,
        ``"autosomes"``, ``"sex_chr"`` or ``"all"``.
    :param padding: bp added to each side of every gene span. Default is
        :data:`TEST_GENE_PADDING`.
    :return: List of (gene_name, interval) tuples, e.g.
        ``("SGCA", "chr17:50164114-50176028")``.
    """
    return [
        (name, f"{contig}:{max(start - padding, 1)}-{end + padding}")
        for name, (contig, start, end) in TEST_GENE_SETS[gene_set].items()
    ]


def _mnv_root_path(
    version: str = CURRENT_VERSION,
    test: bool = False,
    data_type: str = "exomes",
) -> str:
    """Return the root GCS path for MNV resources.

    :param version: gnomAD version string. Default is :data:`CURRENT_VERSION`.
    :param test: Whether to use the testing output bucket. Default is False.
    :param data_type: Data type, either ``"exomes"`` or ``"genomes"``. Default is
        ``"exomes"``.
    :return: Root GCS path string.
    """
    return (
        f"gs://gnomad-tmp/gnomad_v{version}_testing/mnv/{data_type}"
        if test
        else f"gs://gnomad/v{version}/mnv/{data_type}"
    )


def _interval_suffix(interval_names: Optional[List[str]]) -> str:
    """Build a filename suffix labeling the interval(s) a run was scoped to.

    :param interval_names: Optional list of labels (gene names, or a single
        contig name) for the intervals the run was filtered to. If None, empty, or
        all-empty, no suffix is added.
    :return: Filename suffix, e.g. ``".PCSK9_PCNT"``, ``".chr21"``, or ``""``.
    """
    joined = "_".join(interval_names) if interval_names else ""
    return f".{joined}" if joined else ""


def get_discovered_mnvs(
    version: str = CURRENT_VERSION,
    test: bool = False,
    data_type: str = "exomes",
    interval_names: Optional[List[str]] = None,
    ukb_only: bool = False,
    chrom: Optional[str] = None,
    suffix: Optional[str] = None,
) -> TableResource:
    """Get the MNV discovery TableResource.

    Contains raw MNV pair counts (n_hethet, n_homhom, n_hethom, n_total) without
    frequency or VEP annotations.

    :param version: gnomAD version string. Default is :data:`CURRENT_VERSION`.
    :param test: Whether to use the testing output bucket. Default is False.
    :param data_type: Data type, either ``"exomes"`` or ``"genomes"``. Default is
        ``"exomes"``.
    :param interval_names: Optional labels (gene names) for the intervals the run
        was scoped to via ``--intervals`` or ``--test-genes``. Appended to the
        output filename so runs over different intervals don't overwrite each other.
    :param ukb_only: Whether the run used the UKB-only VDS subset. Adds a
        ``.ukb_only`` filename component so UKB output doesn't overwrite
        full-cohort output. Default is False.
    :param chrom: Optional contig name (e.g. ``"chr21"``) for a single-chromosome
        run. Adds a ``.<chrom>`` filename component so each chromosome's output is
        written to a distinct path. Default is None.
    :param suffix: Optional free-form filename suffix (e.g. a username) appended
        last, so a run doesn't overwrite anyone else's output. Default is None.
    :return: TableResource for MNV discovery output.
    """
    if chrom and interval_names:
        raise ValueError(
            "chrom and interval_names both scope the output path; pass only one."
        )

    return TableResource(
        path=(
            f"{_mnv_root_path(version, test, data_type)}"
            f"/gnomad.{data_type}.v{version}.mnv_discovery"
            f"{'.ukb_only' if ukb_only else ''}"
            f"{'.' + chrom if chrom else ''}"
            f"{_interval_suffix(interval_names)}"
            f"{'.' + suffix if suffix else ''}.ht"
        )
    )


def get_annotated_mnvs(
    version: str = CURRENT_VERSION,
    test: bool = False,
    data_type: str = "exomes",
    interval_names: Optional[List[str]] = None,
    ukb_only: bool = False,
    chrom: Optional[str] = None,
    suffix: Optional[str] = None,
) -> TableResource:
    """Get the annotated MNV TableResource.

    Contains MNV pair counts annotated with frequency (AC, AF, filters) and/or VEP
    consequences for both SNVs.

    :param version: gnomAD version string. Default is :data:`CURRENT_VERSION`.
    :param test: Whether to use the testing output bucket. Default is False.
    :param data_type: Data type, either ``"exomes"`` or ``"genomes"``. Default is
        ``"exomes"``.
    :param interval_names: Optional labels (gene names) for the intervals the run
        was scoped to via ``--intervals`` or ``--test-genes``. Appended to the
        output filename so runs over different intervals don't overwrite each other.
    :param ukb_only: Whether the run used the UKB-only VDS subset. Adds a
        ``.ukb_only`` filename component so UKB output doesn't overwrite
        full-cohort output. Default is False.
    :param chrom: Optional contig name (e.g. ``"chr21"``) for a single-chromosome
        run. Adds a ``.<chrom>`` filename component so each chromosome's output is
        written to a distinct path. Default is None.
    :param suffix: Optional free-form filename suffix (e.g. a username) appended
        last, so a run doesn't overwrite anyone else's output. Default is None.
    :return: TableResource for annotated MNV output.
    """
    if chrom and interval_names:
        raise ValueError(
            "chrom and interval_names both scope the output path; pass only one."
        )

    return TableResource(
        path=(
            f"{_mnv_root_path(version, test, data_type)}"
            f"/gnomad.{data_type}.v{version}.mnv_annotated"
            f"{'.ukb_only' if ukb_only else ''}"
            f"{'.' + chrom if chrom else ''}"
            f"{_interval_suffix(interval_names)}"
            f"{'.' + suffix if suffix else ''}.ht"
        )
    )


def get_gnomad_v4_vds(
    filter_intervals: Optional[List[Union[str, hl.tinterval]]] = None,
    high_quality_only: bool = False,
    release_only: bool = False,
) -> hl.vds.VariantDataset:
    """Load the gnomAD v4 VDS (unsplit) with standard sample filtering.

    Returns the VDS without splitting multiallelics. Callers are responsible for
    any pre-filtering and splitting.

    :param filter_intervals: Optional list of intervals to filter the VDS to before
        returning. Accepts strings (e.g., ``"chr1:55039447-55064852"``) or
        ``hl.tinterval`` objects.
    :param high_quality_only: Whether to filter to only high-quality samples. Default
        is False.
    :param release_only: Whether to filter to only release samples. Default is False.
    :return: Unsplit gnomAD v4 VariantDataset.
    """
    from gnomad_qc.v4.resources.basics import get_gnomad_v4_vds as _get_vds

    return _get_vds(
        split=False,
        filter_intervals=filter_intervals,
        high_quality_only=high_quality_only,
        release_only=release_only,
    )


def get_gnomad_v4_ukb_vds(
    filter_intervals: Optional[List[Union[str, hl.tinterval]]] = None,
    high_quality_only: bool = False,
    release_only: bool = False,
) -> hl.vds.VariantDataset:
    """Load the UKB-only gnomAD v4 VDS subset (unsplit) with sample filtering.

    Drop-in replacement for :func:`get_gnomad_v4_vds` that reads the pre-built
    :data:`UKB_VDS` subset instead of the full v4 VDS. Because the subset is a raw
    VDS (not produced by the gnomad_qc loader), interval and sample filtering are
    applied here by hand to mirror what ``gnomad_qc.v4.resources.basics.get_gnomad_v4_vds``
    does internally.

    :param filter_intervals: Optional list of intervals to filter the VDS to before
        returning. Accepts strings (e.g., ``"chr1:55039447-55064852"``) or
        ``hl.tinterval`` objects.
    :param high_quality_only: Whether to filter to only high-quality samples. Default
        is False.
    :param release_only: Whether to filter to only release samples. Default is False.
    :return: Unsplit UKB-only gnomAD v4 VariantDataset.
    """
    from gnomad_qc.v4.resources.meta import meta
    from gnomad_qc.v4.resources.sample_qc import hard_filtered_samples

    vds = UKB_VDS.vds()

    if filter_intervals:
        if isinstance(filter_intervals[0], str):
            filter_intervals = [
                hl.parse_locus_interval(x, reference_genome="GRCh38")
                for x in filter_intervals
            ]
        vds = hl.vds.filter_intervals(
            vds, filter_intervals, split_reference_blocks=True
        )

    if release_only:
        meta_ht = meta().ht()
        vds = hl.vds.filter_samples(vds, meta_ht.filter(meta_ht.release))
    elif high_quality_only:
        meta_ht = meta().ht()
        vds = hl.vds.filter_samples(vds, meta_ht.filter(meta_ht.high_quality))
    else:
        vds = hl.vds.filter_samples(vds, hard_filtered_samples.ht(), keep=False)

    return vds
