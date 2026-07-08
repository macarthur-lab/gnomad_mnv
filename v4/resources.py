"""Resource paths and data loading utilities for gnomAD v4 MNV discovery."""

from typing import List, Optional, Union

import hail as hl
from gnomad.resources.resource_utils import TableResource

CURRENT_VERSION = "4.1"

MNV_ENTRIES_TO_KEEP = ["GT", "PGT", "PID"]
"""Entry fields needed for MNV discovery (global names, post-split)."""


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


def _test_gene_suffix(test: bool, test_genes: Optional[List[str]]) -> str:
    """Build a filename suffix identifying the gene(s) used in a test run.

    :param test: Whether this is a test run. If False, no suffix is added.
    :param test_genes: Optional list of gene names tested. If None or empty while
        ``test`` is True, the suffix is ``.test`` (generic, no genes specified).
    :return: Filename suffix, e.g. ``".test_PCSK9_PCNT"``, ``".test"``, or ``""``.
    """
    if not test:
        return ""
    return f".test_{'_'.join(test_genes)}" if test_genes else ".test"


def mnv_discovery(
    version: str = CURRENT_VERSION,
    test: bool = False,
    data_type: str = "exomes",
    test_genes: Optional[List[str]] = None,
) -> TableResource:
    """Get the MNV discovery TableResource.

    Contains raw MNV pair counts (n_hethet, n_homhom, n_hethom, n_total) without
    frequency or VEP annotations.

    :param version: gnomAD version string. Default is :data:`CURRENT_VERSION`.
    :param test: Whether to use the testing output bucket. Default is False.
    :param data_type: Data type, either ``"exomes"`` or ``"genomes"``. Default is
        ``"exomes"``.
    :param test_genes: Optional list of gene names tested (``test`` must be True).
        Appended to the output filename so that multiple test runs, e.g. for
        different genes, don't overwrite each other's output.
    :return: TableResource for MNV discovery output.
    """
    return TableResource(
        path=(
            f"{_mnv_root_path(version, test, data_type)}"
            f"/gnomad.{data_type}.v{version}.mnv_discovery"
            f"{_test_gene_suffix(test, test_genes)}.ht"
        )
    )


def mnv_annotated(
    version: str = CURRENT_VERSION,
    test: bool = False,
    data_type: str = "exomes",
    test_genes: Optional[List[str]] = None,
) -> TableResource:
    """Get the annotated MNV TableResource.

    Contains MNV pair counts annotated with frequency (AC, AF, filters) and/or VEP
    consequences for both SNVs.

    :param version: gnomAD version string. Default is :data:`CURRENT_VERSION`.
    :param test: Whether to use the testing output bucket. Default is False.
    :param data_type: Data type, either ``"exomes"`` or ``"genomes"``. Default is
        ``"exomes"``.
    :param test_genes: Optional list of gene names tested (``test`` must be True).
        Appended to the output filename so that multiple test runs, e.g. for
        different genes, don't overwrite each other's output.
    :return: TableResource for annotated MNV output.
    """
    return TableResource(
        path=(
            f"{_mnv_root_path(version, test, data_type)}"
            f"/gnomad.{data_type}.v{version}.mnv_annotated"
            f"{_test_gene_suffix(test, test_genes)}.ht"
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
