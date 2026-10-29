"""Shared fixtures for the gnomAD v4 MNV tests.

Tests marked ``gcs`` read the small MNV test dataset (see
``v4/create_test_data.py``) from GCS, so they need gcloud credentials and are
skipped when the ``CI`` environment variable is set. Run only the fast tests with
``pytest -m "not gcs"``.
"""

import os

import hail as hl
import pytest

from v4.discover_mnv import discover_mnv, prepare_genotypes
from v4.resources import TEST_DATA_TABLES, get_test_data

TEST_MAX_DISTANCE = 10
"""Pair distance used for every test run, matching the planned 10 bp default."""

DRIVER_MEMORY = os.environ.get("MNV_TEST_DRIVER_MEMORY", "8g")
"""Local Spark driver heap. The scan holds 730,947-wide entry arrays for the test
data, which overflows Spark's 1 GB default; override with
``MNV_TEST_DRIVER_MEMORY``."""


def pytest_addoption(parser: pytest.Parser) -> None:
    """
    Add ``--update-expected`` to rewrite the stored expected MNV pairs.

    :param parser: pytest argument parser.
    """
    parser.addoption(
        "--update-expected",
        action="store_true",
        help="Rewrite tests/data/expected_mnv_pairs.tsv from the current pipeline.",
    )


def pytest_collection_modifyitems(
    config: pytest.Config, items: list[pytest.Item]
) -> None:
    """
    Skip ``gcs`` tests in CI, where there are no GCS credentials.

    :param config: pytest config.
    :param items: Collected test items.
    """
    if not os.environ.get("CI"):
        return
    skip = pytest.mark.skip(reason="needs GCS access to the MNV test data")
    for item in items:
        if "gcs" in item.keywords:
            item.add_marker(skip)


@pytest.fixture(scope="session")
def hail_session(tmp_path_factory: pytest.TempPathFactory) -> None:
    """
    Start one local Hail session for the whole test run.

    :param tmp_path_factory: pytest temp directory factory.
    """
    tmp = tmp_path_factory.mktemp("hail")
    hl.init(
        quiet=True,
        tmp_dir=str(tmp),
        log=str(tmp / "hail.log"),
        spark_conf={"spark.driver.memory": DRIVER_MEMORY},
    )
    yield
    hl.stop()


@pytest.fixture(scope="session")
def test_data(hail_session: None) -> dict:
    """
    Load the small MNV test dataset from GCS.

    :param hail_session: Hail session fixture.
    :return: Dict with ``"vds"`` and one Table per :data:`TEST_DATA_TABLES` entry.
    """
    data = {name: get_test_data(name).ht() for name in TEST_DATA_TABLES}
    data["vds"] = get_test_data("vds").vds()
    return data


@pytest.fixture(scope="session")
def prepared_mt(test_data: dict) -> hl.MatrixTable:
    """
    Run :func:`prepare_genotypes` on the test dataset once per session.

    :param test_data: Test dataset fixture.
    :return: Adjusted variant data, checkpointed to the session's temp directory.
    """
    mt = prepare_genotypes(
        test_data["vds"], test_data["release"], test_data["info"], test_data["meta"]
    )
    return mt.checkpoint(hl.utils.new_temp_file("prepared", "mt"))


@pytest.fixture(scope="session")
def mnv_ht(test_data: dict) -> hl.Table:
    """
    Run :func:`discover_mnv` on the test dataset once per session.

    :param test_data: Test dataset fixture.
    :return: MNV pair Table, checkpointed to the session's temp directory.
    """
    ht = discover_mnv(
        test_data["vds"],
        test_data["release"],
        test_data["info"],
        test_data["meta"],
        max_distance=TEST_MAX_DISTANCE,
    )
    return ht.checkpoint(hl.utils.new_temp_file("mnv_pairs", "ht"))
