"""Check the pipeline on the small MNV test dataset against the v4.1.1 release.

The test dataset (see ``v4/create_test_data.py``) is a release-only VDS cut down to
a few short windows, so the pipeline's adjusted genotypes should reproduce the
release's adj allele counts exactly. These tests read it from GCS and are marked
``gcs``.
"""

import csv
from pathlib import Path

import hail as hl
import pytest

from tests.conftest import TEST_MAX_DISTANCE
from v4.discover_mnv import _get_carried_alts

pytestmark = pytest.mark.gcs

EXPECTED_PAIRS_PATH = Path(__file__).parent / "data" / "expected_mnv_pairs.tsv"
"""Stored MNV pairs for the test dataset, rewritten with ``--update-expected``."""

PAIR_COUNT_FIELDS = [
    "n_hethet",
    "n_homhom",
    "n_hethom",
    "n_total",
    "n_hethet_adj",
    "n_homhom_adj",
    "n_hethom_adj",
    "n_total_adj",
    "dist",
]
"""Per-pair fields compared against the stored expected pairs."""


def count_alt_alleles(mt: hl.MatrixTable) -> hl.Table:
    """
    Count alt alleles per SNV alt from the pipeline's adjusted genotypes.

    Uses the pipeline's own carried-alt records, so the alleles are min-repped the
    same way the MNV pairs are, and uses each record's ``adj`` when it has one
    (per-alt adj) or the entry's ``adj`` otherwise.

    :param mt: Output of :func:`v4.discover_mnv.prepare_genotypes`.
    :return: Table keyed by ``(locus, alleles)`` with ``adj`` and ``raw`` structs,
        each holding ``AC``, ``homozygote_count`` (diploid only), ``AC_XY``,
        ``n_carriers`` and ``n_hom`` (diploid or haploid).
    """
    mt = mt.annotate_entries(
        _alts=hl.or_missing(
            hl.is_defined(mt.LGT) & mt.LGT.is_non_ref(),
            _get_carried_alts(mt.entry.annotate(_alleles=mt.alleles, _locus=mt.locus)),
        )
    )
    ht = mt.entries()
    ht = ht.filter(hl.is_defined(ht._alts)).explode("_alts")
    alt = ht._alts
    adj = alt.adj if "adj" in alt.dtype.fields else ht.adj
    n_copies = hl.if_else(alt.is_hom, ht.LGT.ploidy, 1)

    def _counts(keep: hl.expr.BooleanExpression) -> hl.expr.StructExpression:
        return hl.struct(
            AC=hl.agg.filter(keep, hl.agg.sum(n_copies)),
            homozygote_count=hl.agg.count_where(
                keep & alt.is_hom & (ht.LGT.ploidy == 2)
            ),
            AC_XY=hl.agg.filter(
                keep & (ht.sex_karyotype == "XY"), hl.agg.sum(n_copies)
            ),
            n_carriers=hl.agg.count_where(keep),
            n_hom=hl.agg.count_where(keep & alt.is_hom),
        )

    return ht.group_by(locus=alt.locus, alleles=alt.alleles).aggregate(
        adj=_counts(adj), raw=_counts(hl.bool(True))
    )


@pytest.fixture(scope="module")
def alt_counts(prepared_mt: hl.MatrixTable) -> hl.Table:
    """
    Per-alt adj and raw counts from the adjusted test genotypes.

    :param prepared_mt: Adjusted test genotypes fixture.
    :return: Output of :func:`count_alt_alleles`, checkpointed.
    """
    return count_alt_alleles(prepared_mt).checkpoint(
        hl.utils.new_temp_file("alt_counts", "ht")
    )


@pytest.fixture(scope="module")
def release_snvs(test_data: dict) -> hl.Table:
    """
    Release SNVs with at least one raw alt allele, with their adj and raw counts.

    :param test_data: Test dataset fixture.
    :return: Release slice keyed by ``(locus, alleles)``.
    """
    release = test_data["release"]
    xy = hl.eval(release.freq_index_dict["XY_adj"])
    release = release.filter(
        hl.is_snp(release.alleles[0], release.alleles[1]) & (release.freq[1].AC > 0)
    )
    return release.select(
        adj=hl.struct(
            AC=release.freq[0].AC,
            homozygote_count=release.freq[0].homozygote_count,
            AC_XY=release.freq[xy].AC,
        ),
        raw=hl.struct(
            AC=release.freq[1].AC, homozygote_count=release.freq[1].homozygote_count
        ),
    )


def _mismatches(
    alt_counts: hl.Table, release_snvs: hl.Table, group: str, fields: list
) -> list:
    """
    Collect alts whose pipeline counts differ from the release's.

    Only alts in both Tables are compared; the coverage tests check the rest.

    :param alt_counts: Output of :func:`count_alt_alleles`.
    :param release_snvs: Output of the ``release_snvs`` fixture.
    :param group: ``"adj"`` or ``"raw"``.
    :param fields: Count fields to compare.
    :return: Rows with both sets of counts, for alts where any field differs.
    """
    ours = alt_counts[release_snvs.key][group]
    theirs = release_snvs[group]
    ht = release_snvs.select(ours=ours.select(*fields), release=theirs.select(*fields))
    ht = ht.filter(
        hl.is_defined(ht.ours) & hl.any([ht.ours[f] != ht.release[f] for f in fields])
    )
    return ht.collect()


def test_every_release_snv_is_counted(
    alt_counts: hl.Table, release_snvs: hl.Table
) -> None:
    """Every release SNV with a raw alt allele shows up in the pipeline's counts."""
    missing = release_snvs.anti_join(alt_counts)
    assert missing.count() == 0, missing.take(5)


def test_extra_alts_are_at_lowqual_sites(
    alt_counts: hl.Table, release_snvs: hl.Table, test_data: dict
) -> None:
    """
    Alts the release doesn't have come from sites with an AS_lowqual alt.

    The pipeline keeps every alt at a site that isn't entirely AS_lowqual; the
    release drops the lowqual ones (see gnomad_production#2333).
    """
    info = test_data["info"]
    lowqual_loci = info.filter(hl.any(info.AS_lowqual)).key_by("locus").select()
    extra = alt_counts.anti_join(release_snvs)
    extra = extra.filter(hl.is_missing(lowqual_loci[extra.locus]))
    assert extra.count() == 0, extra.take(5)


def test_adj_alt_counts_match_release(
    alt_counts: hl.Table, release_snvs: hl.Table
) -> None:
    """
    Every release SNV's adj AC, hom count and XY AC equal the pipeline's.

    The phased alt-first hets here only match on gnomad_methods < 0.9.1 until
    gnomad_production#2300 pins the v4.1 allele-balance reading.
    """
    bad = _mismatches(
        alt_counts, release_snvs, "adj", ["AC", "homozygote_count", "AC_XY"]
    )
    assert not bad, f"{len(bad)} alts differ, e.g. {bad[:5]}"


@pytest.mark.xfail(
    strict=True,
    reason="needs the v4.1-matching hotfix order (gnomad_production#2300)",
)
def test_raw_alt_counts_match_release(
    alt_counts: hl.Table, release_snvs: hl.Table
) -> None:
    """
    Every release SNV's raw AC and hom count equal the pipeline's.

    Raw counts are where the hotfix order shows: v4 exomes only hotfixed calls
    that pass adj, so a high-AB het failing adj stays a het in the release.
    """
    bad = _mismatches(alt_counts, release_snvs, "raw", ["AC", "homozygote_count"])
    assert not bad, f"{len(bad)} alts differ, e.g. {bad[:5]}"


def test_pairs_are_within_max_distance(mnv_ht: hl.Table) -> None:
    """Every pair is 1 to :data:`TEST_MAX_DISTANCE` bp apart."""
    bad = mnv_ht.filter((mnv_ht.dist < 1) | (mnv_ht.dist > TEST_MAX_DISTANCE))
    assert bad.count() == 0


def test_every_pair_class_is_present(mnv_ht: hl.Table) -> None:
    """The test data contains adj het-het, hom-hom and het-hom pairs at 1, 2 and 3+ bp."""
    found = mnv_ht.aggregate(
        hl.struct(
            hethet=hl.agg.any(mnv_ht.n_hethet_adj > 0),
            homhom=hl.agg.any(mnv_ht.n_homhom_adj > 0),
            hethom=hl.agg.any(mnv_ht.n_hethom_adj > 0),
            dist_1=hl.agg.any(mnv_ht.dist == 1),
            dist_2=hl.agg.any(mnv_ht.dist == 2),
            dist_3_plus=hl.agg.any(mnv_ht.dist >= 3),
        )
    )
    assert all(found.values()), found


def test_pair_counts_fit_alt_counts(mnv_ht: hl.Table, alt_counts: hl.Table) -> None:
    """No pair has more adj carriers, or hom-hom carriers, than either of its alts."""
    cur = alt_counts[mnv_ht.locus, mnv_ht.alleles].adj
    prev = alt_counts[mnv_ht.prev_locus, mnv_ht.prev_alleles].adj
    bad = mnv_ht.filter(
        (mnv_ht.n_total_adj > hl.min(cur.n_carriers, prev.n_carriers))
        | (mnv_ht.n_homhom_adj > hl.min(cur.n_hom, prev.n_hom))
    )
    assert bad.count() == 0, bad.take(5)


def test_pairs_match_expected(mnv_ht: hl.Table, request: pytest.FixtureRequest) -> None:
    """The MNV pairs equal the stored expected pairs (``--update-expected`` rewrites)."""
    rows = [
        {
            "locus": str(r.locus),
            "alleles": ",".join(r.alleles),
            "prev_locus": str(r.prev_locus),
            "prev_alleles": ",".join(r.prev_alleles),
            **{f: str(r[f]) for f in PAIR_COUNT_FIELDS},
        }
        for r in mnv_ht.collect()
    ]
    if request.config.getoption("--update-expected"):
        EXPECTED_PAIRS_PATH.parent.mkdir(exist_ok=True)
        with EXPECTED_PAIRS_PATH.open("w", newline="") as f:
            writer = csv.DictWriter(
                f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
        pytest.skip(f"Rewrote {EXPECTED_PAIRS_PATH}.")
    with EXPECTED_PAIRS_PATH.open() as f:
        expected = list(csv.DictReader(f, delimiter="\t"))
    assert rows == expected
