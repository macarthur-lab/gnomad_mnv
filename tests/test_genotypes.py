"""Check genotype adjustment and pair classification on hand-built data.

These tests use only tiny synthetic inputs, so they need a local Hail but no GCS
access. Cases marked ``xfail`` describe the behaviour decided in
broadinstitute/gnomad_production#2300 (match the v4.1 exome release) that the
pipeline doesn't implement yet; ``strict=True`` makes them fail loudly once it
does, so the marker gets removed.
"""

import hail as hl
import pytest

from tests.synthetic import entry, make_inputs
from v4.discover_mnv import _get_carried_alts, discover_mnv, prepare_genotypes

pytestmark = pytest.mark.usefixtures("hail_session")

NEEDS_2300 = pytest.mark.xfail(
    strict=True,
    reason="needs the v4.1-matching genotype adjustments (gnomad_production#2300)",
)

SAMPLES = {
    "xx": {"sex_karyotype": "XX", "fixed_homalt_model": False},
    "xx_fixed": {"sex_karyotype": "XX", "fixed_homalt_model": True},
    "xy": {"sex_karyotype": "XY", "fixed_homalt_model": False},
    "xy2": {"sex_karyotype": "XY", "fixed_homalt_model": False},
}

SITES = [
    # Common site: a high-AB het becomes hom-alt, except in fixed_homalt_model
    # samples; a balanced het is left alone.
    {
        "locus": "chr1:1000",
        "alleles": ["A", "C"],
        "af": [0.05],
        "entries": {
            "xx": entry([0, 1], [1, 19]),
            "xx_fixed": entry([0, 1], [1, 19]),
            "xy": entry([0, 1], [10, 10]),
        },
    },
    # het-non-ref calls are never hotfixed.
    {
        "locus": "chr1:1005",
        "alleles": ["A", "C", "G"],
        "af": [0.05, 0.05],
        "entries": {"xx": entry([1, 2], [0, 19, 1])},
    },
    # Rare site: no hotfix.
    {
        "locus": "chr1:2000",
        "alleles": ["A", "C"],
        "af": [0.001],
        "entries": {"xx": entry([0, 1], [1, 19])},
    },
    # Every alt AS_lowqual: the site is dropped. Mixed: kept.
    {
        "locus": "chr1:3000",
        "alleles": ["A", "C"],
        "af": [0.05],
        "lowqual": [True],
        "entries": {"xx": entry([0, 1], [10, 10])},
    },
    {
        "locus": "chr1:3005",
        "alleles": ["A", "C", "AT"],
        "af": [0.05, 0.05],
        "lowqual": [True, False],
        "entries": {"xx": entry([0, 1], [10, 10, 0])},
    },
    # No SNV alt: the site is dropped.
    {
        "locus": "chr1:4000",
        "alleles": ["A", "AT"],
        "af": [0.05],
        "entries": {"xx": entry([0, 1], [10, 10])},
    },
    # High-AB het that fails adj (GQ 10).
    {
        "locus": "chr1:5000",
        "alleles": ["A", "C"],
        "af": [0.05],
        "entries": {"xx": entry([0, 1], [1, 19], gq=10)},
    },
    # Phased alt-first het: v4.1 read its allele balance from the ref depth
    # (2 / 20 < 0.2), so it fails adj.
    {
        "locus": "chr1:6000",
        "alleles": ["A", "C"],
        "af": [0.001],
        "entries": {"xx_fixed": entry([1, 0], [2, 18], phased=True, pid="p")},
    },
    # Unbalanced 1/2: the release splits first, so C (17 / 20) passes adj and
    # G (3 / 20) fails.
    {
        "locus": "chr1:7000",
        "alleles": ["A", "C", "G"],
        "af": [0.001, 0.001],
        "entries": {"xx": entry([1, 2], [0, 17, 3])},
    },
    # Non-PAR X: an XY het is set missing, an XY hom becomes haploid, XX is
    # unchanged.
    {
        "locus": "chrX:154380000",
        "alleles": ["T", "C"],
        "af": [0.001],
        "entries": {
            "xy": entry([0, 1], [10, 10]),
            "xy2": entry([1, 1], [0, 20]),
            "xx": entry([0, 1], [10, 10]),
        },
    },
    # Non-PAR X high-AB XY het at a common site: v4 exomes applies sex ploidy
    # before the hotfix, so the het is set missing rather than made hom-alt.
    {
        "locus": "chrX:154380010",
        "alleles": ["T", "C"],
        "af": [0.05],
        "entries": {"xy": entry([0, 1], [1, 19])},
    },
    # PAR1: an XY het stays a diploid het.
    {
        "locus": "chrX:630000",
        "alleles": ["T", "C"],
        "af": [0.001],
        "entries": {"xy": entry([0, 1], [10, 10])},
    },
    # Non-PAR Y: XX calls are set missing, XY calls become haploid.
    {
        "locus": "chrY:12905000",
        "alleles": ["T", "C"],
        "af": [0.001],
        "entries": {
            "xx": entry([1, 1], [0, 20]),
            "xy": entry([1, 1], [0, 20]),
        },
    },
]


@pytest.fixture(scope="module")
def genotypes() -> dict:
    """
    Run :func:`prepare_genotypes` on :data:`SITES` and collect the entries.

    :return: Dict keyed by ``(locus string, sample)`` with each entry's adjusted
        ``LGT``, ``adj`` and carried-alt records.
    """
    inputs = make_inputs(SITES, SAMPLES)
    mt = prepare_genotypes(
        inputs.vds, inputs.release_ht, inputs.info_ht, inputs.meta_ht
    )
    mt = mt.annotate_entries(
        _alts=hl.or_missing(
            hl.is_defined(mt.LGT) & mt.LGT.is_non_ref(),
            _get_carried_alts(mt.entry.annotate(_alleles=mt.alleles, _locus=mt.locus)),
        )
    )
    entries = mt.entries()
    entries = entries.filter(hl.is_defined(entries.LA))
    return {
        (str(e.locus), e.s): e for e in entries.select("LGT", "adj", "_alts").collect()
    }


def test_hotfix_makes_common_high_ab_het_hom_alt(genotypes: dict) -> None:
    """A high-AB het at a common site becomes hom-alt and stays adj."""
    e = genotypes[("chr1:1000", "xx")]
    assert e.LGT == hl.Call([1, 1])
    assert e.adj


def test_hotfix_skips_fixed_homalt_model_samples(genotypes: dict) -> None:
    """Samples called with the fixed hom-alt model keep their het call."""
    assert genotypes[("chr1:1000", "xx_fixed")].LGT == hl.Call([0, 1])


def test_hotfix_skips_balanced_hets(genotypes: dict) -> None:
    """A balanced het at a common site is left alone."""
    assert genotypes[("chr1:1000", "xy")].LGT == hl.Call([0, 1])


def test_hotfix_skips_het_non_ref(genotypes: dict) -> None:
    """A 1/2 call is never converted, even with one dominant alt."""
    assert genotypes[("chr1:1005", "xx")].LGT == hl.Call([1, 2])


def test_hotfix_skips_rare_sites(genotypes: dict) -> None:
    """A high-AB het at a site with release AF ≤ 1% is left alone."""
    assert genotypes[("chr1:2000", "xx")].LGT == hl.Call([0, 1])


@NEEDS_2300
def test_hotfix_skips_non_adj_calls(genotypes: dict) -> None:
    """v4 exomes only hotfixed calls that pass adj."""
    assert genotypes[("chr1:5000", "xx")].LGT == hl.Call([0, 1])


def test_all_lowqual_site_is_dropped(genotypes: dict) -> None:
    """A site whose every alt is AS_lowqual is removed; a mixed site is kept."""
    loci = {locus for locus, _ in genotypes}
    assert "chr1:3000" not in loci
    assert "chr1:3005" in loci


def test_site_without_snv_alt_is_dropped(genotypes: dict) -> None:
    """A site with only indel alts can't form an MNV and is removed."""
    assert "chr1:4000" not in {locus for locus, _ in genotypes}


def test_alt_first_phased_het_uses_v4_1_allele_balance(genotypes: dict) -> None:
    """
    A phased alt-first het reads allele balance from the ref depth, as in v4.1.

    Passes on gnomad_methods < 0.9.1; once gnomad_production#2300 lands it must
    pass on any version (``legacy_het_ab=True``).
    """
    assert not genotypes[("chr1:6000", "xx_fixed")].adj


@NEEDS_2300
def test_het_non_ref_adj_is_per_alt(genotypes: dict) -> None:
    """Each alt of a 1/2 call gets its own adj, as the release does after splitting."""
    alts = {r.alleles[1]: r.adj for r in genotypes[("chr1:7000", "xx")]._alts}
    assert alts == {"C": True, "G": False}


def test_xy_non_par_het_is_set_missing(genotypes: dict) -> None:
    """Sex ploidy sets an XY het on non-PAR X to missing."""
    assert genotypes[("chrX:154380000", "xy")].LGT is None


def test_xy_non_par_hom_becomes_haploid(genotypes: dict) -> None:
    """Sex ploidy makes an XY hom-alt on non-PAR X a haploid alt call."""
    assert genotypes[("chrX:154380000", "xy2")].LGT == hl.Call([1])


def test_xx_non_par_het_is_unchanged(genotypes: dict) -> None:
    """XX calls on non-PAR X stay diploid."""
    assert genotypes[("chrX:154380000", "xx")].LGT == hl.Call([0, 1])


@NEEDS_2300
def test_xy_non_par_high_ab_het_is_set_missing(genotypes: dict) -> None:
    """v4 exomes applied sex ploidy before the hotfix, so this het is dropped."""
    assert genotypes[("chrX:154380010", "xy")].LGT is None


def test_par_het_stays_diploid(genotypes: dict) -> None:
    """An XY het in PAR1 is left as a diploid het."""
    assert genotypes[("chrX:630000", "xy")].LGT == hl.Call([0, 1])


def test_chry_calls_follow_sex_karyotype(genotypes: dict) -> None:
    """On non-PAR Y, XX calls are set missing and XY calls become haploid."""
    assert genotypes[("chrY:12905000", "xx")].LGT is None
    assert genotypes[("chrY:12905000", "xy")].LGT == hl.Call([1])


PAIR_SAMPLES = {
    s: {"sex_karyotype": "XX", "fixed_homalt_model": False}
    for s in ["cis", "trans", "unphased", "hom", "het_hom", "far"]
}

PAIR_SITES = [
    {
        "locus": "chr2:100",
        "alleles": ["A", "C"],
        "af": [0.001],
        "entries": {
            "cis": entry([0, 1], [10, 10], phased=True, pid="a"),
            "trans": entry([0, 1], [10, 10], phased=True, pid="b"),
            "unphased": entry([0, 1], [10, 10]),
            "hom": entry([1, 1], [0, 20]),
            "het_hom": entry([0, 1], [10, 10]),
        },
    },
    {
        "locus": "chr2:101",
        "alleles": ["A", "G"],
        "af": [0.001],
        "entries": {
            "cis": entry([0, 1], [10, 10], phased=True, pid="a"),
            "trans": entry([1, 0], [10, 10], phased=True, pid="b"),
            "unphased": entry([0, 1], [10, 10]),
            "hom": entry([1, 1], [0, 20]),
            "het_hom": entry([1, 1], [0, 20]),
        },
    },
    # 10 bp from chr2:100: still a pair.
    {
        "locus": "chr2:110",
        "alleles": ["A", "T"],
        "af": [0.001],
        "entries": {"hom": entry([1, 1], [0, 20]), "far": entry([1, 1], [0, 20])},
    },
    # 11 bp from chr2:110: too far.
    {
        "locus": "chr2:121",
        "alleles": ["A", "C"],
        "af": [0.001],
        "entries": {"far": entry([1, 1], [0, 20])},
    },
]


def test_pair_classification() -> None:
    """
    Count pairs by class, leaving out trans and unphased het pairs.

    The 10 bp limit is inclusive.
    """
    inputs = make_inputs(PAIR_SITES, PAIR_SAMPLES)
    ht = discover_mnv(
        inputs.vds,
        inputs.release_ht,
        inputs.info_ht,
        inputs.meta_ht,
        max_distance=10,
    )
    pairs = {
        (str(r.prev_locus), str(r.locus)): (
            r.n_hethet,
            r.n_homhom,
            r.n_hethom,
            r.n_total_adj,
            r.dist,
        )
        for r in ht.collect()
    }
    assert pairs == {
        ("chr2:100", "chr2:101"): (1, 1, 1, 3, 1),
        ("chr2:100", "chr2:110"): (0, 1, 0, 1, 10),
        ("chr2:101", "chr2:110"): (0, 1, 0, 1, 9),
    }


# KNOWN_ISSUES #1: the scan pairs sites by their site locus, but emits each SNV alt
# at its min-repped locus, which moves forward when the SNP sits after the first
# base of a padded multiallelic ref. A scan of every v4 exome site found no such
# SNV alt, so these cases only exist in hand-built data. The release AF lookups for
# the padded alts miss (their release key is min-repped), which only means no
# hotfix.
MIN_REP_SAMPLES = {"s": {"sex_karyotype": "XX", "fixed_homalt_model": False}}

MIN_REP_SITES = [
    # Clean SNP 1 bp before a padded site whose SNP alt min-reps to chr3:1002, so
    # the pair is formed (site distance 1) but emitted 3 bp apart.
    {
        "locus": "chr3:999",
        "alleles": ["A", "G"],
        "af": [0.001],
        "entries": {"s": entry([1, 1], [0, 20])},
    },
    {
        "locus": "chr3:1000",
        "alleles": ["GCC", "G", "GCT"],
        "af": [0.001, 0.001],
        "entries": {"s": entry([2, 2], [0, 0, 20])},
    },
    # Padded site whose SNP alt min-reps to chr3:2001, 2 bp before a clean SNP at
    # chr3:2003: within 2 bp after min-rep, but 3 bp apart by site locus.
    {
        "locus": "chr3:2000",
        "alleles": ["GCC", "G", "GTC"],
        "af": [0.001, 0.001],
        "entries": {"s": entry([2, 2], [0, 0, 20])},
    },
    {
        "locus": "chr3:2003",
        "alleles": ["A", "G"],
        "af": [0.001],
        "entries": {"s": entry([1, 1], [0, 20])},
    },
]


@pytest.fixture(scope="module")
def min_rep_pairs() -> dict:
    """
    Run :func:`discover_mnv` at ``max_distance=2`` on :data:`MIN_REP_SITES`.

    :return: Dict keyed by ``(prev locus, locus)`` strings, valued by ``dist``.
    """
    inputs = make_inputs(MIN_REP_SITES, MIN_REP_SAMPLES)
    ht = discover_mnv(
        inputs.vds, inputs.release_ht, inputs.info_ht, inputs.meta_ht, max_distance=2
    )
    return {(str(r.prev_locus), str(r.locus)): r.dist for r in ht.collect()}


def test_min_rep_shift_out_of_range_pair_is_emitted(min_rep_pairs: dict) -> None:
    """
    A pair 1 bp apart by site locus but 3 bp apart after min-rep is still emitted.

    ``main`` writes such pairs to ``*.dist_outliers.ht`` (KNOWN_ISSUES #1, defect 2).
    """
    assert min_rep_pairs.get(("chr3:999", "chr3:1002")) == 3


@pytest.mark.xfail(
    strict=True,
    reason="KNOWN_ISSUES #1: the scan pairs by site locus; never seen in v4 exomes",
)
def test_min_rep_shift_pair_within_range_is_found(min_rep_pairs: dict) -> None:
    """A pair 2 bp apart after min-rep but 3 bp apart by site locus is found."""
    assert min_rep_pairs.get(("chr3:2001", "chr3:2003")) == 2
