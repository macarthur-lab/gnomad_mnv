"""Build the small MNV test dataset used by the local test suite.

The test dataset is a release-only gnomAD v4 exomes VDS cut down to a few short
windows on the autosomes, chrX (PAR, the PAR1 boundary and non-PAR) and chrY,
chosen so that together they exercise everything the MNV pipeline has to get
right: MNV pairs of each class, het-non-ref (``1/2``) carriers, alt-first phased
hets, the hom-alt depletion hotfix, AS_lowqual sites and sex ploidy.

Building it takes three steps:

1. ``--compute-site-features`` scores every site in the test genes for those
   cases and writes the scores to the test bucket. The windows are then picked
   from these scores and recorded in ``v4/resources.py``.
2. ``--build-vds`` writes the windowed VDS to a permanent location.
3. ``--build-tables`` writes matching slices of the release, info, VEP and sample
   metadata Tables next to it.

Usage::

    python -m v4.create_test_data --compute-site-features --overwrite
    python -m v4.create_test_data --build-vds --build-tables --overwrite
"""

import argparse
import logging

import hail as hl
from gnomad.resources.grch38.gnomad import public_release
from gnomad_qc.v4.resources.annotations import get_info, get_vep
from gnomad_qc.v4.resources.meta import meta

from v4.resources import (
    TEST_DATA_WINDOWS,
    get_gnomad_v4_vds,
    get_test_data,
    get_test_gene_intervals,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("create_test_data")
logger.setLevel(logging.INFO)

SITE_FEATURES_PATH = (
    "gs://gnomad-tmp/gnomad_v4.1_testing/mnv/exomes/test_data_selection/"
    "gnomad.exomes.v4.1.mnv_test_site_features.ht"
)
"""Per-site test-case scores over the test genes, used to pick the windows."""

HOM_ALT_FIX_AF_CUTOFF = 0.01
"""Release AF above which a high-AB het-ref call is a hom-alt hotfix candidate."""

HOM_ALT_FIX_AB_CUTOFF = 0.9
"""Allele balance above which a het-ref call is a hom-alt hotfix candidate."""


def compute_site_features(
    mt: hl.MatrixTable,
    release_ht: hl.Table,
    info_ht: hl.Table,
    meta_ht: hl.Table,
) -> hl.Table:
    """
    Count, per site, the carriers that exercise each MNV test case.

    Works on raw (unadjusted) local genotypes, so the counts say what a window
    would contain before the pipeline's own adjustments. Adj is approximated by
    the GQ/DP/allele-balance thresholds, which is enough to rank sites.

    :param mt: Unsplit variant data MatrixTable of release samples.
    :param release_ht: Split release sites Table (``freq[0].AF`` is the adj AF).
    :param info_ht: Unsplit info Table with ``AS_lowqual`` (one bool per alt).
    :param meta_ht: Sample metadata Table.
    :return: Table keyed by ``(locus, alleles)`` with one count per test case.
    """
    # Release AF of each alt, looked up by its min-repped biallelic key.
    alts = mt.rows().select()
    alts = alts.annotate(_ai=hl.range(1, hl.len(alts.alleles))).explode("_ai")
    alts = alts.annotate(
        _mr=hl.min_rep(alts.locus, hl.array([alts.alleles[0], alts.alleles[alts._ai]]))
    )
    alts = alts.annotate(
        _af=release_ht[alts._mr.locus, alts._mr.alleles].freq[0].AF,
        _shift=alts._mr.locus.position - alts.locus.position,
        _snv=hl.is_snp(alts.alleles[0], alts.alleles[alts._ai]),
    )
    alts = alts.group_by("locus", "alleles").aggregate(
        alt_af=hl.agg.collect((alts._ai, alts._af)),
        n_snv_alts=hl.agg.count_where(alts._snv),
        min_rep_shift=hl.agg.any(alts._snv & (alts._shift > 0)),
    )
    alts = alts.annotate(
        alt_af=hl.sorted(alts.alt_af, key=lambda t: t[0]).map(lambda t: t[1])
    )

    mt = mt.annotate_rows(
        **alts[mt.row_key],
        as_lowqual=info_ht[mt.row_key].AS_lowqual,
        nonpar_sex=mt.locus.in_x_nonpar() | mt.locus.in_y_nonpar(),
    )
    sample_meta = meta_ht[mt.s]
    mt = mt.annotate_cols(
        xy=sample_meta.sex_imputation.sex_karyotype == "XY",
        fixed_homalt_model=sample_meta.project_meta.fixed_homalt_model,
        ukb=sample_meta.project_meta.ukb_sample,
    )

    gt = mt.LGT
    het_ref = gt.is_het_ref() & (gt.ploidy == 2)
    # Only het-ref calls are scored below. or_missing keeps the indexing off other
    # calls: a haploid call has no gt[1], and a hom-ref call has no alt.
    alt_li = hl.or_missing(het_ref, hl.max(gt[0], gt[1]))
    alt_ab = hl.float64(mt.LAD[alt_li]) / mt.DP
    alt_af = mt.alt_af[mt.LA[alt_li] - 1]
    gq_dp_adj = (mt.GQ >= 20) & (mt.DP >= 10)
    hotfix = (
        het_ref
        & ~mt.fixed_homalt_model
        & (alt_af > HOM_ALT_FIX_AF_CUTOFF)
        & (alt_ab > HOM_ALT_FIX_AB_CUTOFF)
    )
    # Old adj read a phased alt-first het's allele balance from the ref depth.
    alt_first = het_ref & gt.phased & (gt[0] > 0)
    old_ab_pass = hl.or_missing(het_ref, hl.float64(mt.LAD[gt[1]]) / mt.DP >= 0.2)
    new_ab_pass = alt_ab >= 0.2
    phased_het = (
        gt.is_het() & hl.is_defined(mt.LPGT) & mt.LPGT.phased & hl.is_defined(mt.PID)
    )

    return (
        mt.annotate_rows(
            n_carriers=hl.agg.count_where(gt.is_non_ref()),
            n_hom=hl.agg.count_where(gt.is_hom_var()),
            n_het_non_ref=hl.agg.count_where(gt.is_het_non_ref()),
            n_phased_het=hl.agg.count_where(phased_het),
            n_unphased_het=hl.agg.count_where(gt.is_het() & ~phased_het),
            n_alt_first=hl.agg.count_where(alt_first),
            n_alt_first_adj_differs=hl.agg.count_where(
                alt_first & gq_dp_adj & (old_ab_pass != new_ab_pass)
            ),
            n_hotfix_adj=hl.agg.count_where(hotfix & gq_dp_adj),
            n_hotfix_non_adj=hl.agg.count_where(hotfix & ~gq_dp_adj),
            n_xy_nonpar_het=hl.agg.count_where(mt.nonpar_sex & mt.xy & gt.is_het()),
            n_xy_nonpar_hotfix=hl.agg.count_where(mt.nonpar_sex & mt.xy & hotfix),
            n_xx_chry=hl.agg.count_where(
                mt.locus.in_y_nonpar() & ~mt.xy & gt.is_non_ref()
            ),
            n_ukb_carriers=hl.agg.count_where(mt.ukb & gt.is_non_ref()),
            n_fixed_homalt_model_carriers=hl.agg.count_where(
                mt.fixed_homalt_model & gt.is_non_ref()
            ),
        )
        .rows()
        .select(
            "n_snv_alts",
            "min_rep_shift",
            "as_lowqual",
            "alt_af",
            "n_carriers",
            "n_hom",
            "n_het_non_ref",
            "n_phased_het",
            "n_unphased_het",
            "n_alt_first",
            "n_alt_first_adj_differs",
            "n_hotfix_adj",
            "n_hotfix_non_adj",
            "n_xy_nonpar_het",
            "n_xy_nonpar_hotfix",
            "n_xx_chry",
            "n_ukb_carriers",
            "n_fixed_homalt_model_carriers",
        )
    )


def main(args: argparse.Namespace) -> None:
    """
    Run the requested test-data steps.

    :param args: Parsed argparse namespace.
    """
    hl.init(
        log="create_test_data.log",
        tmp_dir="gs://gnomad-tmp-4day/create_test_data",
    )

    if args.compute_site_features:
        gene_intervals = [
            hl.parse_locus_interval(i, reference_genome="GRCh38")
            for _, i in get_test_gene_intervals("all")
        ]
        vds = get_gnomad_v4_vds(
            filter_intervals=[i for _, i in get_test_gene_intervals("all")],
            release_only=True,
        )
        ht = compute_site_features(
            vds.variant_data,
            hl.filter_intervals(public_release("exomes").ht(), gene_intervals),
            hl.filter_intervals(get_info(split=False).ht(), gene_intervals),
            meta().ht(),
        )
        logger.info("Writing site features to %s.", SITE_FEATURES_PATH)
        ht.write(SITE_FEATURES_PATH, overwrite=args.overwrite)

    intervals = [
        f"{contig}:{start}-{end}" for contig, start, end in TEST_DATA_WINDOWS.values()
    ]
    parsed = [hl.parse_locus_interval(i, reference_genome="GRCh38") for i in intervals]

    if args.build_vds:
        # get_gnomad_v4_vds splits reference blocks at the window edges, so the
        # reference data only covers the windows.
        vds = get_gnomad_v4_vds(filter_intervals=intervals, release_only=True)
        # Checkpoint at the filtered VDS's own partitioning before coalescing:
        # coalescing straight to one partition runs the interval and sample
        # filtering, and the whole write, in a single task.
        vds = vds.checkpoint(hl.utils.new_temp_file("mnv_test_vds", "vds"))
        vds = hl.vds.VariantDataset(
            vds.reference_data.naive_coalesce(1), vds.variant_data.naive_coalesce(1)
        )
        logger.info("Writing the test VDS to %s.", get_test_data("vds").path)
        vds.write(get_test_data("vds").path, overwrite=args.overwrite)

    if args.build_tables:
        samples = get_test_data("vds").vds().variant_data.cols()
        meta_ht = meta().ht()
        meta_ht = meta_ht.semi_join(samples)
        release_ht = hl.filter_intervals(public_release("exomes").ht(), parsed)
        info_ht = hl.filter_intervals(get_info(split=False).ht(), parsed)
        tables = {
            "release": release_ht.select("freq", "filters").select_globals(
                "freq_meta", "freq_index_dict", "version"
            ),
            "info": info_ht.select("AS_lowqual").select_globals(),
            "vep": hl.filter_intervals(get_vep(data_type="exomes").ht(), parsed).select(
                "vep"
            ),
            "meta": meta_ht.select(
                release=meta_ht.release,
                sex_imputation=hl.struct(
                    sex_karyotype=meta_ht.sex_imputation.sex_karyotype
                ),
                project_meta=hl.struct(
                    fixed_homalt_model=meta_ht.project_meta.fixed_homalt_model,
                    ukb_sample=meta_ht.project_meta.ukb_sample,
                ),
            ).select_globals(),
        }
        for name, ht in tables.items():
            logger.info(
                "Writing the test %s Table to %s.", name, get_test_data(name).path
            )
            ht.naive_coalesce(1).write(
                get_test_data(name).path, overwrite=args.overwrite
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build the small MNV test dataset used by the local test suite."
    )
    parser.add_argument(
        "--compute-site-features",
        help=(
            "Score every site in the test genes for the MNV test cases (release"
            " samples) and write the scores to the test bucket."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--build-vds",
        help=(
            "Write the test VDS (release samples, TEST_DATA_WINDOWS in"
            " v4/resources.py, one partition) to the permanent test_data/ folder."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--build-tables",
        help=(
            "Write release, info, VEP and metadata slices matching the test VDS to"
            " the permanent test_data/ folder. Needs the test VDS."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--overwrite",
        help="Whether to overwrite existing output files.",
        action="store_true",
    )
    args = parser.parse_args()
    if not (args.compute_site_features or args.build_vds or args.build_tables):
        parser.error(
            "Pass a step to run: --compute-site-features, --build-vds or"
            " --build-tables."
        )
    main(args)
