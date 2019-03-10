# -*- coding: utf-8 -*-
__author__ = 'QingboWang'


import pandas as pd
import time as tm
import hail as hl
import hail.expr.aggregators as agg
hl.init()

#functions for importing gnomAD vcf
#taken from https://github.com/macarthur-lab/gnomad_hail/blob/master/resources/basics.py
from typing import *

CURRENT_HAIL_VERSION = "0.2"
CURRENT_RELEASE = "2.0.2"
CURRENT_GENOME_META = "2018-04-25"  # YYYY-MM-DD
CURRENT_EXOME_META = "2018-04-25"
CURRENT_FAM = '2018-04-12'
CURRENT_DUPS = '2017-10-04'

RELEASES = ["2.0.1", "2.0.2"]

GENOME_POPS = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH']
EXOME_POPS = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS']
EXAC_POPS = ["AFR", "AMR", "EAS", "FIN", "NFE", "OTH", "SAS"]

def public_exomes_mt_path(split=True, version=CURRENT_RELEASE):
    return 'gs://gnomad-public/release/{0}/mt/exomes/gnomad.exomes.r{0}.sites{1}.mt'.format(version, "" if split else ".unsplit")


def public_genomes_mt_path(split=True, version=CURRENT_RELEASE):
    return 'gs://gnomad-public/release/{0}/mt/genomes/gnomad.genomes.r{0}.sites{1}.mt'.format(version, "" if split else ".unsplit")


def get_gnomad_public_data(data_type, split=True, version=CURRENT_RELEASE):
    """
    Wrapper function to get public gnomAD data as VDS.
    :param str data_type: One of `exomes` or `genomes`
    :param bool split: Whether the dataset should be split
    :param str version: One of the RELEASEs
    :return: Chosen VDS
    :rtype: MatrixTable
    """
    return hl.read_matrix_table(get_gnomad_public_data_path(data_type, split=split, version=version))


def get_gnomad_data(data_type: str, adj: bool = False, split: bool = True, raw: bool = False,
                    non_refs_only: bool = False, hail_version: str = CURRENT_HAIL_VERSION,
                    meta_version: str = None, meta_root: Optional[str] = 'meta', full_meta: bool = False,
                    fam_version: str = CURRENT_FAM, fam_root: str = None, duplicate_mapping_root: str = None,
                    release_samples: bool = False, release_annotations: bool = None) -> hl.MatrixTable:
    """
    Wrapper function to get gnomAD data as VDS. By default, returns split hardcalls (with adj annotated but not filtered)
    :param str data_type: One of `exomes` or `genomes`
    :param bool adj: Whether the returned data should be filtered to adj genotypes
    :param bool split: Whether the dataset should be split (only applies to raw=False)
    :param bool raw: Whether to return the raw (10T+) data (not recommended: unsplit, and no special consideration on sex chromosomes)
    :param bool non_refs_only: Whether to return the non-ref-genotype only MT (warning: no special consideration on sex chromosomes)
    :param str hail_version: One of the HAIL_VERSIONs
    :param str meta_version: Version of metadata (None for current)
    :param str meta_root: Where to put metadata. Set to None if no metadata is desired.
    :param str full_meta: Whether to add all metadata (warning: large)
    :param str fam_version: Version of metadata (default to current)
    :param str fam_root: Where to put the pedigree information. Set to None if no pedigree information is desired.
    :param str duplicate_mapping_root: Where to put the duplicate genome/exome samples ID mapping (default is None -- do not annotate)
    :param bool release_samples: When set, filters the data to release samples only
    :param str release_annotations: One of the RELEASES to add variant annotations (into va), or None for no data
    :return: gnomAD hardcalls dataset with chosen annotations
    :rtype: MatrixTable
    """
    #from gnomad_hail.utils import filter_to_adj # なぜかsubmitするときerrorるので...

    if raw and split:
        raise DataException('No split raw data. Use of hardcalls is recommended.')

    if non_refs_only:
        mt = hl.read_matrix_table(get_gnomad_data_path(data_type, split=split, non_refs_only=non_refs_only, hail_version=hail_version))
    else:
        mt = hl.read_matrix_table(get_gnomad_data_path(data_type, hardcalls=not raw, split=split, hail_version=hail_version))

    if adj:
        pass  # なぜかsubmitするときerrorるので...
        #mt = filter_to_adj(mt)

    if meta_root:
        meta_ht = get_gnomad_meta(data_type, meta_version, full_meta=full_meta)
        mt = mt.annotate_cols(**{meta_root: meta_ht[mt.s]})

    if duplicate_mapping_root:
        dup_ht = hl.import_table(genomes_exomes_duplicate_ids_tsv_path, impute=True,
                                 key='exome_id' if data_type == "exomes" else 'genome_id')
        mt = mt.annotate_cols(**{duplicate_mapping_root: dup_ht[mt.s]})

    if fam_root:
        fam_ht = hl.import_fam(fam_path(data_type, fam_version))
        mt = mt.annotate_cols(**{fam_root: fam_ht[mt.s]})

    if release_samples:
        mt = mt.filter_cols(mt.meta.release)

    if release_annotations:
        sites_mt = get_gnomad_public_data(data_type, split, release_annotations)
        mt = mt.select_rows(release=sites_mt[mt.v, :])  # TODO: replace with ** to nuke old annotations

    return mt


def get_gnomad_meta(data_type: str, version: str = None, full_meta: bool = False) -> hl.Table:
    """
    Wrapper function to get gnomAD metadata as Table
    :param str data_type: One of `exomes` or `genomes`
    :param str version: Metadata version (None for current)
    :param bool full_meta: Whether to annotate full metadata (rather than just summarized version)
    :return: Metadata Table
    :rtype: Table
    """
    ht = hl.read_table(get_gnomad_meta_path(data_type, version)).key_by('s')
    if not full_meta:
        columns = ['age', 'sex',
                   'hard_filters', 'perm_filters', 'pop_platform_filters', 'related',
                   'data_type', 'product', 'product_simplified', 'qc_platform',
                   'project_id', 'project_description', 'internal', 'investigator',
                   'known_pop', 'known_subpop', 'pop', 'subpop',
                   'neuro', 'control',
                   'high_quality', 'release']
        if data_type == 'genomes':
            columns.extend(['pcr_free', 'project_name', 'release_2_0_2'])
        else:
            columns.extend(['diabetes', 'exac_joint'])
        ht = ht.select(*columns)
    return ht


def get_gnomad_public_data_path(data_type, split=True, version=CURRENT_RELEASE):
    """
    Wrapper function to get paths to gnomAD data
    :param str data_type: One of `exomes` or `genomes`
    :param bool split: Whether the dataset should be split
    :param str version: One of the RELEASEs
    :return: Path to chosen VDS
    :rtype: str
    """
    if version not in RELEASES:
        return DataException("Select version as one of: {}".format(','.join(RELEASES)))

    if data_type == 'exomes':
        return public_exomes_mt_path(split, version)
    elif data_type == 'genomes':
        return public_genomes_mt_path(split, version)
    return DataException("Select data_type as one of 'genomes' or 'exomes'")


def get_gnomad_data_path(data_type, hardcalls=False, split=True, non_refs_only=False, hail_version=CURRENT_HAIL_VERSION):
    """
    Wrapper function to get paths to gnomAD data
    :param str data_type: One of `exomes` or `genomes`
    :param bool hardcalls: Whether hardcalls should be returned
    :param bool split: Whether the dataset should be split (applies to hardcalls and non_refs_only)
    :param bool non_refs_only: Whether non-ref-genotype only MT should be returned
    :param str hail_version: One of the HAIL_VERSIONs
    :return: Path to chosen VDS
    :rtype: str
    """
    if hardcalls and non_refs_only:
        raise DataException('No dataset with hardcalls and non_refs_only')
    if data_type not in ('exomes', 'genomes'):
        raise DataException("Select data_type as one of 'genomes' or 'exomes'")
    if hardcalls:
        return hardcalls_mt_path(data_type, split, hail_version)
    elif non_refs_only:
        return non_refs_only_mt_path(data_type, split)
    else:
        return raw_exomes_mt_path(hail_version) if data_type == 'exomes' else raw_genomes_mt_path(hail_version)


def get_gnomad_meta_path(data_type, version=None):
    """
    Wrapper function to get paths to gnomAD metadata
    :param str data_type: One of `exomes` or `genomes`
    :param str version: String with version (date) for metadata
    :return: Path to chosen metadata file
    :rtype: str
    """
    if data_type == 'exomes':
        if version:
            return metadata_exomes_ht_path(version)
        return metadata_exomes_ht_path()
    elif data_type == 'genomes':
        if version:
            return metadata_genomes_ht_path(version)
        return metadata_genomes_ht_path()
    return DataException("Select data_type as one of 'genomes' or 'exomes'")


def raw_exomes_mt_path(hail_version=CURRENT_HAIL_VERSION):
    """
    Warning: unsplit and no special consideration on sex chromosomes
    """
    return 'gs://gnomad/raw/hail-{0}/mt/exomes/gnomad.exomes.mt'.format(hail_version)


def raw_genomes_mt_path(hail_version=CURRENT_HAIL_VERSION):
    """
    Warning: unsplit and no special consideration on sex chromosomes
    """
    return 'gs://gnomad/raw/hail-{0}/mt/genomes/gnomad.genomes.mt'.format(hail_version)


def raw_exac_mt_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad/raw/hail-{0}/mt/exac/exac.mt'.format(hail_version)


def exac_release_sites_mt_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad/raw/hail-{}/mt/exac/exac.r1.sites.vep.mt'.format(hail_version)


def hardcalls_mt_path(data_type, split=True, hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad/hardcalls/hail-{0}/mt/{1}/gnomad.{1}{2}.mt'.format(hail_version, data_type,
                                                                           "" if split else ".unsplit")


def non_refs_only_mt_path(data_type, split=True):
    return f'gs://gnomad/non_refs_only/hail-0.2/mt/{data_type}/gnomad.{data_type}{"" if split else ".unsplit"}.mt'


def annotations_ht_path(data_type, annotation_type, hail_version=CURRENT_HAIL_VERSION):
    """
    Get sites-level annotations
    :param str data_type: One of "exomes" or "genomes"
    :param str annotation_type: One of "vep", "qc_stats", "frequencies", "rf", "concordance"
    :param str hail_version: One of the HAIL_VERSIONs
    :return: Path to annotations Table
    :rtype: str
    """
    return 'gs://gnomad/annotations/hail-{0}/ht/{1}/gnomad.{1}.{2}.ht'.format(hail_version, data_type,
                                                                              annotation_type)


def sample_annotations_table_path(data_type, annotation_type, hail_version=CURRENT_HAIL_VERSION):
    """
    Get samples-level annotations
    :param str data_type: One of "exomes" or "genomes"
    :param str annotation_type: One of "family_stats", "downsampling"
    :param str hail_version: One of the HAIL_VERSIONs
    :return: Path to annotations VDS
    :rtype: str
    """
    return 'gs://gnomad/annotations/hail-{0}/sample_tables/{1}/gnomad.{1}.{2}.ht'.format(hail_version, data_type,
                                                                                          annotation_type)

gnomad_pca_mt_path = "gs://gnomad-genomes/sampleqc/gnomad.pca.mt"


def gnomad_public_pca_mt_path(version=CURRENT_RELEASE):
    """
    Returns the path for the public gnomAD VDS containing sites and loadings from the PCA
    :param str version: One of the RELEASEs
    :return: path to gnomAD public PCA VDS
    :rtype: str
    """
    return "gs://gnomad-public/release/{}/pca/gnomad_pca_loadings.mt".format(version)


def metadata_genomes_tsv_path(version=CURRENT_GENOME_META):
    return 'gs://gnomad/metadata/genomes/gnomad.genomes.metadata.{0}.tsv.bgz'.format(version)


def metadata_exomes_tsv_path(version=CURRENT_EXOME_META):
    return 'gs://gnomad/metadata/exomes/gnomad.exomes.metadata.{0}.tsv.bgz'.format(version)


def metadata_genomes_ht_path(version=CURRENT_GENOME_META):
    return 'gs://gnomad/metadata/genomes/gnomad.genomes.metadata.{0}.ht'.format(version)


def metadata_exomes_ht_path(version=CURRENT_EXOME_META):
    return 'gs://gnomad/metadata/exomes/gnomad.exomes.metadata.{0}.ht'.format(version)


def coverage_mt_path(data_type) -> str:
    return f'gs://gnomad/coverage/hail-0.2/coverage/{data_type}/mt/gnomad.{data_type}.coverage.mt'


def coverage_ht_path(data_type, by_population: bool = False, by_platform: bool = False) -> str:
    if by_population and by_population:
        raise DataException('Cannot assess coverage by both population and platform... yet...')
    by = '.population' if by_population else '.platform' if by_platform else ''
    return f'gs://gnomad/coverage/hail-0.2/coverage/{data_type}/ht/gnomad.{data_type}.coverage{by}.summary.ht'


def fam_path(data_type: str, version: str = CURRENT_FAM) -> str:
        return f"gs://gnomad/metadata/{data_type}/gnomad.{data_type}.{version}.fam"


def genomes_exomes_duplicate_ids_tsv_path(version: str = CURRENT_DUPS) -> str:
    return f"gs://gnomad/metadata/join/gnomad.genomes_exomes.{version}.duplicate_ids.tsv"


def omni_mt_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad-public/truth-sets/hail-{0}/1000G_omni2.5.b37.mt'.format(hail_version)


def mills_mt_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad-public/truth-sets/hail-{0}/Mills_and_1000G_gold_standard.indels.b37.mt'.format(hail_version)


def hapmap_mt_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad-public/truth-sets/hail-{0}/hapmap_3.3.b37.mt'.format(hail_version)


def kgp_high_conf_snvs_mt_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad-public/truth-sets/hail-{0}/1000G_phase1.snps.high_confidence.b37.mt'.format(hail_version)


def NA12878_mt_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad-public/truth-sets/hail-{0}/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.mt'.format(hail_version)


def syndip_mt_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad-public/truth-sets/hail-{0}/hybrid.m37m.mt'.format(hail_version)


def cpg_sites_mt_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad-public/resources/hail-{}/cpg.mt'.format(hail_version)


def methylation_sites_mt_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad-resources/methylation/hail-{}/methylation.ht'.format(hail_version)


dbsnp_vcf_path = "gs://gnomad-public/truth-sets/source/All_20160601.vcf.bgz"

NA12878_high_conf_regions_bed_path = "gs://gnomad-public/truth-sets/source/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.bed"
NA12878_high_conf_exome_regions_bed_path = "gs://gnomad-public/truth-sets/source/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.18_2mindatasets_5minYesNoRatio.bed"
syndip_high_conf_regions_bed_path = "gs://gnomad-public/truth-sets/source/hybrid.m37m.bed"
clinvar_tsv_path = "gs://gnomad-resources/clinvar/source/clinvar_alleles.single.b37.tsv.gz"
clinvar_mt_path = "gs://gnomad-resources/clinvar/hail-0.2/clinvar_alleles.single.b37.mt"

# Useful intervals
lcr_intervals_path = "gs://gnomad-public/intervals/LCR.GRCh37_compliant.interval_list"  # "gs://gnomad-public/intervals/LCR.interval_list"
decoy_intervals_path = "gs://gnomad-public/intervals/mm-2-merged.GRCh37_compliant.bed"  # "gs://gnomad-public/intervals/mm-2-merged.bed.gz"
purcell5k_intervals_path = "gs://gnomad-public/intervals/purcell5k.interval_list"
segdup_intervals_path = "gs://gnomad-public/intervals/hg19_self_chain_split_both.bed"

# Exome intervals
exomes_high_conf_regions_intervals_path = "gs://gnomad-public/intervals/exomes_high_coverage.auto.interval_list"
exome_calling_intervals_path = 'gs://gnomad-public/intervals/exome_calling_regions.v1.interval_list'
evaluation_intervals_path = 'gs://gnomad-public/intervals/exome_evaluation_regions.v1.noheader.interval_list'
high_coverage_intervals_path = 'gs://gnomad-public/intervals/high_coverage.auto.interval_list'

# Genome intervals
genome_evaluation_intervals_path = "gs://gnomad-public/intervals/hg19-v0-wgs_evaluation_regions.v1.interval_list"  # from Broad GP
genome_evaluation_intervals_path_hg38 = "gs://gnomad-public/intervals/hg38-v0-wgs_evaluation_regions.hg38.interval_list"
# More can be found at gs://broad-references/hg19

vep_config = "/vep/vep-gcloud.properties"

# Annotations
context_mt_path = 'gs://gnomad-resources/constraint/context_processed.mt'


class DataException(Exception):
    pass



#now testing the new function, which should be faster
def phase_sensitivity_fast(mt, windowsize=1, adj=True):  # trying to make the above faster.
    # takes matrix table that has PID, GT, PBT_GT, calculate the phase sensitivity, sum of all individuals
    # for window size k, get the result of window size=1, 2, ... k
    import pandas as pd
    mt = hl.filter_alleles(mt, lambda allele, i: hl.is_snp(mt.alleles[0], allele))  # currently take only SNP
    mt = mt.select_rows()  # throw away unwanted rows
    mt = mt.filter_entries(mt.GT.is_het())  # throw away unwanted entries (non alt)
    mt = hl.window_by_locus(mt, windowsize)
    mt = mt.filter_entries((hl.is_defined(mt.GT) & (mt.prev_entries.length() > 0)))
    mt = mt.filter_entries(mt.prev_entries.filter(lambda x: x.GT.is_het()).length() > 0)
    et = mt.entries()
    et = et.annotate(indices=hl.range(0, hl.len(et.prev_rows)))
    et = et.explode('indices')
    et = et.transmute(prev_row=et.prev_rows[et.indices], prev_entry=et.prev_entries[et.indices])
    et = et.filter(hl.is_defined(et.prev_entry.GT))  # and remove non-corresponding entries


    if adj:  # restrict to adj pass
        et = et.filter(et.adj & et.prev_entry.adj)
    print("\n et created and filtered. \n Starting to look at phase info \n" + tm.ctime())
    # annotate columns
    et = et.annotate(dist=et.locus.position - et.prev_row.locus.position,
                     pair_phased=(et.GT.phased) & (et.prev_entry.GT.phased),
                     has_PBT=(hl.is_defined(et.PBT_GT)) & (hl.is_defined(et.prev_entry.PBT_GT))
                     )
    et = et.annotate(is_mnv=(et.pair_phased & (et.PID == et.prev_entry.PID) & (et.GT == et.prev_entry.GT)))
    et = et.annotate(flipped_GT=hl.call(et.GT[1], et.GT[0], phased=et.GT.phased),
                     prev_entry_flipped_GT=hl.call(et.prev_entry.GT[1], et.prev_entry.GT[0],
                                                   phased=et.prev_entry.GT.phased))
    et = et.annotate(agrees_PBT=(((et.GT == et.PBT_GT) & (et.prev_entry.GT == et.prev_entry.PBT_GT)) |
                                 ((et.flipped_GT == et.PBT_GT) & (et.prev_entry_flipped_GT == et.prev_entry.PBT_GT))))
    et = et.annotate(agrees_PID=((et.pair_phased)&(et.PID == et.prev_entry.PID) & hl.is_defined(et.PID)))
    #agrees PID only if they are phased at all

    # define each categ
    et_has_PBT = et.filter(et.has_PBT)
    et_agrees_PBT = et.filter(et.agrees_PBT)
    et_phased = et.filter(et.pair_phased)
    et_phased_and_has_PBT = et_phased.filter(et_phased.has_PBT)
    et_phased_and_agrees_PBT = et_phased.filter(et_phased.agrees_PBT)
    et_same_PID = et.filter(et.agrees_PID)
    et_same_PID_and_has_PBT = et_same_PID.filter(et_same_PID.has_PBT)
    et_same_PID_and_agrees_PBT = et_same_PID.filter(et_same_PID.agrees_PBT)
    et_mnv = et.filter(et.is_mnv)
    et_mnv_and_has_PBT = et_mnv.filter(et_mnv.has_PBT)
    et_mnv_and_agrees_PBT = et_mnv.filter(et_mnv.agrees_PBT)
    print("Starting to aggregate \n" + tm.ctime())
    n_all = et.aggregate(hl.agg.counter(et.dist))
    n_has_PBT = et_has_PBT.aggregate(hl.agg.counter(et_has_PBT.dist))
    n_agrees_PBT = et_agrees_PBT.aggregate(hl.agg.counter(et_agrees_PBT.dist))
    n_phased = et_phased.aggregate(hl.agg.counter(et_phased.dist))
    n_phased_and_has_PBT = et_phased_and_has_PBT.aggregate(hl.agg.counter(et_phased_and_has_PBT.dist))
    n_phased_and_agrees_PBT = et_phased_and_agrees_PBT.aggregate(hl.agg.counter(et_phased_and_agrees_PBT.dist))
    n_same_PID = et_same_PID.aggregate(hl.agg.counter(et_same_PID.dist))
    n_same_PID_and_has_PBT = et_same_PID_and_has_PBT.aggregate(hl.agg.counter(et_same_PID_and_has_PBT.dist))
    n_same_PID_and_agrees_PBT = et_same_PID_and_agrees_PBT.aggregate(hl.agg.counter(et_same_PID_and_agrees_PBT.dist))
    n_mnv = et_mnv.aggregate(hl.agg.counter(et_mnv.dist))
    n_mnv_and_has_PBT = et_mnv_and_has_PBT.aggregate(hl.agg.counter(et_mnv_and_has_PBT.dist))
    n_mnv_and_agrees_PBT = et_mnv_and_agrees_PBT.aggregate(hl.agg.counter(et_mnv_and_agrees_PBT.dist))

    #also some missing: same PID and has PBT in general (not restricting to MNV) / those that agrees.
    #no we actually have it.


    print("Done aggregate \n" + tm.ctime())
    # and if we return these we are done
    df = pd.DataFrame(n_all, index=["n_all"])
    df2 = pd.DataFrame(n_has_PBT, index=["n_has_PBT"])
    df3 = pd.DataFrame(n_agrees_PBT, index=["n_agrees_PBT"])
    df4 = pd.DataFrame(n_phased, index=["n_phased"])
    df5 = pd.DataFrame(n_phased_and_has_PBT, index=["n_phased_and_has_PBT"])
    df6 = pd.DataFrame(n_phased_and_agrees_PBT, index=["n_phased_and_agrees_PBT"])
    df7 = pd.DataFrame(n_mnv, index=["n_mnv"])
    df8 = pd.DataFrame(n_mnv_and_has_PBT, index=["n_mnv_and_has_PBT"])
    df9 = pd.DataFrame(n_mnv_and_agrees_PBT, index=["n_mnv_and_agrees_PBT"])
    df10 = pd.DataFrame(n_same_PID, index=["n_same_PID"])
    df11 = pd.DataFrame(n_same_PID_and_has_PBT, index=["n_same_PID_and_has_PBT"])
    df12 = pd.DataFrame(n_same_PID_and_agrees_PBT, index=["n_same_PID_and_agrees_PBT"])
    print(n_all)
    return (pd.concat([df, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12]))
def phase_sensitivity_per_indv(mt, indv, windowsize=1):
    # takes matrix table that has PID, GT, PBT_GT, plus an individual to focus calculate the phase sensitivity, for that specific individual
    sample_to_keep = hl.literal({indv})
    mt = mt.filter_cols(sample_to_keep.contains(mt['s']))  # filter to that specific individual
    return (phase_sensitivity_fast(mt, windowsize=windowsize))  # and then run the same function
def pbt_phased_trios_mt_path(data_type: str, split: bool = True, hail_version: str = CURRENT_HAIL_VERSION):
    return "gs://gnomad/hardcalls/hail-{0}/mt/{1}/gnomad.{1}.trios.pbt_phased{2}.mt".format(hail_version, data_type,
                                                                                            "" if split else ".unsplit")
exomes = hl.read_matrix_table(pbt_phased_trios_mt_path("exomes"))
exomes = exomes.filter_cols(exomes.s == exomes.source_trio.proband.s)
df = phase_sensitivity_fast(exomes, windowsize=100) #should be dealable, for a single individual
print ("per indv exome done" + tm.ctime())
df["categ"] = df.index
hl.Table.from_pandas(df).export("gs://gnomad-qingbowang/MNV/phase_sensitivity_exome_proband_w100.tsv")

genomes = hl.read_matrix_table(pbt_phased_trios_mt_path("genomes"))
fam_ht = hl.import_fam(fam_path("genomes"), delimiter="\t") #for genomes, we need to annotate this
genomes = genomes.annotate_cols(source_trio = fam_ht[genomes.s])
genomes = genomes.filter_cols(hl.len(genomes.source_trio.fam_id)>0) #filtering to child only
df = phase_sensitivity_fast(genomes, windowsize=100) #should be dealable, for a single individual
print ("per indv genome done" + tm.ctime())
df["categ"] = df.index
hl.Table.from_pandas(df).export("gs://gnomad-qingbowang/MNV/phase_sensitivity_genome_proband_w100.tsv")
