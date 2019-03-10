# -*- coding: utf-8 -*-
__author__ = 'QingboWang'


#get the enrichment vs mnv_10

# -*- coding: utf-8 -*-
__author__ = 'QingboWang'

output_path = "gs://gnomad-qingbowang/MNV/wholegenome"
import pandas as pd
import numpy as np
#from matplotlib import pyplot as plt
#import seaborn as sns
import time as tm
from scipy import stats
import hail as hl
import hail.expr.aggregators as agg
from typing import *

CURRENT_HAIL_VERSION = "0.2"
CURRENT_RELEASE = "2.0.2"
CURRENT_GENOME_META = "2018-08-16"  # YYYY-MM-DD
CURRENT_EXOME_META = "2018-08-16"
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
    #from gnomad_hail.utils import filter_to_adj

    if raw and split:
        raise DataException('No split raw data. Use of hardcalls is recommended.')

    if non_refs_only:
        mt = hl.read_matrix_table(get_gnomad_data_path(data_type, split=split, non_refs_only=non_refs_only, hail_version=hail_version))
    else:
        mt = hl.read_matrix_table(get_gnomad_data_path(data_type, hardcalls=not raw, split=split, hail_version=hail_version))

    if adj:
        mt = filter_to_adj(mt)

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
            columns.extend(['diabetes', 'exac_joint', 'tcga'])
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


def pbt_phased_trios_mt_path(data_type: str, split: bool = True, hail_version : str = CURRENT_HAIL_VERSION):
    return "gs://gnomad/hardcalls/hail-{0}/mt/{1}/gnomad.{1}.trios.pbt_phased{2}.mt".format(hail_version, data_type,
                                                                           "" if split else ".unsplit")

def annotations_ht_path(data_type, annotation_type, hail_version=CURRENT_HAIL_VERSION):
    """
    Get sites-level annotations
    :param str data_type: One of "exomes" or "genomes"
    :param str annotation_type: One of "vep", "qc_stats", "frequencies", "rf", "omes_concordance", "NA12878_concordance", "syndip_concordance", "omes_by_platform_concordance"
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
    :param str annotation_type: One of "family_stats", "downsampling", "omes_concordance", "NA12878_concordance", "syndip_concordance"
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


def fam_path(data_type: str, version: str = CURRENT_FAM, true_trios: bool = False) -> str:
    if not true_trios:
        return f"gs://gnomad/metadata/{data_type}/gnomad.{data_type}.{version}.fam"
    else:
        return f"gs://gnomad/metadata/{data_type}/gnomad.{data_type}.{version}.true_trios.fam"


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


def kgp_phase3_genotypes_mt_path(split: bool = True, hail_version=CURRENT_HAIL_VERSION) -> str:
    """
    1000 Genomes Phase 3 with genotypes (b37)
    Imported from: gs://genomics-public-data/1000-genomes-phase-3/vcf-20150220/ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf
    Samples populations from: gs://gnomad-public/truth-sets/hail-0.2/1000G.GRCh38.20130502.phase3.sequence.index
    :param bool split: Whether to load to split or non-split version
    :param str hail_version: Hail version
    :return: Path to 1000 Genomes MT
    :rtype: str
    """
    return 'gs://gnomad-public/truth-sets/hail-{0}/1000Genomes_phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes{1}.mt'.format(hail_version, '.split' if split else '')


def NA12878_mt_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad-public/truth-sets/hail-{0}/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.mt'.format(hail_version)


def syndip_mt_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad-public/truth-sets/hail-{0}/hybrid.m37m.mt'.format(hail_version)


def cpg_sites_mt_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad-public/resources/hail-{}/cpg.mt'.format(hail_version)


def methylation_sites_mt_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad-resources/methylation/hail-{}/methylation.ht'.format(hail_version)


dbsnp_vcf_path = "gs://gnomad-public/truth-sets/source/All_20180423.vcf.bgz"
dbsnp_ht_path = "gs://gnomad-public/truth-sets/source/All_20180423.ht"

NA12878_high_conf_regions_bed_path = "gs://gnomad-public/truth-sets/source/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.bed"
NA12878_high_conf_exome_regions_bed_path = "gs://gnomad-public/truth-sets/source/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.18_2mindatasets_5minYesNoRatio.bed"
syndip_high_conf_regions_bed_path = "gs://gnomad-public/truth-sets/source/hybrid.m37m.bed"
clinvar_tsv_path = "gs://gnomad-resources/clinvar/source/clinvar_alleles.single.b37.tsv.bgz"
clinvar_mt_path = "gs://gnomad-resources/clinvar/hail-0.2/clinvar_alleles.single.b37.vep.mt"

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

vep_config = 'gs://hail-common/vep/vep/vep85-gcloud.json'

# Annotations
context_mt_path = 'gs://gnomad-resources/constraint/context_processed.mt'


# Sample QC files
def qc_mt_path(data_type: str):
    return 'gs://gnomad/sample_qc/mt/gnomad.{}.high_callrate_common_biallelic_snps.mt'.format(data_type)


def qc_ht_path(data_type: str):
    return 'gs://gnomad/sample_qc/ht/gnomad.{}.high_callrate_common_biallelic_snps.ht'.format(data_type)


def qc_temp_data_prefix(data_type: str):
    return 'gs://gnomad/sample_qc/temp/{0}/gnomad.{0}'.format(data_type)


def qc_meta_path(data_type: str):
    if data_type == 'exomes':
        return 'gs://gnomad/sample_qc/input_meta/gnomad.exomes.streamlined_metadata.2018-03-21.txt.bgz'
    else:
        return 'gs://gnomad/sample_qc/input_meta/gnomad.genomes.streamlined_metadata.2018-03-21.txt.bgz'



class DataException(Exception):
    pass


lowcv_reg = "gs://gnomad-qingbowang/MNV/cov_leq15_reg.bed"

categ = ["Coding_UCSC", "DHS_Trynka", "Enhancer_Hoffman", "H3K27ac_PGC2", "H3K4me1_Trynka", "H3K4me3_Trynka",
         "H3K9ac_Trynka", "Intron_UCSC", "TSS_Hoffman",
         "Promoter_UCSC", "Transcribed_Hoffman", "UTR_3_UCSC", "UTR_5_UCSC", "TFBS_ENCODE"]


comp = {}
comp["A"] = "T"
comp["T"] = "A"
comp["G"] = "C"
comp["C"] = "G"
comp["N"] = "N"

def revcomp(seq):
    out = ""
    for i in seq[::-1]:
        out = out + comp[i]
    return (out)

def collapse_crstb_to_revcomp(crstb):
    # collaspse to a table, instead of another matrix
    flt = crstb.stack().reset_index()
    flt.columns = ['refs', 'alts', 'cnt']
    for i in range(flt.shape[0]):
        if i in flt.index:  # if it hasn't been deleted yet
            refs_revcomp = revcomp(flt.refs[i])
            alts_revcomp = revcomp(flt.alts[i])
            ix_revcomp = flt[(flt.refs == refs_revcomp) & (flt.alts == alts_revcomp)].index[0]
            if i != ix_revcomp:  # if revcomp is not yourself
                flt.loc[i, "cnt"] = flt.loc[i, "cnt"] + flt.loc[ix_revcomp, "cnt"]
                flt.drop(ix_revcomp, inplace=True)
    flt = flt[flt.cnt > 0]  # deleting the no dNVs
    flt.reset_index(inplace=True)
    del flt["index"]
    return (flt)

def ht_cnt_mat_to_pd(ht_cnt_mat):
    pds = ht_cnt_mat.to_pandas()
    pds.fillna(value=0, inplace=True)
    pds = pds.iloc[:, 2:]
    pds = pd.pivot_table(pds, index="refs", columns="alts")
    pds.fillna(value=0, inplace=True)
    pds.columns = pds.columns.get_level_values(1)  # to set the columns name properly.
    return (pds) #returning the pandas table

def get_density_enrichment(mnvs_mat, hg19_mat):
    df_pval = pd.DataFrame(index=mnvs_mat.index, columns = mnvs_mat.columns)
    df_ratio = pd.DataFrame(index=mnvs_mat.index, columns = mnvs_mat.columns)
    for i in mnvs_mat.index:
        for j in mnvs_mat.columns:
            if np.isnan(mnvs_mat.loc[i,j]):
                df_pval.loc[i,j] = np.nan
                df_pval.loc[i,j] = np.nan
                print ("{0},{1} done. was nan".format(i,j))
            else:
                refs = i[:2]
                x1 = mnvs_mat.loc[i,j]#case1
                x2 = hg19_mat.loc[refs,str(j)]#cont1
                y1 = mnvs_mat.loc[i,9]#case2
                y2 = hg19_mat.loc[refs,str(9)]#cont2
                oddsratio, pvalue = stats.fisher_exact([[x1, x2], [y1, y2]])
                print (pvalue)
                print (oddsratio)
                print (type(pvalue))
                print (type(oddsratio))
                if isinstance(pvalue, np.float64):
                    df_pval.loc[i, j] = (pvalue).item()  # .item() turns to hail compatible float
                else: df_pval.loc[i, j] = (pvalue)
                if isinstance(oddsratio, np.float64):
                    df_ratio.loc[i,j] = (oddsratio).item() #because this is what hail takes. not 64 but 32..
                else: df_ratio.loc[i, j] = (oddsratio)
                #only for pvalue, it is normal float, when the pvalue is 1. else, everything is np.float64 and is not compatible
                print ("{0},{1} done".format(i,j))
    return ((df_pval, df_ratio))


def hl_strc_to_pd_df(strc,lname):
    #from hail struct, create a dataframe of single line
    dict = {}
    for k in strc.keys():
        dict[k] = strc[k]
    return (pd.DataFrame(dict, index=[lname]))

def read_cnt_matrix(tsv_path):
    with hl.hadoop_open(tsv_path, 'r') as f:
        df = pd.read_table(f)
    df.index = df.columns
    return (df)


#from functions import *
hl.init()

#create a dataframe of row=(refs,alts), col=(Ns in middle)
cnts = pd.DataFrame()
for d in range(1,11):
    c = read_cnt_matrix("{0}/cnt_mat_d{1}.tsv".format(output_path, str(d)))
    c.columns = c.columns.str.replace("N","") #since N is not needed in this context..
    c.index = c.columns
    c = collapse_crstb_to_revcomp(c)
    c.index = c.refs + "->" + c.alts
    del c["refs"]
    del c["alts"]
    c.columns = [d-1]
    if d==1:#first, set colnames etc
        inds = c.index
    cnts = pd.concat([cnts,c.loc[inds,d-1]], axis=1, sort=True) #can be pd.concat([cnts,c.loc[:,"cnt"]]) if large enough file.

#import the denominator
dnm = hl.import_table("gs://gnomad-qingbowang/MNV/background_dist_all.tsv", impute=True).to_pandas()
dnm.index = dnm[""]
del dnm[""]

m = get_density_enrichment(cnts.fillna(value=0), dnm) #fill na for now. with larger set, we don't need to

#and write to table. again, as tsv rather than as ht because they are small enough
vsmnv10_pval = hl.Table.from_pandas(m[0]).export("{0}/vsmnv10_pval.tsv".format(output_path))
vsmnv10_ratio = hl.Table.from_pandas(m[1]).export("{0}/vsmnv10_ratio.tsv".format(output_path))


