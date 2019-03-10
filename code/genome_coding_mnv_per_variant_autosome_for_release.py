# -*- coding: utf-8 -*-
__author__ = 'QingboWang'




#first step -- get the variants
#neglecting the partition problem, first work on unphased problem.
#from functions import *


output_path = "gs://gnomad-qingbowang/MNV/1206_genome"
from typing import *
import pandas as pd
import time as tm

""" not used
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import time as tm
from scipy import stats
"""
import hail as hl
import hail.expr.aggregators as agg
from typing import *

#this is renewed to the newest version

CURRENT_HAIL_VERSION = "0.2"
CURRENT_RELEASE = "2.1"
CURRENT_GENOME_META = "2018-10-11"  # YYYY-MM-DD
CURRENT_EXOME_META = "2018-10-11"
CURRENT_FAM = '2018-04-12'
CURRENT_DUPS = '2017-10-04'

RELEASES = ["2.0.1", "2.0.2", "2.1"]

SUBPOPS = {'NFE': ['BGR', 'EST', 'NWE', 'SEU', 'SWE', 'ONF'],
           'EAS': ['KOR', 'JPN', 'OEA']
           }
GENOME_POPS = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH']
EXOME_POPS = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS']
EXAC_POPS = ["AFR", "AMR", "EAS", "FIN", "NFE", "OTH", "SAS"]


def public_exomes_ht_path(split=True, version=CURRENT_RELEASE):
    if version == '2.1':
        return f'gs://gnomad-public/release/{version}/ht/exomes/gnomad.exomes.r{version}.sites.ht'
    else:
        return 'gs://gnomad-public/release/{0}/vds/exomes/gnomad.exomes.r{0}.sites{1}.vds'.format(version, ".split" if split else "")


def public_genomes_ht_path(split=True, version=CURRENT_RELEASE):
    if version == '2.1':
        return f'gs://gnomad-public/release/{version}/ht/genomes/gnomad.genomes.r{version}.sites.ht'
    else:
        return 'gs://gnomad-public/release/{0}/vds/genomes/gnomad.genomes.r{0}.sites{1}.vds'.format(version, ".split" if split else "")


def get_gnomad_public_data(data_type, split=True, version=CURRENT_RELEASE):
    """
    Wrapper function to get public gnomAD data as VDS.
    :param str data_type: One of `exomes` or `genomes`
    :param bool split: Whether the dataset should be split
    :param str version: One of the RELEASEs
    :return: Chosen VDS
    :rtype: MatrixTable
    """
    return hl.read_table(get_gnomad_public_data_path(data_type, split=split, version=version))


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
        sites_ht = get_gnomad_public_data(data_type, split)
        mt = mt.select_rows(**sites_ht[mt.row_key])

    mt = mt.select_globals()  # Required since a backward-incompatible change in Hail

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
                   'neuro', 'control', 'topmed',
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
        return public_exomes_ht_path(split, version)
    elif data_type == 'genomes':
        return public_genomes_ht_path(split, version)
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


def coverage_mt_path(data_type, by_population: bool = False, by_platform: bool = False) -> str:
    by = '.population' if by_population else ''
    by += '.platform' if by_platform else ''
    return f'gs://gnomad/coverage/hail-0.2/coverage/{data_type}/mt/gnomad.{data_type}.coverage{by}.mt'


def coverage_ht_path(data_type, by_population: bool = False, by_platform: bool = False) -> str:
    by = '.population' if by_population else ''
    by += '.platform' if by_platform else ''
    return f'gs://gnomad/coverage/hail-0.2/coverage/{data_type}/ht/gnomad.{data_type}.coverage{by}.summary.ht'


def fam_path(data_type: str, version: str = CURRENT_FAM, true_trios: bool = False) -> str:
    """
    Returns the path to gnomAD pedigree file.
    :param str data_type: One of 'exomes' or 'genomes'
    :param str version: Version of the fam file to get
    :param bool true_trios: If set, removes all families with more than one offspring
    :return: Path to fam file
    :rtype: str
    """
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
clinvar_loftee_beta_mt_path = "gs://gnomad-resources/clinvar/hail-0.2/clinvar_alleles.single.b37.loftee.beta.vep.mt"

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

vep_config = 'gs://hail-common/vep/vep/vep85-loftee-gcloud.json'

# Annotations
context_mt_path = 'gs://gnomad-resources/context/hail-0.2/context_processed.mt'


# Sample QC files
def qc_mt_path(data_type: str):
    return 'gs://gnomad/sample_qc/mt/gnomad.{}.high_callrate_common_biallelic_snps.mt'.format(data_type)


def qc_ht_path(data_type: str):
    return 'gs://gnomad/sample_qc/ht/gnomad.{}.high_callrate_common_biallelic_snps.ht'.format(data_type)


def qc_temp_data_prefix(data_type: str):
    return 'gs://gnomad/sample_qc/temp/{0}/gnomad.{0}'.format(data_type)


def qc_meta_path(data_type: str):
    if data_type == 'exomes':
        return 'gs://gnomad/sample_qc/input_meta/gnomad.exomes.streamlined_metadata.2018-10-10.txt.bgz'
    else:
        return 'gs://gnomad/sample_qc/input_meta/gnomad.genomes.streamlined_metadata.2018-10-11.txt.bgz'



def filter_to_adj(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Filter genotypes to adj criteria
    """
    if 'adj' not in list(mt.entry):
        mt = annotate_adj(mt)
    mt = mt.filter_entries(mt.adj)
    return mt.drop(mt.adj)

def annotate_adj(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Annotate genotypes with adj criteria (assumes diploid)
    """
    adj_gq = 20
    adj_dp = 10
    adj_ab = 0.2

    return mt.annotate_entries(adj=(mt.GQ >= adj_gq) & (mt.DP >= adj_dp) & (
                                   ~mt.GT.is_het() |
                                   ((mt.GT[0] == 0) & (mt.AD[mt.GT[1]] / mt.DP >= adj_ab)) |
                                   ((mt.GT[0] > 0) & (mt.AD[mt.GT[0]] / mt.DP >= adj_ab) &
                                    (mt.AD[mt.GT[1]] / mt.DP >= adj_ab))
                               ))




class DataException(Exception):
    pass




lowcv_reg = "gs://gnomad-qingbowang/MNV/cov_leq15_reg.bed"



def hl_strc_to_pd_df(strc,lname):
    #from hail struct, create a dataframe of single line
    dict = {}
    for k in strc.keys():
        dict[k] = strc[k]
    return (pd.DataFrame(dict, index=[lname]))




hl.init(tmp_dir="gs://gnomad-qingbowang/tmp")
import hail.expr.aggregators as agg

#get gnomAD genome
mt_all = get_gnomad_data("genomes", release_samples=True, adj=True, release_annotations=True)

exomereg = hl.import_bed("gs://gnomad-qingbowang/MNV/agg_stats/exome_calling_regions.v1.asbed.bed",
                         skip_invalid_intervals=True)
mt_all = mt_all.filter_rows(hl.is_defined(exomereg[mt_all.locus]))

###check whether this works or not


#let's do per chromosome as well.
#also sex chromosome, later.

for chr in range(22,0,-1):
    #filter to that range
    chr = str(chr)
    import time as tm
    print ("starting chr{0}".format(chr))
    print (tm.ctime())
    #repartition -actually not needed. 10000 from the beginning.
    #mt = hl.filter_intervals(mt_all, [hl.parse_locus_interval(chr)])
    mt = hl.filter_intervals(mt_all, [hl.parse_locus_interval(chr)])
    #repartition -- result is subjective to this repartitioning...............
    #so I will start with 200 partitions first. just as inquery.
    #mt = mt.repartition(1000)#let's do no repartitioning here
    #keep also AF etc info
    mt = mt.select_cols()
    mt = mt.select_rows(mt.freq, mt.filters)
    mt = mt.annotate_rows(AC = mt.freq[0].AC, AF = mt.freq[0].AF) #re-annotating the AF/AC
    #and delete the "freq" -> this makes things lighter, hopefully
    mt = mt.select_rows(mt.AC, mt.AF, mt.filters)
    mt = mt.filter_entries(mt.GT.is_non_ref())
    mt = hl.window_by_locus(mt, 2) #partition in window -- only within codon reading frame
    mt = mt.filter_entries((hl.is_defined(mt.GT) & (mt.prev_entries.length() > 0))) #throwing away no MNV SNPs
    mt = mt.filter_entries(mt.prev_entries.filter(lambda x: x.GT.is_non_ref()).length() > 0) #same
    et = mt.key_cols_by().entries() # Matrix with 1000 rows (variant) + 1000 cols (sample)=> 1 million entries
    et = et.annotate(indices = hl.range(0, hl.len(et.prev_rows)))
    et = et.explode('indices')
    et = et.transmute(prev_row = et.prev_rows[et.indices],
                      prev_entry = et.prev_entries[et.indices])    
    et = et.annotate(dist=et.locus.position - et.prev_row.locus.position) #annotating the distance
    #et.cache() #should make everything faster -> no, actually seems like making it slower..
    #and annotate 3 types of MNV here.
    et = et.annotate(homhom = (et.GT.is_hom_var() & (et.prev_entry.GT.is_hom_var())),
                     hethom = ((et.GT.is_hom_var() & et.prev_entry.GT.is_het_ref()) | (et.GT.is_het_ref() & et.prev_entry.GT.is_hom_var())),#including hom-het, just not distinguishing them two.
                     hethet = ( hl.is_defined(et.PID) & hl.is_defined(et.prev_entry.PID) & (et.PID==et.prev_entry.PID) & (et.GT.phased&et.prev_entry.GT.phased) & (et.GT.is_het_ref()&et.prev_entry.GT.is_het_ref()) & (et.GT==et.prev_entry.GT) ),
                     hethet2_cand = ( hl.is_defined(et.PID) & hl.is_defined(et.prev_entry.PID) & (et.PID==et.prev_entry.PID) & (et.GT.phased) & (et.GT.is_het_ref()&et.prev_entry.GT.is_het_ref()) )
                     )#There are also some candidates at this stage (for those PID edge -> prevGT might not be phased / GT!=prevGT possible)
    
    #het x het (excluding het het 2)

    et_het = et.filter(et.hethet)
    #et_het = et_het.repartition(200) #again, no more repartitioning
    et_het.write("{0}/tmp_MNV_genome_chr{1}_et_het.ht".format(output_path, chr))
    et_het = hl.read_table("{0}/tmp_MNV_genome_chr{1}_et_het.ht".format(output_path, chr))
    per_variant_het = et_het.group_by('locus', 'alleles', "prev_row").aggregate(n=hl.agg.count())
    #and we can annotate back AF, AC, filter
    et_het = et_het.key_by("locus", "alleles", "prev_row")
    per_variant_het = per_variant_het.annotate(dist = et_het[per_variant_het.key].dist,
                          AF = et_het[per_variant_het.key].AF,
                          AC = et_het[per_variant_het.key].AC,
                          filters = et_het[per_variant_het.key].filters)
    per_variant_het = per_variant_het.annotate(prev_locus = per_variant_het.prev_row.locus,
                          prev_alleles = per_variant_het.prev_row.alleles,
                          prev_filters = per_variant_het.prev_row.filters,
                          prev_AC = per_variant_het.prev_row.AC,
                          prev_AF = per_variant_het.prev_row.AF)

    #and filter out non SNPs
    per_variant_het = per_variant_het.filter((per_variant_het.alleles[0].length() == 1) & (per_variant_het.alleles[1].length() == 1)\
             & (per_variant_het.prev_alleles[0].length() == 1) & (per_variant_het.prev_alleles[1].length() == 1))

    #and filter out d=0
    per_variant_het = per_variant_het.filter(per_variant_het.dist!=0)

    per_variant_het = per_variant_het.key_by()
    per_variant_het = per_variant_het.drop("prev_row") #dropping off unnecessaries
    import time as tm
    print ("start writing het het")
    print (tm.ctime())
    per_variant_het.write("{0}/MNV_genome_chr{1}_het.ht".format(output_path,chr))
    print ("wrote het het")
    print (tm.ctime())

    #hom x hom
    et_hom_hom = et.filter(et.homhom)
    #et_hom_hom = et_hom_hom.repartition(200)
    et_hom_hom.write("{0}/tmp_MNV_genome_chr{1}_et_hom_hom.ht".format(output_path, chr))
    et_hom_hom = hl.read_table("{0}/tmp_MNV_genome_chr{1}_et_hom_hom.ht".format(output_path, chr))
    per_variant_hom_hom = et_hom_hom.group_by('locus', 'alleles', "prev_row").aggregate(n=hl.agg.count())
    #and we can annotate back AF, AC, filter, rf_filter
    et_hom_hom = et_hom_hom.key_by("locus", "alleles", "prev_row")
    per_variant_hom_hom = per_variant_hom_hom.annotate(dist = et_hom_hom[per_variant_hom_hom.key].dist,
                                                  AF = et_hom_hom[per_variant_hom_hom.key].AF,
                                                  AC = et_hom_hom[per_variant_hom_hom.key].AC,
                                                  filters = et_hom_hom[per_variant_hom_hom.key].filters)
    per_variant_hom_hom = per_variant_hom_hom.annotate(prev_locus = per_variant_hom_hom.prev_row.locus,
                                                  prev_alleles = per_variant_hom_hom.prev_row.alleles,
                                                  prev_filters = per_variant_hom_hom.prev_row.filters,
                                                  prev_AC = per_variant_hom_hom.prev_row.AC,
                                                  prev_AF = per_variant_hom_hom.prev_row.AF)
    #and filter out non SNPs
    per_variant_hom_hom = per_variant_hom_hom.filter((per_variant_hom_hom.alleles[0].length() == 1) & (per_variant_hom_hom.alleles[1].length() == 1)\
                     & (per_variant_hom_hom.prev_alleles[0].length() == 1) & (per_variant_hom_hom.prev_alleles[1].length() == 1))
    per_variant_hom_hom = per_variant_hom_hom.filter(per_variant_hom_hom.dist != 0)

    per_variant_hom_hom = per_variant_hom_hom.key_by()
    per_variant_hom_hom = per_variant_hom_hom.drop("prev_row") #dropping off unnecessaries
    import time as tm
    print ("start writing hom hom")
    print (tm.ctime())
    per_variant_hom_hom.write("{0}/MNV_genome_chr{1}_hom_hom.ht".format(output_path, chr))
    print ("wrote hom hom")
    import time as tm
    print (tm.ctime())

    #het x hom, hom x het
    et_partially_hom = et.filter(et.hethom)

    #et_partially_hom = et_partially_hom.repartition(200)
    et_partially_hom.write("{0}/tmp_MNV_genome_chr{1}_et_partially_hom.ht".format(output_path, chr))
    et_partially_hom = hl.read_table("{0}/tmp_MNV_genome_chr{1}_et_partially_hom.ht".format(output_path, chr))

    per_variant_partially_hom = et_partially_hom.group_by('locus', 'alleles', "prev_row").aggregate(n=hl.agg.count())
    #and we can annotate back AF, AC, filter, rf_filter
    et_partially_hom = et_partially_hom.key_by("locus", "alleles", "prev_row")
    per_variant_partially_hom = per_variant_partially_hom.annotate(dist = et_partially_hom[per_variant_partially_hom.key].dist,
                                                  AF = et_partially_hom[per_variant_partially_hom.key].AF,
                                                  AC = et_partially_hom[per_variant_partially_hom.key].AC,
                                                  filters = et_partially_hom[per_variant_partially_hom.key].filters)
    per_variant_partially_hom = per_variant_partially_hom.annotate(prev_locus = per_variant_partially_hom.prev_row.locus,
                                                  prev_alleles = per_variant_partially_hom.prev_row.alleles,
                                                  prev_filters = per_variant_partially_hom.prev_row.filters,
                                                  prev_AC = per_variant_partially_hom.prev_row.AC,
                                                  prev_AF = per_variant_partially_hom.prev_row.AF)
    per_variant_partially_hom = per_variant_partially_hom.filter((per_variant_partially_hom.alleles[0].length() == 1) & (per_variant_partially_hom.alleles[1].length() == 1)\
                     & (per_variant_partially_hom.prev_alleles[0].length() == 1) & (per_variant_partially_hom.prev_alleles[1].length() == 1))
    per_variant_partially_hom = per_variant_partially_hom.filter(per_variant_partially_hom.dist != 0)
    per_variant_partially_hom = per_variant_partially_hom.key_by()
    per_variant_partially_hom = per_variant_partially_hom.drop("prev_row") #dropping off unnecessaries
    print ("start writing partially hom")
    print (tm.ctime())
    per_variant_partially_hom.write("{0}/MNV_genome_chr{1}_partially_hom.ht".format(output_path,chr))
    print ("wrote partially hom")
    import time as tm
    print (tm.ctime())


    #het het, PID edge unphased case
    et_het2 = et.filter(et.hethet2_cand)
    et_het2 = et_het2.annotate(PID_pivot=et_het2.prev_entry.PID.split("_")[0])
    et_het2 = et_het2.annotate(pos_str=hl.format('%s', et_het2.prev_row.locus.position))
    et_het2 = et_het2.annotate(is_edge=(et_het2.PID_pivot == et_het2.pos_str))
    et_het2 = et_het2.filter(et_het2.is_edge) #filtering to those whose prev_entry is on the edge
    et_het2 = et_het2.annotate(edge_is_phased=et_het2.prev_entry.GT.phased)
    et_het2 = et_het2.filter(~et_het2.edge_is_phased) #filter to those whose edge is unphased (potentially missing ones)
    et_het2 = et_het2.annotate(prev_GT = hl.cond(et_het2.prev_entry.GT.phased, et_het2.prev_entry.GT, hl.call(0,1, phased=True)))
    #if unphased, forcing it to be phased (and it is 0|1, not 1|0 by definition)
    et_het2 = et_het2.filter(et_het2.GT==et_het2.prev_GT)
    #et_het2 = et_het2.repartition(200)
    et_het2.write("{0}/tmp_MNV_genome_chr{1}_et_het2.ht".format(output_path, chr))
    et_het2 = hl.read_table("{0}/tmp_MNV_genome_chr{1}_et_het2.ht".format(output_path, chr))
    per_variant_het2 = et_het2.group_by('locus', 'alleles', "prev_row").aggregate(n=hl.agg.count())
    # and we can annotate back AF, AC, filter
    et_het2 = et_het2.key_by("locus", "alleles", "prev_row")
    per_variant_het2 = per_variant_het2.annotate(dist=et_het2[per_variant_het2.key].dist,
                                               AF=et_het2[per_variant_het2.key].AF,
                                               AC=et_het2[per_variant_het2.key].AC,
                                               filters=et_het2[per_variant_het2.key].filters)
    per_variant_het2 = per_variant_het2.annotate(prev_locus=per_variant_het2.prev_row.locus,
                                               prev_alleles=per_variant_het2.prev_row.alleles,
                                               prev_filters=per_variant_het2.prev_row.filters,
                                               prev_AC=per_variant_het2.prev_row.AC,
                                               prev_AF=per_variant_het2.prev_row.AF)    
    # and filter out non SNPs
    per_variant_het2 = per_variant_het2.filter(
        (per_variant_het2.alleles[0].length() == 1) & (per_variant_het2.alleles[1].length() == 1) \
        & (per_variant_het2.prev_alleles[0].length() == 1) & (per_variant_het2.prev_alleles[1].length() == 1))
    per_variant_het2 = per_variant_het2.filter(per_variant_het2.dist != 0)
    per_variant_het2 = per_variant_het2.key_by()
    per_variant_het2 = per_variant_het2.drop("prev_row")  # dropping off unnecessaries
    import time as tm
    
    print("start writing het2 het2")
    print(tm.ctime())
    per_variant_het2.write("{0}/MNV_genome_chr{1}_het2.ht".format(output_path, chr))
    print("wrote het2 het2")
    print(tm.ctime())

#assembl to a single file, filter to SNP only / filter pass only, and write
for chr in range(22,0,-1):
    chr = str(chr)
    het = hl.read_table("{0}/MNV_genome_chr{1}_het.ht".format(output_path, chr))
    het2 = hl.read_table("{0}/MNV_genome_chr{1}_het2.ht".format(output_path, chr))
    het_hom = hl.read_table("{0}/MNV_genome_chr{1}_partially_hom.ht".format(output_path, chr))
    hom = hl.read_table("{0}/MNV_genome_chr{1}_hom_hom.ht".format(output_path, chr))
    comb = het.key_by("locus", "alleles","prev_locus","prev_alleles","dist","AF","AC","filters","prev_AF","prev_AC","prev_filters") \
        .join(het2.key_by("locus", "alleles", "prev_locus", "prev_alleles", "dist", "AF", "AC", "filters", "prev_AF","prev_AC", "prev_filters"), how='outer') \
        .join(het_hom.key_by("locus", "alleles","prev_locus","prev_alleles","dist","AF","AC","filters","prev_AF","prev_AC","prev_filters"), how='outer')\
        .join(hom.key_by("locus", "alleles","prev_locus","prev_alleles","dist","AF","AC","filters","prev_AF","prev_AC","prev_filters"), how='outer')
    comb = comb.transmute(n_hethet=hl.or_else(comb.n, 0), n_hethet2=hl.or_else(comb.n_1, 0), n_hethom=hl.or_else(comb.n_2, 0), n_homhom=hl.or_else(comb.n_3, 0))
    comb = comb.select("n_hethet","n_hethet2", "n_hethom","n_homhom")

    print ("start writing chr{0}".format(chr))
    print (tm.ctime())
    comb.write("{0}/MNV_genome_chr{1}_combined.ht".format(output_path, chr))
    comb.export("{0}/MNV_genome_chr{1}_combined.tsv".format(output_path, chr))
    print("wrote chr{0}".format(chr))
    print(tm.ctime())



#and still need to annotate the downstream, for excluding no codon ones.
def annotate_vep_mnv(t0, block_size=100,dist=1):
    #t0: per variant table, already properly filtered
    #also filter to SNP / QC pass only -- no that's already done
        t = t0.filter(((t0.locus.position - t0.prev_locus.position) == dist))  # filter to that specific distance
        t0.unpersist()
        #t = t.repartition(16, shuffle=False)  # 4 didn't work, 8 as well
        t = t.repartition(8, shuffle=False)  # hoping that this repartition helps at all -> worked. too large repartition will kill things
        #t = t.repartition(20)  # hoping that this repartition helps at all -> worked. too large repartition will kill things
        #t = t.repartition(40)  # hoping that this repartition helps at all -> worked. too large repartition will kill things
        if dist==1:#looking at the case of dnv
            t = t.annotate(refs=t.prev_alleles[0] +t.alleles[0], alts=t.prev_alleles[1] +t.alleles[1])  # annotate combined refs
            #and make it a single thing
        if dist==2:#looking at the case of gMNV
            t = t.annotate(refs = t.prev_alleles[0] + t.prev_locus.sequence_context(before=-1, after=1) + t.alleles[0],
                           alts = t.prev_alleles[1] + t.prev_locus.sequence_context(before=-1, after=1) + t.alleles[1])
        t = t.annotate(mnv_alleles = [t.refs,t.alts])
        print ("annotation done. d={0}".format(dist))
        #vep the SNP2 (defined as 5'-SNP1-SNP2-3')
        #which is the original locus, alleles
        t = t.key_by("locus","alleles")
        t = hl.vep(t, vep_config, name="snp2_vep", block_size=block_size)
        print ("SNP2 vep done")
        #vep the SNP1
        t = t.key_by()
        t = t.rename({'locus' : 'snp2_locus', 'alleles' : 'snp2_alleles',"prev_locus":"locus", "prev_alleles":"alleles"})
        t = t.key_by('locus', 'alleles') #and re-key
        t = hl.vep(t, vep_config, name="snp1_vep", block_size=block_size)
        print ("SNP1 vep done")
        #vep the MNV
        t = t.key_by()
        t = t.rename({'alleles' : 'snp1_alleles',"mnv_alleles":"alleles"})
        t = t.key_by('locus', 'alleles') #and re-key
        t = hl.vep(t, vep_config, name="mnv_vep", block_size=block_size)
        print ("MNV vep done")
        return (t)
        #return as a first step (will annotate in detail later)
        #ここでtを再利用しているからだめなのか..?


def mnv_category_by_aa_change(snp1_con, snp2_con, mnv_con,aa1,aa2,aa3):
    if (snp1_con, snp2_con, mnv_con)==("synonymous_variant","missense_variant","missense_variant"):
        if aa2==aa3: return ("Unchanged")
        else: return ("Changed missense")
    elif (snp1_con, snp2_con, mnv_con)==("missense_variant","synonymous_variant","missense_variant"):
        if aa1==aa3: return ("Unchanged")
        else: return ("Changed missense")
    elif (snp1_con, snp2_con, mnv_con)==("missense_variant","missense_variant","missense_variant"):
        if ((aa1==aa3) or (aa2==aa3)):
            return ("Partially changed missense")
        else: return ("Changed missense")
    else: return ("something wrong going on")

def cons_term_most_severe(cons_term_array):
    if "start_lost" in cons_term_array: return ("start_lost") #by definition this is sufficient to determine the mnv is uninteresting
    elif "stop_lost" in cons_term_array: return ("stop_lost")
    elif "stop_gained" in cons_term_array: return ("stop_gained")
    elif "missense_variant" in cons_term_array: return ("missense_variant")
    elif "stop_retained_variant" in cons_term_array: return ("stop_retained_variant")
    elif "synonymous_variant" in cons_term_array: return ("synonymous_variant")
    else: return ("Noncoding_or_else")

def mnv_category(snp1_con, snp2_con, mnv_con, aa1, aa2, aa3):
    # return the MNV consequence category, such as gained PTV, unchange, etc.
    # just classify everything of 3*3*3=27 pattern.
    #plus, case where we see start_lost / stop_lost
    # and if undeterminisitic, look at the aa change and determine.
    if snp1_con == "synonymous_variant":
        if snp2_con == "synonymous_variant":
            if mnv_con == "synonymous_variant":
                return ("Unchanged")
            elif mnv_con == "missense_variant":
                return ("Gained missense")
            elif mnv_con == "stop_gained":
                return ("Gained PTV")
            else:
                return ("Noncoding_or_else")
        if snp2_con == "missense_variant":
            if mnv_con == "synonymous_variant":
                return ("Lost missense")
            elif mnv_con == "missense_variant":
                return (mnv_category_by_aa_change(snp1_con, snp2_con, mnv_con, aa1, aa2, aa3))  # go look at aa change.
            elif mnv_con == "stop_gained":
                return ("Gained PTV")
            else:
                return ("Noncoding_or_else")
        if snp2_con == "stop_gained":
            if mnv_con == "synonymous_variant":
                return ("Rescued PTV")
            elif mnv_con == "missense_variant":
                return ("Rescued PTV")
            elif mnv_con == "stop_gained":
                return ("Unchanged")
            else:
                return ("Noncoding_or_else")
    if snp1_con == "missense_variant":
        if snp2_con == "synonymous_variant":
            if mnv_con == "synonymous_variant":
                return ("Lost missense")
            elif mnv_con == "missense_variant":
                return (mnv_category_by_aa_change(snp1_con, snp2_con, mnv_con, aa1, aa2, aa3))  # go look at aa change.
            elif mnv_con == "stop_gained":
                return ("Gained PTV")
            else:
                return ("Noncoding_or_else")
        if snp2_con == "missense_variant":
            if mnv_con == "synonymous_variant":
                return ("Lost missense")
            elif mnv_con == "missense_variant":
                return (mnv_category_by_aa_change(snp1_con, snp2_con, mnv_con, aa1, aa2, aa3))  # go look at aa change.
            elif mnv_con == "stop_gained":
                return ("Gained PTV")
            else:
                return ("Noncoding_or_else")
        if snp2_con == "stop_gained":
            if mnv_con == "synonymous_variant":
                return ("Rescued PTV")
            elif mnv_con == "missense_variant":
                return ("Rescued PTV")
            elif mnv_con == "stop_gained":
                return ("Unchanged")
            else:
                return ("Noncoding_or_else")
        else:
            return ("Noncoding_or_else")
    if snp1_con == "stop_gained":
        if snp2_con == "synonymous_variant":
            if mnv_con == "synonymous_variant":
                return ("Rescued PTV")
            elif mnv_con == "missense_variant":
                return ("Rescued PTV")
            elif mnv_con == "stop_gained":
                return ("Unchanged")
            else:
                return ("Noncoding_or_else")
        if snp2_con == "missense_variant":
            if mnv_con == "synonymous_variant":
                return ("Rescued PTV")
            elif mnv_con == "missense_variant":
                return ("Rescued PTV")
            elif mnv_con == "stop_gained":
                return ("Unchanged")
            else:
                return ("Noncoding_or_else")
        if snp2_con == "stop_gained":
            if mnv_con == "synonymous_variant":
                return ("Rescued PTV")
            elif mnv_con == "missense_variant":
                return ("Rescued PTV")
            elif mnv_con == "stop_gained":
                return ("Unchanged")
            else:
                return ("Noncoding_or_else")
        else:
            return ("Noncoding_or_else")
    #else, involving start_loss etc -> look at mnv cons first.
    elif mnv_con=="start_lost": return "Unchanged" #by definition individual effect is also start loss
    elif mnv_con=="stop_lost":
        if ((snp1_con=="stop_retained_variant") & (snp2_con=="stop_retained_variant")): return "gained_stop_loss"
        else: return ("Unchanged")
    elif mnv_con=="stop_retained_variant":#this case, by definition one of the variant is stop_lost, and the other is stop_retained
        return ("Rescued stop loss")
    else:
        return ("Noncoding_or_else")
def filter_vep_to_canonical_transcripts(mt: Union[hl.MatrixTable, hl.Table],
                                        vep_root: str = 'vep') -> Union[hl.MatrixTable, hl.Table]:
    canonical = mt[vep_root].transcript_consequences.filter(lambda csq: csq.canonical == 1)
    vep_data = mt[vep_root].annotate(transcript_consequences=canonical)
    return mt.annotate_rows(**{vep_root: vep_data}) if isinstance(mt, hl.MatrixTable) else mt.annotate(**{vep_root: vep_data})


grch37 = hl.get_reference('GRCh37')
grch37_fasta = 'gs://hail-common/references/human_g1k_v37.fasta.gz'
grch37_fai = 'gs://hail-common/references/human_g1k_v37.fasta.fai'
grch37.add_sequence(grch37_fasta, grch37_fai)
vep_config = "gs://gnomad-resources/loftee-beta/vep85-loftee-gcloud.json"  # this is the config that actually works!

import time as tm

for chr in range(22, 0, -1):  # start from chr22 to make things easier
    chr = str(chr)
    print("starting chr{0}".format(chr))
    print(tm.ctime())
    t = hl.read_table("{0}/MNV_genome_chr{1}_combined.ht".format(output_path, chr))

    t = t.filter((t.filters.length() == 0) & (t.prev_filters.length() == 0))  # filter pass only

    t = t.filter((t.alleles[0].length() == 1) & (t.alleles[1].length() == 1) \
                 & (t.prev_alleles[0].length() == 1) & (
                 t.prev_alleles[1].length() == 1))  # also filter to snps first, to make things easier

    # and also filter to exome region, beforehand
    # exomereg = hl.import_bed("gs://gnomad-qingbowang/MNV/agg_stats/exome_calling_regions.v1.asbed.bed",
    #                             skip_invalid_intervals=True)
    # t = t.filter(hl.is_defined(exomereg[t.locus]))
    # no, it can be filtered downstream anyways

    print("starting vep")
    print(tm.ctime())
    vepped_d1 = annotate_vep_mnv(t, dist=1)

    print("starting downstream filtering")
    print(tm.ctime())

    # keep only the essential columns
    vepped_d1_essense = vepped_d1.key_by("locus", "refs", "alts")  ###このkey byで変に結合している可能性なくはない..
    vepped_d1_essense = vepped_d1_essense.select("AC", "prev_AC", "prev_AF", "n_hethet", "n_hethet2", "n_hethom",
                                                 "n_homhom", "snp1_vep", "snp2_vep", "mnv_vep")

    # filter to canonicals:
    canon_cons_d1 = filter_vep_to_canonical_transcripts(
        filter_vep_to_canonical_transcripts(filter_vep_to_canonical_transcripts(vepped_d1_essense, vep_root="snp1_vep"),
                                            vep_root="snp2_vep"), vep_root="mnv_vep")

    ###explode by consequence
    canon_cons_d1 = canon_cons_d1.annotate(indices=hl.range(0, hl.min(
        hl.len(canon_cons_d1.snp1_vep.transcript_consequences),
        hl.len(canon_cons_d1.snp2_vep.transcript_consequences))))
    canon_cons_d1 = canon_cons_d1.explode('indices')
    canon_cons_d1 = canon_cons_d1.transmute(
        snp1_transcript_consequences=canon_cons_d1.snp1_vep.transcript_consequences[canon_cons_d1.indices],
        snp2_transcript_consequences=canon_cons_d1.snp2_vep.transcript_consequences[canon_cons_d1.indices],
        mnv_transcript_consequences=canon_cons_d1.mnv_vep.transcript_consequences[canon_cons_d1.indices],
        )
    # keep necessary columns
    canon_cons_d1 = canon_cons_d1.select("AC", "prev_AC", "prev_AF", "n_hethet", "n_hethet2", "n_hethom", "n_homhom",
                                         "snp1_transcript_consequences",
                                         "snp2_transcript_consequences",
                                         "mnv_transcript_consequences")

    # further annotate necessary ones
    # codons, consequence_terms, lof, lof_flags, transcript_id
    canon_cons_d1 = canon_cons_d1.annotate(snp1_cons_term=canon_cons_d1.snp1_transcript_consequences.consequence_terms,
                                           snp2_cons_term=canon_cons_d1.snp2_transcript_consequences.consequence_terms,
                                           mnv_cons_term=canon_cons_d1.mnv_transcript_consequences.consequence_terms,
                                           snp1_codons=canon_cons_d1.snp1_transcript_consequences.codons,
                                           snp2_codons=canon_cons_d1.snp2_transcript_consequences.codons,
                                           mnv_codons=canon_cons_d1.mnv_transcript_consequences.codons,
                                           snp1_amino_acids=canon_cons_d1.snp1_transcript_consequences.amino_acids,
                                           snp2_amino_acids=canon_cons_d1.snp2_transcript_consequences.amino_acids,
                                           mnv_amino_acids=canon_cons_d1.mnv_transcript_consequences.amino_acids,
                                           snp1_lof=canon_cons_d1.snp1_transcript_consequences.lof,
                                           snp2_lof=canon_cons_d1.snp2_transcript_consequences.lof,
                                           mnv_lof=canon_cons_d1.mnv_transcript_consequences.lof,
                                           transcript_id=canon_cons_d1.snp1_transcript_consequences.transcript_id
                                           )

    # filter to those that the codon are changed within a single reading frame
    canon_cons_d1 = canon_cons_d1.filter(
        (canon_cons_d1.snp1_codons.length() == 7) & (canon_cons_d1.snp2_codons.length() == 7) & (
        canon_cons_d1.mnv_codons.length() == 7))

    # (subtle but) annotate the AC total, and AF total -- maybe don't do the AF total (since it seems unstable)
    canon_cons_d1 = canon_cons_d1.annotate(
        AC_mnv=canon_cons_d1.n_hethet + canon_cons_d1.n_hethet2 + canon_cons_d1.n_hethom + canon_cons_d1.n_homhom * 2)
    # canon_cons_d1 = canon_cons_d1.annotate(AF_mnv=canon_cons_d1.AC_mnv / (canon_cons_d1.prev_AC / canon_cons_d1.prev_AF)) #AFいいや
    print("starting turing to pandas")
    print(tm.ctime())
    # and turn to pd
    canon_cons_pd1 = canon_cons_d1.select("snp1_cons_term", "snp2_cons_term", "mnv_cons_term", "snp1_codons",
                                          "snp2_codons", "mnv_codons", "snp1_amino_acids", "snp2_amino_acids",
                                          "mnv_amino_acids", "snp1_lof", "snp2_lof", "mnv_lof", "transcript_id", "AC",
                                          "prev_AC", "AC_mnv", "n_homhom").to_pandas()
    # get the most severe
    canon_cons_pd1["snp1_sev"] = canon_cons_pd1.snp1_cons_term.apply(lambda x: cons_term_most_severe(x))
    canon_cons_pd1["snp2_sev"] = canon_cons_pd1.snp2_cons_term.apply(lambda x: cons_term_most_severe(x))
    canon_cons_pd1["mnv_sev"] = canon_cons_pd1.mnv_cons_term.apply(lambda x: cons_term_most_severe(x))
    # annotate the categ
    canon_cons_pd1["categ"] = canon_cons_pd1.apply(
        lambda x: mnv_category(x["snp1_sev"], x["snp2_sev"], x["mnv_sev"], x["snp1_amino_acids"], x["snp2_amino_acids"],
                               x["mnv_amino_acids"]), axis=1)
    # just in case if lof column is all None (NA) and that causes error
    canon_cons_pd1.snp1_lof = canon_cons_pd1.snp1_lof.astype(str)
    canon_cons_pd1.snp2_lof = canon_cons_pd1.snp2_lof.astype(str)
    canon_cons_pd1.mnv_lof = canon_cons_pd1.mnv_lof.astype(str)
    print("starting turing back to hail table and write")
    print(tm.ctime())
    # turn back to hail table, and write as hail table
    hl.Table.from_pandas(canon_cons_pd1).write("{0}/v2_consequence_genome_chr{1}_d1.ht".format(output_path, chr))
    hl.Table.from_pandas(canon_cons_pd1).export("{0}/v2_consequence_genome_chr{1}_d1.tsv".format(output_path, chr))

    # and d2
    t = hl.read_table("{0}/MNV_genome_chr{1}_combined.ht".format(output_path, chr))
    t = t.filter((t.filters.length() == 0) & (t.prev_filters.length() == 0))  # filter pass only
    t = t.filter((t.alleles[0].length() == 1) & (t.alleles[1].length() == 1) \
                 & (t.prev_alleles[0].length() == 1) & (
                     t.prev_alleles[1].length() == 1))  # also filter to snps first, to make things easier
    # and also filter to exome region, beforehand
    # exomereg = hl.import_bed("gs://gnomad-qingbowang/MNV/agg_stats/exome_calling_regions.v1.asbed.bed",
    #                         skip_invalid_intervals=True)
    # t = t.filter(hl.is_defined(exomereg[t.locus]))
    print("starting vep")
    print(tm.ctime())
    vepped_d2 = annotate_vep_mnv(t, dist=2)
    # write vep so that we can start from the downstream later -- let's not do this, for some unknown error
    # print ("saving vep")
    # print (tm.ctime())
    # vepped[0].write("{0}/MNV_full_combined_vepped_d2.ht".format(output_path))
    print("starting downstream filtering")
    print(tm.ctime())
    # keep only the essential columns
    vepped_d2_essense = vepped_d2.key_by("locus", "refs", "alts")
    vepped_d2_essense = vepped_d2_essense.select("AC", "prev_AC", "prev_AF", "n_hethet", "n_hethet2", "n_hethom",
                                                 "n_homhom", "snp1_vep",
                                                 "snp2_vep", "mnv_vep")

    # filter to canonicals:
    canon_cons_d2 = filter_vep_to_canonical_transcripts(
        filter_vep_to_canonical_transcripts(filter_vep_to_canonical_transcripts(vepped_d2_essense, vep_root="snp1_vep"),
                                            vep_root="snp2_vep"), vep_root="mnv_vep")

    ###explode by consequence
    canon_cons_d2 = canon_cons_d2.annotate(indices=hl.range(0, hl.min(
        hl.len(canon_cons_d2.snp1_vep.transcript_consequences),
        hl.len(canon_cons_d2.snp2_vep.transcript_consequences))))
    canon_cons_d2 = canon_cons_d2.explode('indices')
    canon_cons_d2 = canon_cons_d2.transmute(
        snp1_transcript_consequences=canon_cons_d2.snp1_vep.transcript_consequences[canon_cons_d2.indices],
        snp2_transcript_consequences=canon_cons_d2.snp2_vep.transcript_consequences[canon_cons_d2.indices],
        mnv_transcript_consequences=canon_cons_d2.mnv_vep.transcript_consequences[canon_cons_d2.indices],
        )
    # keep necessary columns
    canon_cons_d2 = canon_cons_d2.select("AC", "prev_AC", "prev_AF", "n_hethet", "n_hethet2", "n_hethom", "n_homhom",
                                         "snp1_transcript_consequences",
                                         "snp2_transcript_consequences",
                                         "mnv_transcript_consequences")

    # further annotate necessary ones
    # codons, consequence_terms, lof, lof_flags, transcript_id
    canon_cons_d2 = canon_cons_d2.annotate(snp1_cons_term=canon_cons_d2.snp1_transcript_consequences.consequence_terms,
                                           snp2_cons_term=canon_cons_d2.snp2_transcript_consequences.consequence_terms,
                                           mnv_cons_term=canon_cons_d2.mnv_transcript_consequences.consequence_terms,
                                           snp1_codons=canon_cons_d2.snp1_transcript_consequences.codons,
                                           snp2_codons=canon_cons_d2.snp2_transcript_consequences.codons,
                                           mnv_codons=canon_cons_d2.mnv_transcript_consequences.codons,
                                           snp1_amino_acids=canon_cons_d2.snp1_transcript_consequences.amino_acids,
                                           snp2_amino_acids=canon_cons_d2.snp2_transcript_consequences.amino_acids,
                                           mnv_amino_acids=canon_cons_d2.mnv_transcript_consequences.amino_acids,
                                           snp1_lof=canon_cons_d2.snp1_transcript_consequences.lof,
                                           snp2_lof=canon_cons_d2.snp2_transcript_consequences.lof,
                                           mnv_lof=canon_cons_d2.mnv_transcript_consequences.lof,
                                           transcript_id=canon_cons_d2.snp1_transcript_consequences.transcript_id
                                           )

    # filter to those that the codon are changed within a single reading frame
    canon_cons_d2 = canon_cons_d2.filter(
        (canon_cons_d2.snp1_codons.length() == 7) & (canon_cons_d2.snp2_codons.length() == 7) & (
            canon_cons_d2.mnv_codons.length() == 7))

    # (subtle but) annotate the AC total, and AF total -- maybe don't do the AF total (since it seems unstable)
    canon_cons_d2 = canon_cons_d2.annotate(
        AC_mnv=canon_cons_d2.n_hethet + canon_cons_d2.n_hethet2 + canon_cons_d2.n_hethom + canon_cons_d2.n_homhom * 2)
    # canon_cons_d2 = canon_cons_d2.annotate(AF_mnv=canon_cons_d2.AC_mnv / (canon_cons_d2.prev_AC / canon_cons_d2.prev_AF))
    print("starting turing to pandas")
    print(tm.ctime())
    # and turn to pd --out of memory here... -> create new cluster!
    canon_cons_pd2 = canon_cons_d2.select("snp1_cons_term", "snp2_cons_term", "mnv_cons_term", "snp1_codons",
                                          "snp2_codons",
                                          "mnv_codons", "snp1_amino_acids", "snp2_amino_acids", "mnv_amino_acids",
                                          "snp1_lof", "snp2_lof", "mnv_lof", "transcript_id", "AC", "prev_AC", "AC_mnv",
                                          "n_homhom").to_pandas()
    # get the most severe
    canon_cons_pd2["snp1_sev"] = canon_cons_pd2.snp1_cons_term.apply(lambda x: cons_term_most_severe(x))
    canon_cons_pd2["snp2_sev"] = canon_cons_pd2.snp2_cons_term.apply(lambda x: cons_term_most_severe(x))
    canon_cons_pd2["mnv_sev"] = canon_cons_pd2.mnv_cons_term.apply(lambda x: cons_term_most_severe(x))
    # annotate the categ
    canon_cons_pd2["categ"] = canon_cons_pd2.apply(
        lambda x: mnv_category(x["snp1_sev"], x["snp2_sev"], x["mnv_sev"], x["snp1_amino_acids"], x["snp2_amino_acids"],
                               x["mnv_amino_acids"]), axis=1)
    # just in case if lof column is all None (NA) and that causes error
    canon_cons_pd2.snp1_lof = canon_cons_pd2.snp1_lof.astype(str)
    canon_cons_pd2.snp2_lof = canon_cons_pd2.snp2_lof.astype(str)
    canon_cons_pd2.mnv_lof = canon_cons_pd2.mnv_lof.astype(str)

    print("starting turing back to hail table and write")
    print(tm.ctime())
    # turn back to hail table, and write as hail table
    hl.Table.from_pandas(canon_cons_pd2).write("{0}/v2_consequence_genome_chr{1}_d2.ht".format(output_path, chr))
    hl.Table.from_pandas(canon_cons_pd2).export("{0}/v2_consequence_genome_chr{1}_d2.tsv".format(output_path, chr))
    del canon_cons_pd2  # to free the memory, hopefully...



