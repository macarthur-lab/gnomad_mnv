# -*- coding: utf-8 -*-
__author__ = 'QingboWang'

import hail as hl
import hail.expr.aggregators as agg
from typing import *
hl.init()

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


# check PGT==NA --  most of them are het het parents pairs
def pbt_phased_trios_mt_path(data_type: str, split: bool = True, hail_version: str = CURRENT_HAIL_VERSION):
    return "gs://gnomad/hardcalls/hail-{0}/mt/{1}/gnomad.{1}.trios.pbt_phased{2}.mt".format(hail_version, data_type,
                                                                                            "" if split else ".unsplit")
exomes = hl.read_matrix_table(pbt_phased_trios_mt_path("exomes"))
exomes_proband = exomes.filter_cols(exomes.s == exomes.source_trio.proband.s)
exomes_fath = exomes.filter_cols(exomes.s == exomes.source_trio.father.s)
exomes_moth = exomes.filter_cols(exomes.s == exomes.source_trio.mother.s)

exomes_proband = exomes_proband.key_cols_by(exomes_proband['source_trio'].fam_id)
exomes_fath = exomes_fath.key_cols_by(exomes_fath['source_trio'].fam_id)
exomes_moth = exomes_moth.key_cols_by(exomes_moth['source_trio'].fam_id)

exomes_proband = exomes_proband.annotate_entries(mother_PBT_GT=exomes_moth[exomes_proband.row_key, exomes_proband.col_key]["PBT_GT"], father_PBT_GT=exomes_fath[exomes_proband.row_key,exomes_proband.col_key]["PBT_GT"])
exomes_proband = exomes_proband.annotate_entries(hethet =  ((exomes_proband.mother_PBT_GT.is_het_ref()) & (exomes_proband.father_PBT_GT.is_het_ref())))

exomes_proband = hl.filter_alleles(exomes_proband, lambda allele, i: hl.is_snp(exomes_proband.alleles[0], allele))  # currently take only SNP
exomes_proband = exomes_proband.filter_entries(exomes_proband.GT.is_het())  # throw away unwanted entries (non alt)


exomes_proband_et = exomes_proband.entries()
exomes_proband_et = exomes_proband_et.filter(exomes_proband_et.adj)
#exomes_proband_et = exomes_proband_et.filter(exomes_proband_et.GT.is_het_ref() & exomes_proband_et.GT.is_diploid())
aggstats = exomes_proband_et.group_by("GT",'PBT_GT', 'mother_PBT_GT', "father_PBT_GT").aggregate(n=hl.agg.count())
aggstats.write("gs://gnomad_qingbowang/MNV/hethet_aggstats_exome_wGT_re.ht")


















