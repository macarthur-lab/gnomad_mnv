# -*- coding: utf-8 -*-
__author__ = 'QingboWang'

#from functions import *

# -*- coding: utf-8 -*-
__author__ = 'QingboWang'


output_path = "gs://gnomad-qingbowang/MNV/wholegenome"
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import time as tm
from scipy import stats
import hail as hl
import hail.expr.aggregators as agg
from typing import *



def ht_cnt_mat_to_pd(ht_cnt_mat):
    pds = ht_cnt_mat.to_pandas()
    pds.fillna(value=0, inplace=True)
    pds = pds.iloc[:, 2:]
    pds = pd.pivot_table(pds, index="refs", columns="alts")
    pds.fillna(value=0, inplace=True)
    pds.columns = pds.columns.get_level_values(1)  # to set the columns name properly.
    return (pds) #returning the pandas table


def hl_strc_to_pd_df(strc,lname):
    #from hail struct, create a dataframe of single line
    dict = {}
    for k in strc.keys():
        dict[k] = strc[k]
    return (pd.DataFrame(dict, index=[lname]))



hl.init(tmp_dir="gs://gnomad-qingbowang/tmp")

grch37 = hl.get_reference('GRCh37')
grch37_fasta = 'gs://hail-common/references/human_g1k_v37.fasta.gz'
grch37_fai = 'gs://hail-common/references/human_g1k_v37.fasta.fai'
grch37.add_sequence(grch37_fasta, grch37_fai)
import time as tm
for chr in range(1,23):
    chr = str(chr)
    print ("starting chr{0}".format(chr))
    print (tm.ctime())
    #get MNV
    mnv0 = hl.read_table("{0}/MNV_chr{1}_combined.ht".format(output_path, chr))
    #for all the d, just get the context. downstream could be done local.
    #filtering beforehand
    mnv0 = mnv0.filter((mnv0.filters.length()==0) & (mnv0.prev_filters.length()==0))
    mnv0 = mnv0.filter((mnv0.alleles[0].length()==1) &
                       (mnv0.alleles[1].length()==1) &
                       (mnv0.prev_alleles[0].length()==1) &
                       (mnv0.prev_alleles[1].length()==1))
    mnv0 = mnv0.annotate(refs = mnv0.prev_alleles[0] + mnv0.alleles[0], alts = mnv0.prev_alleles[1] + mnv0.alleles[1])
    mnv0 = mnv0.transmute(ac1 = mnv0.AC, ac2 = mnv0.prev_AC, ac_mnv = mnv0.n_hethet + mnv0.n_hethom + mnv0.n_homhom * 2)
    mnv0 = mnv0.repartition(40)  # less partition better?


    #also annotate AC_mnv, ac1, ac2

    for d in range(1,11):
            print("chr{0} d{1} export start".format(chr, d))
            print (tm.ctime())
            mnv = mnv0.filter((mnv0.locus.position - mnv0.prev_locus.position)==d)
            #add context
            mnv = mnv.annotate(context_ref = mnv.locus.sequence_context(before=4+d, after=4))
            #actually no need to to_pd. Can just export this as tsv. That gives us more information as well.
            #mnv.export("{0}/MNV_chr{1}_d{2}_context.tsv".format(output_path, chr, d))
            #make it lighter
            #no this repartition kills everything.. and it freezes.. ... どうしようin order to let this work...
            mnv = mnv.key_by("locus","refs","alts") #throwing away unneeded
            mnv = mnv.select("context_ref", "ac1","ac2","ac_mnv") #throwing away unneeded
            #mnv = mnv.repartition(160)  # less partition better?
            #mnv.export("{0}/MNV_chr{1}_d{2}_context.tsv".format(output_path, chr, d), parallel="header_per_shard")
            mnv = mnv.repartition(40) #less partition, as always.
            mnv.export("{0}/MNV_chr{1}_d{2}_context.tsv".format(output_path, chr, d))
            #parallel solves what we want?? hopefully... -> not really
            print ("chr{0} d{1} export done".format(chr, d))
            print (tm.ctime())



