# -*- coding: utf-8 -*-
__author__ = 'QingboWang'
import pandas as pd
#initialize hail
import hail as hl
import hail.expr.aggregators as agg
hl.init()


output_path = "gs://gnomad-qingbowang/MNV/wholegenome"
idx_met = ['Coding_UCSC', 'DHS_Trynka', 'Enhancer_Hoffman', 'H3K4me1_Trynka',
       'H3K4me3_Trynka', 'H3K9ac_Trynka', 'Intron_UCSC', 'TSS_Hoffman',
       'Promoter_UCSC', 'Transcribed_Hoffman', 'UTR_3_UCSC', 'UTR_5_UCSC',
       'TFBS_ENCODE', 'H3K27ac_PGC2', 'all']
idx_refs = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT',
       'TA', 'TC', 'TG', 'TT']

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
                flt = flt.drop(ix_revcomp)
    flt = flt[(flt.refs.str[0]!=flt.alts.str[0])&(flt.refs.str[1]!=flt.alts.str[1])]  # deleting trivial ones
    flt.reset_index(inplace=True)
    del flt["index"]
    return (flt)

#also collapse the reference counts
def collapse_ref_cnts(ref_cnt_matrix):
    bp = ref_cnt_matrix.index
    for i in bp: #for all the 2bp
        if i in ref_cnt_matrix.index: #if hasn't been deleted yet
            i_revcomp = revcomp(i)
            if i != i_revcomp:  # if revcomp is not yourself
                ref_cnt_matrix.loc[i, :] = ref_cnt_matrix.loc[i, :] + ref_cnt_matrix.loc[i_revcomp, :]
                ref_cnt_matrix = ref_cnt_matrix.drop(i_revcomp)
    return (ref_cnt_matrix)

#partition by functional annotation, and look at MNV density etc by
#functional annotation, and also per MNV pattern.

#the category that we are interested in
categ = ["Coding_UCSC", "DHS_Trynka", "Enhancer_Hoffman", "H3K27ac_PGC2", "H3K4me1_Trynka", "H3K4me3_Trynka", "H3K9ac_Trynka", "Intron_UCSC", "TSS_Hoffman",
        "Promoter_UCSC", "Transcribed_Hoffman", "UTR_3_UCSC", "UTR_5_UCSC", "TFBS_ENCODE"]


#methylation level per annotation -- also need to re-do this once we have extra functional category..
m = "gs://gnomad-qingbowang/MNV/agg_stats/met_score_final.tsv"
with hl.hadoop_open(m, 'r') as f:
        met = pd.read_table(f)


for d in range(1,11):
    df = pd.DataFrame()
    m = "gs://gnomad-qingbowang/MNV/agg_stats/background_refcnts_d{0}.tsv".format(d)
    with hl.hadoop_open(m, 'r') as f:
        refcnts = pd.read_table(f) #ref cnts
    refcnts.index = idx_refs
    refcnts = collapse_ref_cnts(refcnts)
    for c in categ:
            m = "gs://gnomad-qingbowang/MNV/wholegenome/cnt_mat_d{0}_{1}.tsv".format(d, c)
            with hl.hadoop_open(m, 'r') as f:
                casecnts = pd.read_table(f) #case cnts
            casecnts.index = casecnts.columns 
            casecnts.index = idx_refs
            casecnts.columns = idx_refs #remove the Ns for now.
            casecnts = collapse_crstb_to_revcomp(casecnts)
            ref_of_case = refcnts.loc[:,c]
            casecnts[c] = casecnts.apply(lambda x: x["cnt"] / ref_of_case[x["refs"]], axis=1)
            df = pd.concat([df,casecnts[c]], axis=1, sort=True)
            print ("{0} done".format(c))
    ix = casecnts.refs.str[0] + "N"*(d-1) + casecnts.refs.str[1] + "->" + casecnts.alts.str[0] + "N"*(d-1) + casecnts.alts.str[1]
    df.index = ix
    print ("d={0} done".format(d))
    hl.Table.from_pandas(df).export("{0}/mnv_density_percateg_d{0}.tsv".format(output_path,d))
