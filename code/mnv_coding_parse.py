# -*- coding: utf-8 -*-
__author__ = 'QingboWang'
import pandas as pd
import numpy as np

#take exome and genome, combine them, and sort out everything for release
#(might want to do this in cloud?)
import hail as hl
import hail.expr.aggregators as agg
from typing import *


#this code works in local

df1s = []
df2s = []
dfg1s = []
dfg2s = []
for chr in list(np.arange(1,23)) + ["X","Y"]:
    chr = str(chr)
    df1s.append(pd.read_csv("~/Downloads/v2_consequence_exome_chr{0}_d1.tsv".format(chr),sep="\t"))
    df2s.append(pd.read_csv("~/Downloads/v2_consequence_exome_chr{0}_d2.tsv".format(chr),sep="\t"))
    if chr!="Y":
        dfg1s.append(pd.read_csv("~/Downloads/v2_consequence_genome_chr{0}_d1.tsv".format(chr),sep="\t"))
        dfg2s.append(pd.read_csv("~/Downloads/v2_consequence_genome_chr{0}_d2.tsv".format(chr), sep="\t"))
df1s = pd.concat(df1s)
df2s = pd.concat(df2s)
df2s.mnv_codons = df2s.mnv_codons.str[0] + df2s.mnv_codons.str[1].apply(lambda x: x.lower()) + df2s.mnv_codons.str[2:5] + df2s.mnv_codons.str[5].apply(lambda x: x.lower()) + df2s.mnv_codons.str[6:]
dfe = pd.concat([df1s, df2s])
dfg1s = pd.concat(dfg1s)
dfg2s = pd.concat(dfg2s)
dfg2s.mnv_codons = dfg2s.mnv_codons.str[0] + dfg2s.mnv_codons.str[1].apply(lambda x: x.lower()) + dfg2s.mnv_codons.str[2:5] + dfg2s.mnv_codons.str[5].apply(lambda x: x.lower()) + dfg2s.mnv_codons.str[6:]
dfg = pd.concat([dfg1s, dfg2s])

#as a result, much more MNV than previously discovered... Why? because of the edge phasing problem?
#should be fine.

dfe["snp1_consequence"] = dfe.snp1_cons_term.apply(lambda x: cons_term_most_severe(x))
dfe["snp2_consequence"] = dfe.snp2_cons_term.apply(lambda x: cons_term_most_severe(x))
dfe["mnv_consequence"] = dfe.mnv_cons_term.apply(lambda x: cons_term_most_severe(x))
# annotate the categ
dfe["categ_re"] = dfe.apply(
    lambda x: mnv_category(x["snp1_consequence"], x["snp2_consequence"], x["mnv_consequence"], x["snp1_amino_acids"], x["snp2_amino_acids"],
                           x["mnv_amino_acids"]), axis=1)
print (dfe.categ.value_counts())
print (dfe.categ_re.value_counts()) #to check that noncoding_or_else is basically gone
dfe[dfe.categ!=dfe.categ_re][["snp1_amino_acids","snp2_amino_acids","mnv_amino_acids","categ","categ_re","snp1_consequence","snp2_consequence","mnv_consequence"]]
dfe.categ = dfe.categ_re
dfe.drop(["categ_re"], axis=1, inplace=True)


#do it for dfg as well
dfg["snp1_consequence"] = dfg.snp1_cons_term.apply(lambda x: cons_term_most_severe(x))
dfg["snp2_consequence"] = dfg.snp2_cons_term.apply(lambda x: cons_term_most_severe(x))
dfg["mnv_consequence"] = dfg.mnv_cons_term.apply(lambda x: cons_term_most_severe(x))
# annotate the categ
dfg["categ"] = dfg.apply(
    lambda x: mnv_category(x["snp1_consequence"], x["snp2_consequence"], x["mnv_consequence"], x["snp1_amino_acids"], x["snp2_amino_acids"],
                           x["mnv_amino_acids"]), axis=1)




#annotate some additional columns
dfe.drop(["snp1_cons_term","snp2_cons_term","mnv_cons_term"], axis=1, inplace=True)
dfe["n_indv_ex"] = dfe.AC_mnv - dfe.n_homhom
dfe.rename(index=str, columns={"AC":"AC_ex","prev_AC":"prev_AC_ex",'AC_mnv':"AC_mnv_ex", 'n_homhom':"n_homhom_ex"}, inplace=True)

#concat with genome
dfe.index = dfe["locus.contig"].astype(str)+"-"+dfe["locus.position"].astype(str)+"-"+dfe.refs+"-"+dfe.alts + "-"+dfe.transcript_id
dfg.index = dfg["locus.contig"].astype(str)+"-"+dfg["locus.position"].astype(str)+"-"+dfg.refs+"-"+dfg.alts + "-"+dfg.transcript_id
dfg.drop(["snp1_cons_term","snp2_cons_term","mnv_cons_term"], axis=1, inplace=True)
dfg["n_indv_gen"] = dfg.AC_mnv - dfg.n_homhom
dfg.rename(index=str, columns={"AC":"AC_gen","prev_AC":"prev_AC_gen",'AC_mnv':"AC_mnv_gen", 'n_homhom':"n_homhom_gen"}, inplace=True)

t = dfe.join(dfg, how="outer", rsuffix="_gen")
t.rename(index=str, columns={"AC_ex":"AC_snp2_ex","prev_AC_ex":"AC_snp1_ex","AC_gen":"AC_snp2_gen","prev_AC_gen":"AC_snp1_gen"}, inplace=True)


#t.refs.isnull().sum() #Those only in genome, not in exome
#len(t.refs) #all
only_gen = t.refs.isnull()
for c in t.columns:
    if c + "_gen" in t.columns: #if it is common annotation (basically all other than AC_mnv etc)
        print (c)#just for test
        t.loc[only_gen, c] = t.loc[only_gen, c+"_gen"]
        del t[c+"_gen"]#remove the _gen part
t[['AC_mnv_ex','n_homhom_ex','n_indv_ex','AC_mnv_gen','n_homhom_gen','n_indv_gen']] = t[['AC_mnv_ex','n_homhom_ex','n_indv_ex','AC_mnv_gen','n_homhom_gen','n_indv_gen']].fillna(value=0)
for cols in ["locus.position", 'AC_mnv_ex', 'n_homhom_ex', #turn to int
             'n_indv_ex', 'AC_mnv_gen', 'n_homhom_gen','n_indv_gen']: #'AC_snp2_ex', 'AC_snp1_ex','AC_snp2_gen', 'AC_snp1_gen' はNAが入っているので保留
    print (cols)
    t[cols] = t[cols].astype(int)


#add some columns for the browser
t["snp1_ref"] = t.refs.str[0]
t["snp2_ref"] = t.refs.str[-1]
t["snp1_alt"] = t.alts.str[0]
t["snp2_alt"] = t.alts.str[-1]
t["n_indv"] = t.n_indv_ex + t.n_indv_gen
t["AC_mnv"] = t.AC_mnv_ex + t.AC_mnv_gen
t["n_homhom"] = t.n_homhom_ex + t.n_homhom_gen

t["snp1"] = t["locus.contig"].astype(str) +"-"+ t["locus.position"].astype(str) +"-"+ t.snp1_ref +"-"+ t.snp1_alt
t["snp2"] = t["locus.contig"].astype(str) +"-"+ (t["locus.position"] + t.refs.apply(lambda x: len(x)-1)).astype(str) +"-"+ t.snp2_ref +"-"+ t.snp2_alt
t["mnv"] = t["locus.contig"].astype(str) +"-"+ t["locus.position"].astype(str) +"-"+ t.refs +"-"+ t.alts

#btw check the example of >1 transcript on a same pos, refs, alts
t.mnv.value_counts() #mostly not that much
#crazy example:2-234675782-GC-TT       9 canonical transcripts.. but it's real.
t[['locus.position',"locus.contig" ,'transcript_id', 'categ',"snp1","snp2","mnv",'snp1_consequence','snp2_consequence', 'mnv_consequence',
    "snp1_codons","snp2_codons","mnv_codons",
    "snp1_amino_acids","snp2_amino_acids","mnv_amino_acids",
    'snp1_lof', 'snp2_lof', 'mnv_lof',
    "n_indv",'AC_mnv', 'n_homhom', "n_indv_ex",'AC_mnv_ex', 'n_homhom_ex',"n_indv_gen",'AC_mnv_gen', 'n_homhom_gen']].to_csv("~/Downloads/mnv_exome_and_genome_alpha_before_ACannot.tsv",sep="\t", index=False)

#also export the exome only, for figure2 of mnv paper (and also filtering to HC)


#->annotated the ac etc in hail (but with the old and wrong "categ"):
#so this is not really a "mnv_final", but just calling it final.
mnv_final = pd.read_csv("~/Downloads/mnv_exome_and_genome_final.tsv",sep="\t")
#NaN -> 0
mnv_final[["AC_snp1_ex",  "AC_snp1_gen",  "AC_snp2_ex",  "AC_snp2_gen"]] = mnv_final[["AC_snp1_ex",  "AC_snp1_gen",  "AC_snp2_ex",  "AC_snp2_gen"]].fillna(value=0).astype(int)
mnv_final["AC_snp1"] = mnv_final.AC_snp1_ex + mnv_final.AC_snp1_gen
mnv_final["AC_snp2"] = mnv_final.AC_snp2_ex + mnv_final.AC_snp2_gen

#X and Y --looking ok
mnv_final[mnv_final["locus.contig"]=="X"].head() #2l looking good.
mnv_final[mnv_final["locus.contig"]=="Y"].head() #1l looking good.
print(t.shape,mnv_final.shape) #make sure they are the same row num
#also same order
(np.array(t.transcript_id) != np.array(mnv_final.transcript_id).sum())
#oh we need to set the order again.. since t is order such that d=1 comes first -> order by...
#contig, position, mnv, transcript_id
t = t.sort_values(by=["locus.contig","locus.position","mnv","transcript_id"])
mnv_final = mnv_final.sort_values(by=["locus.contig","locus.position","mnv","transcript_id"])
(np.array(t.transcript_id) != np.array(mnv_final.transcript_id)).sum() #ok everything matches
t.drop(["AC_snp1_ex",  "AC_snp1_gen",  "AC_snp2_ex",  "AC_snp2_gen"], axis=1, inplace=True)
mnv_final = pd.concat([t.reset_index(),mnv_final[["AC_snp1_ex",  "AC_snp1_gen",  "AC_snp2_ex",  "AC_snp2_gen","AC_snp1","AC_snp2"]]],axis=1)
#final output, for public release
mnv_final[['locus.position',"locus.contig" ,'transcript_id', 'categ',"snp1","snp2","mnv",'snp1_consequence','snp2_consequence', 'mnv_consequence',
    "snp1_codons","snp2_codons","mnv_codons",
    "snp1_amino_acids","snp2_amino_acids","mnv_amino_acids",
    'snp1_lof', 'snp2_lof', 'mnv_lof',
    "n_indv",'AC_mnv', 'n_homhom', "n_indv_ex",'AC_mnv_ex', 'n_homhom_ex',"n_indv_gen",'AC_mnv_gen', 'n_homhom_gen',
   "AC_snp1","AC_snp1_ex","AC_snp1_gen","AC_snp2","AC_snp2_ex","AC_snp2_gen"]].to_csv("~/Downloads/mnv_coding_final_release.tsv",sep="\t", index=False)
#and restrict to exome/genome for subset release
ex_final = mnv_final[mnv_final.AC_mnv_ex>0]
gen_final = mnv_final[mnv_final.AC_mnv_gen>0]
ex_final[['locus.position',"locus.contig" ,'transcript_id', 'categ',"snp1","snp2","mnv",'snp1_consequence','snp2_consequence', 'mnv_consequence',
    "snp1_codons","snp2_codons","mnv_codons",
    "snp1_amino_acids","snp2_amino_acids","mnv_amino_acids",
    'snp1_lof', 'snp2_lof', 'mnv_lof',
    "n_indv_ex",'AC_mnv_ex', 'n_homhom_ex',"AC_snp1_ex","AC_snp2_ex"]].to_csv("~/Downloads/mnv_coding_final_release_exome.tsv",sep="\t", index=False)
gen_final[['locus.position',"locus.contig" ,'transcript_id', 'categ',"snp1","snp2","mnv",'snp1_consequence','snp2_consequence', 'mnv_consequence',
    "snp1_codons","snp2_codons","mnv_codons",
    "snp1_amino_acids","snp2_amino_acids","mnv_amino_acids",
    'snp1_lof', 'snp2_lof', 'mnv_lof',
    "n_indv_gen",'AC_mnv_gen', 'n_homhom_gen',"AC_snp1_gen","AC_snp2_gen"]].to_csv("~/Downloads/mnv_coding_final_release_genome.tsv",sep="\t", index=False)

#and I am going to use this ex_final for the analysis in the paper.
#do it as final version in fig2.py

#.... forgot to filter by RF.... -> renamed to mnv_coding_final_release_nonRF.tsv etc, annotating and filtering in cloud.
#doing in 1217 note - mnv_coding_annotate_filters.pyでやった
t = pd.read_csv("~/downloads/mnv_exome_and_genome_beta_filter_annotated.tsv",sep="\t")
#genでNaNな部分は何なんだろう?
t.filter_snp1_ex.value_counts(dropna=False)
t.filter_snp2_ex.value_counts(dropna=False)
t.filter_snp1_gen.value_counts(dropna=False)
t.filter_snp2_gen.value_counts(dropna=False)
#no variant case check some. l2, 4, 5 for snp1
#yes it is.
# -> for exome(genome), look at exome(genome) filtering
ex_final = t[(t.AC_mnv_ex>0)&(t.filter_snp1_ex=="[]")&(t.filter_snp2_ex=="[]")]
gen_final = t[(t.AC_mnv_gen>0)&(t.filter_snp1_gen=="[]")&(t.filter_snp2_gen=="[]")]
ex_final[['locus.position',"locus.contig" ,'transcript_id', 'categ',"snp1","snp2","mnv",'snp1_consequence','snp2_consequence', 'mnv_consequence',
    "snp1_codons","snp2_codons","mnv_codons",
    "snp1_amino_acids","snp2_amino_acids","mnv_amino_acids",
    'snp1_lof', 'snp2_lof', 'mnv_lof',
    "n_indv_ex",'AC_mnv_ex', 'n_homhom_ex',"AC_snp1_ex","AC_snp2_ex"]].to_csv("~/Downloads/mnv_coding_final_release_exome.tsv",sep="\t", index=False)
gen_final[['locus.position',"locus.contig" ,'transcript_id', 'categ',"snp1","snp2","mnv",'snp1_consequence','snp2_consequence', 'mnv_consequence',
    "snp1_codons","snp2_codons","mnv_codons",
    "snp1_amino_acids","snp2_amino_acids","mnv_amino_acids",
    'snp1_lof', 'snp2_lof', 'mnv_lof',
    "n_indv_gen",'AC_mnv_gen', 'n_homhom_gen',"AC_snp1_gen","AC_snp2_gen"]].to_csv("~/Downloads/mnv_coding_final_release_genome.tsv",sep="\t", index=False)
#for the conbined data, we actually want to keep the NaNs. = Filter those where at least one RF etc is found.
#so let's just get rid of sites where both genome and exome are RF filtered, but
#keep it if at least one of them survived. eg. http://gnomad.broadinstitute.org/variant/1-138593-G-T
#and for them, we can still put the total ac including the ones filtered out (as in the example above)

#checking those that are filtered out by both
t[(t.filter_snp1_ex=='["RF"]') & (t.filter_snp1_gen=='["RF"]')]
#そのまま残しておけばいいや..
t[['locus.position',"locus.contig" ,'transcript_id', 'categ',"snp1","snp2","mnv",'snp1_consequence','snp2_consequence', 'mnv_consequence',
    "snp1_codons","snp2_codons","mnv_codons",
    "snp1_amino_acids","snp2_amino_acids","mnv_amino_acids",
    'snp1_lof', 'snp2_lof', 'mnv_lof',
    "n_indv",'AC_mnv', 'n_homhom', "n_indv_ex",'AC_mnv_ex', 'n_homhom_ex',"n_indv_gen",'AC_mnv_gen', 'n_homhom_gen',
   "AC_snp1","AC_snp1_ex","AC_snp1_gen","AC_snp2","AC_snp2_ex","AC_snp2_gen",
   "filter_snp1_ex","filter_snp1_gen","filter_snp2_ex","filter_snp2_gen"]].to_csv("~/Downloads/mnv_coding_final_release.tsv",sep="\t", index=False, quotechar="?") #dummy quote char


