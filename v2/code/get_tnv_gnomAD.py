# -*- coding: utf-8 -*-
__author__ = 'QingboWang'

import hail as hl
import hail.expr.aggregators as agg
from typing import *


#first filter the et to relevant ones (bed file specifying the regions where MNV exists at all)
output_path = "gs://gnomad-qingbowang/MNV/1206_genome"
tnvbed = hl.import_bed("gs://gnomad-qingbowang/MNV/MNV_coding_exomes_pos_possibletnv.bed")

output_path = "gs://gnomad-qingbowang/MNV/1206_exome"
for chr in range(22,0,-1):
    et_het = hl.read_table("{0}/tmp_MNV_exome_chr{1}_et_het.ht".format(output_path, chr))
    et_het = et_het.annotate(snp1=hl.str(et_het.prev_row.locus.contig)+"-"+hl.str(et_het.prev_row.locus.position)+"-"+ et_het.prev_row.alleles[0] + "-" + et_het.prev_row.alleles[1],
                             snp2=hl.str(et_het.locus.contig)+"-"+hl.str(et_het.locus.position)+"-"+ et_het.alleles[0] + "-" + et_het.alleles[1])
    et_het = et_het.filter((et_het.alleles[0].length() == 1) & (et_het.alleles[1].length() == 1) \
                           & (et_het.prev_row.alleles[0].length() == 1) & (et_het.prev_row.alleles[1].length() == 1))
    et_het = et_het.annotate(dist=et_het.locus.position - et_het.prev_row.locus.position)  # annotating the distance
    et_het = et_het.filter(et_het.dist != 0)
    et_het = et_het.key_by("snp1","snp2","s","locus","prev_row","homhom")
    et_het2 = hl.read_table("{0}/tmp_MNV_exome_chr{1}_et_het2.ht".format(output_path, chr))
    et_het2 = et_het2.annotate(snp1=hl.str(et_het2.prev_row.locus.contig)+"-"+hl.str(et_het2.prev_row.locus.position)+"-"+ et_het2.prev_row.alleles[0] + "-" + et_het2.prev_row.alleles[1],
                             snp2=hl.str(et_het2.locus.contig)+"-"+hl.str(et_het2.locus.position)+"-"+ et_het2.alleles[0] + "-" + et_het2.alleles[1])
    et_het2 = et_het2.key_by("snp1","snp2","s","locus","prev_row","homhom")
    et_het2 = et_het2.filter((et_het2.alleles[0].length() == 1) & (et_het2.alleles[1].length() == 1) \
                           & (et_het2.prev_row.alleles[0].length() == 1) & (et_het2.prev_row.alleles[1].length() == 1))
    et_het2 = et_het2.annotate(dist=et_het2.locus.position - et_het2.prev_row.locus.position)  # annotating the distance
    et_het2 = et_het2.filter(et_het2.dist != 0)

    et_partially_hom = hl.read_table("{0}/tmp_MNV_exome_chr{1}_et_partially_hom.ht".format(output_path, chr))
    et_partially_hom = et_partially_hom.annotate(snp1=hl.str(et_partially_hom.prev_row.locus.contig)+"-"+hl.str(et_partially_hom.prev_row.locus.position)+"-"+ et_partially_hom.prev_row.alleles[0] + "-" + et_partially_hom.prev_row.alleles[1],
                             snp2=hl.str(et_partially_hom.locus.contig)+"-"+hl.str(et_partially_hom.locus.position)+"-"+ et_partially_hom.alleles[0] + "-" + et_partially_hom.alleles[1])
    et_partially_hom = et_partially_hom.key_by("snp1","snp2","s","locus","prev_row","homhom")
    et_partially_hom = et_partially_hom.filter((et_partially_hom.alleles[0].length() == 1) & (et_partially_hom.alleles[1].length() == 1) \
                           & (et_partially_hom.prev_row.alleles[0].length() == 1) & (et_partially_hom.prev_row.alleles[1].length() == 1))
    et_partially_hom = et_partially_hom.annotate(dist=et_partially_hom.locus.position - et_partially_hom.prev_row.locus.position)  # annotating the distance
    et_partially_hom = et_partially_hom.filter(et_partially_hom.dist != 0)

    et_hom_hom = hl.read_table("{0}/tmp_MNV_exome_chr{1}_et_hom_hom.ht".format(output_path, chr))
    et_hom_hom = et_hom_hom.annotate(snp1=hl.str(et_hom_hom.prev_row.locus.contig)+"-"+hl.str(et_hom_hom.prev_row.locus.position)+"-"+ et_hom_hom.prev_row.alleles[0] + "-" + et_hom_hom.prev_row.alleles[1],
                             snp2=hl.str(et_hom_hom.locus.contig)+"-"+hl.str(et_hom_hom.locus.position)+"-"+ et_hom_hom.alleles[0] + "-" + et_hom_hom.alleles[1])
    et_hom_hom = et_hom_hom.key_by("snp1","snp2","s","locus","prev_row","homhom")
    et_hom_hom = et_hom_hom.filter((et_hom_hom.alleles[0].length() == 1) & (et_hom_hom.alleles[1].length() == 1) \
                           & (et_hom_hom.prev_row.alleles[0].length() == 1) & (et_hom_hom.prev_row.alleles[1].length() == 1))
    et_hom_hom = et_hom_hom.annotate(dist=et_hom_hom.locus.position - et_hom_hom.prev_row.locus.position)  # annotating the distance
    et_hom_hom = et_hom_hom.filter(et_hom_hom.dist != 0)

    et_union = et_het.join(et_het2, how='outer').join(et_partially_hom, how="outer").join(et_hom_hom, how="outer")
    if chr==22:
        et_all = et_union
    else:
        et_all = et_all.join(et_union, how="outer")
    print ("done chr {0}".format(chr))
et_all = et_all.filter((hl.is_defined(tnvbed[et_all.locus])) & (hl.is_defined(tnvbed[et_all.prev_row.locus])))
et_all = et_all.key_by()
et_all = et_all.select("snp1","snp2","s","homhom")
et_all.write("gs://gnomad-qingbowang/MNV/exome_tnvcandidates.ht", overwrite=True)


output_path = "gs://gnomad-qingbowang/MNV/1206_genome"
for chr in range(22,0,-1):
    et_het = hl.read_table("{0}/tmp_MNV_genome_chr{1}_et_het.ht".format(output_path, chr))
    et_het = et_het.annotate(snp1=hl.str(et_het.prev_row.locus.contig)+"-"+hl.str(et_het.prev_row.locus.position)+"-"+ et_het.prev_row.alleles[0] + "-" + et_het.prev_row.alleles[1],
                             snp2=hl.str(et_het.locus.contig)+"-"+hl.str(et_het.locus.position)+"-"+ et_het.alleles[0] + "-" + et_het.alleles[1])
    et_het = et_het.filter((et_het.alleles[0].length() == 1) & (et_het.alleles[1].length() == 1) \
                           & (et_het.prev_row.alleles[0].length() == 1) & (et_het.prev_row.alleles[1].length() == 1))
    et_het = et_het.annotate(dist=et_het.locus.position - et_het.prev_row.locus.position)  # annotating the distance
    et_het = et_het.filter(et_het.dist != 0)
    et_het = et_het.key_by("snp1","snp2","s","locus","prev_row","homhom")
    et_het2 = hl.read_table("{0}/tmp_MNV_genome_chr{1}_et_het2.ht".format(output_path, chr))
    et_het2 = et_het2.annotate(snp1=hl.str(et_het2.prev_row.locus.contig)+"-"+hl.str(et_het2.prev_row.locus.position)+"-"+ et_het2.prev_row.alleles[0] + "-" + et_het2.prev_row.alleles[1],
                             snp2=hl.str(et_het2.locus.contig)+"-"+hl.str(et_het2.locus.position)+"-"+ et_het2.alleles[0] + "-" + et_het2.alleles[1])
    et_het2 = et_het2.key_by("snp1","snp2","s","locus","prev_row","homhom")
    et_het2 = et_het2.filter((et_het2.alleles[0].length() == 1) & (et_het2.alleles[1].length() == 1) \
                           & (et_het2.prev_row.alleles[0].length() == 1) & (et_het2.prev_row.alleles[1].length() == 1))
    et_het2 = et_het2.annotate(dist=et_het2.locus.position - et_het2.prev_row.locus.position)  # annotating the distance
    et_het2 = et_het2.filter(et_het2.dist != 0)

    et_partially_hom = hl.read_table("{0}/tmp_MNV_genome_chr{1}_et_partially_hom.ht".format(output_path, chr))
    et_partially_hom = et_partially_hom.annotate(snp1=hl.str(et_partially_hom.prev_row.locus.contig)+"-"+hl.str(et_partially_hom.prev_row.locus.position)+"-"+ et_partially_hom.prev_row.alleles[0] + "-" + et_partially_hom.prev_row.alleles[1],
                             snp2=hl.str(et_partially_hom.locus.contig)+"-"+hl.str(et_partially_hom.locus.position)+"-"+ et_partially_hom.alleles[0] + "-" + et_partially_hom.alleles[1])
    et_partially_hom = et_partially_hom.key_by("snp1","snp2","s","locus","prev_row","homhom")
    et_partially_hom = et_partially_hom.filter((et_partially_hom.alleles[0].length() == 1) & (et_partially_hom.alleles[1].length() == 1) \
                           & (et_partially_hom.prev_row.alleles[0].length() == 1) & (et_partially_hom.prev_row.alleles[1].length() == 1))
    et_partially_hom = et_partially_hom.annotate(dist=et_partially_hom.locus.position - et_partially_hom.prev_row.locus.position)  # annotating the distance
    et_partially_hom = et_partially_hom.filter(et_partially_hom.dist != 0)

    et_hom_hom = hl.read_table("{0}/tmp_MNV_genome_chr{1}_et_hom_hom.ht".format(output_path, chr))
    et_hom_hom = et_hom_hom.annotate(snp1=hl.str(et_hom_hom.prev_row.locus.contig)+"-"+hl.str(et_hom_hom.prev_row.locus.position)+"-"+ et_hom_hom.prev_row.alleles[0] + "-" + et_hom_hom.prev_row.alleles[1],
                             snp2=hl.str(et_hom_hom.locus.contig)+"-"+hl.str(et_hom_hom.locus.position)+"-"+ et_hom_hom.alleles[0] + "-" + et_hom_hom.alleles[1])
    et_hom_hom = et_hom_hom.key_by("snp1","snp2","s","locus","prev_row","homhom")
    et_hom_hom = et_hom_hom.filter((et_hom_hom.alleles[0].length() == 1) & (et_hom_hom.alleles[1].length() == 1) \
                           & (et_hom_hom.prev_row.alleles[0].length() == 1) & (et_hom_hom.prev_row.alleles[1].length() == 1))
    et_hom_hom = et_hom_hom.annotate(dist=et_hom_hom.locus.position - et_hom_hom.prev_row.locus.position)  # annotating the distance
    et_hom_hom = et_hom_hom.filter(et_hom_hom.dist != 0)

    et_union = et_het.join(et_het2, how='outer').join(et_partially_hom, how="outer").join(et_hom_hom, how="outer")
    if chr==22:
        et_all = et_union
    else:
        et_all = et_all.join(et_union, how="outer")
    print ("done chr {0}".format(chr))
et_all = et_all.filter((hl.is_defined(tnvbed[et_all.locus])) & (hl.is_defined(tnvbed[et_all.prev_row.locus])))
et_all = et_all.key_by()
et_all = et_all.select("snp1","snp2","s","homhom") #あとはあとでannotateすればよさげ.
et_all.write("gs://gnomad-qingbowang/MNV/genome_tnvcandidates.ht", overwrite=True)


#now filter these tnv candidate further,
tnvcand_ex = hl.read_table("gs://gnomad-qingbowang/MNV/exome_tnvcandidates.ht")
tnvcand_ex = tnvcand_ex.annotate(dist = hl.int(tnvcand_ex.snp2.split("-")[1])- hl.int(tnvcand_ex.snp1.split("-")[1]))
#start from the necessary condition: snp1 is tagged with more than one downstream snp for a individual (let the position be x)
neces = tnvcand_ex.group_by("snp1","s").aggregate(n=agg.count())
neces = neces.filter(neces.n>1)
#and if and only if a snp2, in position x+2, is also tagged with more than one upstream snp, this is MNV
neces2 = tnvcand_ex.group_by("snp2","s").aggregate(n=agg.count())
neces2 = neces2.filter(neces2.n>1)
#we are going to match them by snp2 position (x+2)
neces = neces.annotate(snp2_chrpos = neces.snp1.split("-")[0] + "-" + hl.str(hl.int(neces.snp1.split("-")[1])+2))
neces2 = neces2.annotate(snp2_chrpos = neces2.snp2.split("-")[0] + "-" + neces2.snp2.split("-")[1])
joined = neces.key_by("snp2_chrpos","s").join(neces2.key_by("snp2_chrpos","s"), how="inner")
#and pull back the snp in between (position x+1)
d1 = tnvcand_ex.filter(tnvcand_ex.dist==1)
d1 = d1.key_by("snp1","s")
joined = joined.key_by("snp1","s")
joined = joined.annotate(snp1_2 = d1[joined.key].snp2)
#aggregate to get the count
agg = joined.group_by("snp1", "snp1_2", "snp2").aggregate(n_s = agg.count())
agg.write("gs://gnomad-qingbowang/MNV/exome_tnv_agg_all.ht", overwrite=True)

#a bit less efficient, but let's do this for hom combination as well, same calculation again, to get the hom count
tnvcand_exh = hl.read_table("gs://gnomad-qingbowang/MNV/exome_tnvcandidates.ht")
tnvcand_exh = tnvcand_exh.annotate(dist = hl.int(tnvcand_exh.snp2.split("-")[1])- hl.int(tnvcand_exh.snp1.split("-")[1]))
tnvcand_exh = tnvcand_exh.filter(tnvcand_exh.homhom)
necesh = tnvcand_exh.group_by("snp1","s").aggregate(n=agg.count())
necesh = necesh.filter(necesh.n>1)
#and if and only if a snp2, in position x+2, is also tagged with more than one upstream snp, this is MNV
necesh2 = tnvcand_exh.group_by("snp2","s").aggregate(n=agg.count())
necesh2 = necesh2.filter(necesh2.n>1)
#we are going to match them by snp2 position (x+2)
necesh = necesh.annotate(snp2_chrpos = necesh.snp1.split("-")[0] + "-" + hl.str(hl.int(necesh.snp1.split("-")[1])+2))
necesh2 = necesh2.annotate(snp2_chrpos = necesh2.snp2.split("-")[0] + "-" + necesh2.snp2.split("-")[1])
joined = necesh.key_by("snp2_chrpos","s").join(necesh2.key_by("snp2_chrpos","s"), how="inner")
#and pull back the snp in between (position x+1)
d1 = tnvcand_exh.filter(tnvcand_exh.dist==1)
d1 = d1.key_by("snp1","s")
joined = joined.key_by("snp1","s")
joined = joined.annotate(snp1_2 = d1[joined.key].snp2)
#aggregate to get the count
agghom = joined.group_by("snp1", "snp1_2", "snp2").aggregate(n_s = agg.count())
agghom.write("gs://gnomad-qingbowang/MNV/exome_tnv_agg_hom.ht", overwrite=True)


#do the same thing for genome:
tnvcand_gen = hl.read_table("gs://gnomad-qingbowang/MNV/genome_tnvcandidates.ht")
tnvcand_gen = tnvcand_gen.annotate(dist = hl.int(tnvcand_gen.snp2.split("-")[1])- hl.int(tnvcand_gen.snp1.split("-")[1]))
#start from the necessary condition: snp1 is tagged with more than one downstream snp for a individual (let the position be x)
neces = tnvcand_gen.group_by("snp1","s").aggregate(n=agg.count())
neces = neces.filter(neces.n>1)
#and if and only if a snp2, in position x+2, is also tagged with more than one upstream snp, this is MNV
neces2 = tnvcand_gen.group_by("snp2","s").aggregate(n=agg.count())
neces2 = neces2.filter(neces2.n>1)
#we are going to match them by snp2 position (x+2)
neces = neces.annotate(snp2_chrpos = neces.snp1.split("-")[0] + "-" + hl.str(hl.int(neces.snp1.split("-")[1])+2))
neces2 = neces2.annotate(snp2_chrpos = neces2.snp2.split("-")[0] + "-" + neces2.snp2.split("-")[1])
joined = neces.key_by("snp2_chrpos","s").join(neces2.key_by("snp2_chrpos","s"), how="inner")
#and pull back the snp in between (position x+1)
d1 = tnvcand_gen.filter(tnvcand_gen.dist==1)
d1 = d1.key_by("snp1","s")
joined = joined.key_by("snp1","s")
joined = joined.annotate(snp1_2 = d1[joined.key].snp2)
#aggregate to get the count
agg = joined.group_by("snp1", "snp1_2", "snp2").aggregate(n_s = agg.count())
agg.write("gs://gnomad-qingbowang/MNV/genome_tnv_agg_all.ht", overwrite=True)

#a bit less efficient, but let's do this for hom combination as well, same calculation again, to get the hom count
tnvcand_genh = hl.read_table("gs://gnomad-qingbowang/MNV/genome_tnvcandidates.ht")
tnvcand_genh = tnvcand_genh.annotate(dist = hl.int(tnvcand_genh.snp2.split("-")[1])- hl.int(tnvcand_genh.snp1.split("-")[1]))
tnvcand_genh = tnvcand_genh.filter(tnvcand_genh.homhom)
necesh = tnvcand_genh.group_by("snp1","s").aggregate(n=agg.count())
necesh = necesh.filter(necesh.n>1)
#and if and only if a snp2, in position x+2, is also tagged with more than one upstream snp, this is MNV
necesh2 = tnvcand_genh.group_by("snp2","s").aggregate(n=agg.count())
necesh2 = necesh2.filter(necesh2.n>1)
#we are going to match them by snp2 position (x+2)
necesh = necesh.annotate(snp2_chrpos = necesh.snp1.split("-")[0] + "-" + hl.str(hl.int(necesh.snp1.split("-")[1])+2))
necesh2 = necesh2.annotate(snp2_chrpos = necesh2.snp2.split("-")[0] + "-" + necesh2.snp2.split("-")[1])
joined = necesh.key_by("snp2_chrpos","s").join(necesh2.key_by("snp2_chrpos","s"), how="inner")
#and pull back the snp in between (position x+1)
d1 = tnvcand_genh.filter(tnvcand_genh.dist==1)
d1 = d1.key_by("snp1","s")
joined = joined.key_by("snp1","s")
joined = joined.annotate(snp1_2 = d1[joined.key].snp2)
#aggregate to get the count
agghom = joined.group_by("snp1", "snp1_2", "snp2").aggregate(n_s = agg.count())
agghom.write("gs://gnomad-qingbowang/MNV/genome_tnv_agg_hom.ht", overwrite=True)


#でいろいろannotateして終わらせる.

#First, concat these two (hom)
#gen:
agghom = hl.read_table("gs://gnomad-qingbowang/MNV/genome_tnv_agg_hom.ht")
aggall = hl.read_table("gs://gnomad-qingbowang/MNV/genome_tnv_agg_all.ht")
aggall_j = aggall.key_by("snp1","snp1_2","snp2").join(agghom.key_by("snp1","snp1_2","snp2"), how="left")
aggall_j = aggall_j.transmute(n=aggall_j.n_s,n_hom=hl.or_else(aggall_j.n_s_1,0))
aggall_j = aggall_j.annotate(AC_tnv = aggall_j.n + aggall_j.n_hom)
aggall_j.write("gs://gnomad-qingbowang/MNV/genome_tnv_agg.ht")

#ex:
agghom = hl.read_table("gs://gnomad-qingbowang/MNV/exome_tnv_agg_hom.ht")
aggall = hl.read_table("gs://gnomad-qingbowang/MNV/exome_tnv_agg_all.ht")
aggall_j = aggall.key_by("snp1","snp1_2","snp2").join(agghom.key_by("snp1","snp1_2","snp2"), how="left")
aggall_j = aggall_j.transmute(n=aggall_j.n_s,n_hom=hl.or_else(aggall_j.n_s_1,0))
aggall_j = aggall_j.annotate(AC_tnv = aggall_j.n + aggall_j.n_hom)
aggall_j.write("gs://gnomad-qingbowang/MNV/exomes_tnv_agg.ht")

#combine those two
ex = hl.read_table("gs://gnomad-qingbowang/MNV/exomes_tnv_agg.ht")
gen = hl.read_table("gs://gnomad-qingbowang/MNV/genome_tnv_agg.ht")
exgen = ex.key_by("snp1","snp1_2","snp2").join(gen.key_by("snp1","snp1_2","snp2"), how="outer")
exgen = exgen.transmute(n_tnv_ex=hl.or_else(exgen.n,0), n_tnv_hom_ex = hl.or_else(exgen.n_hom,0), AC_tnv_ex = hl.or_else(exgen.AC_tnv,0),
                        n_tnv_gen=hl.or_else(exgen.n_1,0), n_tnv_hom_gen = hl.or_else(exgen.n_hom_1,0), AC_tnv_gen = hl.or_else(exgen.AC_tnv_1,0))
exgen.write("gs://gnomad-qingbowang/MNV/gnomAD_tNV_before_vep.ht")

#also add the component SNPs info
tnv = hl.read_table("gs://gnomad-qingbowang/MNV/gnomAD_tNV_before_vep.ht")
final = hl.import_table("gs://gnomad-public/release/2.1/mnv/gnomad_mnv_coding.tsv")
final = final.key_by("snp1","snp2")
tnv = tnv.key_by("snp1","snp1_2")
tnv = tnv.annotate(snp12_mnv_categ=final[tnv.key].categ,
                   snp12_mnv_n_indv_ex=hl.int(final[tnv.key].n_indv_ex),
                   snp12_mnv_n_indv_gen=hl.int(final[tnv.key].n_indv_gen),
                   snp1_cons = final[tnv.key].snp1_consequence,
                   snp2_cons = final[tnv.key].snp2_consequence)
tnv = tnv.key_by("snp1","snp2")
tnv = tnv.annotate(snp13_mnv_categ=final[tnv.key].categ,
                   snp13_mnv_n_indv_ex=hl.int(final[tnv.key].n_indv_ex),
                   snp13_mnv_n_indv_gen=hl.int(final[tnv.key].n_indv_gen),
                   snp3_cons = final[tnv.key].snp2_consequence)
tnv = tnv.key_by("snp1_2","snp2")
tnv = tnv.annotate(snp23_mnv_categ=final[tnv.key].categ,
                   snp23_mnv_n_indv_ex=hl.int(final[tnv.key].n_indv_ex),
                   snp23_mnv_n_indv_gen=hl.int(final[tnv.key].n_indv_gen))


#then, vep them

tnv = tnv.annotate(tnv = tnv.snp1.split("-")[0]+"-"+tnv.snp1.split("-")[1]+"-"
                   +tnv.snp1.split("-")[2]+tnv.snp1_2.split("-")[2]+tnv.snp2.split("-")[2]+"-"
                   +tnv.snp1.split("-")[3]+tnv.snp1_2.split("-")[3]+tnv.snp2.split("-")[3])
tnv = tnv.annotate(alleles = [tnv.tnv.split("-")[2], tnv.tnv.split("-")[3]],
                  locus = hl.parse_locus(tnv.tnv.split("-")[0]+":"+tnv.tnv.split("-")[1]))
vep_config = 'gs://hail-common/vep/vep/vep85-loftee-gcloud.json'
tnv = tnv.key_by("locus", "alleles")
tnv = hl.vep(tnv, vep_config)

tnv = tnv.explode(tnv.vep.transcript_consequences)
tnv = tnv.filter(tnv.vep.transcript_consequences.canonical == 1) #per transcript, and only canonical ones
tnv = tnv.filter(tnv.vep.transcript_consequences.codons.length() == 7)#also only the one that exactly changes the codon consequence
tnv = tnv.annotate(tnv_cons=tnv.vep.most_severe_consequence, transcript_id = tnv.vep.transcript_consequences.transcript_id)
tnv = tnv.transmute(snp3 = tnv.snp2)
tnv = tnv.transmute(snp2 = tnv.snp1_2)
tnv.select("snp1","snp1_2","snp2","tnv","n_tnv_ex", "n_tnv_hom_ex", "AC_tnv_ex",
                "n_tnv_gen", "n_tnv_hom_gen", "AC_tnv_gen",
                "snp12_mnv_categ","snp13_mnv_categ","snp23_mnv_categ",
          "snp1_cons","snp2_cons","snp3_cons","snp12_mnv_n_indv_ex","snp23_mnv_n_indv_ex","snp13_mnv_n_indv_ex",
                "snp12_mnv_n_indv_gen","snp23_mnv_n_indv_gen","snp13_mnv_n_indv_gen",
          "transcript_id","tnv_cons").export("gs://gnomad-qingbowang/MNV/gnomad_mnv_coding_3bp.tsv")


