# -*- coding: utf-8 -*-
__author__ = 'QingboWang'

import hail as hl
import hail.expr.aggregators as agg
from typing import *


output_path = "gs://gnomad-qingbowang/MNV/1206_exome"
for chr in range(22,0,-1):
    et_het = hl.read_table("{0}/tmp_MNV_exome_chr{1}_et_het.ht".format(output_path, chr))
    et_het = et_het.annotate(snp1=hl.str(et_het.prev_row.locus.contig)+"-"+hl.str(et_het.prev_row.locus.position)+"-"+ et_het.prev_row.alleles[0] + "-" + et_het.prev_row.alleles[1],
                             snp2=hl.str(et_het.locus.contig)+"-"+hl.str(et_het.locus.position)+"-"+ et_het.alleles[0] + "-" + et_het.alleles[1])
    et_het = et_het.key_by("snp1","snp2")
    et_het2 = hl.read_table("{0}/tmp_MNV_exome_chr{1}_et_het2.ht".format(output_path, chr))
    et_het2 = et_het2.annotate(snp1=hl.str(et_het2.prev_row.locus.contig)+"-"+hl.str(et_het2.prev_row.locus.position)+"-"+ et_het2.prev_row.alleles[0] + "-" + et_het2.prev_row.alleles[1],
                             snp2=hl.str(et_het2.locus.contig)+"-"+hl.str(et_het2.locus.position)+"-"+ et_het2.alleles[0] + "-" + et_het2.alleles[1])
    et_het2 = et_het2.key_by("snp1","snp2")
    et_partially_hom = hl.read_table("{0}/tmp_MNV_exome_chr{1}_et_partially_hom.ht".format(output_path, chr))
    et_partially_hom = et_partially_hom.annotate(snp1=hl.str(et_partially_hom.prev_row.locus.contig)+"-"+hl.str(et_partially_hom.prev_row.locus.position)+"-"+ et_partially_hom.prev_row.alleles[0] + "-" + et_partially_hom.prev_row.alleles[1],
                             snp2=hl.str(et_partially_hom.locus.contig)+"-"+hl.str(et_partially_hom.locus.position)+"-"+ et_partially_hom.alleles[0] + "-" + et_partially_hom.alleles[1])
    et_partially_hom = et_partially_hom.key_by("snp1","snp2")
    et_hom_hom = hl.read_table("{0}/tmp_MNV_exome_chr{1}_et_hom_hom.ht".format(output_path, chr))
    et_hom_hom = et_hom_hom.annotate(snp1=hl.str(et_hom_hom.prev_row.locus.contig)+"-"+hl.str(et_hom_hom.prev_row.locus.position)+"-"+ et_hom_hom.prev_row.alleles[0] + "-" + et_hom_hom.prev_row.alleles[1],
                             snp2=hl.str(et_hom_hom.locus.contig)+"-"+hl.str(et_hom_hom.locus.position)+"-"+ et_hom_hom.alleles[0] + "-" + et_hom_hom.alleles[1])
    et_hom_hom = et_hom_hom.key_by("snp1","snp2")
    et_het = et_het.key_by("snp1","snp2","s")
    et_het2 = et_het2.key_by("snp1","snp2","s")
    et_partially_hom = et_partially_hom.key_by("snp1","snp2","s")
    et_hom_hom = et_hom_hom.key_by("snp1","snp2","s") #s is actually not needed cuz mutually exclusive? no it's needed. 同じmnvをたくさんの人が持ってるから.
    et_union = et_het.join(et_het2, how='outer').join(et_partially_hom, how="outer").join(et_hom_hom, how="outer")
    if chr==22:
        et_all = et_union
    else:
        et_all = et_all.join(et_union, how="outer")
    print ("done chr {0}".format(chr))
et_all.select("s").write("gs://gnomad-qingbowang/MNV/tmp_exome_snp_pairs.ht", overwrite=True)
final = hl.import_table("gs://gnomad-public/release/2.1/mnv/gnomad_mnv_coding.tsv", types={'n_indv_ex': hl.tint32})
final = final.filter(final.n_indv_ex>0)

#et_all.show()
et_all= hl.read_table("gs://gnomad-qingbowang/MNV/tmp_exome_snp_pairs.ht")
final = final.key_by()
final = final.key_by("snp1","snp2")
et_all = et_all.key_by()
et_all = et_all.key_by("snp1","snp2")
et_all = et_all.annotate(categ = final[et_all.key].categ)

cnt = et_all.group_by("s","categ").aggregate(n = agg.count())
sums = cnt.group_by("categ").aggregate(sums=hl.agg.sum(cnt.n))

sums.show()
#output:
"""
| "Changed missense"           |  5499378 |
| "Gained PTV"                 |     8079 |
| "Gained missense"            |      202 |
| "Lost missense"              |    10434 |
| "Partially changed missense" |   843375 |
| "Rescued PTV"                |   556142 |
| "Rescued stop loss"          |        1 |
| "Unchanged"                  |  4953574 |
| "gained_stop_loss"           |    22762 |
| NA                           | 57642684 |
"""
#print per indv
import numpy as np
n = 125748.
v = np.array([5499378,8079,202,10434, 843375, 556142, 4953574, 22762])#, 57642684])
print (v/n)
print ((sum(v)-4953574)/n) #Unchangedを差し引いた