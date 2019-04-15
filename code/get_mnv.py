# -*- coding: utf-8 -*-
__author__ = 'QingboWang'


import hail as hl
import hail.expr.aggregators as agg
from typing import *
import sys

vcf = hl.import_vcf(sys.argv[0], call_fields=["GT"]) #change the call_fields according to the hail documentation
vcf = hl.split_multi_hts(vcf)
vcf.write(sys.argv[0] + ".mt", overwrite=True)
mt = hl.read_matrix_table(sys.argv[0] + ".mt")

#calling
mt = mt.select_cols() #dropping unneeded  columns makes things faster
mt = mt.annotate_rows(AC = mt.info.AC[mt.a_index-1], AF = mt.info.AF[mt.a_index-1]) #for case of multiallelic
mt = mt.select_rows(mt.filters, mt.AC, mt.AF) #or any rows that you want to store for future investigation
mt = mt.filter_entries(hl.is_defined(mt.GT) & mt.GT.is_non_ref()) #interested in non-ref only.
mt = hl.window_by_locus(mt, 2) #partition in window -- we only care within codon reading frame, so the max distance is set to be 2
mt = mt.filter_entries((mt.prev_entries.length() > 0))
mt = mt.filter_entries((hl.is_defined(mt.GT) & (mt.prev_entries.length() > 0))) #throwing away no MNV SNPs
mt = mt.filter_entries(mt.prev_entries.filter(lambda x: x.GT.is_non_ref()).length() > 0) #same

et = mt.key_cols_by().entries() # Matrix with 1000 rows (variant) + 1000 cols (sample)=> 1 million entries
et = et.annotate(indices = hl.range(0, hl.len(et.prev_rows)))
et = et.explode('indices') #for the case where there are more than one prev_row for a variant
et = et.transmute(prev_row = et.prev_rows[et.indices],
                      prev_entry = et.prev_entries[et.indices])#tracking back the one with matched indices
et = et.annotate(dist=et.locus.position - et.prev_row.locus.position) #annotating the distance

#filtering
et = et.filter(et.dist>0) #distance=0 is just multiallelic
et = et.filter( (et.alleles[0].length()==1) & (et.alleles[1].length()==1) \
                 & (et.prev_row.alleles[0].length()==1) & (et.prev_row.alleles[1].length()==1) )#interested in SNP only
et = et.filter((et.filters.length()==0) & (et.prev_row.filters.length()==0)) #if you are interested in filter pass variants only

et = et.annotate(hom = ((et.GT.is_diploid()) & (et.prev_entry.GT.is_diploid()) & et.GT.is_hom_var() & (et.prev_entry.GT.is_hom_var())),
                 hethom = ((et.GT.is_diploid()) & (et.prev_entry.GT.is_diploid()) & (et.GT.is_hom_var() & et.prev_entry.GT.is_het_ref()) | (et.GT.is_het_ref() & et.prev_entry.GT.is_hom_var())),#including hom-het, just not distinguishing them two.
                 hethet = ((et.GT.is_diploid()) & (et.prev_entry.GT.is_diploid()) & (et.GT.phased&et.prev_entry.GT.phased) & (et.GT.is_het_ref()&et.prev_entry.GT.is_het_ref()) & (et.GT==et.prev_entry.GT) ))
#write
et_het = et.filter(et.hethet)
et_het.write(sys.argv[0] + "et_het.ht", overwrite=True) #and we are going to save them
et_hethom = et.filter(et.hethom)
et_hethom.write(sys.argv[0] + "et_hethom.ht", overwrite=True)
et_hom = et.filter(et.hom)
et_hom.write(sys.argv[0] + "et_hom.ht", overwrite=True)

#aggregate
per_variant_het = et_het.group_by('locus', 'alleles', "prev_row").aggregate(n=hl.agg.count())
per_variant_hom = et_hom.group_by('locus', 'alleles', "prev_row").aggregate(n=hl.agg.count())
per_variant_hethom = et_hethom.group_by('locus', 'alleles', "prev_row").aggregate(n=hl.agg.count())

#annotate some column names
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
et_hom = et_hom.key_by("locus", "alleles", "prev_row")
per_variant_hom = per_variant_hom.annotate(dist = et_hom[per_variant_hom.key].dist,
                                                  AF = et_hom[per_variant_hom.key].AF,
                                                  AC = et_hom[per_variant_hom.key].AC,
                                                  filters = et_hom[per_variant_hom.key].filters)
per_variant_hom = per_variant_hom.annotate(prev_locus = per_variant_hom.prev_row.locus,
                                                  prev_alleles = per_variant_hom.prev_row.alleles,
                                                  prev_filters = per_variant_hom.prev_row.filters,
                                                  prev_AC = per_variant_hom.prev_row.AC,
                                                  prev_AF = per_variant_hom.prev_row.AF)
et_hethom = et_hethom.key_by("locus", "alleles", "prev_row")
per_variant_hethom = per_variant_hethom.annotate(dist = et_hethom[per_variant_hethom.key].dist,
                                                  AF = et_hethom[per_variant_hethom.key].AF,
                                                  AC = et_hethom[per_variant_hethom.key].AC,
                                                  filters = et_hethom[per_variant_hethom.key].filters)
per_variant_hethom = per_variant_hethom.annotate(prev_locus = per_variant_hethom.prev_row.locus,
                                                  prev_alleles = per_variant_hethom.prev_row.alleles,
                                                  prev_filters = per_variant_hethom.prev_row.filters,
                                                  prev_AC = per_variant_hethom.prev_row.AC,
                                                  prev_AF = per_variant_hethom.prev_row.AF)
#combine het, hom, hethom
comb = per_variant_het.key_by("locus", "alleles","prev_locus","prev_alleles","dist","AF","AC","filters","prev_AF","prev_AC","prev_filters") \
            .join(per_variant_hom.key_by("locus", "alleles", "prev_locus", "prev_alleles", "dist", "AF", "AC", "filters", "prev_AF","prev_AC", "prev_filters"), how='outer') \
            .join(per_variant_hethom.key_by("locus", "alleles","prev_locus","prev_alleles","dist","AF","AC","filters","prev_AF","prev_AC","prev_filters"), how='outer')
comb = comb.transmute(n_hethet=hl.or_else(comb.n, 0), n_hethom=hl.or_else(comb.n_1, 0), n_homhom=hl.or_else(comb.n_2, 0))
comb = comb.select("n_hethet","n_hethom","n_homhom")
comb = comb.annotate(n_total = comb.n_hethet + comb.n_hethom + comb.n_homhom)
#write
comb.write(sys.argv[0] + "mnv_combined.ht", overwrite=True)






