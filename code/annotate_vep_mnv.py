# -*- coding: utf-8 -*-
__author__ = 'QingboWang'


import hail as hl
import hail.expr.aggregators as agg
from typing import *
vep_config = "gs://gnomad-resources/loftee-beta/vep85-loftee-gcloud.json"  # or any config that works for users
import sys
import pandas as pd
import numpy as np

def filter_vep_to_canonical_transcripts(mt: Union[hl.MatrixTable, hl.Table],
                                        vep_root: str = 'vep') -> Union[hl.MatrixTable, hl.Table]:
    canonical = mt[vep_root].transcript_consequences.filter(lambda csq: csq.canonical == 1)
    vep_data = mt[vep_root].annotate(transcript_consequences=canonical)
    return mt.annotate_rows(**{vep_root: vep_data}) if isinstance(mt, hl.MatrixTable) else mt.annotate(**{vep_root: vep_data})

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


#read MNV
mnv = hl.read_table(sys.argv[1])

#annotate snv effects
mnv = mnv.key_by("locus","alleles")
mnv = hl.vep(mnv, vep_config, name="snp2_vep")
mnv = mnv.key_by() #unkey first
mnv = mnv.rename({'locus' : 'snp2_locus', 'alleles' : 'snp2_alleles',"prev_locus":"locus", "prev_alleles":"alleles"}) #rename the snp1 locus and alleles as locus, allele
mnv = mnv.key_by('locus', 'alleles') #and re-key
mnv = hl.vep(mnv, vep_config, name="snp1_vep") #and vep

#annotate MNV effects, specified by the distance
if sys.argv[2]==1:
    t = mnv.filter(mnv.dist==1)
    t = t.annotate(refs=t.alleles[0] +t.snp2_alleles[0], alts=t.alleles[1] +t.snp2_alleles[1])# annotate combined refs, and alts
    t = t.annotate(mnv_alleles = [t.refs,t.alts]) #and let it be mnv_alleles
    t = t.key_by()#to re-key
    t = t.rename({'alleles' : 'snp1_alleles',"mnv_alleles":"alleles"}) #now the alleles is that of mnv (note. locus does not change. same as snp1)
    t = t.key_by('locus', 'alleles') #and re-key
    t = hl.vep(t, vep_config, name="mnv_vep")


if sys.argv[2]==2:
    grch37 = hl.get_reference('GRCh37')
    grch37_fasta = 'gs://hail-common/references/human_g1k_v37.fasta.gz'
    grch37_fai = 'gs://hail-common/references/human_g1k_v37.fasta.fai'
    grch37.add_sequence(grch37_fasta, grch37_fai)
    t = mnv.filter(mnv.dist == 2)
    t = t.annotate(refs=t.alleles[0] + t.locus.sequence_context(before=-1, after=1) + t.snp2_alleles[0],
                     alts=t.alleles[1] + t.locus.sequence_context(before=-1, after=1) + t.snp2_alleles[1])
    t = t.annotate(mnv_alleles=[t.refs, t.alts])  # and let it be mnv_alleles
    t = t.key_by()  # to re-key
    t = t.rename({'alleles': 'snp1_alleles',
                    "mnv_alleles": "alleles"})  # now the alleles is that of mnv (note. locus does not change. same as snp1)
    t = t.key_by('locus', 'alleles')  # and re-key
    t = hl.vep(t, vep_config, name="mnv_vep")




vepped_d_essense = t.key_by("locus", "refs", "alts") 
vepped_d_essense = vepped_d_essense.select("AC", "prev_AC", "prev_AF", "n_hethet", "n_hethom",
                                                     "n_homhom", "snp1_vep", "snp2_vep", "mnv_vep")
# filter to canonicals:
canon_cons_d = filter_vep_to_canonical_transcripts(filter_vep_to_canonical_transcripts(
filter_vep_to_canonical_transcripts(vepped_d_essense, vep_root="snp1_vep"), vep_root="snp2_vep"),
                                                            vep_root="mnv_vep")

###explode by consequence, since there might be more than 1 canonical transcript for a MNV
canon_cons_d = canon_cons_d.annotate(indices=hl.range(0, hl.min(
            hl.len(canon_cons_d.snp1_vep.transcript_consequences), hl.len(canon_cons_d.snp2_vep.transcript_consequences))))
canon_cons_d = canon_cons_d.explode('indices')
canon_cons_d = canon_cons_d.transmute(
            snp1_transcript_consequences=canon_cons_d.snp1_vep.transcript_consequences[canon_cons_d.indices],
            snp2_transcript_consequences=canon_cons_d.snp2_vep.transcript_consequences[canon_cons_d.indices],
            mnv_transcript_consequences=canon_cons_d.mnv_vep.transcript_consequences[canon_cons_d.indices],
            )
# keep necessary columns
canon_cons_d = canon_cons_d.select("AC", "prev_AC", "prev_AF", "n_hethet", "n_hethom",
                                             "n_homhom",
                                             "snp1_transcript_consequences",
                                             "snp2_transcript_consequences",
                                             "mnv_transcript_consequences")

# further annotate necessary ones
# codons, consequence_terms, lof, lof_flags, transcript_id
canon_cons_d = canon_cons_d.annotate(
            snp1_cons_term=canon_cons_d.snp1_transcript_consequences.consequence_terms,
            snp2_cons_term=canon_cons_d.snp2_transcript_consequences.consequence_terms,
            mnv_cons_term=canon_cons_d.mnv_transcript_consequences.consequence_terms,
            snp1_codons=canon_cons_d.snp1_transcript_consequences.codons,
            snp2_codons=canon_cons_d.snp2_transcript_consequences.codons,
            mnv_codons=canon_cons_d.mnv_transcript_consequences.codons,
            snp1_amino_acids=canon_cons_d.snp1_transcript_consequences.amino_acids,
            snp2_amino_acids=canon_cons_d.snp2_transcript_consequences.amino_acids,
            mnv_amino_acids=canon_cons_d.mnv_transcript_consequences.amino_acids,
            snp1_lof=canon_cons_d.snp1_transcript_consequences.lof,
            snp2_lof=canon_cons_d.snp2_transcript_consequences.lof,
            mnv_lof=canon_cons_d.mnv_transcript_consequences.lof,
            transcript_id=canon_cons_d.snp1_transcript_consequences.transcript_id
            )

# filter to those that the codon are changed within a single reading frame
canon_cons_d = canon_cons_d.filter(
            (canon_cons_d.snp1_codons.length() == 7) & (canon_cons_d.snp2_codons.length() == 7) & (
            canon_cons_d.mnv_codons.length() == 7))

canon_cons_pd = canon_cons_d.select("snp1_cons_term", "snp2_cons_term", "mnv_cons_term", "snp1_codons",
                                              "snp2_codons", "mnv_codons", "snp1_amino_acids", "snp2_amino_acids",
                                              "mnv_amino_acids", "snp1_lof", "snp2_lof", "mnv_lof", "transcript_id",
                                              "AC", "prev_AC", "n_hethet", "n_hethom", "n_homhom").to_pandas()
# get the most severe consequence
canon_cons_pd["snp1_sev"] = canon_cons_pd.snp1_cons_term.apply(lambda x: cons_term_most_severe(x))
canon_cons_pd["snp2_sev"] = canon_cons_pd.snp2_cons_term.apply(lambda x: cons_term_most_severe(x))
canon_cons_pd["mnv_sev"] = canon_cons_pd.mnv_cons_term.apply(lambda x: cons_term_most_severe(x))
# annotate the category
canon_cons_pd["categ"] = canon_cons_pd.apply(
            lambda x: mnv_category(x["snp1_sev"], x["snp2_sev"], x["mnv_sev"], x["snp1_amino_acids"],
                                   x["snp2_amino_acids"],
                                   x["mnv_amino_acids"]), axis=1)
#rewrite lof flags just in case if lof column is all None (NA) and that causes error
canon_cons_pd.snp1_lof = canon_cons_pd.snp1_lof.astype(str)
canon_cons_pd.snp2_lof = canon_cons_pd.snp2_lof.astype(str)
canon_cons_pd.mnv_lof = canon_cons_pd.mnv_lof.astype(str)


#and export
hl.Table.from_pandas(canon_cons_pd).write(sys.argv[1] + "_vepped_d{0}.ht".format(sys.argv[2]), overwrite=True)
