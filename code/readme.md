This directory stores the codes used for the [gnomAD MNV preprint](https://www.biorxiv.org/content/10.1101/573378v2), as well as codes for users can run for their own analysis.
(Other than the `Codes for users to run`, simple copy&paste of the code in your local will mostly not work, due to limited access to the individual level data.)

## Codes for users to run

`get_mnv.py` can be used to identify MNVs in your dataset

`annotate_vep_mnv.py` can be used to annotate the functional consequences of MNVs 

(Tweek the `get_tnv_gnomAD.py` below to look into MNVs consisting of 3 SNVs.)

## Codes related to MNV discovery and functional impact annotation:

`get_mnv_per_variant.py` was used to identify the MNVs in gnomAD genome, autosome region

`exome_mnv_per_variant_autosome_for_release.py` was used to identify and annotate the MNVs in gnomAD exome, autosome region

`exome_mnv_per_variant_sexchr_for_release.py` was used to identify and annotate the MNVs in gnomAD exome, sex chromosome (including pseudo-autosomal) region

`genome_coding_mnv_per_variant_autosome_for_release.py` was used to identify and annotate the MNVs in the coding region of gnomAD genome, autosome region

`get_tnv_gnomAD.py` was used to identify and annotate the MNVs in gnomAD exome, specifically those consisting of 3 SNVs

`mnv_coding_parse.py` was used to parse the MNV list in coding region and construct a dataframe for final release


## Codes related to phasing sensitivity and specificity evaluation:

`hethet_trio_check.py` was used to investigate the MNVs unphased by trio based phasing

`inspect_unphased.py` was used to investigate the MNVs unphased by read based phasing

`sensitivity_proband_full.py` was used to compare the phase sensitivity of read and trio based phasing


## Codes related to MNV mechanism exploration:

`annotate_context.py` was used to annotate the local context of MNV

`classify_onestep.py` was used to extract the one-step MNVs (MNVs with AC1==AC2) 

`density_per_func_annot.py` was used to calculate the MNV density per functional annotation

`get_cnt_matrix.py` was used to generate the overall count matrix of MNVs

`get_cnt_matrix_hom.py` was used to generate the count matrix of homozygote MNVs

`get_cnt_matrix_nonpass.py` was used to generate the count matrix of MNVs that did not pass QC (=filtered out)

`get_cnt_matrix_per_annot.py` was used to generate the count matrix of MNVs for different functional annotations

`get_coverage_8bp.py` was used to calculate the average coverage (read depth) in gnomAD (which was used to construct the base line mutation rate)

`per_sample_stats.py` was used to calculate the (aggreagated stats of) number of MNVs per sample

`vs_mnv10_enrichment.py` was used to quantify and compare the density enrichment of MNVs across distance

`MNV_logoanalysis.R` was used to create the bits representation of MNV contexts



See also `tutorials` for related codes/figures.



