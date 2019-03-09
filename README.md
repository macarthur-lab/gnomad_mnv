# Multi-nucleotide variants (MNVs) in gnomAD 2.1

Welcome to the github page for the "Landscape of multi-nucleotide variants in 125,748 human exomes and 15,708 genomes". 

The `tutorials` directory contains all the information necessary to apply MNV discovery and annotation pipeline on your own dataset.
It also allows user to reproduce the main figure and most of the supplementary figure in the gnomAD MNV paper. 
Specifically, the tutorials directory consists of xx Jupyter notebooks:
1. `identify_mnv.ipynb` explains how to extract MNV from a vcf (or a matrix table).
2. `annotate_mnv.ipynb` explains how to annotate the functional consequences and category of MNVs
3. `functional_impact.ipynb` explains how to generate per gene (or per individual) statistics 
4. `global_mechanisms.ipynb` explains how to analyize the generation mechanisms of MNVs, genome wide
5. `per_region_mechanisms.ipynb` explains how to partition the genome into different functional 
6. `phase_sensitivity.ipynb` explains how the phasing sensitivity analysis could be performed
category, and analyze the MNVs across categories

Which figure and table in the paper is generated in which notebook is listed below:

| notebook  | main figure  | supplementary figure  | supplementary table  |   supplementary file  |
|---|---|---|---|---|
|identify_mnv.ipynb   |   | 11  |   |   |
|annotate_mnv.ipynb   | 2a  |   |   | 1  |
|functional_impact.ipynb  | 2  | 2  |   | 1  |
|global_mechanisms.ipynb  | 3, 4a  | 3-10,12-19  |   | 3  |
|per_region_mechanisms.ipynb  | 4b-d  | 20  | 2  |   |
|phase_sensitivity.ipynb  |   | 1,  | 1  |   |

However, since most of the analysis was performed in Hail, we recommend users who are not familiat with Hail to visit the [Hail tutorial page](https://hail.is/docs/0.2/tutorials-landing.html).


All the scripts used in the gnomAD MNV paper are stored in the `code` directory. 
However, note that due to the gnomAD sample data as well as the exome data of rare disease families being not publicly available, 
the scripts cannot be simply run in your local.

`util` contains functions and directories used in the analysis.
