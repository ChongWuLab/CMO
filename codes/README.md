# Code

Code for replicating all the results presented in the manuscript

## Real Data analyses

* step1_prepare_weights.R: process the data for some further analysis

* step2a_generate_weights_gene.R and step2b_generate_weights2_gene.R: mapping imputed DNAm to gene body regions

* step3a_generate_weights_enhancer and step3b_generate_weights2_enhancer.R: mapping imputed DNAm to enhancer regions

* step4_CMO_main.R: conduct CMO test for AD GWAS summary results

* step5_CMO_main_UKBiobankImage.R: conduct CMO test for UK Biobank Imaging related GWAS summary results

* /resAna: Generate all the figures and tables in the manuscript

## Support function

* dist_support.R: hypothesis testing support function
* generate_pbs_files.R: submit jobs more efficiently
* mwas_support.R: CMO support function
* prepare_summary_stat.R: prepare GWAS summary results


