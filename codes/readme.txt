#####################################
# Useful link
#####################################
CrossMap: Transform enhancer bed HG38 to HG19: http://useast.ensembl.org/Homo_sapiens/Tools/AssemblyConverter?db=core

http://gwas-browser.nygenome.org/downloads/gwas-browser/

Gene enhancer: https://genecards.weizmann.ac.il/geneloc_prev/index.shtml
 4.7

Gene ID: ensembl ID, http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=722057481_agxhIRPO1gqvYO93g0N56ItH8Yca&clade=mammal&org=Human&db=hg38&hgta_group=genes&hgta_track=refSeqComposite&hgta_table=0&hgta_regionType=genome&position=chr1%3A11%2C102%2C837-11%2C267%2C747&hgta_outputType=primaryTable&hgta_outFileName=

mQTL weights: http://bbmri.researchlumc.nl/atlas/#data
gwax_AD: http://gwas-browser.nygenome.org/downloads/gwas-browser/


###################################
# code explanation:
###################################

step1_prepare_weights.R: process the data for some further analysis

step2a_generate_weights_gene.R and step2b_generate_weights2_gene.R: mapping imputed DNAm to gene body regions

step3a_generate_weights_enhancer and step3b_generate_weights2_enhancer.R: mapping imputed DNAm to enhancer regions

# conduct test
step4_CMO_main.R: conduct CMO test for AD GWAS summary results
step5_CMO_main_UKBiobankImage.R: conduct CMO test for UK Biobank Imaging related GWAS summary results


########################
## Supporting function
########################
