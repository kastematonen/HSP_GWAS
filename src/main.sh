# HSP GWAS 

# merging files (HSP, IBD and HSCT)
./src/merging_files.sh

# cleaning HSP, IBD and HSCT data, removing bad samples based on FIMM's genotyping reports from HSP, IBD and HSCT data
./src/remove_bad_samples.sh

# genotype lift-over from hg19 to hg38 for merged HSP, IBD and HSCT data 
# changing SNP names in the reference file to match names in HSP, IBD and HSCT data to lift-over more variants
./src/lift-over_snp_names.R
# scripts based on protocol Genotyping chip data lift-over to reference genome build GRCh38/hg38 V.2
./src/genotype_lift-over_hg38.sh

# merging blood donor controls with HSP, IBD and HSCT data
# unifying SNP names of blood donor data for the merge
./src/unify_snp_names.R
# cleaning blood donor data, checking sex information and merging data
./src/merge_controls.sh

# LD pruning, PCA, and IBD
./src/pca_and_ibd.sh # includes removing relatives based on IBD analysis
# drawing plots of PCA results, checking IBD results
./src/pca_plots_and_ibd_list.R

# platform bias analysis by GWAS (HSCT against blood donors)
# creating dosage files
./src/dosage_file_for_platform_bias_gwas.sh #(needs snp lists made with ./src/platform_bias_gwas.R)
# gwas
./src/platform_bias_gwas.R #(needs dosage files made with ./src/dosage_file_for_platform_bias_gwas.sh)

# HLA imputation
# MCH region into a data set
./src/data_for_hla_imputation.sh
# imputing HLA alleles and their protein sequences
./src/hla_imputation.R

# GWAS (SPAtest)
# dosage file for HSP and IBD GWAS
./src/dosage_file_for_spatest.sh
# HSP and IBD GWAS  + HLA GWAS + Manhattan plots
./src/gwas.R

# HSP - IBD genetic correlation analysis
# common SNPs enriched between HSP and IBD GWAS
./src/joint_gwas.R
# annotating common enriched SNPs
#./src/snp_annotation.R

# creating tables from GWAS results 
./src/gwas_table.R
# checking posterior probabilities for HLA imputation
./scr/hla_postprob.R
# Manhattan plot for the HLA imputations
./scr/hla_manhattan.R


