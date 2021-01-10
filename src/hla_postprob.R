library(data.table)
library(tidyverse)

# information for the rows
pca <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/pca/final_pca_ld_pruned_cleaned_merged_1_23.eigenvec', data.table=F) # all samples
hsct.samples <- fread('./data/hsct_register_genotypes/VPU_ILLUMINA_AUG_2018.fam', data.table=F)[, 2]
bd.samples <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_controls/cleaned_merged_1_23.fam', data.table=F)[, 2]
ibd.samples <- fread('./data/IBD.sample.list', data.table=F)[, 2]
hsp.samples <- fread('./data/HSP.sample.list', data.table=F)[, 2]

# samples for all
rows.all <- rep(NA, length(pca$IID))
rows.all[pca$IID %in% hsp.samples] <- TRUE
rows.all[pca$IID %in% hsct.samples] <- TRUE
rows.all[pca$IID %in% bd.samples] <- TRUE
rows.all[pca$IID %in% ibd.samples] <- FALSE
# samples for controls
rows.controls <- rep(NA, length(pca$IID))
rows.controls[pca$IID %in% hsp.samples] <- FALSE
rows.controls[pca$IID %in% hsct.samples] <- TRUE
rows.controls[pca$IID %in% bd.samples] <- TRUE
rows.controls[pca$IID %in% ibd.samples] <- FALSE
# samples for HSP
rows.hsp <- rep(NA, length(pca$IID))
rows.hsp[pca$IID %in% hsp.samples] <- TRUE
rows.hsp[pca$IID %in% hsct.samples] <- FALSE
rows.hsp[pca$IID %in% bd.samples] <- FALSE
rows.hsp[pca$IID %in% ibd.samples] <- FALSE

# table for the results
results <- data.frame(hla.allele = rep(NA,3), mean.hsp = rep(NA,3), mean.controls = rep(NA,3), sd.hsp = rep(NA,3), sd.controls = rep(NA,3))

# info for the genes and alleles
files <- c("/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/hla_imputation//hla_geno_DQB1", "/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/hla_imputation//hla_geno_DQA1", "/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/hla_imputation//hla_geno_DRB1")
genes <- c("DQB1", "DQA1", "DRB1")
alleles <- c("05:01", "01:01", "01:01")

results$hla.allele <- paste(genes, alleles, sep="*")

for (i in 1:length(genes)) {
  hla.geno <- fread(files[i], data.table=F)
  
  # leave only correct people:
  geno.hsp <- hla.geno[rows.hsp,]
  geno.controls <- hla.geno[rows.controls,]
  
  # find rows for the right allele 
  hsp.postprob.rows <- geno.hsp$allele1 == alleles[i] | geno.hsp$allele2 == alleles[i]
  controls.postprob.rows <- geno.controls$allele1 == alleles[i] | geno.controls$allele2 == alleles[i]
  
  # leave rows for the right allele for imputation results
  geno.hsp.allele <- geno.hsp[hsp.postprob.rows,]
  geno.controls.allele <- geno.controls[controls.postprob.rows,]
  
  # mean for post prob values
  mean.hsp <- mean(geno.hsp.allele$prob)
  mean.controls <- mean(geno.controls.allele$prob)
  # standard deviation post prob values
  sd.hsp <- sd(geno.hsp.allele$prob)
  sd.controls <- sd(geno.controls.allele$prob)
  
  # results <- rbind(results, c(alleles[i], mean.hsp, mean.controls, sd.hsp, sd.controls))
  results[i,2:5] <- c(mean.hsp, mean.controls, sd.hsp, sd.controls)
}

write.table(results, "./tmp/hla-postprob.txt", sep=" ", quote=F, row.names=F, col.names=T)








