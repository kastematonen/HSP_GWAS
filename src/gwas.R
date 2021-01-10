library(data.table)
library(tidyverse)
library(parallel)
library(gridExtra)
library(ggpubr)
library(qqman)
library(SPAtest)

# GWAS for one chromosome at a time in a loop

pca <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/pca/final_pca_ld_pruned_cleaned_merged_1_23.eigenvec', data.table=F)
all.samples <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_all/final.fam', data.table=F)[, c(2,5)]
hsct.samples <- fread('./data/hsct_register_genotypes/VPU_ILLUMINA_AUG_2018.fam', data.table=F)[, 2]
bd.samples <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_controls/cleaned_merged_1_23.fam', data.table=F)[, 2]
ibd.samples <- fread('./data/IBD.sample.list', data.table=F)[, 2]
hsp.samples <- fread('./data/HSP.sample.list', data.table=F)[, 2]

# phenotypes for HSP and IBD
phenotypes <- rep(NA, length(pca$IID))
phenotypes[pca$IID %in% hsp.samples] <- 1
phenotypes[pca$IID %in% hsct.samples] <- 0
phenotypes[pca$IID %in% bd.samples] <- 0
phenotypes[pca$IID %in% ibd.samples] <- 1
phenotypes.table <- data.frame(phenotypes)

# samples for HSP
rows.hsp <- rep(NA, length(pca$IID))
rows.hsp[pca$IID %in% hsp.samples] <- TRUE
rows.hsp[pca$IID %in% hsct.samples] <- TRUE
rows.hsp[pca$IID %in% bd.samples] <- TRUE
rows.hsp[pca$IID %in% ibd.samples] <- FALSE
# samples for IBD
rows.ibd <- rep(NA, length(pca$IID))
rows.ibd[pca$IID %in% hsp.samples] <- FALSE
rows.ibd[pca$IID %in% hsct.samples] <- TRUE
rows.ibd[pca$IID %in% bd.samples] <- TRUE
rows.ibd[pca$IID %in% ibd.samples] <- TRUE

# covariates: PC:s, sex, platform
colnames(all.samples) <- c("IID", "sex")
pca$row.numbers <- 1:nrow(pca)
pca.sex <- merge(pca, all.samples, by = "IID", all = TRUE)
pca.sex <- pca.sex[order(pca.sex$row.numbers),]

platform <- rep(0, length(pca$IID))
platform[pca$IID %in% bd.samples] <- 1
pca.sex$platform <- platform

# leave only principal components and sex as covariates
pca.no.id <- pca.sex[,c(paste("PC", 1:3, sep=""), "sex", "platform")]

# function for GWAS
association_test <- function(x, rows, type){
  print(x)
  if(x < 8){
    
    dosage.table <- fread(paste("/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/dosage_chr_",x,"_beginning.raw", sep=""), data.table=F)
    dosage.table <- dosage.table[,7:ncol(dosage.table)]
    print("tiedosto 1 luettu")
    spa.gwas <- ScoreTest_SPA(t(dosage.table[rows,]), phenotypes.table[rows,], pca.no.id[rows,], beta.out=T)
    
    spa.gwas <- data.frame(spa.gwas$p.value, colnames(dosage.table), p.value.NA = spa.gwas$p.value.NA, is.converge = spa.gwas$Is.converge, beta = spa.gwas$beta, SEbeta = spa.gwas$SEbeta)
  
    write.table(spa.gwas, paste("/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/",type,"_results_chr_", x, "_beginning.txt", sep=""), sep=" ", quote=F, row.names=F, col.names=T)
    
    dosage.table <- fread(paste("/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/dosage_chr_", x, "_end.raw", sep=""), data.table=F)
    dosage.table <- dosage.table[,7:ncol(dosage.table)]
    print("tiedosto 2 luettu")
    spa.gwas <- ScoreTest_SPA(t(dosage.table[rows,]), phenotypes.table[rows,], pca.no.id[rows,], beta.out=T)
    spa.gwas <- data.frame(spa.gwas$p.value, colnames(dosage.table), p.value.NA = spa.gwas$p.value.NA, is.converge = spa.gwas$Is.converge, beta = spa.gwas$beta, SEbeta = spa.gwas$SEbeta)
    write.table(spa.gwas, paste("/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/",type,"_results_chr_", x, "_end.txt", sep=""), sep=" ", quote=F, row.names=F, col.names=T)
    
  } else {
    
    dosage.table <- fread(paste("/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/dosage_chr_", x, ".raw", sep=""), data.table=F)
    dosage.table <- dosage.table[,7:ncol(dosage.table)]
    print("tiedosto luettu")
    spa.gwas <- ScoreTest_SPA(t(dosage.table[rows,]), phenotypes.table[rows,], pca.no.id[rows,], beta.out=T)
    spa.gwas <- data.frame(spa.gwas$p.value, colnames(dosage.table), p.value.NA = spa.gwas$p.value.NA, is.converge = spa.gwas$Is.converge, beta = spa.gwas$beta, SEbeta = spa.gwas$SEbeta)
    write.table(spa.gwas, paste("/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/",type,"_results_chr_", x, ".txt", sep=""), sep=" ", quote=F, row.names=F, col.names=T)
    
  }
  return(x)
}

# HSP GWAS
for (x in 1:23) {
  association_test(x, rows.hsp, "hsp")
}
# IBD GWAS
for (x in 1:23) {
  association_test(x, rows.ibd, "ibd")
}

# reading saved GWAS results to one table
combine.results <- function(type){
  all_snps <- data.frame(p.value = double(), snp = character(), stringsAsFactors=F)
  for (i in 1:23) {
      print(i)
      if(i < 8){
      
      results <- fread(paste("/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/",type,"_results_chr_", i, "_beginning.txt", sep=""), data.table=F)
      colnames(results) <- c("p.value", "snp", "p.na", "is.converge", "beta", "SEbeta")
      all_snps <- rbind(all_snps, results)
      results <- fread(paste("/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/",type,"_results_chr_", i, "_end.txt", sep=""), data.table=F)
      colnames(results) <- c("p.value", "snp", "p.na", "is.converge", "beta", "SEbeta")
      all_snps <- rbind(all_snps, results)
      
    } else {
      
      results <- fread(paste("/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/",type,"_results_chr_", i, ".txt", sep=""), data.table=F)
      colnames(results) <- c("p.value", "snp", "p.na", "is.converge", "beta", "SEbeta")
      all_snps <- rbind(all_snps, results)
      
    }
  }
  return(all_snps)
}

# combining HSP results
hsp.results <- combine.results("hsp")
# combining IBD results
ibd.results <- combine.results("ibd")

# Manhattan plots:
# other information for the plots
plots.info <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_all/final.bim', data.table=F)

# information for HSP
hsp.results$CHR <- plots.info$V1
hsp.results$SNP <- plots.info$V2
hsp.results$BP <- plots.info$V4
# Manhattan plot for HSP
jpeg("./results/images/koko_datan/manhattan_plot_hsp.jpeg", height= 15, width= 35, units="cm", res=400)
manhattan(hsp.results, chr="CHR", bp="BP", p="p.value", snp="SNP", ylim= c(0,10), col=c("blue", "brown"), cex = 0.6, cex.axis=0.8, las=2)
dev.off()

# Q-Q plot for HSP
jpeg("./results/images/koko_datan/qq_plot_hsp.jpeg", height= 15, width= 20, units="cm", res=300)
qq(hsp.results$p.value)
dev.off()

# information for IBD
ibd.results$CHR <- plots.info$V1
ibd.results$SNP <- plots.info$V2
ibd.results$BP <- plots.info$V4
# Manhattan plot for IBD
jpeg("./results/images/koko_datan/manhattan_plot_ibd.jpeg", height= 15, width= 35, units="cm", res=400)
manhattan(ibd.results, chr="CHR", bp="BP", p="p.value", snp="SNP", ylim= c(0,10), col=c("blue", "brown"), cex = 0.6, cex.axis=0.8, las=2)
dev.off()
# Q-Q plot for IBD
jpeg("./results/images/koko_datan/qq_plot_ibd.jpeg", height= 15, width= 20, units="cm", res=300)
qq(ibd.results$p.value)
dev.off()

# save GWAS results
write.table(hsp.results, "/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/hsp_results.txt", sep=" ", quote=F, row.names=F, col.names=T)
write.table(ibd.results, "/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/ibd_results.txt", sep=" ", quote=F, row.names=F, col.names=T)
###############################################################################################################
# HLA imputation: association analysis

# loading allelle dosage table
hla.dosage <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/hla_imputation/hla_dosage.txt', data.table=F)[,-1]

# association analysis for HSP
spa.gwas.hla.hsp <- ScoreTest_SPA(t(hla.dosage[rows.hsp,]), phenotypes.table[rows.hsp,], pca.no.id[rows.hsp,], beta.out=T, beta.Cutoff=5*10^-2)
# association analysis for IBD
spa.gwas.hla.ibd <- ScoreTest_SPA(t(hla.dosage[rows.ibd,]), phenotypes.table[rows.ibd,], pca.no.id[rows.ibd,], beta.out=T, beta.Cutoff=5*10^-2)

# saving the results
spa.gwas.hsp.save <- data.frame(name =  colnames(hla.dosage), p.value = spa.gwas.hla.hsp$p.value, p.value.NA = spa.gwas.hla.hsp$p.value.NA, is.converge = spa.gwas.hla.hsp$Is.converge, beta = spa.gwas.hla.hsp$beta, SEbeta = spa.gwas.hla.hsp$SEbeta)
spa.gwas.ibd.save <- data.frame(name =  colnames(hla.dosage), p.value = spa.gwas.hla.ibd$p.value, p.value.NA = spa.gwas.hla.ibd$p.value.NA, is.converge = spa.gwas.hla.ibd$Is.converge, beta = spa.gwas.hla.ibd$beta, SEbeta = spa.gwas.hla.ibd$SEbeta)
write.table(spa.gwas.hsp.save, "/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/hla_hsp_results.txt", sep=" ", quote=F, row.names=F, col.names=T)
write.table(spa.gwas.ibd.save, "/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/hla_ibd_results.txt", sep=" ", quote=F, row.names=F, col.names=T)

# Q-Q plot for HSP
jpeg("./results/images/koko_datan/qq_plot_hla_hsp.jpeg", height= 15, width= 20, units="cm", res=300)
qq(spa.gwas.hsp.save$p.value)
dev.off()
# Q-Q plot for IBD
jpeg("./results/images/koko_datan/qq_plot_hla_ibd.jpeg", height= 15, width= 20, units="cm", res=300)
qq(spa.gwas.ibd.save$p.value)
dev.off()
#####################################################################################################################################
# HLA protein sequences: association analysis

# loading amino acid matrix
hla.aa.dosage <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/hla_imputation/hla_prot_dosage.txt', data.table=F)[,-1]
# spatest for HSP
spa.gwas.hla.hsp.aa <- ScoreTest_SPA(t(hla.aa.dosage[rows.hsp,]), phenotypes.table[rows.hsp,], pca.no.id[rows.hsp,], beta.out=T)
# spatest for IBD
spa.gwas.hla.ibd.aa <- ScoreTest_SPA(t(hla.aa.dosage[rows.ibd,]), phenotypes.table[rows.ibd,], pca.no.id[rows.ibd,], beta.out=T)

# saving the results
spa.gwas.hla.hsp.aa.save <- data.frame(name = colnames(hla.aa.dosage), p.value = spa.gwas.hla.hsp.aa$p.value, p.value.NA = spa.gwas.hla.hsp.aa$p.value.NA, is.converge = spa.gwas.hla.hsp.aa$Is.converge, beta = spa.gwas.hla.hsp.aa$beta, SEbeta = spa.gwas.hla.hsp.aa$SEbeta)
spa.gwas.hla.ibd.aa.save <- data.frame(name = colnames(hla.aa.dosage), p.value = spa.gwas.hla.ibd.aa$p.value, p.value.NA = spa.gwas.hla.ibd.aa$p.value.NA)
write.table(spa.gwas.hla.hsp.aa.save, "/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/hla_prot_hsp_results.txt", sep=" ", quote=F, row.names=F, col.names=T)
write.table(spa.gwas.hla.ibd.aa.save, "/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/hla_prot_ibd_results.txt", sep=" ", quote=F, row.names=F, col.names=T)

# Q-Q plot for HSP
jpeg("./results/images/koko_datan/qq_plot_hla_prot_hsp.jpeg", height= 15, width= 20, units="cm", res=300)
qq(spa.gwas.hla.hsp.aa.save$p.value)
dev.off()
# Q-Q plot for IBD
jpeg("./results/images/koko_datan/qq_plot_hla_prot_ibd.jpeg", height= 15, width= 20, units="cm", res=300)
qq(spa.gwas.hla.ibd.aa.save$p.value)
dev.off()

