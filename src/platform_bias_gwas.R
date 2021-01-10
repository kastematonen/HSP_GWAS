library(data.table)
library(tidyverse)
library(qqman)
library(SPAtest)

# divide data by chromosomes for smaller files: list SNPs for extraction
all_snips <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_all/final_cleaned_merged_1_23.bim', data.table=F)[,1:2]
snips_for_extracting <- function(x){
  snips <- all_snips$V2[all_snips$V1 == x]
  if (x < 8) {
    write.table(snips[1:(as.integer(length(snips)/2))], file=paste("/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/snips_chr_", x, "_beginning.txt", sep = ""), sep="\t", quote=F, row.names=F, col.names=F)
    write.table(snips[(as.integer(length(snips)/2)+1):length(snips)], file=paste("/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/snips_chr_", x, "_end.txt", sep = ""), sep="\t", quote=F, row.names=F, col.names=F)
  } else {
    write.table(snips, file=paste("/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/snips_chr_", x, ".txt", sep = ""), sep="\t", quote=F, row.names=F, col.names=F)
  }
  
  return(x)
}
sapply(unique(all_snips$V1), snips_for_extracting)

# -> making dosage tables with ./src/dosage_file_for_platform_bias_gwas.sh

#############################################################################################################################################

# information for GWAS: which samples to include, phenotypes, covariates

blood.donor.samples <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_controls/cleaned_merged_1_23.fam', data.table=F)[, c(2,5)]
hsct.samples <- fread('./data/hsct_register_genotypes/VPU_ILLUMINA_AUG_2018.fam', data.table=F)
hsct.samples <- hsct.samples[hsct.samples$V6 == 1,c(2,5)]
pca <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/pca/final_pca_ld_pruned_cleaned_merged_1_23.eigenvec', data.table=F)

# phenotypes: 
phenotypes <- rep(NA, length(pca$IID))
phenotypes[pca$IID %in% hsct.samples$V2] <- 1
phenotypes[pca$IID %in% blood.donor.samples$V2] <- 0
phenotypes.table <- data.frame(phenotypes)

# samples to be included
rows <- rep(FALSE, length(pca$IID))
rows[pca$IID %in% hsct.samples$V2] <- TRUE
rows[pca$IID %in% blood.donor.samples$V2] <- TRUE

# covariates: PC:s, sex 
sex.blood.donors <- blood.donor.samples[blood.donor.samples$V2 %in% pca$IID,]
sex.hsct <- hsct.samples[hsct.samples$V2 %in% pca$IID,]
colnames(sex.blood.donors) <- c("IID", "sex")
colnames(sex.hsct) <- c("IID", "sex")

pca$row.numbers <- 1:nrow(pca)

pca.sex <- merge(pca, sex.hsct, by = "IID", all = TRUE)
pca.sex <- merge(pca.sex, sex.blood.donors, by = "IID", all = TRUE)
pca.sex$sex <- rowSums(pca.sex[,c("sex.x", "sex.y")], na.rm=TRUE)

pca.sex <- pca.sex[order(pca.sex$row.numbers),]

# leave only principal components and sex as covariates
pca.no.id <- pca.sex[,c(paste("PC", 1:5, sep=""), "sex")]

#############################################################################################################################################
# platform-bias GWAS

association_test <- function(x){
  print(x)
  if(x < 8){
    dosage.table <- fread(paste("/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/dosage_chr_", x, "_beginning.raw", sep=""), data.table=F)
    dosage.table <- dosage.table[,7:ncol(dosage.table)]
    print("tiedosto 1 luettu")
    spa.gwas <- ScoreTest_SPA(t(dosage.table[rows,]), phenotypes.table[rows,], pca.no.id[rows,], beta.out=T)
    spa.gwas <- data.frame(spa.gwas$p.value, colnames(dosage.table))
    write.table(spa.gwas, paste("/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/results_chr_", x, "_beginning.txt", sep=""), sep=" ", quote=F, row.names=F, col.names=T)
    
    dosage.table <- fread(paste("/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/dosage_chr_", x, "_end.raw", sep=""), data.table=F)
    dosage.table <- dosage.table[,7:ncol(dosage.table)]
    print("tiedosto 2 luettu")
    spa.gwas <- ScoreTest_SPA(t(dosage.table[rows,]), phenotypes.table[rows,], pca.no.id[rows,], beta.out=T)
    spa.gwas <- data.frame(spa.gwas$p.value, colnames(dosage.table))
    write.table(spa.gwas, paste("/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/results_chr_", x, "_end.txt", sep=""), sep=" ", quote=F, row.names=F, col.names=T)
  } else {
    dosage.table <- fread(paste("/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/dosage_chr_", x, ".raw", sep=""), data.table=F)
    dosage.table <- dosage.table[,7:ncol(dosage.table)]
    print("tiedosto luettu")
    spa.gwas <- ScoreTest_SPA(t(dosage.table[rows,]), phenotypes.table[rows,], pca.no.id[rows,], beta.out=T)
    spa.gwas <- data.frame(spa.gwas$p.value, colnames(dosage.table))
    write.table(spa.gwas, paste("/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/results_chr_", x, ".txt", sep=""), sep=" ", quote=F, row.names=F, col.names=T)
  }
  return(x)
}

sapply(1:23, association_test)

# GWAS results to a single file

all_snps <- data.frame(p.value = double(), snp = character(), stringsAsFactors=F)

for (i in 1:23) {
  print(i)
  if(i < 8){

    results <- fread(paste("/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/results_chr_", i, "_beginning.txt", sep=""), data.table=F)
    colnames(results) <- c("p.value", "snp")
    all_snps <- rbind(all_snps, results)
    results <- fread(paste("/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/results_chr_", i, "_end.txt", sep=""), data.table=F)
    colnames(results) <- c("p.value", "snp")
    all_snps <- rbind(all_snps, results)

  } else {

    results <- fread(paste("/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/results_chr_", i, ".txt", sep=""), data.table=F)
    colnames(results) <- c("p.value", "snp")
    all_snps <- rbind(all_snps, results)
    
  }
}

# Manhattan plots of GWAS results

# information for the plot
plot.info <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_all/final_cleaned_merged_1_23.bim', data.table=F)
all_snps$CHR <- plot.info$V1
all_snps$SNP <- plot.info$V2
all_snps$BP <- plot.info$V4

# Manhattan plot 
# jpeg("./results/images/manhattanplotti_platform_bias.jpg", height= 20, width= 35, units="cm", res=300)
jpeg("./results/images/koko_datan/platform_bias.jpg", height= 15, width= 35, units="cm", res=400)
manhattan(all_snps, chr="CHR", bp="BP", p="p.value", snp="SNP", ylim= c(0,10), col=c("blue", "brown"),  cex = 0.6, cex.axis=0.8, las=2)
abline(h = - log10(1 * (10 ^ -4)), col = "black")
dev.off()

#-------------------------------------------------------------------------------------------------------
# Manhattan plot with & without ylim 
jpeg("./results/images/koko_datan/platform_bias_both.jpeg", height= 35, width= 35, units="cm", res=400)
par(mfrow=c(2,1))
manhattan(all_snps, chr="CHR", bp="BP", p="p.value", snp="SNP", ylim= c(0,10), col=c("blue", "brown"), cex = 0.6, cex.axis=0.8, las=2)
abline(h = - log10(1 * (10 ^ -4)), col = "black")
mtext("A)", side=3, line=1, cex=2, adj=-0.05)
manhattan(all_snps, chr="CHR", bp="BP", p="p.value", snp="SNP", col=c("blue", "brown"), cex = 0.6, cex.axis=0.8, las=2)
mtext("B)", side=3, line=1, cex=2, adj=-0.05)
abline(h = - log10(1 * (10 ^ -4)), col = "black")
dev.off()

par(mfrow=c(1,1))
#----------------------------------------------------------------------------------------------------------

# Q-Q plot
jpeg("./results/images/koko_datan/qq_plot_platform_bias.jpeg", height= 15, width= 20, units="cm", res=300)
qq(all_snps$p.value)
dev.off()
  
# removing SNPs associated with platform:
for_removal <- all_snps[all_snps$p.value < 0.0001,] # 251 

write.table(for_removal$SNP, file=paste("/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/platform_bias/snps_for_removal.txt", sep = ""), sep="\t", quote=F, row.names=F, col.names=F)

# -> remove associated SNPs from data with ./src/dosage_file_for_platform_bias_gwas.sh
