library(data.table)
library(tidyverse)
library(parallel)
library(gridExtra)
library(ggpubr)
library(qqman)
library(SPAtest)

spa.gwas.hsp <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/hsp_results.txt', data.table=F)[,1:2]
spa.gwas.ibd <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/ibd_results.txt', data.table=F)[,1:2]

hla.hsp <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/hla_hsp_results.txt', data.table=F)[,1:3]
hla.ibd <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/hla_ibd_results.txt', data.table=F)[,1:3]

hla.hsp.prot <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/hla_prot_hsp_results.txt', data.table=F)[,1:3]
hla.ibd.prot <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/hla_prot_ibd_results.txt', data.table=F)[,1:3]

# order the SNPs by their p value
spa.gwas.hsp <- spa.gwas.hsp[order(spa.gwas.hsp[,1]),]
spa.gwas.ibd <- spa.gwas.ibd[order(spa.gwas.ibd[,1]),]
colnames(spa.gwas.hsp) <- c("p-value", "snp")
colnames(spa.gwas.ibd) <- c("p-value", "snp")

hla.hsp <- hla.hsp[order(hla.hsp[,2]),]
hla.ibd <- hla.ibd[order(hla.ibd[,2]),]
colnames(hla.hsp) <- c("snp", "p-value", "p-value.na")
colnames(hla.ibd) <- c("snp", "p-value", "p-value.na")

hla.hsp.prot <- hla.hsp.prot[order(hla.hsp.prot[,2]),]
hla.ibd.prot <- hla.ibd.prot[order(hla.ibd.prot[,2]),]
colnames(hla.hsp.prot) <- c("snp", "p-value", "p-value.na")
colnames(hla.ibd.prot) <- c("snp", "p-value", "p-value.na")

###############################################################################################
# function for determining enrichment of common SNPs for gwas results at different sizes  of results data
find.enrichment.limit <- function(gwas.results.hsp, gwas.results.ibd, by){
  
  limit.values <- seq(0, nrow(gwas.results.hsp), by = by)
  results <- data.frame("M" = limit.values, "p.value" = rep(NA, length(limit.values)))
  
  for (i in 1:length(limit.values)) {
    hsp.snp <- gwas.results.hsp[1:limit.values[i],"snp"]
    ibd.snp <- gwas.results.ibd[1:limit.values[i],"snp"]
    
    common.snps <- hsp.snp %in% ibd.snp
    
    common.expected <- results$M[i]
    divergent.expected <- nrow(gwas.results.hsp) - common.expected
    
    significance <- phyper(sum(common.snps, na.rm = TRUE)-1, common.expected, divergent.expected, limit.values[i], lower.tail= F)
    results[i,2] <- significance
  }
  
  results$min.p <- log10(results$p.value)*(-1)
  
  return(results)
}
###############################################################################################
# finding enrichment scores for normal gwas
gwas.all <- find.enrichment.limit(spa.gwas.hsp, spa.gwas.ibd, 100)

# how many common SNPs by the enrichment limit
gwas.all.no.na <- gwas.all[gwas.all$p.value != 0,]
m.value <- gwas.all.no.na[gwas.all.no.na$p.value == min(gwas.all.no.na$p.value),1]
hsp.snp <- spa.gwas.hsp[1:m.value,4]
ibd.snp <- spa.gwas.ibd[1:m.value,4]
common.snps <- hsp.snp %in% ibd.snp
snps <- hsp.snp[common.snps]
sum(common.snps, na.rm = TRUE) # 619 common SNPs
write.table(snps, "/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/joint_gwas/enriched_commom_snps.txt", sep=" ", quote=F, row.names=F, col.names=F)

######################################################################################################
# for imputed HLA alleles
hla.hsp.no.na <- hla.hsp[!is.na(hla.hsp$`p-value`),]
hla.ibd.no.na <- hla.ibd[!is.na(hla.ibd$`p-value`),]

gwas.hla <- find.enrichment.limit(hla.hsp.no.na, hla.ibd.no.na, 1)

# how many common SNPs by the enrichment limit
gwas.hla.no.na <- gwas.hla[gwas.hla$p.value != 0,]
m.value <- gwas.hla.no.na[gwas.hla.no.na$p.value == min(gwas.hla.no.na$p.value),1]
hsp.snp <- hla.hsp.no.na[1:m.value,"snp"]
ibd.snp <- hla.ibd.no.na[1:m.value,"snp"]
common.snps <- hsp.snp %in% ibd.snp
snps <- hsp.snp[common.snps]
sum(common.snps, na.rm = TRUE) # 112
write.table(snps, "/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/joint_gwas/hla_enriched_commom_snps.txt", sep=" ", quote=F, row.names=F, col.names=F)

############################################################################################################
# joint-GWAS for HLA protein sequences

hla.hsp.prot.no.na <- hla.hsp.prot[!is.na(hla.hsp.prot$`p-value`),]
hla.ibd.prot.no.na <- hla.ibd.prot[!is.na(hla.ibd.prot$`p-value`),]

gwas.hla.prot <- find.enrichment.limit(hla.hsp.prot.no.na, hla.ibd.prot.no.na, 1)

# how many common SNPs by the enrichment limit
gwas.hla.prot.no.na <- gwas.hla.prot[gwas.hla.prot$p.value != 0,]
m.value <- gwas.hla.prot.no.na[gwas.hla.prot.no.na$p.value == min(gwas.hla.prot.no.na$p.value),1]
hsp.snp <- hla.hsp.prot.no.na[1:m.value,"snp"]
ibd.snp <- hla.ibd.prot.no.na[1:m.value,"snp"]
common.snps <- hsp.snp %in% ibd.snp
snps <- hsp.snp[common.snps]
sum(common.snps, na.rm = TRUE) # 381
write.table(snps, "/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/joint_gwas/hla_prot_enriched_commom_snps.txt", sep=" ", quote=F, row.names=F, col.names=F)

####################################################################################################################
# enrichment with GWAS results in a random order

# all SNPs
spa.gwas.hsp.random <- spa.gwas.hsp[sample(nrow(spa.gwas.hsp)),]
spa.gwas.ibd.random <- spa.gwas.ibd[sample(nrow(spa.gwas.ibd)),]

gwas.all.random <- find.enrichment.limit(spa.gwas.hsp.random, spa.gwas.ibd.random, 100)

# HLA alleles
hla.hsp.no.na.random <- hla.hsp.no.na[sample(nrow(hla.hsp.no.na)),]
hla.ibd.no.na.random <- hla.ibd.no.na[sample(nrow(hla.ibd.no.na)),]
hla.random <- find.enrichment.limit(hla.hsp.no.na.random, hla.ibd.no.na.random, 1)

# HLA amino acids
hla.hsp.prot.no.na.random <- hla.hsp.prot.no.na[sample(nrow(hla.hsp.prot.no.na)),]
hla.ibd.prot.no.na.random <- hla.ibd.prot.no.na[sample(nrow(hla.ibd.prot.no.na)),]

prot.random <- find.enrichment.limit(hla.hsp.prot.no.na.random, hla.ibd.prot.no.na.random, 1)

#################################################################################################################
# plots

gwas.all$data.type <- "ordered"
gwas.all.random$data.type <- "unordered"
gwas.all.combined <- rbind(gwas.all, gwas.all.random)

gwas.hla$data.type <- "ordered"
hla.random$data.type <- "unordered"
hla.combined <- rbind(gwas.hla, hla.random)

gwas.hla.prot$data.type <- "ordered"
prot.random$data.type <- "unordered"
prot.combined <- rbind(gwas.hla.prot, prot.random)

all.plot <- ggplot(gwas.all.combined, aes(x = M, y = min.p)) + geom_line(aes(color = data.type)) + labs(title= "Enrichment of shared SNPs", x = "Length of SNP subset", y= "Enrichment, -log10(p)", tag = "A)", colour= "GWAS results") + theme_minimal() + scale_color_manual(values = c("blue", "brown"))

hla.plot <- ggplot(hla.combined, aes(x = M, y = min.p)) + geom_line(aes(color = data.type)) + labs(title= "Enrichment of shared HLA alleles", x = "Length of HLA allele subset", y= "Enrichment, -log10(p)", tag = "B)") + theme_minimal() + scale_color_manual(values = c("blue", "brown"))

prot.plot <- ggplot(prot.combined, aes(x = M, y = min.p)) + geom_line(aes(color = data.type)) + labs(title= "Enrichment of shared polymorphic amino acids", x = "Length of  amino acid subset", y= "Enrichment, -log10(p)", tag = "C)") + theme_minimal() + scale_color_manual(values = c("blue", "brown"))

ggarrange(all.plot, hla.plot, prot.plot, ncol=1, nrow=3, common.legend=T, legend="bottom")
ggsave("./results/images/koko_datan/joint_gwas_all.jpeg", device="jpeg", height= 35, width= 25, units="cm")







