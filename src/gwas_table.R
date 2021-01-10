library(data.table)
library(tidyverse)
library(parallel)
library(gridExtra)
library(ggpubr)

#####################################################################################################
# table of GWAS results: SNPs

gwas.results <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/hsp_results.txt', data.table=F)
gwas.results <- gwas.results[order(gwas.results[,1]),]
gwas.results <- gwas.results[gwas.results$p.value < 5 * 10 ** -8,]


# SNP names to rs codes
change.names <- gwas.results[,8:9]
for (i in 1:nrow(change.names)) {
  if(change.names[i,1] %like% "-"){
    change.names[i,1] <- str_split_fixed(change.names[i,1], "-", n=2)[2]
  }
}

change.names$SNP <- tolower(change.names$SNP)
colnames(change.names) <- c("snp_rs", "BP")

gwas.results <- merge(gwas.results, change.names, by = "BP")

# minor allele
minor.allele <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_all/final.bim', data.table=F)
in.results <- minor.allele$V2 %in% gwas.results$SNP
minor.allele <- minor.allele[in.results,4:5]
colnames(minor.allele) <- c("BP", "minor.allele")

gwas.results <- merge(gwas.results, minor.allele, by = "BP")

# allele frequencies
dosage.table.beginning <- fread(paste("/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/dosage_chr_6_beginning.raw", sep=""), data.table=F)
rows <- colnames(dosage.table.beginning) %in% gwas.results$snp
dosage.table.beginning <- dosage.table.beginning[,rows]

people <- fread("/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_all/final_chr_6_beginning.fam", data.table=F)

hsct.samples <- fread('./data/hsct_register_genotypes/VPU_ILLUMINA_AUG_2018.fam', data.table=F)[, 2]
bd.samples <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_controls/cleaned_merged_1_23.fam', data.table=F)[, 2]
hsp.samples <- fread('./data/HSP.sample.list', data.table=F)[, 2]
ibd.samples <- fread('./data/IBD.sample.list', data.table=F)[, 2]

rows.hsp <- rep(FALSE, length(people$V2))
rows.hsp[people$V2 %in% hsp.samples] <- TRUE

rows.ibd <- rep(FALSE, length(people$V2))
rows.ibd[people$V2 %in% ibd.samples] <- TRUE

rows.controls <- rep(FALSE, length(people$V2))
rows.controls[people$V2 %in% hsct.samples] <- TRUE
rows.controls[people$V2 %in% bd.samples] <- TRUE

dosage.hsp <- dosage.table.beginning[rows.hsp,]
dosage.controls <- dosage.table.beginning[rows.controls,]

names(gwas.results)[names(gwas.results) == "snp"] <- "name"

#  function for calculating allele frequencies:
allele.freq <- function(dosage.hsp, dosage.controls, results.table){
  
  snp.sum <- colSums(dosage.hsp, na.rm=T)
  people.na <- apply(dosage.hsp, 2, function(x) sum(is.na(x)))
  freq.hsp <- data.table("name" = names(people.na))
  freq.hsp$freq.hsp <- snp.sum / (2*(nrow(dosage.hsp) - people.na))
  results.table <- merge(results.table, freq.hsp, by = "name")
  
  snp.sum <- colSums(dosage.controls, na.rm=T)
  people.na <- apply(dosage.controls, 2, function(x) sum(is.na(x)))
  freq.controls <- data.table("name" = names(people.na))
  freq.controls$freq.controls <- snp.sum / (2*(nrow(dosage.controls) - people.na))
  results.table <- merge(results.table, freq.controls, by = "name")
  
  return(results.table)
}

# allele frequencies
gwas.results <- allele.freq(dosage.hsp, dosage.controls, gwas.results)

# function for calculating odds ratio
odds.ratio <- function(results.table){
  results.table$odds.ratio <- exp(results.table$beta)
  results.table$CI.upper <- exp(results.table$beta + results.table$SEbeta)
  results.table$CI.lower <- exp(results.table$beta - results.table$SEbeta)
  
  return(results.table)
}

# odds ratios
gwas.results <- odds.ratio(gwas.results)

gwas.results.final <- gwas.results[,c(10, 2, 11, 3, 14, 15, 16, 12, 13)]

# function for rounding the results
round.results <- function(results.table){
  results.table$odds.ratio <- round(results.table$odds.ratio,digits=2)
  results.table$p.value <- signif(results.table$p.value,digits=3)
  results.table$CI.upper <- round(results.table$CI.upper,digits=2)
  results.table$CI.lower <- round(results.table$CI.lower,digits=2)
  results.table$freq.controls <- round(results.table$freq.controls,digits=2)
  results.table$freq.hsp <- round(results.table$freq.hsp,digits=2)
  results.table$or <- paste(results.table$odds.ratio, " [", results.table$CI.lower, "-", results.table$CI.upper, "]", sep="")
  return(results.table)
}

# rounding results
gwas.results.final <- round.results(gwas.results.final)

gwas.results.final <- gwas.results.final[,c(1,2,3,4,10,8,9)]
colnames(gwas.results.final) <- c("SNP", "position", "minor allele", "P", "OR[CI 95%]", "freq.hsp", "freq.controls")
gwas.results.final <- gwas.results.final[order(gwas.results.final$P),]
write.table(gwas.results.final, "./results/images/koko_datan/HLA_tables/hsp_table.txt", sep=" ", quote=F, row.names=F, col.names=T)
#####################################################################################################
# hla allele table
hla <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/hla_hsp_results.txt', data.table=F)

hla <- hla[order(hla[,2]),]
hla <- hla[hla$p.value < 5 * 10 ** -8,]
hla <- hla[!(is.na(hla$name)),]

# odds ratios
hla <- odds.ratio(hla)

dosage <- fread("/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/hla_imputation/hla_dosage.txt", data.table=F)
dosage.hsp <- dosage[rows.hsp,2:ncol(dosage)]
dosage.controls <- dosage[rows.controls,2:ncol(dosage)]

# allele frequencies
hla <- allele.freq(dosage.hsp, dosage.controls, hla)

# locations
location.chr6 <- data.table("name" = hla$name)
location.chr6$location <- c("32,637,396-32,654,774", "32,659,464-32,666,689", "32,578,769-32,589,836")
hla <- merge(hla, location.chr6, by = "name")

# gene names
for (i in 1:length(hla$name)) {
  split <- str_split_fixed(hla$name[i], "_", n=2)
  hla$name[i] <- paste("HLA-", toupper(split[1]), "*", split[2], sep="")
}

hla.final <- hla[,c(1,12,2,7,8,9, 11, 10)]

# rounding results
hla.final <- round.results(hla.final)

hla.final <- hla.final[,c(1,2,3,9,8,7)]
colnames(hla.final) <- c("HLA allele", "position", "P", "OR[CI 95%]", "freq.hsp", "freq.controls")
hla.final <- hla.final[order(hla.final$P),]

write.table(hla.final, "./results/images/koko_datan/HLA_tables/hla_hsp_table.txt", sep=" ", quote=F, row.names=F, col.names=T)
#####################################################################################################
# hla allele table for HSP (p < 0.05)
hla <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/hla_hsp_results.txt', data.table=F)

hla <- hla[order(hla[,2]),]
hla <- hla[hla$p.value < 5 * 10 ** -2,]
hla <- hla[!(is.na(hla$name)),]

# odds ratios
hla <- odds.ratio(hla)

dosage <- fread("/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/hla_imputation/hla_dosage.txt", data.table=F)
dosage.hsp <- dosage[rows.hsp,2:ncol(dosage)]
dosage.controls <- dosage[rows.controls,2:ncol(dosage)]

# allele frequencies
hla <- allele.freq(dosage.hsp, dosage.controls, hla)

# gene names
for (i in 1:length(hla$name)) {
  split <- str_split_fixed(hla$name[i], "_", n=2)
  hla$name[i] <- paste("HLA-", toupper(split[1]), "*", split[2], sep="")
}

hla.final <- hla[,c(1,2,7,8,9,11,10)]

# rounding results
hla.final <- round.results(hla.final)

hla.final <- hla.final[,c(1,2,8,7,6)]
colnames(hla.final) <- c("HLA allele", "P", "OR[CI 95%]", "freq.hsp", "freq.controls")
hla.final <- hla.final[order(hla.final$P),]

write.table(hla.final, "./results/images/koko_datan/HLA_tables/hla_hsp_table_0.05.txt", sep=" ", quote=F, row.names=F, col.names=T)
#################################################################################################################
# hla allele table for IBD (p < 0.05)
hla <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/hla_ibd_results.txt', data.table=F)

hla <- hla[order(hla[,2]),]
hla <- hla[hla$p.value < 5 * 10 ** -2,]
hla <- hla[!(is.na(hla$name)),]

# odds ratios
hla <- odds.ratio(hla)

dosage <- fread("/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/hla_imputation/hla_dosage.txt", data.table=F)
dosage.ibd <- dosage[rows.ibd,2:ncol(dosage)]
dosage.controls <- dosage[rows.controls,2:ncol(dosage)]

# allele frequencies
hla <- allele.freq(dosage.ibd, dosage.controls, hla)

# gene names
for (i in 1:length(hla$name)) {
  split <- str_split_fixed(hla$name[i], "_", n=2)
  hla$name[i] <- paste("HLA-", toupper(split[1]), "*", split[2], sep="")
}

hla.final <- hla[,c(1,2,7,8,9, 11, 10)]

# rounding results
hla.final <- round.results(hla.final)

hla.final <- hla.final[,c(1,2,8,7,6)]
colnames(hla.final) <- c("HLA allele", "P", "OR[CI 95%]", "freq.ibd", "freq.controls")
hla.final <- hla.final[order(hla.final$P),]

write.table(hla.final, "./results/images/koko_datan/HLA_tables/hla_ibd_table_0.05.txt", sep=" ", quote=F, row.names=F, col.names=T)
#####################################################################################################
# hla prot seq table

hla.prot <- fread(paste("/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/hla_prot_hsp_results.txt", sep=""), data.table=F)
hla.prot <- hla.prot[order(hla.prot[,2]),]
hla.prot <- hla.prot[hla.prot$p.value < 5 * 10 ** -8,]

# odds ratios
hla.prot <- odds.ratio(hla.prot)

# genes, amino acids and amino acid positions
change.names <- hla.prot[,1]
genes <- c() 
positions <- c()
aa <- c()

for (i in 1:length(change.names)) {
  split <- str_split_fixed(change.names[i], "_", n=2)
  genes <- c(genes, split[1])
  
  positions <- c(positions, substr(split[2], 1, 2))
  
  aa <- c(aa, substr(split[2], 3, 3))
}

hla.prot$gene <- genes
hla.prot$position <- positions
hla.prot$amino.acid <- aa

# frequencies:
hla.prot <- allele.freq(dosage.hsp, dosage.controls, hla.prot)
# prot dosage:
dosage <- fread(paste("/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/hla_imputation/hla_prot_dosage.txt", sep=""), data.table=F)

dosage.hsp <- dosage[rows.hsp,2:ncol(dosage)]
dosage.controls <- dosage[rows.controls,2:ncol(dosage)]

prot.final <- hla.prot[,c(10,11,12,2,7,8,9, 13, 14)]

prot.final <- round.results(prot.final)

prot.final <- prot.final[,c(1,2,3,4,10,8,9)]
colnames(prot.final) <- c("HLA molecule", "position", "amino.acid", "P", "OR[CI 95%]", "freq.hsp", "freq.controls")
prot.final <- prot.final[order(prot.final$P),]

write.table(prot.final, "./results/images/koko_datan/HLA_tables/prot_hsp_table.txt", sep=" ", quote=F, row.names=F, col.names=T)

#################################################################################################################
# add results by others (not included in finished thesis as a table)

# variants - signals
others <- fread('./data/others.results.txt', data.table=F)
colnames(others) <- c("SNP", "position", "minor allele", "P", "OR[CI 95%]")

# unify posotion format
split.pos.1 <- substr(gwas.results.final$position, 1, 2)
split.pos.2 <- substr(gwas.results.final$position, 3, 5)
split.pos.3 <- substr(gwas.results.final$position, 6,8)

gwas.results.final$position <- paste(split.pos.1, split.pos.2, split.pos.3, sep=".")

add.col <- rep(NA, 42)
others$freq.hsp <- add.col
others$freq.controls <- add.col
gwas.results.all <- rbind(gwas.results.final, others)

# all results
tbody.style = tbody_style(color = "black", fill = c(rep("grey95",length.out=48), rep("grey90",length.out=42)), hjust=0, x=0.1)
ggtexttable(gwas.results.all, rows = NULL, theme = ttheme(tbody.style = tbody.style))
ggsave("./results/images/koko_datan/HLA_tables/hsp_table_all.jpeg", device="jpeg", height= 70, width= 20, units="cm")

# our results
tbody.style = tbody_style(color = "black", fill = c("grey95", "grey90"), hjust=0, x=0.1)
ggtexttable(gwas.results.final, rows = NULL, theme = ttheme(tbody.style = tbody.style))
ggsave("./results/images/koko_datan/HLA_tables/hsp_table.jpeg", device="jpeg", height= 40, width= 20, units="cm")

write.table(gwas.results.all, "./results/images/koko_datan/HLA_tables/hsp_table_ALL.txt", sep=" ", quote=F, row.names=F, col.names=T)

# HLA alleles

tbody.style = tbody_style(color = "black", fill = c("grey95", "grey90"), hjust=0, x=0.1)
ggtexttable(hla.final, rows = NULL, theme = ttheme(tbody.style = tbody.style))
ggsave("./results/images/koko_datan/HLA_tables/hla_hsp_table.jpeg", device="jpeg", height= 10, width= 30, units="cm")

# HLA prot seq - potential  signals
molecule <- c(rep("DRB1", 2))
position <- c(13,11)
missing <- c(rep("", 2))
p <- c(6.67E-05, 1.88E-05)

others.prot <- data.table(molecule, position, missing, p, missing, missing, missing)
colnames(others.prot) <- colnames(prot.final)
prot.all <- rbind(prot.final, others.prot)

tbody.style = tbody_style(color = "black", fill = c(rep("grey95",length.out=5), rep("grey90",length.out=2)), hjust=0, x=0.1)
ggtexttable(prot.all, rows = NULL, theme = ttheme(tbody.style = tbody.style))

ggsave("./results/images/koko_datan/HLA_tables/prot_table_all.jpeg", device="jpeg", height= 10, width= 20, units="cm")

write.table(prot.all, "./results/images/koko_datan/HLA_tables/prot_hsp_table_ALL.txt", sep=" ", quote=F, row.names=F, col.names=T)

  


