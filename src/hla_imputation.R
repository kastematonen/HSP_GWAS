library(data.table)
library(tidyverse)
library(HIBAG)
library(parallel)
library(gridExtra)
library(ggpubr)
library(qqman)
library(SAIGE)
library(SPAtest)
library(HIBAG)

# load PLINK files (from which duplicate variants have been removed)  
hlageno <- hlaBED2Geno(bed.fn="/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/hla_imputation/hla_data.bed", 
                       fam.fn="/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/hla_imputation/hla_data.fam", 
                       bim.fn="/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/hla_imputation/hla_data.bim", assembly='hg38')

# load models for imputation
model.a    <- hlaModelFromObj(get(load('../HLA-imputation/models/hg38/Fin_hg38_model_A.RData'))[['A']])
model.b    <- hlaModelFromObj(get(load('../HLA-imputation/models/hg38/Fin_hg38_model_B.RData'))[['B']])
model.c    <- hlaModelFromObj(get(load('../HLA-imputation/models/hg38/Fin_hg38_model_C.RData'))[['C']])
model.drb1 <- hlaModelFromObj(get(load('../HLA-imputation/models/hg38/Fin_hg38_model_DRB1.RData'))[['DRB1']])
model.dqa1 <- hlaModelFromObj(get(load('../HLA-imputation/models/hg38/Fin_hg38_model_DQA1.RData'))[['DQA1']])
model.dqb1 <- hlaModelFromObj(get(load('../HLA-imputation/models/hg38/Fin_hg38_model_DQB1.RData'))[['DQB1']])
model.dpb1 <- hlaModelFromObj(get(load('../HLA-imputation/models/hg38/Fin_hg38_model_DPB1.RData'))[['DPB1']])

# impute HLA alleles
hlageno.a    <- predict(model.a,    hlaGenoSubset(hlageno), type="response+prob", match.type="Position")
hlageno.b    <- predict(model.b,    hlaGenoSubset(hlageno), type="response+prob", match.type="Position")
hlageno.c    <- predict(model.c,    hlaGenoSubset(hlageno), type="response+prob", match.type="Position")
hlageno.drb1 <- predict(model.drb1, hlaGenoSubset(hlageno), type="response+prob", match.type="Position")
hlageno.dqa1 <- predict(model.dqa1, hlaGenoSubset(hlageno), type="response+prob", match.type="Position")
hlageno.dqb1 <- predict(model.dqb1, hlaGenoSubset(hlageno), type="response+prob", match.type="Position")
hlageno.dpb1 <- predict(model.dpb1, hlaGenoSubset(hlageno), type="response+prob", match.type="Position")

hlageno.a <- hlageno.a$value
hlageno.b <- hlageno.b$value
hlageno.c <- hlageno.c$value
hlageno.drb1 <- hlageno.drb1$value
hlageno.dqa1 <- hlageno.dqa1$value
hlageno.dqb1 <- hlageno.dqb1$value
hlageno.dpb1 <- hlageno.dpb1$value

##############################################################################################################
# create a dosage file for spatest
dosage.table.alleles <- data.frame(samplle.id = hlageno.a$sample.id)

# function to count the number of an allele in individuals
fill.table <- function(i, dosage.table, hlageno.table, allele.list){
  
  allele1.count <- rep(NA, length(dosage.table$samplle.id))
  allele2.count <- rep(NA, length(dosage.table$samplle.id))
  
  allele1.count[hlageno.table$allele1 %in% allele.list[i]] <- 1
  allele1.count[!(hlageno.table$allele1 %in% allele.list[i])] <- 0
  allele2.count[hlageno.table$allele2 %in% allele.list[i]] <- 1
  allele2.count[!(hlageno.table$allele2 %in% allele.list[i])] <- 0
  
  allele.count.total <- allele1.count + allele2.count
  return(allele.count.total)
}
# unique alleles for each locus
alleles.a <- unique(c(hlageno.a$allele1, hlageno.a$allele2))
alleles.b <- unique(c(hlageno.b$allele1, hlageno.b$allele2))
alleles.c <- unique(c(hlageno.c$allele1, hlageno.c$allele2))
alleles.drb1 <- unique(c(hlageno.drb1$allele1, hlageno.drb1$allele2))
alleles.dqa1 <- unique(c(hlageno.dqa1$allele1, hlageno.dqa1$allele2))
alleles.dqb1 <- unique(c(hlageno.dqb1$allele1, hlageno.dqb1$allele2))
alleles.dpb1 <- unique(c(hlageno.dpb1$allele1, hlageno.dpb1$allele2))

# filling the dosage table 
# some genes have same allele names so the table is filled in parts + a prefix added to the allele names
for (i in 1:length(alleles.a)) {
  dosage.table.alleles[, paste("a_", alleles.a[i], sep="")] <- fill.table(i, dosage.table.alleles, hlageno.a, alleles.a)
  print(i)
}
for (i in 1:length(alleles.b)) {
  dosage.table.alleles[, paste("b_", alleles.b[i], sep="")] <- fill.table(i, dosage.table.alleles, hlageno.b, alleles.b)
  print(i)
}
for (i in 1:length(alleles.c)) {
  dosage.table.alleles[, paste("c_", alleles.c[i], sep="")] <- fill.table(i, dosage.table.alleles, hlageno.c, alleles.c)
  print(i)
}
for (i in 1:length(alleles.drb1)) {
  dosage.table.alleles[, paste("drb1_", alleles.drb1[i], sep="")] <- fill.table(i, dosage.table.alleles, hlageno.drb1, alleles.drb1)
  print(i)
}
for (i in 1:length(alleles.dqa1)) {
  dosage.table.alleles[, paste("dqa1_", alleles.dqa1[i], sep="")] <- fill.table(i, dosage.table.alleles, hlageno.dqa1, alleles.dqa1)
  print(i)
}
for (i in 1:length(alleles.dqb1)) {
  dosage.table.alleles[, paste("dqb1_", alleles.dqb1[i], sep="")] <- fill.table(i, dosage.table.alleles, hlageno.dqb1, alleles.dqb1)
  print(i)
}
for (i in 1:length(alleles.dpb1)) {
  dosage.table.alleles[, paste("dpb1_", alleles.dpb1[i], sep="")] <- fill.table(i, dosage.table.alleles, hlageno.dpb1, alleles.dpb1)
  print(i)
}

##########################################################################################################

# association analysis with spatest: save dosage file

write.table(dosage.table.alleles, "/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/hla_imputation/hla_dosage.txt", sep=" ", quote=F, row.names=F)

# association analysis in ./src/gwas.R

############################################################################################################
# protein sequences for imputed HLA alleles

seq.a <- hlaConvSequence(hla=alleles.a, locus="A",
                               method=c("protein"),
                               code="P.code.merge",
                               region="auto")
seq.b <- hlaConvSequence(hla=alleles.b, locus="B", method=c("protein"),
                               code="P.code.merge",
                               region="auto")
seq.c <- hlaConvSequence(hla=alleles.c, locus="C", method=c("protein"),
                               code="P.code.merge",
                               region="auto")
seq.drb1 <- hlaConvSequence(hla=alleles.drb1, locus="DRB1", method=c("protein"),
                               code="P.code.merge",
                               region="auto")
seq.dqa1 <- hlaConvSequence(hla=alleles.dqa1, locus="DQA1", method=c("protein"),
                                  code="P.code.merge",
                                  region="auto")
seq.dqb1 <- hlaConvSequence(hla=alleles.dqb1, locus="DQB1", method=c("protein"),
                                  code="P.code.merge",
                                  region="auto")
seq.dpb1 <- hlaConvSequence(hla=alleles.dpb1, locus="DPB1", method=c("protein"),
                                  code="P.code.merge",
                                  region="auto")

#################################################################################################################
# dosage file for amino acids 

hla.dosage <- dosage.table.alleles
# hla.dosage <- fread('./data/cleaned_genotypes/hla_dosage.raw', data.table=F)

fill.table.aa <- function(table, sequences, count){
  for (i in 1:length(sequences)) {
    seq <- unlist(str_split(sequences[i], ""))
    numbers <- c(1:length(seq))
    
    # reference "-" will not be included:
    leave.these <- seq != "-"
    seq <- seq[leave.these]
    numbers <- numbers[leave.these]
    
    for (j in 1:length(seq)) {
      if(!(paste(numbers[j], seq[j], sep="") %in% colnames(table))){
        table[,paste(numbers[j], seq[j], sep="")] <- hla.dosage[, i+count]
      } else {
        table[,paste(numbers[j], seq[j], sep="")] <- table[,paste(numbers[j], seq[j], sep="")] + hla.dosage[, i+count]
      }
    }
  }
  return(table)
}

allele.count.a <- 1
allele.count.b <- allele.count.a + length(alleles.a)
allele.count.c <- allele.count.b + length(alleles.b)
allele.count.drb1 <- allele.count.c + length(alleles.c)
allele.count.dqa1 <- allele.count.drb1 + length(alleles.drb1)
allele.count.dqb1 <- allele.count.dqa1 + length(alleles.dqa1)
allele.count.dpb1 <- allele.count.dqb1 + length(alleles.dqb1)

# filling a table per locus
dosage.table.a <- data.frame(samplle.id = hlageno.a$sample.id)
dosage.table.a <- fill.table.aa(dosage.table.a, seq.a, allele.count.a)
colnames(dosage.table.a) <- paste("A", colnames(dosage.table.a), sep = "_")

dosage.table.b <- data.frame(samplle.id = hlageno.b$sample.id)
dosage.table.b <- fill.table.aa(dosage.table.b, seq.b, allele.count.b)
colnames(dosage.table.b) <- paste("B", colnames(dosage.table.b), sep = "_")

dosage.table.c <- data.frame(samplle.id = hlageno.b$sample.id)
dosage.table.c <- fill.table.aa(dosage.table.c, seq.c, allele.count.c)
colnames(dosage.table.c) <- paste("C", colnames(dosage.table.c), sep = "_")

dosage.table.dpb1 <- data.frame(samplle.id = hlageno.dpb1$sample.id)
dosage.table.dpb1 <- fill.table.aa(dosage.table.dpb1, seq.dpb1, allele.count.dpb1)
colnames(dosage.table.dpb1) <- paste("DPB1", colnames(dosage.table.dpb1), sep = "_")

dosage.table.dqa1 <- data.frame(samplle.id = hlageno.dqa1$sample.id)
dosage.table.dqa1 <- fill.table.aa(dosage.table.dqa1, seq.dqa1, allele.count.dqa1)
colnames(dosage.table.dqa1) <- paste("DQA1", colnames(dosage.table.dqa1), sep = "_")

dosage.table.dqb1 <- data.frame(samplle.id = hlageno.dqb1$sample.id)
dosage.table.dqb1 <- fill.table.aa(dosage.table.dqb1, seq.dqb1, allele.count.dqb1)
colnames(dosage.table.dqb1) <- paste("DQB1", colnames(dosage.table.dqb1), sep = "_")

dosage.table.drb1 <- data.frame(samplle.id = hlageno.drb1$sample.id)
dosage.table.drb1 <- fill.table.aa(dosage.table.drb1, seq.drb1, allele.count.drb1)
colnames(dosage.table.drb1) <- paste("DRB1", colnames(dosage.table.drb1), sep = "_")

# combining individual tables
dosage.table.all <- cbind(dosage.table.a, dosage.table.b[,2:ncol(dosage.table.b)], dosage.table.c[,2:ncol(dosage.table.c)], dosage.table.dpb1[2:ncol(dosage.table.dpb1)], dosage.table.dqa1[,2:ncol(dosage.table.dqa1)], dosage.table.dqb1[,2:ncol(dosage.table.dqb1)], dosage.table.drb1[2:ncol(dosage.table.drb1)])

# remove columns made from NA sequences
dosage.table.all <- dosage.table.all[!grepl("NA", colnames(dosage.table.all), fixed = TRUE)]

# excluding columns whose amino acid substitutions are found in all samples
# is an amino acid substitution in all samples:
in.all <- sapply(dosage.table.all[,2:ncol(dosage.table.all)], min) > 0 # 2
# include sample id in list
in.all <- c(F, in.all)
# exclude columns (amino acid substitutions) 
dosage.table.all <- dosage.table.all[,!in.all]

#########################################################################################################
# saving dosage file for association test 
# -> ./src/gwas.R

write.table(dosage.table.all, "/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/hla_imputation/hla_prot_dosage.txt", sep=" ", quote=F, row.names=F, col.names=T)
