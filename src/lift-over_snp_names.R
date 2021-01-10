library(data.table)
library(tidyverse)

# strand file used in genotype lift-over, SNPs matched to data by SNP names
strand.file <- fread('./GSA-24v2-0_A1-b38.strand', data.table=F)

# proportion of SNP names in common between data and strand file before unifying SNP names
data <- fread('./data/cleaned_genotypes/merged_cleaned_no_duplicates.bim', data.table=F)
shared <- data$V2 %in% strand.file$V1
sum(shared) / length(shared) # 0.5496903

# changing SNP names in strand file from "GSA-rs123" to "GSA-RS123"
change.names <- strand.file[strand.file$V1 %like% "GSA",]
change.names$V1 <- toupper(change.names$V1)
strand.file[rownames(change.names),1] <- change.names$V1

# changing SNP names in strand file from "exm123" to "EXM123"
change.names <- strand.file[strand.file$V1 %like% "exm",]
change.names$V1 <- toupper(change.names$V1)
strand.file[rownames(change.names),1] <- change.names$V1

# changing SNP names in strand file from "seq-rs123" to "SEQ-RS123"
change.names <- strand.file[strand.file$V1 %like% "seq-rs",]
change.names$V1 <- toupper(change.names$V1)
strand.file[rownames(change.names),1] <- change.names$V1

# proportion of SNP names in common between data and strand file after unifying SNP names
shared <- data$V2 %in% strand.file$V1
sum(shared) / length(shared) # 0.9736173 (was 0.5496903 before)

missing <- data[!shared,]
length(missing$V1) # 12204 SNP names in data not matchning names in strand file

write.table(strand.file, file="./GSA-24v2-0_A1-b38_mod.strand", sep="\t", quote=F, row.names=F, col.names=F)


