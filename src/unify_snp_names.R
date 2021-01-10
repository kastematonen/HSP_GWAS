library(data.table)
library(tidyverse)

# Changing SNP names from chromosome:location to SNP names in HSP, IBD and HSCT data
path <- "/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/bim/"
files <- list.files(path="/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/bim/", pattern="*bim")
files <- sort(files)
data <- fread('./data/lift-over/merged_cleaned_no_duplicates_hg38_mod_nonreference_flipped_ambiguous_flipped.bim', data.table=F)

chromosomes <- c(1,10,11,12,13,14,15,16,17,18,19,2,20,21,22,3,4,5,6,7,8,9,23)

# function for fetching correct SNP names (SNPs with the same position after lift-over -> first SNP name used)
fetch.names <- function(x){
  if (in.data[x]) {
    names.new <- names.data[names.data$V4 == chr[x,4],2]
    return(names.new[1])
  } else {
    return(names.orig[x])
  }
}
common.snps.all <- data.frame(V1 = integer(), V2 = character(), V3 = integer(), V4 = integer(), V5 = character(), V6 = character(), stringsAsFactors=F) 
for (i in 1:length(chromosomes)) {
  
  chr <- fread(paste(path,files[i],sep=""), data.table=F)
  names.orig <- chr$V2
  names.data <- data[data$V1 == chromosomes[i],]
  in.data <- chr$V4 %in% names.data$V4
  
  names.new.all <- sapply(1:length(in.data), fetch.names)
  
  chr$V2 <- unlist(names.new.all)
  write.table(chr, file=paste(path,"/uudet/",files[i],sep=""), sep="\t", quote=F, row.names=F, col.names=F)
  
  # assembling a list of common SNPs pre chromosome
  common.snps.chr <- inner_join(data, chr)
  # assembling a list of all common SNPs between HSP, IBD and HSCT data and blood donor data
  common.snps.all <- rbind(common.snps.all, common.snps.chr)
  
  write.table(common.snps.chr$V2, file=paste("./data/merge_controls/common_snps_chr",chromosomes[i],".txt",sep=""), sep="\t", quote=F, row.names=F, col.names=F)
  
}
write.table(common.snps.all$V2, file=paste("./data/merge_controls/common_snps_all.txt",sep=""), sep="\t", quote=F, row.names=F, col.names=F)


