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

variants <- fread('./data/cleaned_genotypes/merged_cleaned_data.bim', data.table=F)

# variants with no unique name or location:
names <- variants[duplicated(variants$V2),] # 0
locations <- variants[duplicated(variants[c("V1", "V4")]),] # 1421 

write.table(locations$V2, file="./data/cleaned_genotypes/duplicates_cleaned_merged.txt", sep=" ", quote=F, row.names=F, col.names=F)
