#!/bin/env Rscript --no-save

# Required packages
library(data.table)

# Input variables
args <- commandArgs(TRUE)
indataset <- args[1]
infile <- args[2]
ref_dataset <- args[3]

# Read in the frequency files
chip <- fread(infile, header = T)
ref_data <- fread(ref_dataset, header =T)

# Take an intersection of the reference and chip data
# based on SNP column (in format CHR_POS_REF_ALT)
isec <- merge(ref_data, chip, by = "SNP")

# Exclude if AF is not within the range
exclude <- !abs(isec$AF.x - isec$AF.y) < 0.1
discrepant <- isec[exclude]

# Add column with a flipped AF value for the chip data
discrepant$AF_flip <- 1-discrepant$AF.y

# Test if the flipped AF matches to the reference data
flippable <- abs(discrepant$AF.x - discrepant$AF_flip) < 0.1

# Keep only ambiguous alleles
ambiguous <- (discrepant$REF.y %in% c("A", "T") 
             & discrepant$ALT.y %in% c("A", "T")) |
         (discrepant$REF.y %in% c("C", "G") 
             & discrepant$ALT.y %in% c("C", "G"))

# Generate output in format CHR \t POS
output <- strsplit(discrepant[flippable & ambiguous]$SNP, "_")
output <- do.call(rbind, output)[,1:2]

# Store the output
write.table(output, paste0(indataset, "_flippable_AF.txt"),
quote = F, row.names = F, col.names = F, sep = "\t")
