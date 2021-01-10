#!/bin/env Rscript --no-save

# Required packages
library(data.table) # For fast fread()

# Input variables
args <- commandArgs(TRUE)
indataset <- args[1]
infile <- args[2]
ref_dataset <- args[3]
refname <- args[4]
af_diff_limit <- as.numeric(args[5])

# Read in the frequency files
chip <- fread(infile, header = T)
ref_data <- fread(ref_dataset, header = T)

# Take an intersection of the reference and chip data 
# based on SNP column (in format CHR_POS_REF_ALT)
isec <- merge(ref_data, chip, by = "SNP")

# Exclude if AF is not within the range
exclude <- !abs(isec$AF.x - isec$AF.y) < af_diff_limit

# Non-ref data variants
nonref <- chip[!(chip$SNP) %in% (isec$SNP)]

# Save the plot as jpg
png(paste0(indataset, "_", refname, "_AF.png"), 
    width = 600, height = 600)
# Plot first all and then excludable variants
plot(isec$AF.x, isec$AF.y, col=1, pch=20,
    main=paste0(indataset, " vs. ", refname, " AF"),
    xlab="Reference data AF",
    ylab="Chip data AF")
points(isec[exclude]$AF.x, isec[exclude]$AF.y, 
    col=2, pch=20)
# Draw a legend
legend("topleft", legend=c(
    paste0("Concordant AF, n = ", nrow(isec[!exclude])),
    paste0("High AF difference, n = ", nrow(isec[exclude])),
    paste0("Non-ref variants, n= ", nrow(nonref))),
    col=c("black", "red", "white"), pch=20, cex=1.2)
dev.off()

# Store the high AF difference variants
output <- rbind(c("chip variants", nrow(chip)),
                c("intersection", nrow(isec)),
                c("intersection/chip variants", 
                     (nrow(isec))/nrow(chip)),
                c("non-ref variants", 
                      nrow(chip) - nrow(isec)),
                c("high AF difference", 
                      nrow(isec[exclude])))

write.table(output, 
    paste0(indataset, "_", refname, "_comparison.txt"),
    quote=F, row.names=F, col.names=F, sep = "\t")
