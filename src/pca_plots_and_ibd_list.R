library(data.table)
library(tidyverse)
library(parallel)
library(gridExtra)
library(ggpubr)

pca <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/pca/pca_ld_pruned_cleaned_merged_1_23.eigenvec', data.table=F)

hsp.samples <- fread('./data/HSP.sample.list', data.table=F)[, 2]
hsct.samples <- fread('./data/hsct_register_genotypes/VPU_ILLUMINA_AUG_2018.fam', data.table=F)[, 2]
ibd.samples <- fread('./data/IBD.sample.list', data.table=F)[, 2]
blood.donor.samples <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_controls/merged_1_23.fam', data.table=F)[, 2]

sample.group <- rep(NA, length(pca$IID))
sample.group[pca$IID %in% hsp.samples] <- 'HSP' 
sample.group[pca$IID %in% hsct.samples] <- 'CONTROL HSCT' 
sample.group[pca$IID %in% ibd.samples] <- 'IBD' 
sample.group[pca$IID %in% blood.donor.samples] <- 'CONTROL BLOOD DONOR' 

#############################################################################################################################################
# plot of the principal components
pca$sample.group <- sample.group
# ordering pca table by sample group -> HSP and IBD are drawn last and thus more visible  on the plots
pca <- pca[order(sample.group),]

# columns 3-7 of the PCA table -> principal components 1-5 to be plotted
column.combinations <- combn(c(3:7), 2)

plot.function <- function(i){
  x.axis <- column.combinations[1,i]
  y.axis <- column.combinations[2,i]
  x.label <- paste("Principal Component ", column.combinations[1,i]-2)
  y.label <- paste("Principal Component ", column.combinations[2,i]-2)
  
  plot <- ggplot(pca, aes(x=pca[,x.axis], y= pca[,y.axis])) + geom_point(aes(colour= sample.group), size = 2) + labs(x = x.label, y= y.label, colour= "Sample") + theme_minimal() + theme(axis.text.y = element_text(angle = 90, size = 6), axis.text.x = element_text(size = 6), axis.title=element_text(size=13)) + scale_color_manual(values=c("#33638DFF", "#3CBB75FF", "#440154FF", "#FDE725FF"))
  
  return(plot)
}

plot.list <- lapply(1:ncol(column.combinations), FUN= plot.function)
ggarrange(plotlist=plot.list, common.legend=T, ncol=2, nrow=5, legend="bottom", labels = paste(LETTERS[1:length(plot.list)], ")", sep=""))

ggsave("./results/images/koko_datan/pca_outliers.jpeg", device="jpeg", height= 40, width= 30, units="cm")

#############################################################################################################################################
# IBD analysis
relatives <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/pca/ibd_ld_pruned_cleaned_merged_1_23.genome', data.table=F)
# first degree relatives for removal: PI_HAT value cut-off 0.4
first_degree <- relatives[relatives$PI_HAT >= 0.375,]

first_subject <- first_degree[,1:2]
second_subject <- first_degree[,3:4]
colnames(second_subject) <- colnames(first_subject)

all_related <- rbind(first_subject, second_subject)
write.table(all_related, file="/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/pca/relatives.txt", sep="\t", quote=F, row.names=F, col.names=F)

# missingness for individuals with ./src/pca_and_ibd.sh
missingness <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/pca/relatives_missingness.imiss', data.table=F)

individuals_for_removal <- function(x){
  individual_1 <- missingness[missingness$IID == first_degree[x,2],]
  individual_2 <- missingness[missingness$IID == first_degree[x,4],]
  if(individual_2$F_MISS > individual_1$F_MISS){
    return(individual_2$IID)
  } else {
    return(individual_1$IID)
  }
}
remove_these <- sapply(1:nrow(first_degree), individuals_for_removal)
remove_these <- unique(remove_these)

rows <- pca$IID %in% remove_these
list_for_removal <- pca[rows,1:2]

write.table(list_for_removal, file="/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/pca/remove_relatives.txt", sep="\t", quote=F, row.names=F, col.names=F)

############################################################################################################################################
# barplot of the eigenvalues
# from final PCA:
variance <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/pca/final_pca_ld_pruned_cleaned_merged_1_23.eigenval', data.table=F, header=F)

variance$plotted <- 100 * (variance$V1/sum(variance$V1))
variance$PC <- 1:20

ggplot(variance, aes(y=plotted, x=PC)) + geom_bar(stat="identity") + labs(x = "Principal Component", y= "% Variance explained") + scale_x_continuous(breaks = 1:20, labels=1:20) + theme_classic2() + theme(axis.line.x = element_line(colour = "white"), axis.ticks.x = element_blank(), axis.text = element_text(size=12), axis.title=element_text(size=14))
ggsave("./results/images/koko_datan/pca_barplot.jpeg", device="jpeg", height= 20, width= 30, units="cm")




