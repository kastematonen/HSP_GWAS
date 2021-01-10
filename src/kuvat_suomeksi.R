# Kuvat suomeksi
# koodit pääasiassa samat kuin alkuperäisissä toteutuksissa, mutta niiden tekstit vaihdettu suomeksi
# tähän tiedostoon poimittu vain kuvien kohdat muista tiedostoista 

sample.group <- rep(NA, length(pca$IID))
sample.group[pca$IID %in% hsp.samples] <- 'HSP' 
sample.group[pca$IID %in% hsct.samples] <- 'HSCT' 
sample.group[pca$IID %in% ibd.samples] <- 'IBD' 
sample.group[pca$IID %in% blood.donor.samples] <- 'Verenluovuttaja' 

############################################################################################################################################
pca <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/pca/pca_ld_pruned_cleaned_merged_1_23.eigenvec', data.table=F)

hsp.samples <- fread('./data/HSP.sample.list', data.table=F)[, 2]
hsct.samples <- fread('./data/hsct_register_genotypes/VPU_ILLUMINA_AUG_2018.fam', data.table=F)[, 2]
ibd.samples <- fread('./data/IBD.sample.list', data.table=F)[, 2]
blood.donor.samples <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/merged_controls/merged_1_23.fam', data.table=F)[, 2]

sample.group.1 <- rep(NA, length(pca$IID))
sample.group.1[pca$IID %in% hsp.samples] <- 'HSP' 
sample.group.1[pca$IID %in% hsct.samples] <- 'CONTROL HSCT' 
sample.group.1[pca$IID %in% ibd.samples] <- 'IBD' 
sample.group.1[pca$IID %in% blood.donor.samples] <- 'CONTROL BLOOD DONOR' 

#############################################################################################################################################
# plot of the principal components
pca$sample.group.1 <- sample.group.1

sample.group <- rep(NA, length(pca$IID))
sample.group[pca$IID %in% hsp.samples] <- 'HSP' 
sample.group[pca$IID %in% hsct.samples] <- 'HSCT' 
sample.group[pca$IID %in% ibd.samples] <- 'IBD' 
sample.group[pca$IID %in% blood.donor.samples] <- 'Verenluovuttaja' 
pca$sample.group <- sample.group
# ordering pca table by sample group -> HSP and IBD are drawn last and thus more visible  on the plots
pca <- pca[order(sample.group.1),]

# columns 3-7 of the PCA table -> principal components 1-5 to be plotted
column.combinations <- combn(c(3:7), 2)


# PCA scatter plot:
plot.function <- function(i){
  x.axis <- column.combinations[1,i]
  y.axis <- column.combinations[2,i]
  x.label <- paste("Pääkomponentti ", column.combinations[1,i]-2)
  y.label <- paste("Pääkomponentti ", column.combinations[2,i]-2)
  
  plot <- ggplot(pca, aes(x=pca[,x.axis], y= pca[,y.axis])) + geom_point(aes(colour= sample.group), size = 1.5) + labs(x = x.label, y= y.label, colour= "Näyte") + theme_minimal() + theme(axis.text.y = element_text(angle = 90, size = 6), axis.text.x = element_text(size = 6), axis.title=element_text(size=13)) + scale_color_manual(values=c("#3CBB75FF", "#440154FF", "#FDE725FF", "#33638DFF"))
  
  return(plot)
}

plot.list <- lapply(1:ncol(column.combinations), FUN= plot.function)
ggarrange(plotlist=plot.list, common.legend=T, ncol=2, nrow=5, legend="bottom", labels = paste(LETTERS[1:length(plot.list)], ")", sep=""))
ggsave("./tmp/kuvia_suomeksi/pca_outliers_suomi.jpeg", device="jpeg", height= 40, width= 30, units="cm")


### alkuperäinen toteutus englanniksi (./src/pca_plots_and_ibd_list.R)
# plot.function <- function(i){
#   x.axis <- column.combinations[1,i]
#   y.axis <- column.combinations[2,i]
#   header <- paste("PC", column.combinations[1,i]-2, " ja ", "PC", column.combinations[2,i]-2)
#   x.label <- paste("Principal component ", column.combinations[1,i]-2)
#   y.label <- paste("Principal component ", column.combinations[2,i]-2)
#   
#   plot <- ggplot(pca, aes(x=pca[,x.axis], y= pca[,y.axis])) + geom_point(aes(colour= sample.group)) + labs(title= header, x = x.label, y= y.label, colour= "näyte") + theme((axis.text.y = element_text(angle = 90))) 
#   return(plot)
# }


# PCA barplot:
ggplot(variance, aes(y=plotted, x=PC)) + geom_bar(stat="identity") + labs(x = "Pääkomponentti", y= "Selittää varianssin suuruudesta (%)") + scale_x_continuous(breaks = 1:20, labels=1:20) + theme_classic2() + theme(axis.line.x = element_line(colour = "white"), axis.ticks.x = element_blank(), axis.text = element_text(size=12), axis.title=element_text(size=14))
ggsave("./tmp/kuvia_suomeksi/pca_barplot.jpeg", device="jpeg", height= 20, width= 30, units="cm")


# alkuperäinen toteutus (./src/pca_plots_and_ibd_list.R):
#ggplot(variance, aes(y=plotted, x=PC)) + geom_bar(stat="identity") + theme_light() + labs(title= "PCA:n pääkomponenttien suhteelliset suuruudet", x = "Pääkomponentti", y= "Suhteellinen suuruus") + scale_x_continuous(breaks = 1:20, labels=1:20)

#ggplot(variance, aes(y=plotted, x=PC)) + geom_bar(stat="identity") + labs(x = "Principal Component", y= "% Variance explained") + scale_x_continuous(breaks = 1:20, labels=1:20) + theme_classic2() + theme(axis.line.x = element_line(colour = "white"), axis.ticks.x = element_blank(), axis.text = element_text(size=12), axis.title=element_text(size=14))
#ggsave("./results/images/koko_datan/pca_barplot.jpeg", device="jpeg", height= 20, width= 30, units="cm")

#############################################################################################################################################
# platform-bias GWAS (tekstit englanniksi)

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
######################################################################################################
# joit-gwas, geneettinen samankaltaisuus

# alkuperäinen toteutus (./src/joint_gwas.R):
#all.plot <- ggplot(gwas.all.combined, aes(x = M, y = min.p)) + geom_line(aes(color = data.type)) + labs(title= "Enrichment of shared SNPs", x = "Length of SNP subset", y= "Enrichment, -log10(p)", tag = "A)", colour= "GWAS results") + theme_minimal() + scale_color_manual(values = c("blue", "brown"))

#hla.plot <- ggplot(hla.combined, aes(x = M, y = min.p)) + geom_line(aes(color = data.type)) + labs(title= "Enrichment of shared HLA alleles", x = "Length of HLA allele subset", y= "Enrichment, -log10(p)", tag = "B)") + theme_minimal() + scale_color_manual(values = c("blue", "brown"))

#prot.plot <- ggplot(prot.combined, aes(x = M, y = min.p)) + geom_line(aes(color = data.type)) + labs(title= "Enrichment of shared polymorphic amino acids", x = "Length of  amino acid subset", y= "Enrichment, -log10(p)", tag = "C)") + theme_minimal() + scale_color_manual(values = c("blue", "brown"))

#ggarrange(all.plot, hla.plot, prot.plot, ncol=1, nrow=3, common.legend=T, legend="bottom")
#ggsave("./results/images/koko_datan/joint_gwas_all.jpeg", device="jpeg", height= 35, width= 25, units="cm")



# tekstits suomeksi:

all.plot <- ggplot(gwas.all.combined, aes(x = M, y = min.p)) + geom_line(aes(color = data.type)) + labs(title= "Yhteisten SNP:iden rikastuminen", x = "SNP-osajoukon pituus", y= "Rikastuminen, -log10(p)", tag = "A)") + scale_color_manual(values = c("blue", "brown"), name="GWAS-tulokset", labels = c("Järjestetty", "Satunnainen järjestys"))  + theme_minimal()

hla.plot <- ggplot(hla.combined, aes(x = M, y = min.p)) + geom_line(aes(color = data.type)) + labs(title= "Yhteisten HLA-alleelien rikastuminen", x = "HLA-alleelien osajoukon pituus", y= "Rikastuminen, -log10(p)", tag = "B)") + theme_minimal() + scale_color_manual(values = c("blue", "brown"))

prot.plot <- ggplot(prot.combined, aes(x = M, y = min.p)) + geom_line(aes(color = data.type)) + labs(title= "Yhteisten aminohappojen rikastuminen", x = "Aminohappojen osajoukon pituus", y= "Rikastuminen, -log10(p)", tag = "C)") + theme_minimal() + scale_color_manual(values = c("blue", "brown"))

ggarrange(all.plot, hla.plot, prot.plot, ncol=1, nrow=3, common.legend=T, legend="bottom")
ggsave("./tmp/kuvia_suomeksi/joint_gwas_all.jpeg", device="jpeg", height= 35, width= 25, units="cm")

