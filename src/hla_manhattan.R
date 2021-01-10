# HLA-geenien Manhattan-kuvaaja
library(ggrepel)

results.hsp <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/hla_hsp_results.txt', data.table=F)[,1:2]
results.ibd <- fread('/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/gwas/hla_ibd_results.txt', data.table=F)[,1:2]

gene <- c()
genes <- strsplit(results.hsp$name, "_")
for (i in 1:length(genes)) {
  element <- genes[i]
  element <- element[[1]][1]
  gene <- c(gene, element)
}
gene <- toupper(gene)
results.hsp$gene <- gene
results.hsp$log.p <- -log10(results.hsp$p.value)
results.hsp$name <- toupper(results.hsp$name)
results.hsp$name <- gsub("_", "*", results.hsp$name)

gene <- c()
genes <- strsplit(results.ibd$name, "_")
for (i in 1:length(genes)) {
  element <- genes[i]
  element <- element[[1]][1]
  gene <- c(gene, element)
}
gene <- toupper(gene)
results.ibd$gene <- gene
results.ibd$log.p <- -log10(results.ibd$p.value)

# plot for HSP:
ggplot(results.hsp, aes(x=gene, y=log.p)) +
  geom_point(size=1, alpha=0.7) + 
  xlab("HLA-geeni") +
  ylab("-log10(P-arvo)") +
  theme_classic2() +
  geom_hline(yintercept=5, size=0.5, col="blue") +
  geom_hline(yintercept=-log10(5e-8), size=0.5, col="red") +
  theme(axis.text=element_text(size=10)) +
  theme(axis.title=element_text(size=11)) +
  
  geom_label_repel(aes(gene, log.p, label=name), data=arrange(results.hsp, desc(log.p))[1:3, ], segment.size=0.3, nudge_x=0.3, nudge_y=-0.3, size=3.2, min.segment.length=0)

ggsave("./tmp/kuvia_suomeksi/hla_manhattan_hsp.jpeg", device="jpeg", height= 10, width= 15, units="cm")
  
# plot for IBD:
ggplot(results.ibd, aes(x=gene, y=log.p)) +
  geom_point(size=1, alpha=0.7) + 
  xlab("HLA-geeni") +
  ylab("-log10(P-arvo)") +
  theme_classic2() +
  geom_hline(yintercept=5, size=0.5, col="blue") +
  geom_hline(yintercept=-log10(5e-8), size=0.5, col="red") +
  theme(axis.text=element_text(size=10)) +
  theme(axis.title=element_text(size=11))
  
ggsave("./tmp/kuvia_suomeksi/hla_manhattan_ibd.jpeg", device="jpeg", height= 10, width= 15, units="cm")


