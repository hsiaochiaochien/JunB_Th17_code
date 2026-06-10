setwd("/Users/chien/Library/CloudStorage/OneDrive-OIST/data/00_All_data_and_script_for_manuscript/06_figure6/figure6e/Rscript/data")
tpm <-read.table("salmon.merged.gene_tpm.txt",header =T)

# only select 6 to 72 hr no restimulate
tpm <- tpm[, c(2,9:(ncol(tpm)-6))]

tpm <-tpm[rowSums(tpm[ , c(2:ncol(tpm))])>0,]

#-----------------------------------------------------------------------
#plot RNAseq time serious
for (gene in c("Nfkb1","Nfkb2","Nr4a1","Nr4a2","Nr4a3","Runx1","Runx2","Runx3","Ets1","Ets2")) {
subgene <-tpm[grepl(gene,tpm$gene_name),2:ncol(tpm)]
subgene <-subgene[1,1:ncol(subgene)]
#if there are two sample using subgene[2,1:ncol(subgene)]
rownames(subgene) <- "tpm"

df<- data.frame(genotype = rep(rep(c("WT", "KO"), each=3),times =6),
               time = rep(c("06hr","12hr","24hr","36hr","48hr","72hr"), each=6),
               t(subgene[,1:ncol(subgene)]))
df$genotype <- factor(df$genotype, levels = c("WT", "KO"))

library(ggplot2)
p<-ggplot(df, aes(x=time, y=tpm, group=genotype)) +
  geom_point(size=2,aes(shape=genotype))+
  scale_shape_manual(values = c(16, 1)) +
  theme_classic()+
  geom_smooth(aes(x=time, y=tpm, color = genotype, fill=genotype)) + 
  labs(title = gene, x = "time", y = "TPM")+
  theme(plot.title = element_text(hjust = 0.5,size=20))+
  scale_color_manual(values=c("blue", "red"))+
  scale_fill_manual(values=c("blue", "red"))
p


savedir_png <- "/Users/chien/Library/CloudStorage/OneDrive-OIST/data/00_All_data_and_script_for_manuscript/06_figure6/figure6e/Rscript/plot"
ggsave(file.path(savedir_png,paste0(gene,"_RNAseq_time_serious_tpm",".png",sep="")),dpi = 300,width = 8,height = 5)

}
#-----------------------------------------------------------------------
