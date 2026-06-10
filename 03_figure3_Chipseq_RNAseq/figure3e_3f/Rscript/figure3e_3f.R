setwd("/Users/chien/Library/CloudStorage/OneDrive-OIST/data/00_All_data_and_script_for_manuscript/03_figure3_Chipseq_RNAseq/figure3e_3f/data")
tpm <-read.table("salmon.merged.gene_tpm.txt",header =T)

# only select 6 to 72 hr no restimulate
tpm <- tpm[, c(2,9:(ncol(tpm)-6))]

tpm <-tpm[rowSums(tpm[ , c(2:ncol(tpm))])>0,]

#-----------------------------------------------------------------------
#plot RNAseq time serious
for (gene in c("Il17a","Il23r","Il21","Ifng","Irf8","Stat5a","Bcl2l11")) {
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


savedir_png <- "/Users/chien/Library/CloudStorage/OneDrive-OIST/data/00_All_data_and_script_for_manuscript/03_figure3_Chipseq_RNAseq/figure3e_3f/Rscript/plot"
ggsave(file.path(savedir_png,paste0(gene,"_RNAseq_time_serious_tpm",".png",sep="")),dpi = 300,width = 8,height = 5)

}
#-----------------------------------------------------------------------
#plot RNAseq time serious no label for illustrator
for (gene in c("Il17a","Il23r","Il21","Ifng","Irf8","Stat5a","Bcl2l11")) {
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
    geom_point(size=3,aes(shape=genotype))+
    scale_shape_manual(values = c(16, 1)) +
    theme_classic()+
    geom_smooth(aes(x=time, y=tpm, color = genotype, fill=genotype)) + 
    theme (legend.position = "none",
           plot.title = element_blank(),  
           axis.title.x=element_blank(),
           axis.title.y=element_blank(),
           axis.text.x=element_blank(),
           axis.text.y=element_blank(),
           axis.line = element_line(colour = "black", size = 1, linetype = "solid"))+
    scale_color_manual(values=c("blue", "red"))+
    scale_fill_manual(values=c("blue", "red"))
  p
  
  
  savedir_png <- "/Users/chien/Library/CloudStorage/OneDrive-OIST/data/00_All_data_and_script_for_manuscript/03_figure3_Chipseq_RNAseq/figure3e_3f/Rscript/plot"
  ggsave(file.path(savedir_png,paste0(gene,"_RNAseq_time_serious_tpm","_no_label",".png",sep="")),dpi = 300,width = 8,height = 5)
  
}


#-----------------------------------------------------------------------



