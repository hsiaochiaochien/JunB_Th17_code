#!/usr/bin/Rscript

setwd("/Users/chien/Library/CloudStorage/OneDrive-OIST/data/00_All_data_and_script_for_manuscript/02_figure1_atacseq_invivo/Rscript/")

#generate data for DESeq2
tag_file_annotated <- read.table("figure1b_invivo_ATACseq_Day3_data.txt",sep="\t",skip=1,quote = "")

title <- c("Row.names","Chr","Start","End","Gene Name","Annotation","Distance_to_nearest_RefSeq_TSS","OT2-WT1","OT2-WT2","OT2-WT3","OT2-KO1","OT2-KO2","OT2-KO3")
names(tag_file_annotated)[1:13] <- title

sampleNames <- c("OT2-WT1","OT2-WT2","OT2-WT3","OT2-KO1","OT2-KO2","OT2-KO3")

countData <- as.matrix(tag_file_annotated[8:13])
mode(countData) <- "integer"
rownames(countData) <- tag_file_annotated$Row.names
database <- data.frame(name=sampleNames,condition=c("WT", "WT","WT","KO","KO","KO"))
rownames(database) <- sampleNames


rm (sampleNames)

#run DESeq2
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData, colData=database, design= ~ condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
res <- results(dds,contrast=c("condition","KO","WT"))
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)

#assing 1 to p value or FDR
resdata[["padj"]][is.na(resdata[["padj"]])] <- 1
resdata[["pvalue"]][is.na(resdata[["pvalue"]])] <- 1

rm(countData)
rm(database)
rm(res)

#annotation

onlyname <-subset(tag_file_annotated,select=c("Row.names","Chr","Start","End","Gene Name","Annotation","Distance_to_nearest_RefSeq_TSS"))

finaltable <- merge(resdata,onlyname,by="Row.names",sort=FALSE)

rm(tag_file_annotated)
rm(onlyname)
rm(resdata)


#scatter plot tag data only by log2FC

finaltable$WT_mean <- rowMeans(finaltable[,c("OT2-WT1","OT2-WT2","OT2-WT3")], na.rm=TRUE)
finaltable$KO_mean <- rowMeans(finaltable[,c("OT2-KO1","OT2-KO2","OT2-KO3")], na.rm=TRUE)
finaltable$WT_mean_log2 <- log(finaltable[,c('WT_mean')], 2)
finaltable$KO_mean_log2 <- log(finaltable[,c('KO_mean')], 2)

#-------------------------------------------------------------------------------------------------
#pvalue<0.05,log2FC1
#add color by using tag only data
p_value = 0.05
finaltable$color <- ifelse(finaltable$pvalue < p_value & abs(finaltable$log2FoldChange) >= 1,ifelse(finaltable$log2FoldChange > 1 ,'red','blue'),'gray')

#count the log2FC
library(dplyr)
df_ko <- dplyr::filter(finaltable, pvalue < p_value & log2FoldChange > 1)
df_wt <- dplyr::filter(finaltable, pvalue < p_value & log2FoldChange < -1)


paste("ko=",nrow(df_ko))
paste("wt=",nrow(df_wt))
df <- dplyr::filter(finaltable, pvalue < p_value )
nodifference <-nrow(finaltable)- nrow(df_ko)-nrow(df_wt)
paste(nodifference)
#rm(finaltable)

color<- c(red = "red", gray = "gray", blue ="blue")

library(ggplot2)
p <- ggplot(finaltable, aes(WT_mean_log2, KO_mean_log2, col= color)) +
  geom_point(size = 0.8) +
  theme_bw() +
  scale_color_manual(values = color) +
  labs(x="WT tags (log2)",y="JunB-KO tags (log2)")+
  ggtitle("Day3 activated OTII invivo") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        panel.grid=element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))+
  annotate("text", x=Inf,y=Inf, hjust=1,vjust=1,label = nodifference ,size= 8,col= "gray" )+
  scale_x_continuous(limits = c(0, 12), breaks = seq(0, 12, by = 4))+
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, by = 4))+
  annotate("text", x=Inf,y=-Inf, hjust=1,vjust=0,label =nrow(df_wt),size= 8,col= "blue" )+
  annotate("text", x=-Inf,y=Inf, hjust=0,vjust=1,label =nrow(df_ko),size= 8,col= "red")
p

# #==============================================================================================
#lable base genes you want
finaltable$`Gene Name`<-as.character(finaltable$`Gene Name`)
finaltable$genelabels <- ifelse(finaltable$`Gene Name` == "Il17a" & finaltable$color == 'blue'|
                         finaltable$`Gene Name` == "Il23r" & finaltable$color == 'blue'|
                        finaltable$`Gene Name` == "Ifng" & finaltable$color == 'red'|
                        finaltable$`Gene Name` == "Il4" & finaltable$color == 'red'|
                        finaltable$`Gene Name` == "Il10" & finaltable$color == 'red'|
                        finaltable$`Gene Name` == "Rorc" & finaltable$color == 'blue'                          
                        , finaltable$`Gene Name`,"")

head(finaltable)


#for specific gene name
library(ggrepel)
options(ggrepel.max.overlaps = 100000)
q <- p + geom_text_repel(label = finaltable$genelabels, colour='black',size = 7, box.padding = 1,label.padding =10 )
q

