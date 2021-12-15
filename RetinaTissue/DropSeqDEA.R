# Load libraries
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("dplyr")
library("ggrepel")

ColorMap <- brewer.pal(11, "PiYG")

# Load Data
read.counts<-read.table("DropSeq/DropSeqLayerFrequencies.csv", header=TRUE, sep=',', stringsAsFactors=FALSE)
row.names(read.counts)<-read.counts[,1]
read.counts<-read.counts[, -c(0:1)]
sample.info<-data.frame(Mouse=c("1", "2",  "3", "4", "5", "6", "1", "2", "3", "4", "5", "6"), condition=c(rep("ONL",6), rep("BP",6)), row.names=names(read.counts))

# Differential Expression analysis
# Create DESeq object
DESeq.ds<-DESeqDataSetFromMatrix(countData=read.counts, colData=sample.info, design = ~ condition)
DESeq.ds <- DESeq(DESeq.ds)

results.BP.ONL <- results(DESeq.ds, pAdjustMethod="BH", contrast = c("condition", "BP","ONL"))
DGE.results <- c(results.BP.ONL)

table(results.BP.ONL$padj<0.05)

# Sort and obtain differentially expressed genes in a csv file
results.BP.ONL.sorted <- results.BP.ONL[order(results.BP.ONL$padj),]

DGEgenes.BP.ONL <- rownames(subset(results.BP.ONL.sorted, padj<0.05))
All.DGEgenes <- c(DGEgenes.BP.ONL)

DE_genes.BP.ONL <- as.data.frame(results.BP.ONL.sorted)
write.csv(DE_genes.BP.ONL, "DropSeq/BPvsONL.csv")

#DE genes that are specific to each population
ONL.Pos <- rownames(subset(results.BP.ONL.sorted, log2FoldChange<0&padj<0.05))
BP.Pos <- rownames(subset(results.BP.ONL.sorted, log2FoldChange>0&padj<0.05))

# Heatmap plot of top 20 differentially expressed genes for each population
DGE_Top<-ONL.Pos[1:20]
hm.mat_DGEgenes<-log.norm.counts[DGE_Top,]
hm.mat_DGEgenes<-hm.mat_DGEgenes[,c("ONL_1", "ONL_2", "ONL_3", "ONL_4", "ONL_5", "BP_1", "BP_2", "BP_3", "BP_4", "BP_5")]
pdf(file="DropSeq/DropSeq_Top_DEGenes_ONL_Pos_Heatmap.pdf", onefile=FALSE)
pheatmap(hm.mat_DGEgenes, cluster_rows=F, cluster_cols=F, scale="row", col = ColorMap)
dev.off()

DGE_Top<-BP.Pos[1:20]
hm.mat_DGEgenes<-log.norm.counts[DGE_Top,]
hm.mat_DGEgenes<-hm.mat_DGEgenes[,c("ONL_1", "ONL_2", "ONL_3", "ONL_4", "ONL_5", "BP_1", "BP_2", "BP_3", "BP_4", "BP_5")]
pdf(file="DropSeq/DropSeq_Top_DEGenes_BP_Pos_Heatmap.pdf", onefile=FALSE)
pheatmap(hm.mat_DGEgenes, cluster_rows=F, cluster_cols=F, scale="row", col = ColorMap)
dev.off()

# Heatmap plots of cell type markers
DGE_Top<-c(ONL.Pos, BP.Pos)
hm.mat_DGEgenes<-log.norm.counts[DGE_Top,]
hm.mat_DGEgenes<-hm.mat_DGEgenes[,c("ONL_1", "ONL_2", "ONL_3", "ONL_4", "ONL_5", "BP_1", "BP_2", "BP_3", "BP_4", "BP_5")]
pdf(file="DropSeq/DropSeq_Top_DEGenes_AllLayers_Pos_Heatmap.pdf", onefile=FALSE)
pheatmap(hm.mat_DGEgenes, cluster_rows=F, cluster_cols=F, scale="row", col = ColorMap)
dev.off()
