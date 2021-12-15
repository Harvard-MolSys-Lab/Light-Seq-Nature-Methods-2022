# Load libraries
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("dplyr")
library("ggrepel")

ColorMap <- brewer.pal(11, "PiYG")

# Load Data
read.counts<-read.table("LightSeq/ReorderedLightSeq.csv", header=TRUE, sep=',', stringsAsFactors=FALSE)
row.names(read.counts)<-read.counts[,1]
read.counts<-read.counts[, -c(0:1)]
sample.info<-data.frame(Mouse=c("1","2",  "3", "4", "1", "2", "3", "4", "1", "2", "3", "4"), condition=c(rep("RGC",4), rep("BP",4), rep("ONL",4)), row.names=names(read.counts))

# Differential Expression analysis
# Create DESeq object
DESeq.ds<-DESeqDataSetFromMatrix(countData=read.counts, colData=sample.info, design = ~ condition)
DESeq.ds <- DESeq(DESeq.ds)

results.RGC.BP <- results(DESeq.ds, pAdjustMethod="BH", contrast = c("condition", "RGC","BP"))
results.BP.ONL <- results(DESeq.ds, pAdjustMethod="BH", contrast = c("condition", "BP","ONL"))
results.ONL.RGC <- results(DESeq.ds, pAdjustMethod="BH", contrast = c("condition", "ONL","RGC"))
DGE.results <- c(results.RGC.BP, results.BP.ONL, results.ONL.RGC)

table(results.RGC.BP$padj<0.05)
table(results.BP.ONL$padj<0.05)
table(results.ONL.RGC$padj<0.05)

# Sort and obtain differentially expressed genes in a csv file
results.RGC.BP.sorted <- results.RGC.BP[order(results.RGC.BP$padj),]
results.BP.ONL.sorted <- results.BP.ONL[order(results.BP.ONL$padj),]
results.ONL.RGC.sorted <- results.ONL.RGC[order(results.ONL.RGC$padj),]

DGEgenes.RGC.BP <- rownames(subset(results.RGC.BP.sorted, padj<0.05))
DGEgenes.BP.ONL <- rownames(subset(results.BP.ONL.sorted, padj<0.05))
DGEgenes.ONL.RGC <- rownames(subset(results.ONL.RGC.sorted, padj<0.05))
All.DGEgenes <- c(DGEgenes.ONL.RGC, DGEgenes.RGC.BP, DGEgenes.BP.ONL)

DE_genes.RGC.BP <- as.data.frame(results.RGC.BP.sorted)
DE_genes.BP.ONL <- as.data.frame(results.BP.ONL.sorted)
DE_genes.ONL.RGC <- as.data.frame(results.ONL.RGC.sorted)
write.csv(DE_genes.RGC.BP, "LightSeq/RGCvsBP.csv")
write.csv(DE_genes.BP.ONL, "LightSeq/BPvsONL.csv")
write.csv(DE_genes.ONL.RGC, "LightSeq/ONLvsRGC.csv")

#DE genes that are specific to each population
ONL.Pos <- intersect(rownames(subset(results.BP.ONL.sorted, log2FoldChange<0&padj<0.05)), rownames(subset(results.ONL.RGC.sorted, log2FoldChange>0&padj<0.05)))
ONL.Neg <- intersect(rownames(subset(results.BP.ONL.sorted, log2FoldChange>0&padj<0.05)), rownames(subset(results.ONL.RGC.sorted, log2FoldChange<0&padj<0.05)))

RGC.Pos <- intersect(rownames(subset(results.ONL.RGC.sorted, log2FoldChange<0&padj<0.05)), rownames(subset(results.RGC.BP.sorted, log2FoldChange>0&padj<0.05)))
RGC.Neg <- intersect(rownames(subset(results.ONL.RGC.sorted, log2FoldChange>0&padj<0.05)), rownames(subset(results.RGC.BP.sorted, log2FoldChange<0&padj<0.05)))

BP.Pos <- intersect(rownames(subset(results.BP.ONL.sorted, log2FoldChange>0&padj<0.05)), rownames(subset(results.RGC.BP.sorted, log2FoldChange<0&padj<0.05)))
BP.Neg <- intersect(rownames(subset(results.BP.ONL.sorted, log2FoldChange<0&padj<0.05)), rownames(subset(results.RGC.BP.sorted, log2FoldChange>0&padj<0.05)))

write.csv(as.data.frame(ONL.Pos), "LightSeq/ONLmarkers.csv")
write.csv(as.data.frame(BP.Pos), "LightSeq/BPmarkers.csv")
write.csv(as.data.frame(RGC.Pos), "LightSeq/RGCmarkers.csv")

# Heatmap plot of top 20 differentially expressed genes for each population
DGE_Top<-ONL.Pos[1:20]
hm.mat_DGEgenes<-log.norm.counts[DGE_Top,]
hm.mat_DGEgenes<-hm.mat_DGEgenes[,c("ONL_1", "ONL_2", "ONL_3", "ONL_4", "BP_1", "BP_2", "BP_3", "BP_4", "RGC_1", "RGC_2", "RGC_3", "RGC_4")]
pdf(file="LightSeq/Top_DEGenes_ONL_Pos_Heatmap.pdf", onefile=FALSE)
pheatmap(hm.mat_DGEgenes, cluster_rows=F, cluster_cols=F, scale="row", col = ColorMap)
dev.off()

DGE_Top<-ONL.Neg[1:20]
hm.mat_DGEgenes<-log.norm.counts[DGE_Top,]
hm.mat_DGEgenes<-hm.mat_DGEgenes[,c("ONL_1", "ONL_2", "ONL_3", "ONL_4", "BP_1", "BP_2", "BP_3", "BP_4", "RGC_1", "RGC_2", "RGC_3", "RGC_4")]
pdf(file="LightSeq/Top_DEGenes_ONL_Neg_Heatmap.pdf", onefile=FALSE)
pheatmap(hm.mat_DGEgenes, cluster_rows=F, cluster_cols=F, scale="row", col = ColorMap)
dev.off()

DGE_Top<-BP.Pos[1:20]
hm.mat_DGEgenes<-log.norm.counts[DGE_Top,]
hm.mat_DGEgenes<-hm.mat_DGEgenes[,c("ONL_1", "ONL_2", "ONL_3", "ONL_4", "BP_1", "BP_2", "BP_3", "BP_4", "RGC_1", "RGC_2", "RGC_3", "RGC_4")]
pdf(file="LightSeq/Top_DEGenes_BP_Pos_Heatmap.pdf", onefile=FALSE)
pheatmap(hm.mat_DGEgenes, cluster_rows=F, cluster_cols=F, scale="row", col = ColorMap)
dev.off()

DGE_Top<-BP.Neg[1:20]
hm.mat_DGEgenes<-log.norm.counts[DGE_Top,]
hm.mat_DGEgenes<-hm.mat_DGEgenes[,c("ONL_1", "ONL_2", "ONL_3", "ONL_4", "BP_1", "BP_2", "BP_3", "BP_4", "RGC_1", "RGC_2", "RGC_3", "RGC_4")]
pdf(file="LightSeq/Top_DEGenes_BP_Neg_Heatmap.pdf", onefile=FALSE)
pheatmap(hm.mat_DGEgenes, cluster_rows=F, cluster_cols=F, scale="row", col = ColorMap)
dev.off()

DGE_Top<-RGC.Pos[1:20]
hm.mat_DGEgenes<-log.norm.counts[DGE_Top,]
hm.mat_DGEgenes<-hm.mat_DGEgenes[,c("ONL_1", "ONL_2", "ONL_3", "ONL_4", "BP_1", "BP_2", "BP_3", "BP_4", "RGC_1", "RGC_2", "RGC_3", "RGC_4")]
pdf(file="LightSeq/Top_DEGenes_RGC_Pos_Heatmap.pdf", onefile=FALSE)
pheatmap(hm.mat_DGEgenes, cluster_rows=F, cluster_cols=F, scale="row", col = ColorMap)
dev.off()

DGE_Top<-RGC.Neg[1:20]
hm.mat_DGEgenes<-log.norm.counts[DGE_Top,]
hm.mat_DGEgenes<-hm.mat_DGEgenes[,c("ONL_1", "ONL_2", "ONL_3", "ONL_4", "BP_1", "BP_2", "BP_3", "BP_4", "RGC_1", "RGC_2", "RGC_3", "RGC_4")]
pdf(file="LightSeq/Top_DEGenes_RGC_Neg_Heatmap.pdf", onefile=FALSE)
pheatmap(hm.mat_DGEgenes, cluster_rows=F, cluster_cols=F, scale="row", col = ColorMap)
dev.off()

# Heatmap plots of cell type markers
DGE_Top<-c(ONL.Pos, BP.Pos, RGC.Pos)
hm.mat_DGEgenes<-log.norm.counts[DGE_Top,]
hm.mat_DGEgenes<-hm.mat_DGEgenes[,c("ONL_1", "ONL_2", "ONL_3", "ONL_4", "BP_1", "BP_2", "BP_3", "BP_4", "RGC_1", "RGC_2", "RGC_3", "RGC_4")]
pdf(file="LightSeq/Top_DEGenes_AllLayers_Pos_Heatmap.pdf", onefile=FALSE)
pheatmap(hm.mat_DGEgenes, cluster_rows=F, cluster_cols=F, scale="row", col = ColorMap)
dev.off()

DGE_Top<-c(ONL.Pos[1:10], BP.Pos[1:10], RGC.Pos[1:10])
hm.mat_DGEgenes<-log.norm.counts[DGE_Top,]
hm.mat_DGEgenes<-hm.mat_DGEgenes[,c("ONL_1", "ONL_2", "ONL_3", "ONL_4", "BP_1", "BP_2", "BP_3", "BP_4", "RGC_1", "RGC_2", "RGC_3", "RGC_4")]
pdf(file="LightSeq/Top_DEGenes_AllLayers_JustTop10_Heatmap.pdf", onefile=FALSE)
pheatmap(hm.mat_DGEgenes, cluster_rows=F, cluster_cols=F, scale="row", col = ColorMap)
dev.off()
