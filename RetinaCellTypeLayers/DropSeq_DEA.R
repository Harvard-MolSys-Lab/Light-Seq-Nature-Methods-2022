# Load libraries
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("dplyr")
library("ggrepel")

ColorMap <- brewer.pal(11, "PiYG")

# Load Data
read.counts<-read.table("./DropSeq/DropSeqLayerFrequencies.csv", header=TRUE, sep=',', stringsAsFactors=FALSE)
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
write.csv(DE_genes.BP.ONL, "./DropSeq/BPvsONL.csv")

#DE genes that are specific to each population in Drop-Seq data
ONL.Pos <- rownames(subset(results.BP.ONL.sorted, log2FoldChange<0&padj<0.05))
BP.Pos <- rownames(subset(results.BP.ONL.sorted, log2FoldChange>0&padj<0.05))
