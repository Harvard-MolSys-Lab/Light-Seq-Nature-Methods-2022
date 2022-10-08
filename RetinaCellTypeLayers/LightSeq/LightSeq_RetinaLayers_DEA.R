# Load libraries
library("DESeq2")
library("pheatmap")
library("ggplot2")
library("dplyr")
library("ggrepel")
library("RColorBrewer")

# Load count data 
countData <- read.csv('../LightSeq/ReorderedLightSeq.csv', header = TRUE, sep = ",")
row.names(countData)<-countData[,1]
countData<-countData[, -c(0:1)]

# Record sample information
sample.info<-data.frame(Mouse=c("1","2",  "3", "4", "1", "2", "3", "4", "1", "2", "3", "4"), condition=c(rep("RGC",4), rep("BP",4), rep("ONL",4)), row.names=names(countData))

# Differential expression analysis
# Create DESeq object, comparing according to "condition", or retinal layer
DESeq.ds<-DESeqDataSetFromMatrix(countData=countData, colData=sample.info, design = ~ condition)

DESeq.ds <- DESeq.ds[rowSums(counts(DESeq.ds)) >0,]
DESeq.ds <- estimateSizeFactors(DESeq.ds)
counts.sf_normalized<-counts(DESeq.ds, normalized=TRUE)
log.norm.counts <- log2(counts.sf_normalized +1)
str(colData(DESeq.ds)$condition)

# Run DESeq
DESeq.ds <- DESeq(DESeq.ds)
res <- results(DESeq.ds)
head(results(DESeq.ds, tidy=TRUE))
summary(res)
res <- res[order(res$padj),]
head(res)

#Vst function will perform variance stabilizing transformation on raw count data
vsdata <- vst(DESeq.ds, blind=FALSE)

# Plot PCA for Supplemental Figure
cbPalette <- c( "#0EE5E8","#65E80E","#E82BA6") # Set colormap to correspond to barcodes

pcaData <- plotPCA(vsdata, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
# Generate PCA plot for export 
p<-ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=10,alpha=.7) + scale_color_manual(values=cbPalette)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +theme_bw()
p

# Pairwuse comparison between cell layers
results.RGC.BP <- results(DESeq.ds,  pAdjustMethod="BH", contrast=c("condition","BP","RGC"))
results.RGC.BP <- results.RGC.BP[order(results.RGC.BP$padj),]

results.BP.ONL <- results(DESeq.ds,  pAdjustMethod="BH", contrast=c("condition","BP","ONL"))
results.BP.ONL <- results.BP.ONL[order(results.BP.ONL$padj),]

results.ONL.RGC <- results(DESeq.ds,  pAdjustMethod="BH", contrast=c("condition","RGC","ONL"))
results.ONL.RGC <- results.ONL.RGC[order(results.ONL.RGC$padj),]

pdf(file="Volcano_Plot_AllLayers.pdf", onefile=FALSE)

## Volcano Plots
par(mfrow=c(1,3))
# Make a basic volcano plot
# List of well-known retinal genes to label on plots
sigr=c("Rho","Rp1","Arr3","Rs1","Pde6a","Pde6g","Opn1sw","Arr3","Msi1","Stk35","Ctbp2","Crx","Rp1l1","Dpysl3","Kcnv2","Opn1mw","Cnga1")
sigb=c("Grik1","Grm6","Vsx2","Car8","Trpm1","Neurod4","Prkca","Prox1","Gpr179","Abca8a","Aqp4","Bhlhe23","Igfn1","Neto1","Vsx1")
sigg=c("Tubb3","Cplx2","Syn1","Nefm","Pou4f1","Cplx2","Dpysl2","Calm3","Nefh","Opn4","Calb1","Calb2","Rbpms","Pou4f2","Pou4f3","Tfap2a")

# BPL vs ONL
with(results.BP.ONL, plot(log2FoldChange, -log10(pvalue), pch=1, lwd=0.2, main="ONL vs. BPL", xlim=c(-6.3,6.3), col=alpha(rgb(1,1,1),0.4),yaxt='n'))
with(subset(results.BP.ONL, results.BP.ONL$padj>0.05 ),points(log2FoldChange, -log10(pvalue), pch=16, col=alpha(rgb(0.2,0.2,0.2),0.4)))
axis(2, at = seq(0, 140, 50))
with(subset(results.BP.ONL, results.BP.ONL$log2FoldChange>0 & results.BP.ONL$padj<0.05 ), points(log2FoldChange, -log10(pvalue), pch=16, col=alpha(rgb(0,.9,.9),0.4)))
with(subset(results.BP.ONL, results.BP.ONL$log2FoldChange<(0) & results.BP.ONL$padj<0.05), points(log2FoldChange, -log10(pvalue), pch=16, col=alpha(rgb(0.9,0,0.9),0.4)))
with(subset(results.BP.ONL, rownames(results.BP.ONL) %in% sigr), points(log2FoldChange, -log10(pvalue), pch=16, col=rgb(.5,.2,.2)))
with(subset(results.BP.ONL, rownames(results.BP.ONL) %in% sigb), points(log2FoldChange, -log10(pvalue), pch=16, col=rgb(.2,.2,.5)))
sign.genes1=which((results.BP.ONL$log2FoldChange)>0 & results.BP.ONL$padj<0.05)
sign.genes2=which((results.BP.ONL$log2FoldChange)< (0) & results.BP.ONL$padj<0.05)
sigrows=which(rownames(results.BP.ONL) %in% sigr)
sigrowsb=which(rownames(results.BP.ONL) %in% sigb)
# Label genes of interest
text(x=results.BP.ONL$log2FoldChange[sigrows] , y=-log10(results.BP.ONL$pvalue[sigrows]), label=row.names(results.BP.ONL)[sigrows], cex=1)
text(x=results.BP.ONL$log2FoldChange[sigrowsb] , y=-log10(results.BP.ONL$pvalue[sigrowsb]), label=row.names(results.BP.ONL)[sigrowsb], cex=1)

# ONL vs GCL
with(results.ONL.RGC, plot(log2FoldChange, -log10(pvalue), pch=1, lwd=0.2, main="ONL vs. RGC",xlim=c(-6.3,6.3) ,col=alpha(rgb(1,1,1),0.4),yaxt='n'))
with(subset(results.ONL.RGC, results.ONL.RGC$padj>0.05 ),points(log2FoldChange, -log10(pvalue), pch=16, col=alpha(rgb(0.2,0.2,0.2),0.4)))
axis(2, at = seq(0, 250, 50))
with(subset(results.ONL.RGC, results.ONL.RGC$padj>0.05), points(log2FoldChange, -log10(pvalue), pch=16, col=alpha(rgb(.3,.3,.3),0.4)))
with(subset(results.ONL.RGC, results.ONL.RGC$log2FoldChange>0 & results.ONL.RGC$padj<0.05 ), points(log2FoldChange, -log10(pvalue), pch=16, col=alpha(rgb(0,.9,0),0.4)))
with(subset(results.ONL.RGC, results.ONL.RGC$log2FoldChange<(0) & results.ONL.RGC$padj<0.05), points(log2FoldChange, -log10(pvalue), pch=16, col=alpha(rgb(.9,0,.9),0.4)))
with(subset(results.ONL.RGC, rownames(results.ONL.RGC) %in% sigr), points(log2FoldChange, -log10(pvalue), pch=16, col=rgb(.5,.2,.2)))
with(subset(results.ONL.RGC, rownames(results.ONL.RGC) %in% sigg), points(log2FoldChange, -log10(pvalue), pch=16, col=rgb(.2,.5,.2)))
sign.genes1=which((results.ONL.RGC$log2FoldChange)>0 & results.ONL.RGC$padj<0.05)
sign.genes2=which((results.ONL.RGC$log2FoldChange)<0 & results.ONL.RGC$padj<0.05)
sigrows=which(rownames(results.ONL.RGC) %in% sigr)
sigrowsg=which(rownames(results.ONL.RGC) %in% sigg)
# Label genes of interest
text(x=results.ONL.RGC$log2FoldChange[sigrows], y=-log10(results.ONL.RGC$pvalue[sigrows]), label=row.names(results.ONL.RGC)[sigrows], cex=1)
text(x=results.ONL.RGC$log2FoldChange[sigrowsg] , y=-log10(results.ONL.RGC$pvalue[sigrowsg]), label=row.names(results.ONL.RGC)[sigrowsg], cex=1)

# BPL vs GCL
with(results.RGC.BP, plot(log2FoldChange, -log10(pvalue), pch=1, lwd=0.2, main="BPL vs. GCL", xlim=c(-6.7,6.7) , ylim=c(0,160),col=alpha(rgb(1,1,1),0.4),yaxt='n'))
with(subset(results.RGC.BP, results.RGC.BP$padj>0.05 ),points(log2FoldChange, -log10(pvalue),pch=16, col=alpha(rgb(0.2,0.2,0.2),0.4)))
axis(2, at = seq(0, 150, 50))
with(subset(results.RGC.BP, results.RGC.BP$padj>0.05 ),points(log2FoldChange, -log10(pvalue),pch=16, col=alpha(rgb(0.2,0.2,0.2),0.4)))
with(subset(results.RGC.BP, results.RGC.BP$log2FoldChange>0 & results.RGC.BP$padj<0.05 ),points(log2FoldChange, -log10(pvalue), pch=16, col=alpha(rgb(0,0.9,0.9),0.4)))
with(subset(results.RGC.BP, results.RGC.BP$log2FoldChange<(0) & results.RGC.BP$padj<0.05), points(log2FoldChange, -log10(pvalue), pch=16, col=alpha(rgb(0,0.9,0),0.4)))
with(subset(results.RGC.BP, rownames(results.RGC.BP) %in% sigg), points(log2FoldChange, -log10(pvalue), pch=16, col=rgb(.2,.5,.2)))
with(subset(results.RGC.BP, rownames(results.RGC.BP) %in% sigb), points(log2FoldChange, -log10(pvalue), pch=16, col=rgb(.2,.2,.5)))
sign.genes1=which(results.RGC.BP$padj<0.05)
sign.genes2=which( results.RGC.BP$padj<0.05)
sigrows=which(rownames(results.RGC.BP) %in% sigg)
sigrowsb=which(rownames(results.RGC.BP) %in% sigb)
# Label genes of interest
text(x=results.RGC.BP$log2FoldChange[sigrows] , y=-log10(results.RGC.BP$pvalue[sigrows]), label=row.names(results.RGC.BP)[sigrows], cex=1)
text(x=results.RGC.BP$log2FoldChange[sigrowsb] , y=-log10(results.RGC.BP$pvalue[sigrowsb]), label=row.names(results.RGC.BP)[sigrowsb], cex=1)

dev.off()

# Create heatmap in FIGURE
ColorMap <- brewer.pal(11, "PiYG")
DGE.results <- c(results.RGC.BP, results.BP.ONL, results.ONL.RGC)

# How many pairwise-significantly differential genes per layer?
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
write.csv(DE_genes.RGC.BP, "RGCvsBP.csv")
write.csv(DE_genes.BP.ONL, "BPvsONL.csv")
write.csv(DE_genes.ONL.RGC, "ONLvsRGC.csv")

#DE genes that are specific to each population
ONL.Pos <- intersect(rownames(subset(results.BP.ONL.sorted, log2FoldChange<0&padj<0.05)), rownames(subset(results.ONL.RGC.sorted, log2FoldChange<0&padj<0.05)))
ONL.Neg <- intersect(rownames(subset(results.BP.ONL.sorted, log2FoldChange>0&padj<0.05)), rownames(subset(results.ONL.RGC.sorted, log2FoldChange>0&padj<0.05)))

RGC.Pos <- intersect(rownames(subset(results.ONL.RGC.sorted, log2FoldChange>0&padj<0.05)), rownames(subset(results.RGC.BP.sorted, log2FoldChange<0&padj<0.05)))
RGC.Neg <- intersect(rownames(subset(results.ONL.RGC.sorted, log2FoldChange<0&padj<0.05)), rownames(subset(results.RGC.BP.sorted, log2FoldChange>0&padj<0.05)))

BP.Pos <- intersect(rownames(subset(results.BP.ONL.sorted, log2FoldChange>0&padj<0.05)), rownames(subset(results.RGC.BP.sorted, log2FoldChange>0&padj<0.05)))
BP.Neg <- intersect(rownames(subset(results.BP.ONL.sorted, log2FoldChange<0&padj<0.05)), rownames(subset(results.RGC.BP.sorted, log2FoldChange<0&padj<0.05)))

write.csv(as.data.frame(ONL.Pos), "ONLmarkers.csv")
write.csv(as.data.frame(BP.Pos), "BPmarkers.csv")
write.csv(as.data.frame(RGC.Pos), "RGCmarkers.csv")

# Heatmap plot of top 20 differentially expressed genes for each population
DGE_Top<-ONL.Pos[1:20]
hm.mat_DGEgenes<-log.norm.counts[DGE_Top,]
hm.mat_DGEgenes<-hm.mat_DGEgenes[,c("ONL_1", "ONL_2", "ONL_3", "ONL_4", "BP_1", "BP_2", "BP_3", "BP_4", "RGC_1", "RGC_2", "RGC_3", "RGC_4")]
pdf(file="Top_DEGenes_ONL_Pos_Heatmap.pdf", onefile=FALSE)
pheatmap(hm.mat_DGEgenes, cluster_rows=F, cluster_cols=F, scale="row", col = ColorMap)
dev.off()

DGE_Top<-BP.Pos[1:20]
hm.mat_DGEgenes<-log.norm.counts[DGE_Top,]
hm.mat_DGEgenes<-hm.mat_DGEgenes[,c("ONL_1", "ONL_2", "ONL_3", "ONL_4", "BP_1", "BP_2", "BP_3", "BP_4", "RGC_1", "RGC_2", "RGC_3", "RGC_4")]
pdf(file="Top_DEGenes_BP_Pos_Heatmap.pdf", onefile=FALSE)
pheatmap(hm.mat_DGEgenes, cluster_rows=F, cluster_cols=F, scale="row", col = ColorMap)
dev.off()

DGE_Top<-RGC.Pos[1:20]
hm.mat_DGEgenes<-log.norm.counts[DGE_Top,]
hm.mat_DGEgenes<-hm.mat_DGEgenes[,c("ONL_1", "ONL_2", "ONL_3", "ONL_4", "BP_1", "BP_2", "BP_3", "BP_4", "RGC_1", "RGC_2", "RGC_3", "RGC_4")]
pdf(file="Top_DEGenes_RGC_Pos_Heatmap.pdf", onefile=FALSE)
pheatmap(hm.mat_DGEgenes, cluster_rows=F, cluster_cols=F, scale="row", col = ColorMap)
dev.off()

# Heatmap plots of cell type markers
DGE_Top<-c(ONL.Pos, BP.Pos, RGC.Pos)
hm.mat_DGEgenes<-log.norm.counts[DGE_Top,]
hm.mat_DGEgenes<-hm.mat_DGEgenes[,c("ONL_1", "ONL_2", "ONL_3", "ONL_4", "BP_1", "BP_2", "BP_3", "BP_4", "RGC_1", "RGC_2", "RGC_3", "RGC_4")]
pdf(file="Top_DEGenes_AllLayers_Pos_Heatmap.pdf", onefile=FALSE)
pheatmap(hm.mat_DGEgenes, cluster_rows=F, cluster_cols=F, scale="row", col = ColorMap)
dev.off()
