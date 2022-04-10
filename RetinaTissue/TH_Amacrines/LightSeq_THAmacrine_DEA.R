# Load libraries
library("DESeq2")
library("pheatmap")
library("dplyr")

# Import Counts
countData <- read.csv('./TH_Amacrines/ReorderedLightSeq.csv', header = TRUE, sep = ",")
head(countData)
mydata2 = select(countData, -1)
row.names(mydata2) <- countData$Gene
sample.info<-data.frame(Replicate=c("1", "2",  "3", "4", "5", "1", "2", "3", "4", "5"), condition=c(rep("TH",5), rep("THneg",5)), row.names=names(mydata2))

# Perform differential expression analysis
DESeq.ds <- DESeqDataSetFromMatrix(countData=mydata2, 
                                   colData=sample.info, 
                                   design=~condition, tidy = FALSE)

DESeq.ds <- DESeq.ds[rowSums(counts(DESeq.ds)) >0,]
DESeq.ds <- estimateSizeFactors(DESeq.ds)
counts.sf_normalized<-counts(DESeq.ds, normalized=TRUE)
log.norm.counts <- log2(counts.sf_normalized +1)
DESeq.ds <- DESeq(DESeq.ds)
res <- results(DESeq.ds)
head(results(DESeq.ds, tidy=TRUE))
summary(res)
res <- res[order(res$padj),]
head(res)

# Sort by adjusted p-value
res.sorted <- res[order(res$padj),]
#DE genes that are specific to TH+ population
Pos <- (rownames(subset(res.sorted, padj<(.05) & log2FoldChange<1)))

# Heatmap plot of top significantly enriched genes
DGE_Top<-c(Pos)
hm.mat_DGEgenes<-log.norm.counts[DGE_Top,]
hm.mat_DGEgenes<-hm.mat_DGEgenes[]
pheatmap(hm.mat_DGEgenes, cluster_rows=F, cluster_cols=F, scale="row")

# Pull out top marker list for validation with padj<0.05 and log2FoldChange>3 for RNA FISH validation
# Note: Sstr2 was excluded because the RNA FISH probe preparation failed and therefore it was not applied to tissue for validation
Markers=res[Pos,]
TopMarkers=subset(Markers,log2FoldChange>3)

DE_genes.TH.v.noTH <- as.data.frame(res.sorted)
write.csv(DE_genes.TH.v.noTH, "./TH_Amacrines/TH_Amacrine_DEGs.csv")
