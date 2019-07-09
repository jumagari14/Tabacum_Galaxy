## Author: Juan Manuel Garcia
## Based on the workflow from https://gist.github.com/stephenturner/f60c1934405c127f09a6

setwd("/media/jumagari/JUANMA/Stage/")

## Install libraries and load them
##If necessary, uncomment these section to install the packages: 
##  install.packages("gplots")
##  source("https://bioconductor.org/biocLite.R")
##    biocLite("DESeq2")
require(DESeq2)
require(RColorBrewer)
require(gplots)
library(getopt, quietly=TRUE, warn.conflicts=FALSE)

## Manage input data 
# Collect arguments from command line
args <- commandArgs(trailingOnly=TRUE)

# Get options, using the spec as defined by the enclosed list.
# Read the options from the default: commandArgs(TRUE).
spec <- matrix(c(
  "outPath", "o", 1, "character",
  "matrixPath", "m", 2, "character",
  "factInput", "i", 2, "character",
  "minReadCount", "r", 1, "integer"),
  byrow=TRUE, ncol=4)
opt <- getopt(spec)

## Remove erroneous data from file obtained from featureCounts
fc_res <- read.table(opt$matrixPath, header = TRUE)
fc_res[rowSums(fc_res[,-1])>opt$minReadCount,] -> fc_res
fc_res[!duplicated(as.list(fc_res))]-> fc_res
row.names(fc_res) <- fc_res$Geneid
fc_res<- fc_res[,-1]

Factor_gr<- unlist(strsplit(opt$factInput,","))

colData<- data.frame(row.names = colnames(fc_res),Factor_gr)

# Build model matrix 
Model<- DESeqDataSetFromMatrix(fc_res, colData=colData, design=~Factor_gr, tidy = FALSE)

## Run DESeq2 pipeline 

deseq2<- DESeq(Model)

# Get DE genes results 
res <- results(deseq2)


## Plot dispersions
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(deseq2, main="Dispersion plot")

# R-log transformation for heatmaps 
r_log <- rlogTransformation(deseq2)

mycols <- brewer.pal(8, "Dark2")[1:length(unique(Factor_gr))]
sampleDists <- as.matrix(dist(t(assay(r_log))))
png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[Factor_gr], RowSideColors=mycols[Factor_gr],
          margin=c(10, 10), main="Sample Distance Matrix")

# PCA analysis 
png("qc-pca_pack.png", 1000, 1000, pointsize=30)
plotPCA(r_log,intgroup="Factor_gr")

# Get differential expression results
res <- results(deseq2,alpha = 0.05)
table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(deseq2, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
## Write results
write.csv(resdata, file="diffexpr-results.csv")

## Examine plot of p-values
png("hist_pvalues.png",500,500)
hist(res$pvalue, breaks=50, col="grey")

# MA plot 
maplot <- function (res, thresh=0.05, labelsig=FALSE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
png("diffexpr-maplot1.png", 1500, 1000, pointsize=20)
maplot(resdata, main="MA Plot")
dev.off()

## Volcano plot 
png("volcano_plot.png", 500, 500)

with(resdata, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-10,10)))
# Add colors
with(subset(resdata, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(resdata, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(resdata, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
# Add labels to significant genes 
library(calibrate)
with(subset(resdata, padj<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=0.8))

dev.off()

## Second volcano plot
volcanoplot <- function (res, lfcthresh, sigthresh, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-10, 10))
dev.off()