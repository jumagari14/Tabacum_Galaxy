## Author: Juan Manuel Garcia
## Based on the workflow from https://gist.github.com/stephenturner/f60c1934405c127f09a6

setwd("/media/jumagari/JUANMA/Stage/")

## Install libraries and load them
##If necessary, uncomment these seciton to install the packages: 
##  install.packages("gplots")
##  source("https://bioconductor.org/biocLite.R")
##    biocLite("DESeq2")
require(DESeq2)
require(data.table)
require(RColorBrewer)

## Remove erroneous data from file obtained from featureCounts
fc_res <- read.table("Matrix_data.tabular", header = TRUE)
fc_res[!duplicated(as.list(fc_res))]-> fc_res
row.names(fc_res) <- fc_res$Geneid
fc_res<- fc_res[,-c(1)]

Factor_gr<- c(rep("NonStr",3),rep("Str",3))
colData<- data.frame(row.names = colnames(fc_res),Factor_gr)

# Build model matrix 
Model<- DESeqDataSetFromMatrix(fc_res, colData=colData, design=~Factor_gr, tidy = FALSE)

## Run DESeq2 pipeline 

deseq2<- DESeq(Model)

# Get DE genes results 
res <- results(deseq2)


## Plot dispersions 
plotDispEsts(deseq2, main="Dispersion plot")

# R-log transformation for heatmaps 
r_log <- rlogTransformation(deseq2)

mycols <- brewer.pal(8, "Dark2")[1:length(unique(Factor_gr))]
sampleDists <- as.matrix(dist(t(assay(r_log))))
require(gplots)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[Factor_gr], RowSideColors=mycols[Factor_gr],
          margin=c(10, 10), main="Sample Distance Matrix")

