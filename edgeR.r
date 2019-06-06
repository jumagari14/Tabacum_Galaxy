## Author: Juan Manuel Garcia

setwd("/media/jumagari/JUANMA/Stage/")
## Install libraries and load them
require(edgeR)


## Remove erroneous data from count matrix
fc_res <- read.table("Matrix_data.tabular", header = TRUE)
fc_res[!duplicated(as.list(fc_res))]-> fc_res

row.names(fc_res) <- fc_res$Geneid

fc_res<- fc_res[,-c(1)]

Factor_gr<- c(rep("NonStr",3),rep("Str",3))

## Build model matrix

d<- DGEList(counts = fc_res,group = factor(Factor_gr))

myCPM<- cpm(fc_res)

first_count<- log2(d$counts+1)

plotMD(first_count)

## Estimate dispersion 

d1<- DGEList(d$counts[apply(d$counts,1,sum)!=0,],group = d$samples$group)

dFull <- calcNormFactors(d1, method="TMM")

final_deg<- estimateGLMCommonDisp(dFull)

final_deg<- estimateGLMTrendedDisp(final_deg)

final_deg<- estimateGLMTagwiseDisp(final_deg)

## Differential expression analysis 
deg_test<- exactTest(final_deg)

results<- topTags(deg_test,n=nrow(deg_test$table))

results$table$logFC[results$table$FDR<0.01]