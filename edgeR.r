setwd("/media/jumagari/JUANMA/Stage")
library(edgeR)
fc_res <- read.table("Count_matrix.tabular", header = TRUE)
ext_fc<- read.table("Matrix_data", header = T)
fc_res<- cbind(Geneid=fc_res$Geneid,ext_fc)

row.names(fc_res) <- fc_res$Geneid

fc_res<- fc_res[,-c(1)]

Factor_gr<- c(rep("NonStr",3),rep("Str",3))

d<- DGEList(counts = fc_res,group = factor(Factor_gr))

## Build model matrix 

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