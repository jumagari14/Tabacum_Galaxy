# Author: Juan Manuel Garcia 

setwd("/media/jumagari/JUANMA/Stage/Python/Python-scripts/")
require(RColorBrewer)

#gff_data<- read.csv("test.gff", header = T,sep = '\t')
gff_data<-gff_data[gff_data[,"Class"]!="SRR",]

length_chr<-unique(gff_data[,which(colnames(gff_data)=="Chromosome"):which(colnames(gff_data)=="Length_Chr")])

name_TE<- unique(gff_data[,"TE_Name"])

per_chr<-data.frame(matrix(nrow = length(name_TE),ncol = nrow(length_chr)))
colnames(per_chr)<-length_chr[,1]

for (i in 1:nrow(length_chr)){
    name<-length_chr[i,1]
    subdata<- gff_data[gff_data[,"Chromosome"]==length_chr[i,1],]
    coverage<-list()
    
    for (j in 1:length(name_TE)){
      subdata_2=subdata[subdata[,"TE_Name"]==name_TE[j],]
      percent=round(Reduce("+",subdata_2[,"Length"])/unique(subdata_2[,"Length_Chr"]),4)*100
      if (nrow(subdata_2)==0){
        percent=0
      }
      coverage<- append(coverage,percent)
    }
  per_chr[[name]]<- coverage
}

per_chr<-as.matrix(sapply(per_chr,as.numeric))



n <- length(name_TE)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

colors<-sample(col_vector,n)
barplot(per_chr,col =colors,ylim = c(0,100),cex.names = 0.4,xlab = "Chromosome",ylab = "Coverage (%)")
legend("topleft",legend = name_TE,fill =colors,ncol=3,cex = 0.4)
 

## 
length_chr<-length_chr[,c(2,1)]

length_chr_2<-as.matrix(sapply(length_chr[,1],as.numeric))

matrix<-matrix<-per_chr %*% length_chr_2
matrix<-matrix/sum(length_chr_2)
matrix<-t(matrix)
colnames(matrix)<-name_TE

barplot(matrix,xlab = "TE Name",ylab = "Coverage (%)")

