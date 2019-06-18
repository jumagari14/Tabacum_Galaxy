  # Author: Juan Manuel Garcia 
  
  setwd("/media/jumagari/JUANMA/Stage/Python/Python-scripts/")
  library(plyr)
  library(ggplot2)
  
  gff_data<- read.csv("test.gff", header = T,sep = '\t')
  
  ## Histograms 
  
  TE_Chr<- count(gff_data,"Chromosome")
  names(TE_Chr)[2]<- "Occurences"
  TEName<- count(gff_data,"TE_Name")
  names(TEName)[2]<- "Occurences"
  TEName_strand<- count(gff_data,c("TE_Name","Strand"))
  names(TEName_strand)[3]<- "Occurences"
  TEName_chr<- count(gff_data,c("TE_Name","Chromosome"))
  names(TEName_chr)[3]<- "Occurences"
  TEClass<- count(gff_data,"Class")
  names(TEClass)[2]<- "Occurences"
  TEClass_chr<- count(gff_data,c("Class","Chromosome"))
  names(TEClass_chr)[3]<- "Occurences"

  SRR<- count(gff_data[gff_data[,"Code"]=="SRR",],"Code")
  names(SRR)[2]<- "Occurences"
  SRR_chr<- count(gff_data[gff_data[,"Code"]=="SRR",],"Chromosome")
  names(SRR_chr)[2]<- "Occurences"
  
  save(list = ls(.GlobalEnv,pattern = "[TE.*|SRR.*]"),file = "gff.Rdata")
  
  png("TE_Name_count.png",1000,1000)
  g1<-ggplot(gff_data,aes(x=TE_Name),stat="count")+ 
    geom_bar(fill="blue")
  g1+ylab("Occurences")
  dev.off()
  
  png("TE_Name_count_noSRR.png",1000,1000)
  g1_noSSR<-ggplot(gff_data[gff_data[,"Class"]!="SRR",],aes(x=TE_Name),stat="count")+
    geom_bar()
  g1_noSSR+ylab("Occurences")
  dev.off()
  
  png("TE_Name&Strand_count.png",1000,1000)
  g2<-ggplot(gff_data,aes(x=TE_Name, fill=Strand),stat="count")+
    geom_bar()
  g2+ylab("Occurences")
  dev.off()
  
  png("TE_Name&Strand_noSRR.png",1000,1000)
  g2_noSSR<-ggplot(gff_data[gff_data[,"Class"]!="SRR",],aes(x=TE_Name,fill=Strand),stat="count")+
    geom_bar()
  g2_noSSR+ylab("Occurences")
  dev.off()
  
  png("TE_Name&Chromosome_count.png",1000,1000)
  g3<-ggplot(gff_data,aes(x=TE_Name,fill=Chromosome),stat="count")+
    geom_bar()
  g3+ylab("Occurences")
  dev.off()
  
  png("TE_Name&Chromosome_noSRR.png",1000,1000)
  g3_noSSR<-ggplot(gff_data[gff_data[,"Class"]!="SRR",],aes(x=TE_Name,fill=Chromosome),stat="count")+
    geom_bar()
  g3_noSSR+ylab("Occurences")
  dev.off()
  
  png("TE_Class_count.png",1000,1000)
  g4<-ggplot(gff_data,aes(x=Class),stat="count")+
    geom_bar()
  g4+ylab("Occurences")
  dev.off()
  
  png("TE_Class_noSRR.png",1000,1000)
  g4_noSSR<-ggplot(gff_data[gff_data[,"Class"]!="SRR",],aes(x=Class),stat="count")+
    geom_bar()
  g4_noSSR+ylab("Occurences")
  dev.off()
  
  
  ## Get distribution across chromosomes
  length_chr<-unique(gff_data[,which(colnames(gff_data)=="Chromosome"):which(colnames(gff_data)=="Length_Chr")])
  
  chr_dis<-setNames(data.frame(matrix(ncol = nrow(length_chr),nrow = 0)),length_chr$Chromosome)
  
  library("karyoploteR")
  library("regioneR")
  
  for (i in 1:nrow(length_chr)){
   subdata<- gff_data[gff_data[,"Chromosome"]==length_chr[i,1],]
   png(paste0("TE_Density_",length_chr[i,1],".png"),1000,1000)
   TEDensity<-ggplot(subdata) +
     geom_histogram(aes(x=Start),binwidth=3000) + 
     geom_density(aes(x=Start,y=3000*..count..))+
     ggtitle("Density of TEs") +
     xlab(paste("Start position in the genome ",length_chr[i,1])) +
     ylab("TE density") +
     theme_bw() 
  dev.off()
  # interval_length<- max(subdata$Length)
  # bins<-round(length_chr[i,2]/interval_length)
  # intervals<-seq(from=1,to=length_chr[i,2],by=interval_length)
  # intervals<- append(intervals,length_chr[i,2])

  # count_vect<-vector()

  # for (j in 1:length(intervals)-1){
  #   count_sam<-subdata[subdata[,"Start"]>=intervals[j] & subdata[,"End"]<=intervals[j+1],]
  #   count_vect<-append(count_vect,nrow(count_sam))
  # }
  # distr_bar<- ggplot(data.frame(x=intervals,counts=count_vect))+geom_bar(aes(x=x, y=counts), stat="identity")
  # distr_bar+xlab(length_chr[i,1])+
  # ylab("Occurences")
  # print(distr_bar)
  # # histo<- ggplot(count_sam)+geom_histogram(aes(x=Start),binwidth = )
  # # print(histo)
  #  #chr_dis$length=count_vect

  ## Add rainfall plot to show density 
  test<- toGRanges(gff_data[gff_data[,"Chromosome"]==length_chr[i,1],c("Chromosome","Start","End")])
  custom.genome<- toGRanges(data.frame(Chromosome=length_chr[i,1],start=1,end=length_chr[i,2]))
  
  png(paste0("Rainfall_",length_chr[i,1],".png"),1000,1000)
  kp <- plotKaryotype(plot.type=1,genome = custom.genome)
  kpAddMainTitle(kp, main=paste0("Density of all TE elements (SRR included) across ",length_chr[i,1]), cex=1)
  kpPlotRainfall(kp, data = test)
  kpPlotDensity(kp, data = test, window.size = 3000, r0=0.62, r1=0.8)
  kpAxis(kp, ymax = 7, tick.pos = 1:7, r0=0, r1=0.6)
  kpAddLabels(kp, labels = c("Distance between TEs (log10)"), srt=90, pos=1, label.margin = 0.04, r0=0, r1=0.6)
  kpAddLabels(kp, labels = "Density", srt=90, pos=1, label.margin = 0.04, r0=0.62, r1=1)
  dev.off()
  # ANalysis of SRR 
  test<- toGRanges(gff_data[gff_data[,"Chromosome"]==length_chr[i,1] & gff_data[,"Code"]=="SRR",c("Chromosome","Start","End")])
  custom.genome<- toGRanges(data.frame(Chromosome=length_chr[i,1],start=1,end=length_chr[i,2]))
  
  png(paste0("Rainfall_SRR",length_chr[i,1],".png"),1000,1000)
  kp <- plotKaryotype(plot.type=1,genome = custom.genome)
  kpAddMainTitle(kp, main=paste0("Density of SRR across ",length_chr[i,1]), cex=1)
  kpPlotRainfall(kp, data = test)
  kpPlotDensity(kp, data = test, window.size = 3000, r0=0.62, r1=0.8)
  kpAxis(kp, ymax = 7, tick.pos = 1:7, r0=0, r1=0.6)
  kpAddLabels(kp, labels = c("Distance between SRR (log10)"), srt=90, pos=1, label.margin = 0.04, r0=0, r1=0.6)
  kpAddLabels(kp, labels = "Density", srt=90, pos=1, label.margin = 0.04, r0=0.62, r1=1)
  dev.off()
  }

  zip("png_figures",list.files(pattern = "\\.png$"))  
  unlink(list.files(pattern = "\\.png$"),recursive = TRUE)
  
