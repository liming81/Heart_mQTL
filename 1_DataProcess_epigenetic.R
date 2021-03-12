
# source("https://bioconductor.org/biocLite.R")
# biocLite("minfi")
# biocLite("minfiData")
# install.packages('stringi')
# biocLite("IlluminaHumanMethylationEPICmanifest")
# biocLite("bumphunter")
# biocLite("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")


rm(list=ls())
library(minfi)
setwd('d:/MyWork/RockyAim3/');
set.seed(20181001)

baseDir.baylor <- "./Data/Epigenetic/Baylor"
targets.baylor <- read.metharray.sheet(baseDir.baylor)
RGSet.baylor <- read.metharray.exp(targets = targets.baylor)
RGSet.baylor
manifest.baylor <- getManifest(RGSet.baylor)
manifest.baylor

baseDir.col <- "./Data/Epigenetic/Columbia"
targets.col <- read.metharray.sheet(baseDir.col)
info<-read.table('./Data/SampleInfo.txt',sep='\t',stringsAsFactors = F,header=T)
info<-info[info$type!='Thymus',]
targets.col<-targets.col[targets.col$Sample_Name%in%info$Sample_Name,]
RGset.col <- read.metharray.exp(targets = targets.col)
RGset.col
manifest.col <- getManifest(RGset.col)
manifest.col

RGSet<-combineArrays(RGset.col,RGSet.baylor)
RGSet
MSet<-preprocessRaw(RGSet)
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
GRSet <- mapToGenome(RSet)
snps <- getSnpInfo(GRSet)
detP <- detectionP(RGSet);

QC<-minfiQC(MSet)
plotQC(QC$qc);
plotSex(QC$object);
pSex<-QC$qc[,'predictedSex',drop=F];
qcReport(RGSet)

manifest <- getManifest(RGSet)
manifest
head(getProbeInfo(manifest))

ratioSet <- preprocessFunnorm(RGSet)


gset <- mapToGenome(ratioSet)
gset <- addSnpInfo(gset)
gset <- dropLociWithSnps(gset)

beta <- getBeta(gset)
mval<-getM(gset)
annotation<-getAnnotation(gset)
detP<-detP[rownames(beta),colnames(beta)]
beta0.01<-beta;
beta0.01[detP>0.01]<-NA;
mval0.01<-mval;
mval0.01[detP>0.01]<-NA;



ID1<-info[,1];
ID1<-gsub(" ", "", ID1)
ID2<-paste(info[,2],info[,3],sep='_');
ID<-cbind(ID1,ID2);
ID<-ID[order(ID2),];
ID1<-ID[,1];
ID2<-ID[,2];
ID2[1:50]<-paste(substring(ID2[1:50],1,nchar(ID2[1:50])-2),'01',sep='')
rownames(ID)<-ID2;
colnames(beta)<-ID[colnames(beta),1]
colnames(detP)<-ID[colnames(detP),1]
colnames(beta0.01)<-ID[colnames(beta0.01),1]
colnames(mval)<-ID[colnames(mval),1]
colnames(mval0.01)<-ID[colnames(mval0.01),1]


write.csv(cbind(site=rownames(beta),beta),file='./Data/Beta.csv',row.names=F,quote=F)
write.csv(cbind(site=rownames(mval),mval),file='./Data/Mval.csv',row.names=F,quote=F)
write.csv(cbind(site=rownames(detP),detP),file='./Data/detP.csv',row.names=F,quote=F)
write.csv(cbind(site=rownames(annotation),annotation),file='./Data/annotation.csv',row.names=F,quote=F)

write.csv(cbind(site=rownames(beta0.01),beta0.01),file='./Data/Beta_detP_0.01.csv',row.names=F,quote=F)
write.csv(cbind(site=rownames(mval0.01),mval0.01),file='./Data/Mval_detP_0.01.csv',row.names=F,quote=F)

write.csv(cbind(site=rownames(pSex),ID=ID[rownames(pSex),1],pSex),file='./Data/predictedSex_epigenetic.csv',row.names=F,quote=F)
write.table(ID,file='./Data/ID_epigenetic.txt',row.names=F,col.names=F,quote=F,sep='\t')




