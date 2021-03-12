rm(list=ls());
setwd('D:/MyWork/RockyAim3/Data/AfterQC/');


####Beta value

Betaval<-read.csv('./Beta_detP_0.01_miss_CpG_0.05_miss_sub_0.05_20181010.csv',stringsAsFactors = F)
info<-read.csv('./annotation_for_beta_20181010.csv',stringsAsFactors=F)
sum(info$Name!=Betaval$site);
CpgInfo<-info[,c('Name','chr','pos')]
range<-cbind(CpgInfo[,3]-75000,CpgInfo[,3]+75000);
CpgInfo<-cbind(CpgInfo,range);

library(intervals)
RG<-NULL;
for(ch in 1:22){
  pos<-CpgInfo[,2]==paste('chr',ch,sep='');
  pos[is.na(pos)]<-FALSE;
  temp <- Intervals(CpgInfo[pos,4:5],closed = c( TRUE, TRUE ),type = "Z")
  IV<-interval_union(temp)
  IV<-cbind(ch,IV);
  RG<-rbind(RG,IV);
}

RG<-cbind(RG,paste('R',1:nrow(RG),sep=''))
RG[as.numeric(RG[,2])<0,2]<-0;
#write.table(RG,file='./CHR/Range_beta.txt',row.names=F,col.names=F,quote=F,sep='\t');

sn<-order(CpgInfo[,2],CpgInfo[,3]);
CpgInfo<-CpgInfo[sn,];
Betaval<-Betaval[sn,]
#write.table(CpgInfo,file='./CHR/SiteInfo_beta.txt',row.names=F,col.names=T,quote=F,sep='\t')
#write.table(Betaval,file='./CHR/Beta.txt',row.names=F,col.names=T,quote=F,sep='\t')

for(ch in 1:22){
  pos<-CpgInfo[,2]==paste('chr',ch,sep='');
  Betatmp<-Betaval[pos,];
  InfoTemp<-CpgInfo[pos,];
  write.table(InfoTemp,file=paste('./CHR/Info_beta_CHR',ch,'.txt',sep=''),row.names=F,col.names=T,quote=F,sep='\t')
  write.table(Betatmp,file=paste('./CHR/Beta_CHR',ch,'.txt',sep=''),row.names=F,col.names=T,quote=F,sep='\t')
}



