rm(list=ls());


setwd('./MQTL/');

genofile<-paste('./CHR/GWAS_CHR_',1:22,'.tped',sep='');
famfile<-paste('./CHR/GWAS_CHR_',1:22,'.tfam',sep='')
Infofile<-paste('./CHR/Info_beta_CHR',1:22,'.txt',sep='');
Mfile<-paste('./CHR/Beta_CHR',1:22,'.txt',sep='');
COVfile<-'./CHR/COV_beta.csv';
outfile<-paste('./OUT/Beta_CHR_',1:22,'.txt',sep='');

ch=6;
GT<-read.table(genofile[ch],stringsAsFactors=F);
SNP<-GT[,1:4];
GT<-GT[,-(1:4)];
sn1<-seq(1,ncol(GT),by=2);
sn2<-seq(2,ncol(GT),by=2);
GT<-as.matrix(2-GT[,sn1])+as.matrix(2-GT[,sn2]);
GT[GT==4]<-NA;
TAB<-cbind(rowSums(GT==0,na.rm=T),rowSums(GT==1,na.rm=T),rowSums(GT==2,na.rm=T));
pos<-rowSums(TAB<3)==0	
TAB<-TAB[pos,];
GT<-GT[pos,];
SNP<-SNP[pos,];

info<-read.table(Infofile[ch],header=T,stringsAsFactors=F);
MDt<-read.table(Mfile[ch],header=T);		
pos<-!is.na(info[,2]);
info<-info[pos,]
MDt<-MDt[pos,];
COV<-read.csv(COVfile,stringsAsFactors = F);
rownames(COV)<-gsub('Hobbs','HOBBS',COV[,2]);
fam<-read.table(famfile[ch],header=F,stringsAsFactors = F)
COV<-COV[fam[,1],]


Mval<-MDt
names(Mval)<-gsub('Hobbs','HOBBS',names(Mval))
names(Mval)<-gsub('X','',names(Mval))
Mval<-Mval[,fam[,1]]

PC.geno<-read.table('./CHR/plink.eigenvec',stringsAsFactors = F)

for(i in 1:nrow(info)){
  pos<-abs(SNP[,4]-info[i,3])<75000 & abs(SNP[,4]-info[i,3])>10;
  if(sum(pos)>0){
    y<-as.numeric(Mval[i,]);
    gt<-GT[pos,,drop=F];
    for(j in 1:nrow(gt)){
      snp<-gt[j,];
      fit<-lm(y ~ snp + factor(COV[,3]) + factor(COV[,4]) + as.matrix(COV[,5:9]) + as.matrix(PC.geno[,3:7]));		
      out<-c(unlist(info[i,1:3]),unlist(SNP[pos,c(2,1,4),drop=F][j,]),TAB[pos,,drop=F][j,],summary(fit)$coef[2,],summary(fit)$r.squared,summary(fit)$adj.r.squared)			
      write.table(matrix(out,nrow=1),file=outfile[ch],row.names=F,col.names=F,quote=F,sep='\t',append=T);
    }
  }
}

