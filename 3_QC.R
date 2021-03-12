rm(list=ls());
setwd('D:/MyWork/RockyAim3/Data');

#system('plink --bfile GWAS_Heart_20180930 --check-sex --noweb');
#system('plink --bfile GWAS_Heart_20180930 --missing --noweb');



sex.genetic<-read.table('./plink.sexcheck',header=T,stringsAsFactors = F)
sex.epigenetic<-read.csv('./predictedSex_epigenetic.csv',header=T)
rownames(sex.genetic)<-sex.genetic[,1];
rownames(sex.epigenetic)<-sex.epigenetic[,2];
sex.genetic<-sex.genetic[order(sex.genetic[,1]),]
sex.epigenetic<-sex.epigenetic[order(sex.epigenetic[,2]),]
table(sex.genetic[,4],sex.epigenetic[,3])


imiss<-read.table('./plink.imiss',header=T,stringsAsFactors = F)
lmiss<-read.table('./plink.lmiss',header=T,stringsAsFactors = F)
hist(imiss$F_MISS,xlab='Missing Rate',main='Histogram by subjects');
hist(lmiss$F_MISS,xlab='Missing Rate',main='Histogram by SNPs');


#system('plink --bfile GWAS_Heart_20180930 --mind 0.05 --geno 0.05 --make-bed --out GWAS_Heart_miss_sub_0.05_miss_snp_0.05_20181001 --noweb');
sub.rm<-imiss[imiss$F_MISS>0.05,1]
sex.genetic<-sex.genetic[!sex.genetic$FID%in%sub.rm,]
sex.epigenetic<-sex.epigenetic[!sex.epigenetic$ID%in%sub.rm,]
table(sex.genetic[,4],sex.epigenetic[,3])


####Checking beta values


betaval<-read.csv('./Beta_detP_0.01.csv',header=T,stringsAsFactors = F);
rownames(betaval)<-betaval[,1];
betaval<-betaval[,-1]


mis.sub<-colMeans(is.na(betaval))
mis.CpG<-rowMeans(is.na(betaval))
hist(mis.CpG,xlab='Missing Rate',main='Histogram by CpG');
hist(mis.sub,xlab='Missing Rate',main='Histogram by subjects');


imiss<-read.table('./plink.imiss',header=T,stringsAsFactors = F)
sub.rm<-imiss[imiss$F_MISS>0.05,1]
betaval2<-betaval[mis.CpG<0.05,!names(betaval)%in%paste('X',sub.rm,sep='')]


betaval3<-betaval2;
betaave<-rowMeans(betaval2,na.rm=T)
for(i in 1:ncol(betaval3)){
  pos<-is.na(betaval3[,i]);
  betaval3[pos,i]<-betaave[pos];
}

betaval.pca<-prcomp(t(betaval3), center = TRUE,scale. = TRUE)

id<-names(betaval3);
id<-gsub("X",'',id)


library(ggfortify)
autoplot(betaval.pca,shape=F)

pca_test=betaval.pca
pca_test$x=pca_test$x[,c(3,4)]
colnames(pca_test$x)=c("PC1","PC2")
pca_test$rotation=pca_test$rotation[,c(3,4)]
colnames(pca_test$rotation)=c("PC1","PC2")
autoplot(pca_test,shape=F)+xlab("PC3")+ylab('PC4')

pca_test=betaval.pca
pca_test$x=pca_test$x[,c(5,6)]
colnames(pca_test$x)=c("PC1","PC2")
pca_test$rotation=pca_test$rotation[,c(5,6)]
colnames(pca_test$rotation)=c("PC1","PC2")
autoplot(pca_test,shape=F)+xlab("PC5")+ylab('PC6');

betaval4<-betaval2[,names(betaval2)!='X2560'];
#write.csv(cbind(site=rownames(betaval4),betaval4),file='./AfterQC/Beta_detP_0.01_miss_CpG_0.05_miss_sub_0.05_20181010.csv',row.names=F,quote=F)


betaval5<-betaval4;
betaave<-rowMeans(betaval5,na.rm=T)
for(i in 1:ncol(betaval5)){
  pos<-is.na(betaval5[,i]);
  betaval5[pos,i]<-betaave[pos];
}

betaval.pca2<-prcomp(t(betaval5), center = TRUE,scale. = TRUE)
autoplot(betaval.pca2,shape=F,xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))
PC<-predict(betaval.pca2)
rownames(PC)<-sub('X','',rownames(PC))

fam<-read.table('./GWAS_Heart_miss_sub_0.05_miss_snp_0.05_sub83_20181001.fam',stringsAsFactors = F, header=F)
rownames(fam)<-gsub('HOBBS','Hobbs',fam[,1])
sex.epigenetic<-read.csv('./predictedSex_epigenetic.csv',header=T)
rownames(sex.epigenetic)<-sex.epigenetic[,2];

COV<-cbind(sex.epigenetic[rownames(PC),],disease=fam[rownames(PC),6],PC)
#write.csv(COV,file='./AfterQC/COV_beta.csv',row.names=F,quote=F)


annotation<-read.csv('./annotation.csv',header=T,stringsAsFactors = F);
rownames(annotation)<-annotation[,1];
annotation2<-annotation[rownames(betaval4),]
#write.csv(annotation2,file='./AfterQC/annotation_for_beta_20181010.csv',row.names=F,quote=F)
