rm(list=ls())

textfiles<-dir('d:/MyWork/BirthDefect/BaylorSampleGenotype/Baylor/','.txt',full.n=T)
textfiles<-textfiles[1:51];
manifest<-read.csv('I:/RockyData/Aim3/Baylor Genotyping Data 5-30-2018/InfiniumOmni5Exome-4v1-3_SampleSheet_Baylor Heart 1-50_5-30-2018.csv',skip=10,header=T)
rownames(manifest)<-paste(manifest[,2],manifest[,3],sep='_');
manifest$ID<-paste('HOBBS',manifest$Sample_ID,sep='')
  
i=1


for(i in 1:length(textfiles)){
  raw<-read.table(textfiles[i],header=T,stringsAsFactors = F,sep='\t',skip=10)
  ID<-manifest[raw$Sample.ID[1],]$ID
  LGEN<-cbind(ID,ID,raw[,c('SNP.Name',"Allele1...Top",'Allele2...Top')]);
  fam<-c(ID,ID,0,0,0,2);
  map<-cbind(raw[,c('Chr','SNP.Name')],0,raw[,'Position'])
  write.table(LGEN,file=paste('D:/MyWork/BirthDefect/BaylorSampleGenotype/LGEN/',ID,'.lgen',sep=''),row.names=F,col.names=F,quote=F,sep=' ')
  write.table(matrix(fam,nrow=1),file=paste('D:/MyWork/BirthDefect/BaylorSampleGenotype/LGEN/',ID,'.fam',sep=''),row.names=F,col.names=F,quote=F,sep=' ')
  write.table(map,file=paste('D:/MyWork/BirthDefect/BaylorSampleGenotype/LGEN/',ID,'.map',sep=''),row.names=F,col.names=F,quote=F,sep=' ')
}



setwd('D:/MyWork/BirthDefect/BaylorSampleGenotype/LGEN/');
for(i in 1:50){
  cmd<-paste('plink --lgen HOBBS',i,'.lgen --fam HOBBS',i,'.fam --map HOBBS',i,'.map --out HOBBS',i, ' --make-bed --noweb',sep='');
  system(cmd)
}

MergeList<-paste('./HOBBS',2:50,'.bed',' ./HOBBS',2:50,'.bim',' ./HOBBS',2:50,'.fam',sep='')
write.table(MergeList,file='MergeList.txt',row.names=F,col.names=F,quote=F,sep='\n')

system('plink --bfile HOBBS1 --merge-list MergeList.txt --make-bed --out GWAS_Baylor_20180929 --noweb')

bim1<-read.table('HOBBS1.bim',header=F,stringsAsFactors = F)
bim2<-read.table('HOBBS2.bim',header=F,stringsAsFactors = F)

missnp<-read.table('GWAS_Baylor_20180929.missnp',header=F,stringsAsFactors = F)
tmp1<-bim1[bim1[,2]%in%missnp[,1],]
tmp2<-bim2[bim2[,2]%in%missnp[,1],]


bimfiles<-dir(,'bim',full.names = T)


for(i in 1:length(bimfiles)){
  bim<-read.table(bimfiles[i],header=F,stringsAsFactors = F);
  bim[bim=='-']<-'0';
  write.table(bim,file=bimfiles[i],row.names=F,col.names=F,quote=F,sep=' ');
}



####genomic position between two datasets

bim1<-read.table('d:/MyWork/RockyAim3/Data/Genetic/Baylor/GWAS_Baylor_20180929.bim',header=F,stringsAsFactors = F)
bim2<-read.table('d:/MyWork/RockyAim3/Data/Genetic/Columbia/GWAS20151119.bim',header=F,stringsAsFactors = F)

rownames(bim1)<-bim1[,2];
rownames(bim2)<-bim2[,2];
overlapsnp<-intersect(bim1[,2],bim2[,2])
tmp1<-bim1[overlapsnp,];
tmp2<-bim2[overlapsnp,]


pos<-tmp1[,4]!=tmp2[,4] & tmp1[,4]!=0 & tmp2[,4]!=0
sum(pos);
tmp<-cbind(tmp1[pos,],tmp2[pos,])

pos<-tmp1[,1]!=tmp2[,1] & tmp1[,1]!=0 & tmp2[,1]!=0
tmp<-cbind(tmp1[pos,],tmp2[pos,])

bim2[overlapsnp,1:4]<-bim1[overlapsnp,1:4]
write.table(bim2,file='d:/MyWork/RockyAim3/Data/Genetic/Columbia/GWAS_Columbia_20151119.bim',row.names=F,col.names=F,quote=F,sep=' ');


