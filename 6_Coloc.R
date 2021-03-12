###########################################################################
# Colocalization
###########################################################################
rm(list = ls())
# set workspace
setwd("D:/MyWork/RockyAim3/Colocalization/Data")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install('snpStats')

#install.packages('coloc');
#install.packages('snpStats');
#install.packages('readxl');

library(coloc)
library(readxl)

# GWAS
gwas <- read.table("./GWAS_summary_statistics.txt", stringsAsFactors=F, header = T)
gwas$CHR <- paste("chr",gwas$CHR, sep = '')

# eQTL
eqtl <- read_excel("./eQTL_summary_statistics.xlsx", sheet = 1)
eqtl <- eqtl[,-c(10:11)]



# coloc

# coloc
coloc <- function(chr, study=c("GWAS1","GWAS2"), sample=c("Artery Aorta","Artery Coronary","Artery Tibial","Heart Atrial Appendage","Heart Left Ventricle"),type=c("GWAS","eQTL")){
#  sink("/dev/null")
  
  # GGRF data
  mqtl <- read.table(paste('D:/MyWork/RockyAim3/Colocalization/Data20/Mqtl_snp_level_',chr,'.txt',sep=''), stringsAsFactors=F, header = T)
  region <- unique(mqtl$SNP.region)
  
  eqtl_chr <- eqtl[eqtl$chr == chr,]
  eqtl_chr <- eqtl_chr[eqtl_chr$sample == sample,]
  eqtl_chr <- eqtl_chr[order(eqtl_chr$variant_pos,eqtl_chr$pval_nominal),]      
  eqtl_chr <- eqtl_chr[!duplicated(eqtl_chr$variant_pos),];
  eqtl_chr$type <- "quant"
  eqtl_chr <- eqtl_chr[,c(1,7,13,6,8:9,14)]
        
  gwas_chr <- gwas[gwas$CHR == chr,]
  gwas_chr <- gwas_chr[gwas_chr$study == study,]
  gwas_chr <- gwas_chr[order(gwas_chr$BP),]      
  gwas_chr$type <- "cc"
  colnames(gwas_chr)[11] <- "s"
  gwas_chr <- gwas_chr[,c(1,7,10,9,12,11)]
        
  coloc <- NULL
  for (i in 1:length(region)){
    re <- region[i]
    mqtl_re <- mqtl[mqtl$SNP.region == re,]
      
    if (type == "GWAS"){
      list1 <- mqtl_re$ID[mqtl_re$ID %in% gwas_chr$SNP];
    }
    if (type == "eQTL"){
        list1 <- mqtl_re$ID[mqtl_re$ID %in% eqtl_chr$rs_id_dbSNP151_GRCh38p7]
    }
   
    if (length(list1) != 0){
      # coloc input format
      coloc_mqtl_re <- mqtl_re[,c(1,2,6,3,4,10,5)]
	coloc_mqtl_re <- coloc_mqtl_re[which(coloc_mqtl_re$ID %in% list1),]
      coloc_mqtl_re <- coloc_mqtl_re[,c(4:7,3)]
      colnames(coloc_mqtl_re) <- c("pvalues","N","MAF","type","snp")
      coloc_mqtl_re <- as.list(coloc_mqtl_re)
     
      if (type == "GWAS"){
        # coloc format
        coloc_gwas_re <- gwas_chr[which(gwas_chr$SNP %in% list1),]
        coloc_gwas_re <- coloc_gwas_re[,c(2:6,1)]
        colnames(coloc_gwas_re) <- c("pvalues","N","MAF","type","s","snp")
        coloc_gwas_re <- as.list(coloc_gwas_re)
        
        out <- coloc.abf(coloc_mqtl_re, coloc_gwas_re, p1 = 1e-4, p2 = 0.05, p12 = 1e-5)
        out <- as.data.frame(out$summary)
        colnames(out) <- "value"
        p.temp <- c(chr, re, type, study, sample, out$value[1], out$value[2], out$value[3], out$value[4], out$value[5], out$value[6])
        coloc <- rbind(coloc, p.temp)
      }
      
      if (type == "eQTL"){
        # coloc format
        eqtl_re <- eqtl_chr[which(eqtl_chr$rs_id_dbSNP151_GRCh38p7 %in% list1),]
        eqtl_re <- eqtl_re[,c(2:7,1)]
        colnames(eqtl_re) <- c("pvalues","N","MAF","beta","varbeta","type","snp")
        eqtl_re <- as.list(eqtl_re)

        out <- coloc.abf(coloc_mqtl_re, eqtl_re, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
        out <- as.data.frame(out$summary)
        colnames(out) <- "value"
        p.temp <- c(chr, re, type, study, sample,out$value[1], out$value[2], out$value[3], out$value[4], out$value[5], out$value[6])
        coloc <- rbind(coloc, p.temp)
#        closeAllConnections()
      }
    }
    else{print(paste("Warning: no overlapped SNP on ",chr,",",re," between mQTLs and ",type,"(", study,",", sample,")",sep = ""))}
  }
  colnames(coloc) <- c("Chr", "Region", "Type", "Gwas_Study", "Sample","nsnps","PP0","PP1","PP2","PP3","PP4")
  return(coloc)
}


# colocalization
# run on chr1-22
chromo <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15",
            "chr16","chr17","chr18","chr19","chr20","chr21","chr22")
for (i in 1:length(chromo)){
  chro <- chromo[i]
  gwas1 <- coloc(chr=chro,study="GWAS1",sample="Artery Aorta",type="GWAS")
  write.table(gwas1,file= paste("D:/MyWork/RockyAim3/Colocalization/OUT/Coloc_",chro,"_mQTLs_vs_GWAS1.txt",sep=''),row.names=F,col.names=T,quote=F,sep='\t')
  gwas2 <- coloc(chr=chro,study="GWAS2",sample="Artery Aorta",type="GWAS")
  write.table(gwas2,file= paste("D:/MyWork/RockyAim3/Colocalization/OUT/Coloc_",chro,"_mQTLs_vs_GWAS2.txt",sep=''),row.names=F,col.names=T,quote=F,sep='\t')
  eqtl_aa <- coloc(chr=chro,study="GWAS1",sample="Artery Aorta",type="eQTL")
  write.table(eqtl_aa,file= paste("D:/MyWork/RockyAim3/Colocalization/OUT/Coloc_",chro,"_mQTLs_vs_eqtl_aa.txt",sep=''),row.names=F,col.names=T,quote=F,sep='\t')
  eqtl_ac <- coloc(chr=chro,study="GWAS1",sample="Artery Coronary",type="eQTL")
  write.table(eqtl_ac,file= paste("D:/MyWork/RockyAim3/Colocalization/OUT/Coloc_",chro,"_mQTLs_vs_eqtl_ac.txt",sep=''),row.names=F,col.names=T,quote=F,sep='\t')
  eqtl_at <- coloc(chr=chro,study="GWAS1",sample="Artery Tibial",type="eQTL")
  write.table(eqtl_at,file= paste("D:/MyWork/RockyAim3/Colocalization/OUT/Coloc_",chro,"_mQTLs_vs_eqtl_at.txt",sep=''),row.names=F,col.names=T,quote=F,sep='\t')
  eqtl_ha <- coloc(chr=chro,study="GWAS1",sample="Heart Atrial Appendage",type="eQTL")
  write.table(eqtl_ha,file= paste("D:/MyWork/RockyAim3/Colocalization/OUT/Coloc_",chro,"_mQTLs_vs_eqtl_ha.txt",sep=''),row.names=F,col.names=T,quote=F,sep='\t')
  eqtl_hlv <- coloc(chr=chro,study="GWAS1",sample="Heart Left Ventricle",type="eQTL")
  write.table(eqtl_hlv,file= paste("D:/MyWork/RockyAim3/Colocalization/OUT/Coloc_",chro,"_mQTLs_vs_eqtl_hlv.txt",sep=''),row.names=F,col.names=T,quote=F,sep='\t')
}








# summarize results
# GWAS1
chromo <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15",
            "chr16","chr17","chr18","chr19","chr20","chr21","chr22")
Coloc_GWAS1 <- NULL
for (i in 1:length(chromo)){
  chro <- chromo[i]
  gwas1 <- read.table(paste("D:/MyWork/RockyAim3/Colocalization/OUT/Coloc_",chro,"_mQTLs_vs_GWAS1.txt",sep=''),header = T,sep='\t',stringsAsFactors=F)
  Coloc_GWAS1 <- rbind(Coloc_GWAS1, gwas1)
}
write.table(Coloc_GWAS1,"D:/MyWork/RockyAim3/Colocalization/OUT/Coloc_mQTLs_vs_GWAS1.txt",row.names=F,col.names=T,quote=F,sep='\t')

# GWAS2
Coloc_GWAS2 <- NULL
for (i in 1:length(chromo)){
  chro <- chromo[i]
  gwas2 <- read.table(paste("D:/MyWork/RockyAim3/Colocalization/OUT/Coloc_",chro,"_mQTLs_vs_GWAS2.txt",sep=''),header = T,sep='\t',stringsAsFactors=F)
  Coloc_GWAS2 <- rbind(Coloc_GWAS2, gwas2)
}
write.table(Coloc_GWAS2,"D:/MyWork/RockyAim3/Colocalization/OUT/Coloc_mQTLs_vs_GWAS2.txt",row.names=F,col.names=T,quote=F,sep='\t')

# eQTLs_aa
Coloc_eqtl_aa <- NULL
for (i in 1:length(chromo)){
  chro <- chromo[i]
  eqtls_aa <- read.table(paste("D:/MyWork/RockyAim3/Colocalization/OUT/Coloc_",chro,"_mQTLs_vs_eqtl_aa.txt",sep=''),header = T,sep='\t',stringsAsFactors=F)
  Coloc_eqtl_aa <- rbind(Coloc_eqtl_aa, eqtls_aa)
}
write.table(Coloc_eqtl_aa,"D:/MyWork/RockyAim3/Colocalization/OUT/Coloc_mQTLs_vs_eqtl_aa.txt",row.names=F,col.names=T,quote=F,sep='\t')

# eQTLs_ac
Coloc_eqtl_ac <- NULL
for (i in 1:length(chromo)){
  chro <- chromo[i]
  eqtls_ac <- read.table(paste("D:/MyWork/RockyAim3/Colocalization/OUT/Coloc_",chro,"_mQTLs_vs_eqtl_ac.txt",sep=''),header = T,sep='\t',stringsAsFactors=F)
  Coloc_eqtl_ac <- rbind(Coloc_eqtl_ac, eqtls_ac)
}
write.table(Coloc_eqtl_ac,"D:/MyWork/RockyAim3/Colocalization/OUT/Coloc_mQTLs_vs_eqtl_ac.txt",row.names=F,col.names=T,quote=F,sep='\t')

# eQTLs_at
Coloc_eqtl_at <- NULL
for (i in 1:length(chromo)){
  chro <- chromo[i]
  eqtls_at <- read.table(paste("D:/MyWork/RockyAim3/Colocalization/OUT/Coloc_",chro,"_mQTLs_vs_eqtl_at.txt",sep=''),header = T,sep='\t',stringsAsFactors=F)
  Coloc_eqtl_at <- rbind(Coloc_eqtl_at, eqtls_at)
}
write.table(Coloc_eqtl_at,"D:/MyWork/RockyAim3/Colocalization/OUT/Coloc_mQTLs_vs_eqtl_at.txt",row.names=F,col.names=T,quote=F,sep='\t')

# eQTLs_ha
Coloc_eqtl_ha <- NULL
for (i in 1:length(chromo)){
  chro <- chromo[i]
  eqtls_ha <- read.table(paste("D:/MyWork/RockyAim3/Colocalization/OUT/Coloc_",chro,"_mQTLs_vs_eqtl_ha.txt",sep=''),header = T,sep='\t',stringsAsFactors=F)
  Coloc_eqtl_ha <- rbind(Coloc_eqtl_ha, eqtls_ha)
}
write.table(Coloc_eqtl_ha,"D:/MyWork/RockyAim3/Colocalization/OUT/Coloc_mQTLs_vs_eqtl_ha.txt",row.names=F,col.names=T,quote=F,sep='\t')

# eQTLs_hlv
Coloc_eqtl_hlv <- NULL
for (i in 1:length(chromo)){
  chro <- chromo[i]
  eqtls_hlv <- read.table(paste("D:/MyWork/RockyAim3/Colocalization/OUT/Coloc_",chro,"_mQTLs_vs_eqtl_hlv.txt",sep=''),header = T,sep='\t',stringsAsFactors=F)
  Coloc_eqtl_hlv <- rbind(Coloc_eqtl_hlv, eqtls_hlv)
}
write.table(Coloc_eqtl_hlv,"D:/MyWork/RockyAim3/Colocalization/OUT/Coloc_mQTLs_vs_eqtl_hlv.txt",row.names=F,col.names=T,quote=F,sep='\t')

##########################################################################
# Selected results
##########################################################################
rm(list = ls())
# set workspace

GWAS1 <- read.table("D:/MyWork/RockyAim3/Colocalization/OUT/Coloc_mQTLs_vs_GWAS1.txt",header = T, sep='\t',stringsAsFactors=F)
GWAS2 <- read.table("D:/MyWork/RockyAim3/Colocalization/OUT/Coloc_mQTLs_vs_GWAS2.txt",header = T, sep='\t',stringsAsFactors=F)
eqtls_aa <- read.table("D:/MyWork/RockyAim3/Colocalization/OUT/Coloc_mQTLs_vs_eqtl_aa.txt",header = T, sep='\t',stringsAsFactors=F)
eqtls_ac <- read.table("D:/MyWork/RockyAim3/Colocalization/OUT/Coloc_mQTLs_vs_eqtl_ac.txt",header = T, sep='\t',stringsAsFactors=F)
eqtls_at <- read.table("D:/MyWork/RockyAim3/Colocalization/OUT/Coloc_mQTLs_vs_eqtl_at.txt",header = T, sep='\t',stringsAsFactors=F)
eqtls_ha <- read.table("D:/MyWork/RockyAim3/Colocalization/OUT/Coloc_mQTLs_vs_eqtl_ha.txt",header = T, sep='\t',stringsAsFactors=F)
eqtls_hlv <- read.table("D:/MyWork/RockyAim3/Colocalization/OUT/Coloc_mQTLs_vs_eqtl_hlv.txt",header = T, sep='\t',stringsAsFactors=F)

sig_GWAS1 <- GWAS1[GWAS1$PP4 >=0.95,]
sig_GWAS2 <- GWAS2[GWAS2$PP4 >=0.95,]
sig_eqtls_aa <- eqtls_aa[eqtls_aa$PP4 >=0.95,]
sig_eqtls_at <- eqtls_at[eqtls_at$PP4 >=0.95,]
sig_eqtls_ha <- eqtls_ha[eqtls_ha$PP4 >=0.95,]
sig_eqtls_hlv <- eqtls_hlv[eqtls_hlv$PP4 >=0.95,]

sig_coloc <- rbind(sig_GWAS1,sig_GWAS2,sig_eqtls_aa,sig_eqtls_at,sig_eqtls_ha,sig_eqtls_hlv)
write.table(sig_coloc,"D:/MyWork/RockyAim3/Colocalization/OUT/Coloc_mQTLs_selected_results.txt",row.names=F,col.names=T,quote=F,sep='\t')
