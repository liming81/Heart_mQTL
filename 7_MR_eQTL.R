rm(list = ls())
# set workspace
setwd("/Users/chenlyu/Documents/LC/GGRF/Data")

library(readxl)
library(lettercase)
library(utils)
library(KRIS)
library(tools)
library(dplyr)

# MR-eQTL
# data input
eqtl <- read_excel("./colocalization/eQTL_summary_statistics_updated2.xlsx", sheet = 1)
eqtl <- eqtl[,-c(14:15)]
colnames(eqtl)[5] <-"SNP"

gwas <- read.table("./colocalization/GWAS_summary_statistics.txt", header = T)
gwas$CHR <- paste("chr",gwas$CHR, sep = '')
GWAS1 <- gwas[gwas$study == "GWAS1",]
GWAS2 <- gwas[gwas$study == "GWAS2",]

# overlapped SNPs
snp_eQTL <- as.character(unique(eqtl$SNP))
# 102358
snp_overlap <- snp_eQTL[snp_eQTL %in% GWAS1$SNP & (snp_eQTL %in% GWAS2$SNP)]
# 8283
write.table(as.matrix(snp_overlap,ncol=1), file= "./MR/plink/snp_overlap_eQTL.txt",row.names=F,col.names=F,quote=F,sep='\t')

eQTL <- eqtl[eqtl$SNP %in% snp_overlap,]
# 116339
write.table(eQTL, file= "./MR/eQTLs_overlap.txt",row.names=F,col.names=T,quote=T,sep='\t')
GWAS1 <- GWAS1[GWAS1$SNP %in% snp_overlap,]
# 8283
write.table(GWAS1, file= "./MR/GWAS1_overlap_eQTL.txt",row.names=F,col.names=T,quote=T,sep='\t')
GWAS2 <- GWAS2[GWAS2$SNP %in% snp_overlap,]
#8283
write.table(GWAS2, file= "./MR/GWAS2_overlap_eQTL.txt",row.names=F,col.names=T,quote=T,sep='\t')

# get GWAS bfile with overlapped snps
# plink
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile ./GWAS_new/MRGWAS_Heart_miss_sub_0.05_miss_snp_0.05_sub83_20181001')
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile ./GWAS_new/GWAS_Heart_miss_sub_0.05_miss_snp_0.05_sub83_20181001 --make-bed --extract ./MR/plink/snp_overlap_eQTL.txt --allow-no-sex --out ./MR/plink/GWAS_overlap_eQTL')

# double check
bed <- "./MR/plink/GWAS_overlap_eQTL.bed"
fam <- "./MR/plink/GWAS_overlap_eQTL.fam"
bim <- "./MR/plink/GWAS_overlap_eQTL.bim"
chr_snp <- read.bed(bed, bim, fam)
# 83 obs * 8274 snps

# make haploview format
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile ./MR/plink/GWAS_overlap_eQTL --recodeHV --chr 1 --allow-no-sex --out ./MR/HaploView/eQTL_haploview_1')
# 938
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile ./MR/plink/GWAS_overlap_eQTL --recodeHV --chr 2 --allow-no-sex --out ./MR/HaploView/eQTL_haploview_2')
# 558
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile ./MR/plink/GWAS_overlap_eQTL --recodeHV --chr 3 --allow-no-sex --out ./MR/HaploView/eQTL_haploview_3')
# 508
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile ./MR/plink/GWAS_overlap_eQTL --recodeHV --chr 4 --allow-no-sex --out ./MR/HaploView/eQTL_haploview_4')
# 296
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile ./MR/plink/GWAS_overlap_eQTL --recodeHV --chr 5 --allow-no-sex --out ./MR/HaploView/eQTL_haploview_5')
# 452
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile ./MR/plink/GWAS_overlap_eQTL --recodeHV --chr 6 --allow-no-sex --out ./MR/HaploView/eQTL_haploview_6')
# 431
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile ./MR/plink/GWAS_overlap_eQTL --recodeHV --chr 7 --allow-no-sex --out ./MR/HaploView/eQTL_haploview_7')
# 348
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile ./MR/plink/GWAS_overlap_eQTL --recodeHV --chr 8 --allow-no-sex --out ./MR/HaploView/eQTL_haploview_8')
# 319
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile ./MR/plink/GWAS_overlap_eQTL --recodeHV --chr 9 --allow-no-sex --out ./MR/HaploView/eQTL_haploview_9')
# 304
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile ./MR/plink/GWAS_overlap_eQTL --recodeHV --chr 10 --allow-no-sex --out ./MR/HaploView/eQTL_haploview_10')
# 355
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile ./MR/plink/GWAS_overlap_eQTL --recodeHV --chr 11 --allow-no-sex --out ./MR/HaploView/eQTL_haploview_11')
# 467
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile ./MR/plink/GWAS_overlap_eQTL --recodeHV --chr 12 --allow-no-sex --out ./MR/HaploView/eQTL_haploview_12')
# 525
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile ./MR/plink/GWAS_overlap_eQTL --recodeHV --chr 13 --allow-no-sex --out ./MR/HaploView/eQTL_haploview_13')
# 179
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile ./MR/plink/GWAS_overlap_eQTL --recodeHV --chr 14 --allow-no-sex --out ./MR/HaploView/eQTL_haploview_14')
# 277
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile ./MR/plink/GWAS_overlap_eQTL --recodeHV --chr 15 --allow-no-sex --out ./MR/HaploView/eQTL_haploview_15')
# 287
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile ./MR/plink/GWAS_overlap_eQTL --recodeHV --chr 16 --allow-no-sex --out ./MR/HaploView/eQTL_haploview_16')
# 288
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile ./MR/plink/GWAS_overlap_eQTL --recodeHV --chr 17 --allow-no-sex --out ./MR/HaploView/eQTL_haploview_17')
# 430
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile ./MR/plink/GWAS_overlap_eQTL --recodeHV --chr 18 --allow-no-sex --out ./MR/HaploView/eQTL_haploview_18')
# 164
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile ./MR/plink/GWAS_overlap_eQTL --recodeHV --chr 19 --allow-no-sex --out ./MR/HaploView/eQTL_haploview_19')
# 403
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile ./MR/plink/GWAS_overlap_eQTL --recodeHV --chr 20 --allow-no-sex --out ./MR/HaploView/eQTL_haploview_20')
# 202
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile ./MR/plink/GWAS_overlap_eQTL --recodeHV --chr 21 --allow-no-sex --out ./MR/HaploView/eQTL_haploview_21')
# 104
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile ./MR/plink/GWAS_overlap_eQTL --recodeHV --chr 22 --allow-no-sex --out ./MR/HaploView/eQTL_haploview_22')
# 169

# check LD r^2 to make sure independent SNPs for each CpG
# get tag SNP for each CpG

# run on terminal
# cd "Documents/LC/Horton/HaploView"
# java -jar Haploview.jar
# tagger; pairwise tagging; r^2 threshold = 0.8; LOD threshold = 3

# get tag SNP list
# 7731 for aggressive tagging (2-3 markers)
tagSNP <- NULL
for (i in 1:22){
  list <- read.table(paste("./MR/tagSNP/tagSNP_eQTL_chr", i, sep = ''))
  tagSNP <- rbind(tagSNP, list)
}
write.table(tagSNP, file= "./MR/tagSNP/tagSNPlist_eQTL.txt",row.names=F,col.names=F,quote=T,sep='\t')

# get MR datasets
tagSNP<- read.table("./MR/tagSNP/tagSNPlist_eQTL.txt",header = F)
tagSNP <- as.character(tagSNP[,1])

GWAS1 <- read.table("./MR/GWAS1_overlap_eQTL.txt", header = T)
GWAS1 <- GWAS1[GWAS1$SNP %in% tagSNP,]
write.table(GWAS1, file= "./MR/GWAS1_eQTL_MR.txt",row.names=F,col.names=T,quote=T,sep='\t')

GWAS2 <- read.table("./MR/GWAS2_overlap_eQTL.txt", header = T)
GWAS2 <- GWAS2[GWAS2$SNP %in% tagSNP,]
write.table(GWAS2, file= "./MR/GWAS2_eQTL_MR.txt",row.names=F,col.names=T,quote=T,sep='\t')

eQTL <- read.table("./MR/eQTLs_overlap.txt", header = T)
eQTL_MR <- eQTL[eQTL$SNP %in% tagSNP,]
# 7731 SNPs, 8909 observations, 6942 Genes
write.table(eQTL_MR, file= "./MR/eQTLs_MR.txt",row.names=F,col.names=T,quote=T,sep='\t')

################################################################
# Mendelian Randomization
################################################################
rm(list = ls())
# set workspace
setwd("/Users/chenlyu/Documents/LC/GGRF/Data")

library(readxl)
library(lettercase)
library(utils)
library(KRIS)
library(tools)
library(dplyr)
library(TwoSampleMR)
library(MRInstruments)
library(ggplot2)
library(ggrepel)

# data formatting
# load MR datasets
eQTL_MR <- read.table("./MR/eQTLs_MR.txt", header = T)
GWAS1_MR <- read.table("./MR/GWAS1_eQTL_MR.txt", header = T)
GWAS2_MR <- read.table("./MR/GWAS2_eQTL_MR.txt", header = T)

# exposure datasets
eQTL_MR <- eQTL_MR[,c(5,2,9,8,10:13,16:17)]
colnames(eQTL_MR)[c(2:8,10)] <- c("exposure","effect_allele","other_allele","eaf","pval","beta","se","samplesize")
eQTL_MR$unit <- "beta value"
colnames(eQTL_MR)[c(3:8,10:11)] <- paste(colnames(eQTL_MR)[c(3:8,10:11)] , ".exposure", sep = '')
eQTL_MR_aa <- eQTL_MR[eQTL_MR$sample == "Artery Aorta",c(1:8,10:11)]
eQTL_MR_ac <- eQTL_MR[eQTL_MR$sample == "Artery Coronary",c(1:8,10:11)]
eQTL_MR_at <- eQTL_MR[eQTL_MR$sample == "Artery Tibial",c(1:8,10:11)]
eQTL_MR_ha <- eQTL_MR[eQTL_MR$sample == "Heart Atrial Appendage",c(1:8,10:11)]
eQTL_MR_hlv <- eQTL_MR[eQTL_MR$sample == "Heart Left Ventricle",c(1:8,10:11)]
write.table(eQTL_MR, file= "./MR/eQTLs_MR_format.txt",row.names=F,col.names=T,quote=T,sep='\t')
write.table(eQTL_MR_aa, file= "./MR/eQTLs_MR_aa_format.txt",row.names=F,col.names=T,quote=T,sep='\t')
write.table(eQTL_MR_ac, file= "./MR/eQTLs_MR_ac_format.txt",row.names=F,col.names=T,quote=T,sep='\t')
write.table(eQTL_MR_at, file= "./MR/eQTLs_MR_at_format.txt",row.names=F,col.names=T,quote=T,sep='\t')
write.table(eQTL_MR_ha, file= "./MR/eQTLs_MR_ha_format.txt",row.names=F,col.names=T,quote=T,sep='\t')
write.table(eQTL_MR_hlv, file= "./MR/eQTLs_MR_hlv_format.txt",row.names=F,col.names=T,quote=T,sep='\t')

# outcome datasets
GWAS1_MR$beta <- log(GWAS1_MR$OR)
GWAS1_MR$unit <- "log odds"
GWAS1_MR$zscore <- qnorm(1-GWAS1_MR$P/2)
GWAS1_MR$se <- abs(GWAS1_MR$beta/GWAS1_MR$zscore)
colnames(GWAS1_MR)[c(4:5,7,9:10)] <- c("effect_allele","other_allele","pval","eaf","samplesize")
GWAS1_MR$ncase <- 2470
GWAS1_MR$ncontrol <- 4020 - 2470
GWAS1_MR <- GWAS1_MR[,c(1,12,15,7,13,4:5,9:10,16:17)]
colnames(GWAS1_MR)[c(2:11)] <- paste(colnames(GWAS1_MR)[c(2:11)] , ".outcome", sep = '')
GWAS1_MR$outcome <- "CHD"
write.table(GWAS1_MR, file= "./MR/GWAS1_eQTL_MR_format.txt",row.names=F,col.names=T,quote=T,sep='\t')

GWAS2_MR$beta <- log(GWAS2_MR$OR)
GWAS2_MR$unit <- "log odds"
GWAS2_MR$zscore <- qnorm(1-GWAS2_MR$P/2)
GWAS2_MR$se <- abs(GWAS2_MR$beta/GWAS2_MR$zscore)
colnames(GWAS2_MR)[c(4:5,7,9:10)] <- c("effect_allele","other_allele","pval","eaf","samplesize")
GWAS2_MR$ncase <- 1629
GWAS2_MR$ncontrol <- 2692 - 1629
GWAS2_MR <- GWAS2_MR[,c(1,12,15,7,13,4:5,9:10,16:17)]
colnames(GWAS2_MR)[c(2:11)] <- paste(colnames(GWAS2_MR)[c(2:11)] , ".outcome", sep = '')
GWAS2_MR$outcome <- "CHD"
write.table(GWAS2_MR, file= "./MR/GWAS2_eQTL_MR_format.txt",row.names=F,col.names=T,quote=T,sep='\t')

# MR
# load MR datasets
eQTL_MR_aa <- read_exposure_data("./MR/eQTLs_MR_aa_format.txt", sep='\t', phenotype_col = "exposure",
                              snp_col = "SNP",
                              beta_col = "beta.exposure",
                              se_col = "se.exposure",
                              eaf_col = "eaf.exposure",
                              effect_allele_col = "effect_allele.exposure",
                              other_allele_col = "other_allele.exposure",
                              pval_col = "pval.exposure",
                              units_col = "units.exposure",
                              samplesize_col = "samplesize")
eQTL_MR_ac <- read_exposure_data("./MR/eQTLs_MR_ac_format.txt", sep='\t', phenotype_col = "exposure",
                                 snp_col = "SNP",
                                 beta_col = "beta.exposure",
                                 se_col = "se.exposure",
                                 eaf_col = "eaf.exposure",
                                 effect_allele_col = "effect_allele.exposure",
                                 other_allele_col = "other_allele.exposure",
                                 pval_col = "pval.exposure",
                                 units_col = "units.exposure",
                                 samplesize_col = "samplesize")
eQTL_MR_at <- read_exposure_data("./MR/eQTLs_MR_at_format.txt", sep='\t', phenotype_col = "exposure",
                                 snp_col = "SNP",
                                 beta_col = "beta.exposure",
                                 se_col = "se.exposure",
                                 eaf_col = "eaf.exposure",
                                 effect_allele_col = "effect_allele.exposure",
                                 other_allele_col = "other_allele.exposure",
                                 pval_col = "pval.exposure",
                                 units_col = "units.exposure",
                                 samplesize_col = "samplesize")
eQTL_MR_ha <- read_exposure_data("./MR/eQTLs_MR_ha_format.txt", sep='\t', phenotype_col = "exposure",
                                 snp_col = "SNP",
                                 beta_col = "beta.exposure",
                                 se_col = "se.exposure",
                                 eaf_col = "eaf.exposure",
                                 effect_allele_col = "effect_allele.exposure",
                                 other_allele_col = "other_allele.exposure",
                                 pval_col = "pval.exposure",
                                 units_col = "units.exposure",
                                 samplesize_col = "samplesize")
eQTL_MR_hlv <- read_exposure_data("./MR/eQTLs_MR_hlv_format.txt", sep='\t', phenotype_col = "exposure",
                                 snp_col = "SNP",
                                 beta_col = "beta.exposure",
                                 se_col = "se.exposure",
                                 eaf_col = "eaf.exposure",
                                 effect_allele_col = "effect_allele.exposure",
                                 other_allele_col = "other_allele.exposure",
                                 pval_col = "pval.exposure",
                                 units_col = "units.exposure",
                                 samplesize_col = "samplesize")
GWAS1_MR <- read_outcome_data("./MR/GWAS1_eQTL_MR_format.txt", sep='\t', phenotype_col = "outcome",
                              snp_col = "SNP",
                              beta_col = "beta.outcome",
                              se_col = "se.outcome",
                              eaf_col = "eaf.outcome",
                              effect_allele_col = "effect_allele.outcome",
                              other_allele_col = "other_allele.outcome",
                              pval_col = "pval.outcome",
                              units_col = "units.outcome",
                              ncase_col = "ncase.outcome",
                              ncontrol_col = "ncontrol.outcome",
                              samplesize_col = "samplesize.outcome")
GWAS2_MR <- read_outcome_data("./MR/GWAS2_eQTL_MR_format.txt", sep='\t', phenotype_col = "outcome",
                              snp_col = "SNP",
                              beta_col = "beta.outcome",
                              se_col = "se.outcome",
                              eaf_col = "eaf.outcome",
                              effect_allele_col = "effect_allele.outcome",
                              other_allele_col = "other_allele.outcome",
                              pval_col = "pval.outcome",
                              units_col = "units.outcome",
                              ncase_col = "ncase.outcome",
                              ncontrol_col = "ncontrol.outcome",
                              samplesize_col = "samplesize.outcome")

MR1_aa <- harmonise_data(eQTL_MR_aa, GWAS1_MR, action = 1)
MR2_aa <- harmonise_data(eQTL_MR_aa, GWAS2_MR, action = 1)

MR1_ac <- harmonise_data(eQTL_MR_ac, GWAS1_MR, action = 1)
MR2_ac <- harmonise_data(eQTL_MR_ac, GWAS2_MR, action = 1)

MR1_at <- harmonise_data(eQTL_MR_at, GWAS1_MR, action = 1)
MR2_at <- harmonise_data(eQTL_MR_at, GWAS2_MR, action = 1)

MR1_ha <- harmonise_data(eQTL_MR_ha, GWAS1_MR, action = 1)
MR2_ha <- harmonise_data(eQTL_MR_ha, GWAS2_MR, action = 1)

MR1_hlv <- harmonise_data(eQTL_MR_hlv, GWAS1_MR, action = 1)
MR2_hlv <- harmonise_data(eQTL_MR_hlv, GWAS2_MR, action = 1)

mr_results1_aa <- mr(MR1_aa, method_list=c("mr_wald_ratio","mr_ivw","mr_egger_regression"))
mr_results2_aa <- mr(MR2_aa, method_list=c("mr_wald_ratio","mr_ivw","mr_egger_regression"))

mr_results1_ac <- mr(MR1_ac, method_list=c("mr_wald_ratio","mr_ivw","mr_egger_regression"))
mr_results2_ac <- mr(MR2_ac, method_list=c("mr_wald_ratio","mr_ivw","mr_egger_regression"))

mr_results1_at <- mr(MR1_at, method_list=c("mr_wald_ratio","mr_ivw","mr_egger_regression"))
mr_results2_at <- mr(MR2_at, method_list=c("mr_wald_ratio","mr_ivw","mr_egger_regression"))

mr_results1_ha <- mr(MR1_ha, method_list=c("mr_wald_ratio","mr_ivw","mr_egger_regression"))
mr_results2_ha <- mr(MR2_ha, method_list=c("mr_wald_ratio","mr_ivw","mr_egger_regression"))

mr_results1_hlv <- mr(MR1_hlv, method_list=c("mr_wald_ratio","mr_ivw","mr_egger_regression"))
mr_results2_hlv <- mr(MR2_hlv, method_list=c("mr_wald_ratio","mr_ivw","mr_egger_regression"))

results1_aa <- cbind.data.frame(mr_results1_aa$exposure,mr_results1_aa$outcome,mr_results1_aa$nsnp,mr_results1_aa$method,mr_results1_aa$b,mr_results1_aa$se,mr_results1_aa$pval)
results2_aa <- cbind.data.frame(mr_results2_aa$exposure,mr_results2_aa$outcome,mr_results2_aa$nsnp,mr_results2_aa$method,mr_results2_aa$b,mr_results2_aa$se,mr_results2_aa$pval)

results1_ac <- cbind.data.frame(mr_results1_ac$exposure,mr_results1_ac$outcome,mr_results1_ac$nsnp,mr_results1_ac$method,mr_results1_ac$b,mr_results1_ac$se,mr_results1_ac$pval)
results2_ac <- cbind.data.frame(mr_results2_ac$exposure,mr_results2_ac$outcome,mr_results2_ac$nsnp,mr_results2_ac$method,mr_results2_ac$b,mr_results2_ac$se,mr_results2_ac$pval)

results1_at <- cbind.data.frame(mr_results1_at$exposure,mr_results1_at$outcome,mr_results1_at$nsnp,mr_results1_at$method,mr_results1_at$b,mr_results1_at$se,mr_results1_at$pval)
results2_at <- cbind.data.frame(mr_results2_at$exposure,mr_results2_at$outcome,mr_results2_at$nsnp,mr_results2_at$method,mr_results2_at$b,mr_results2_at$se,mr_results2_at$pval)

results1_ha <- cbind.data.frame(mr_results1_ha$exposure,mr_results1_ha$outcome,mr_results1_ha$nsnp,mr_results1_ha$method,mr_results1_ha$b,mr_results1_ha$se,mr_results1_ha$pval)
results2_ha <- cbind.data.frame(mr_results2_ha$exposure,mr_results2_ha$outcome,mr_results2_ha$nsnp,mr_results2_ha$method,mr_results2_ha$b,mr_results2_ha$se,mr_results2_ha$pval)

results1_hlv <- cbind.data.frame(mr_results1_hlv$exposure,mr_results1_hlv$outcome,mr_results1_hlv$nsnp,mr_results1_hlv$method,mr_results1_hlv$b,mr_results1_hlv$se,mr_results1_hlv$pval)
results2_hlv <- cbind.data.frame(mr_results2_hlv$exposure,mr_results2_hlv$outcome,mr_results2_hlv$nsnp,mr_results2_hlv$method,mr_results2_hlv$b,mr_results2_hlv$se,mr_results2_hlv$pval)

write.csv(results1_aa,"./MR/MRresults/eQTL_aa_GWAS1_MRresults.csv")
write.csv(results2_aa,"./MR/MRresults/eQTL_aa_GWAS2_MRresults.csv")

write.csv(results1_ac,"./MR/MRresults/eQTL_ac_GWAS1_MRresults.csv")
write.csv(results2_ac,"./MR/MRresults/eQTL_ac_GWAS2_MRresults.csv")

write.csv(results1_at,"./MR/MRresults/eQTL_at_GWAS1_MRresults.csv")
write.csv(results2_at,"./MR/MRresults/eQTL_at_GWAS2_MRresults.csv")

write.csv(results1_ha,"./MR/MRresults/eQTL_ha_GWAS1_MRresults.csv")
write.csv(results2_ha,"./MR/MRresults/eQTL_ha_GWAS2_MRresults.csv")

write.csv(results1_hlv,"./MR/MRresults/eQTL_hlv_GWAS1_MRresults.csv")
write.csv(results2_hlv,"./MR/MRresults/eQTL_hlv_GWAS2_MRresults.csv")

# sensitivity analysis
# aa
mr_het1_aa <- mr_heterogeneity(MR1_aa)
nrow(mr_het1_aa[mr_het1_aa$Q_pval <= 0.01,])
# 0
mr_het2_aa <- mr_heterogeneity(MR2_aa)
nrow(mr_het2_aa[mr_het2_aa$Q_pval <= 0.01,])
# 0
# Not enough SNPs available for heterogeneity analysis
mr_pt1_aa <- mr_pleiotropy_test(MR1_aa)
nrow(mr_pt1_aa[mr_pt1_aa$pval <= 0.01,])
# 0
mr_pt2_aa <-mr_pleiotropy_test(MR2_aa)
nrow(mr_pt2_aa[mr_pt2_aa$pval <= 0.01,])
# 0
# Not enough SNPs available for pleiotropy analysis
res_single1_aa <- mr_singlesnp(MR1_aa)
res_single2_aa <- mr_singlesnp(MR2_aa)
write.csv(res_single1_aa,"./MR/MRresults/eQTL_aa_GWAS1_MRsinglesnp.csv")
write.csv(res_single2_aa,"./MR/MRresults/eQTL_aa_GWAS2_MRsinglesnp.csv")

# ac
mr_het1_ac <- mr_heterogeneity(MR1_ac)
nrow(mr_het1_ac[mr_het1_ac$Q_pval <= 0.01,])
# 0
mr_het2_ac <- mr_heterogeneity(MR2_ac)
nrow(mr_het2_ac[mr_het2_ac$Q_pval <= 0.01,])
# 0
mr_pt1_ac <- mr_pleiotropy_test(MR1_ac)
nrow(mr_pt1_ac[mr_pt1_ac$pval <= 0.01,])
# 0
mr_pt2_ac <-mr_pleiotropy_test(MR2_ac)
nrow(mr_pt2_ac[mr_pt2_ac$pval <= 0.01,])
# 0
res_single1_ac <- mr_singlesnp(MR1_ac)
res_single2_ac <- mr_singlesnp(MR2_ac)
write.csv(res_single1_ac,"./MR/MRresults/eQTL_ac_GWAS1_MRsinglesnp.csv")
write.csv(res_single2_ac,"./MR/MRresults/eQTL_ac_GWAS2_MRsinglesnp.csv")

# at
mr_het1_at <- mr_heterogeneity(MR1_at)
nrow(mr_het1_at[mr_het1_at$Q_pval <= 0.01,])
# 0
mr_het2_at <- mr_heterogeneity(MR2_at)
nrow(mr_het2_at[mr_het2_at$Q_pval <= 0.01,])
# 0
mr_pt1_at <- mr_pleiotropy_test(MR1_at)
nrow(mr_pt1_at[mr_pt1_at$pval <= 0.01,])
# 0
mr_pt2_at <-mr_pleiotropy_test(MR2_at)
nrow(mr_pt2_at[mr_pt2_at$pval <= 0.01,])
# 0
res_single1_at <- mr_singlesnp(MR1_at)
res_single2_at <- mr_singlesnp(MR2_at)
write.csv(res_single1_at,"./MR/MRresults/eQTL_at_GWAS1_MRsinglesnp.csv")
write.csv(res_single2_at,"./MR/MRresults/eQTL_at_GWAS2_MRsinglesnp.csv")

# ha
mr_het1_ha <- mr_heterogeneity(MR1_ha)
nrow(mr_het1_ha[mr_het1_ha$Q_pval <= 0.01,])
# 0
mr_het2_ha <- mr_heterogeneity(MR2_ha)
nrow(mr_het2_ha[mr_het2_ha$Q_pval <= 0.01,])
# 0
mr_pt1_ha <- mr_pleiotropy_test(MR1_ha)
nrow(mr_pt1_ha[mr_pt1_ha$pval <= 0.01,])
# 0
mr_pt2_ha <-mr_pleiotropy_test(MR2_ha)
nrow(mr_pt2_ha[mr_pt2_ha$pval <= 0.01,])
# 0
res_single1_ha <- mr_singlesnp(MR1_ha)
res_single2_ha <- mr_singlesnp(MR2_ha)
write.csv(res_single1_ha,"./MR/MRresults/eQTL_ha_GWAS1_MRsinglesnp.csv")
write.csv(res_single2_ha,"./MR/MRresults/eQTL_ha_GWAS2_MRsinglesnp.csv")

# hlv
mr_het1_hlv <- mr_heterogeneity(MR1_hlv)
nrow(mr_het1_hlv[mr_het1_hlv$Q_pval <= 0.01,])
# 0
mr_het2_hlv <- mr_heterogeneity(MR2_hlv)
nrow(mr_het2_hlv[mr_het2_hlv$Q_pval <= 0.01,])
# 0
mr_pt1_hlv <- mr_pleiotropy_test(MR1_hlv)
nrow(mr_pt1_hlv[mr_pt1_hlv$pval <= 0.01,])
# 0
mr_pt2_hlv <-mr_pleiotropy_test(MR2_hlv)
nrow(mr_pt2_hlv[mr_pt2_hlv$pval <= 0.01,])
# 0
res_single1_hlv <- mr_singlesnp(MR1_hlv)
res_single2_hlv <- mr_singlesnp(MR2_hlv)
write.csv(res_single1_hlv,"./MR/MRresults/eQTL_hlv_GWAS1_MRsinglesnp.csv")
write.csv(res_single2_hlv,"./MR/MRresults/eQTL_hlv_GWAS2_MRsinglesnp.csv")

# We can also create a volcano plot for multiple MR results 
# aa, GWAS1
fres1_aa <- mr_results1_aa
fres1_aa$effect <- exp(fres1_aa$b)
fres1_aa$up_ci <- exp(fres1_aa$b+(1.96*fres1_aa$se))
fres1_aa$lo_ci <- exp(fres1_aa$b-(1.96*fres1_aa$se))
fres1_aa<-fres1_aa[!is.na(fres1_aa$effect),]
fres1_aa$bon <- p.adjust(fres1_aa$pval, method="bonferroni")
fres1_aa$sig <- fres1_aa$bon < 0.05
fres1_aa$category <-"<0.05"
fres1_aa$sample <- "Artery Aorta"
fres1_aa <- merge(fres1_aa, MR1_aa[,c(1,23)], by="exposure", all.x = T)
write.csv(fres1_aa,"./MR/MRresults/eQTL_aa_GWAS1_MRresults_adjusted.csv")

ggplot(fres1_aa, aes(x=effect, y=-log10(pval))) +
  geom_vline(xintercept=1, linetype="dotted") +
  geom_point(data=subset(fres1_aa, !sig)) +
  geom_point(data=subset(fres1_aa, sig), aes(colour=category, size=bon < 0.05)) +
  facet_grid(. ~ outcome, scale="free") +
  geom_label_repel(data=subset(fres1_aa, bon<0.05), aes(label=exposure, fill=category), colour="white", segment.colour="black", point.padding = unit(0.7, "lines"), box.padding = unit(0.7, "lines"), segment.size=0.5, force=2, max.iter=3e3) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.x=element_text(size=20)
  ) +
  scale_colour_brewer(type="qual", palette="Dark2") +
  scale_fill_brewer(type="qual", palette="Dark2") +
  labs(x="OR for CHD per SD increase in Gene Expression", size="Bonferroni",y="P value (-log10)") +
  xlim(c(0.8, 1.1))+
  theme(axis.title.y=element_text(size=18),axis.title.x=element_text(size=18),axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),legend.position="none")
ggsave("./MR/MRresults/volcanoplot_eQTL_aa_CHD1.png", width=8, height=8)

# aa, GWAS2
fres2_aa <- mr_results2_aa
fres2_aa$effect <- exp(fres2_aa$b)
fres2_aa$up_ci <- exp(fres2_aa$b+(1.96*fres2_aa$se))
fres2_aa$lo_ci <- exp(fres2_aa$b-(1.96*fres2_aa$se))
fres2_aa<-fres2_aa[!is.na(fres2_aa$effect),]
fres2_aa$bon <- p.adjust(fres2_aa$pval, method="bonferroni")
fres2_aa$sig <- fres2_aa$bon < 0.05
fres2_aa$category <-"<0.05"
fres2_aa$sample <- "Artery Aorta"
fres2_aa <- merge(fres2_aa, MR2_aa[,c(1,23)], by="exposure", all.x = T)
write.csv(fres2_aa,"./MR/MRresults/eQTL_aa_GWAS2_MRresults_adjusted.csv")

ggplot(fres2_aa, aes(x=effect, y=-log10(pval))) +
  geom_vline(xintercept=1, linetype="dotted") +
  geom_point(data=subset(fres2_aa, !sig)) +
  geom_point(data=subset(fres2_aa, sig), aes(colour=category, size=bon < 0.05)) +
  facet_grid(. ~ outcome, scale="free") +
  geom_label_repel(data=subset(fres2_aa, bon<0.05), aes(label=exposure, fill=category), colour="white", segment.colour="black", point.padding = unit(0.7, "lines"), box.padding = unit(0.7, "lines"), segment.size=0.5, force=2, max.iter=3e3) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.x=element_text(size=20)
  ) +
  scale_colour_brewer(type="qual", palette="Dark2") +
  scale_fill_brewer(type="qual", palette="Dark2") +
  labs(x="OR for CHD per SD increase in Gene Expression", size="Bonferroni",y="P value (-log10)") +
  xlim(c(0.8, 1.1))+
  theme(axis.title.y=element_text(size=18),axis.title.x=element_text(size=18),axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),legend.position="none")
ggsave("./MR/MRresults/volcanoplot_eQTL_aa_CHD2.png", width=8, height=8)

# ac, GWAS1
fres1_ac <- mr_results1_ac
fres1_ac$effect <- exp(fres1_ac$b)
fres1_ac$up_ci <- exp(fres1_ac$b+(1.96*fres1_ac$se))
fres1_ac$lo_ci <- exp(fres1_ac$b-(1.96*fres1_ac$se))
fres1_ac<-fres1_ac[!is.na(fres1_ac$effect),]
fres1_ac$bon <- p.adjust(fres1_ac$pval, method="bonferroni")
fres1_ac$sig <- fres1_ac$bon < 0.05
fres1_ac$category <-"<0.05"
fres1_ac$sample <- "Artery Coronary"
fres1_ac <- merge(fres1_ac, MR1_ac[,c(1,23)], by="exposure", all.x = T)
write.csv(fres1_ac,"./MR/MRresults/eQTL_ac_GWAS1_MRresults_adjusted.csv")

ggplot(fres1_ac, aes(x=effect, y=-log10(pval))) +
  geom_vline(xintercept=1, linetype="dotted") +
  geom_point(data=subset(fres1_ac, !sig)) +
  geom_point(data=subset(fres1_ac, sig), aes(colour=category, size=bon < 0.05)) +
  facet_grid(. ~ outcome, scale="free") +
  geom_label_repel(data=subset(fres1_ac, bon<0.05), aes(label=exposure, fill=category), colour="white", segment.colour="black", point.padding = unit(0.7, "lines"), box.padding = unit(0.7, "lines"), segment.size=0.5, force=2, max.iter=3e3) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.x=element_text(size=20)
  ) +
  scale_colour_brewer(type="qual", palette="Dark2") +
  scale_fill_brewer(type="qual", palette="Dark2") +
  labs(x="OR for CHD per SD increase in Gene Expression", size="Bonferroni",y="P value (-log10)") +
  xlim(c(0.8, 1.1))+
  theme(axis.title.y=element_text(size=18),axis.title.x=element_text(size=18),axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),legend.position="none")
ggsave("./MR/MRresults/volcanoplot_eQTL_ac_CHD1.png", width=8, height=8)

# ac, GWAS2
fres2_ac <- mr_results2_ac
fres2_ac$effect <- exp(fres2_ac$b)
fres2_ac$up_ci <- exp(fres2_ac$b+(1.96*fres2_ac$se))
fres2_ac$lo_ci <- exp(fres2_ac$b-(1.96*fres2_ac$se))
fres2_ac<-fres2_ac[!is.na(fres2_ac$effect),]
fres2_ac$bon <- p.adjust(fres2_ac$pval, method="bonferroni")
fres2_ac$sig <- fres2_ac$bon < 0.05
fres2_ac$category <-"<0.05"
fres2_ac$sample <- "Artery Coronary"
fres2_ac <- merge(fres2_ac, MR2_ac[,c(1,23)], by="exposure", all.x = T)
write.csv(fres2_ac,"./MR/MRresults/eQTL_ac_GWAS2_MRresults_adjusted.csv")

ggplot(fres2_ac, aes(x=effect, y=-log10(pval))) +
  geom_vline(xintercept=1, linetype="dotted") +
  geom_point(data=subset(fres2_ac, !sig)) +
  geom_point(data=subset(fres2_ac, sig), aes(colour=category, size=bon < 0.05)) +
  facet_grid(. ~ outcome, scale="free") +
  geom_label_repel(data=subset(fres2_ac, bon<0.05), aes(label=exposure, fill=category), colour="white", segment.colour="black", point.padding = unit(0.7, "lines"), box.padding = unit(0.7, "lines"), segment.size=0.5, force=2, max.iter=3e3) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.x=element_text(size=20)
  ) +
  scale_colour_brewer(type="qual", palette="Dark2") +
  scale_fill_brewer(type="qual", palette="Dark2") +
  labs(x="OR for CHD per SD increase in Gene Expression", size="Bonferroni",y="P value (-log10)") +
  xlim(c(0.8, 1.1))+
  theme(axis.title.y=element_text(size=18),axis.title.x=element_text(size=18),axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),legend.position="none")
ggsave("./MR/MRresults/volcanoplot_eQTL_ac_CHD2.png", width=8, height=8)

# at, GWAS1
fres1_at <- mr_results1_at
fres1_at$effect <- exp(fres1_at$b)
fres1_at$up_ci <- exp(fres1_at$b+(1.96*fres1_at$se))
fres1_at$lo_ci <- exp(fres1_at$b-(1.96*fres1_at$se))
fres1_at<-fres1_at[!is.na(fres1_at$effect),]
fres1_at$bon <- p.adjust(fres1_at$pval, method="bonferroni")
fres1_at$sig <- fres1_at$bon < 0.05
fres1_at$category <-"<0.05"
fres1_at$sample <- "Artery Tibial"
fres1_at <- merge(fres1_at, MR1_at[,c(1,23)], by="exposure", all.x = T)
write.csv(fres1_at,"./MR/MRresults/eQTL_at_GWAS1_MRresults_adjusted.csv")

ggplot(fres1_at, aes(x=effect, y=-log10(pval))) +
  geom_vline(xintercept=1, linetype="dotted") +
  geom_point(data=subset(fres1_at, !sig)) +
  geom_point(data=subset(fres1_at, sig), aes(colour=category, size=bon < 0.05)) +
  facet_grid(. ~ outcome, scale="free") +
  geom_label_repel(data=subset(fres1_at, bon<0.05), aes(label=exposure, fill=category), colour="white", segment.colour="black", point.padding = unit(0.7, "lines"), box.padding = unit(0.7, "lines"), segment.size=0.5, force=2, max.iter=3e3) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.x=element_text(size=20)
  ) +
  scale_colour_brewer(type="qual", palette="Dark2") +
  scale_fill_brewer(type="qual", palette="Dark2") +
  labs(x="OR for CHD per SD increase in Gene Expression", size="Bonferroni",y="P value (-log10)") +
  xlim(c(0.8, 1.1))+
  theme(axis.title.y=element_text(size=18),axis.title.x=element_text(size=18),axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),legend.position="none")
ggsave("./MR/MRresults/volcanoplot_eQTL_at_CHD1.png", width=8, height=8)

# at, GWAS2
fres2_at <- mr_results2_at
fres2_at$effect <- exp(fres2_at$b)
fres2_at$up_ci <- exp(fres2_at$b+(1.96*fres2_at$se))
fres2_at$lo_ci <- exp(fres2_at$b-(1.96*fres2_at$se))
fres2_at<-fres2_at[!is.na(fres2_at$effect),]
fres2_at$bon <- p.adjust(fres2_at$pval, method="bonferroni")
fres2_at$sig <- fres2_at$bon < 0.05
fres2_at$category <-"<0.05"
fres2_at$sample <- "Artery Tibial"
fres2_at <- merge(fres2_at, MR2_at[,c(1,23)], by="exposure", all.x = T)
write.csv(fres2_at,"./MR/MRresults/eQTL_at_GWAS2_MRresults_adjusted.csv")

ggplot(fres2_at, aes(x=effect, y=-log10(pval))) +
  geom_vline(xintercept=1, linetype="dotted") +
  geom_point(data=subset(fres2_at, !sig)) +
  geom_point(data=subset(fres2_at, sig), aes(colour=category, size=bon < 0.05)) +
  facet_grid(. ~ outcome, scale="free") +
  geom_label_repel(data=subset(fres2_at, bon<0.05), aes(label=exposure, fill=category), colour="white", segment.colour="black", point.padding = unit(0.7, "lines"), box.padding = unit(0.7, "lines"), segment.size=0.5, force=2, max.iter=3e3) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.x=element_text(size=20)
  ) +
  scale_colour_brewer(type="qual", palette="Dark2") +
  scale_fill_brewer(type="qual", palette="Dark2") +
  labs(x="OR for CHD per SD increase in Gene Expression", size="Bonferroni",y="P value (-log10)") +
  xlim(c(0.8, 1.1))+
  theme(axis.title.y=element_text(size=18),axis.title.x=element_text(size=18),axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),legend.position="none")
ggsave("./MR/MRresults/volcanoplot_eQTL_at_CHD2.png", width=8, height=8)

# ha, GWAS1
fres1_ha <- mr_results1_ha
fres1_ha$effect <- exp(fres1_ha$b)
fres1_ha$up_ci <- exp(fres1_ha$b+(1.96*fres1_ha$se))
fres1_ha$lo_ci <- exp(fres1_ha$b-(1.96*fres1_ha$se))
fres1_ha<-fres1_ha[!is.na(fres1_ha$effect),]
fres1_ha$bon <- p.adjust(fres1_ha$pval, method="bonferroni")
fres1_ha$sig <- fres1_ha$bon < 0.05
fres1_ha$category <-"<0.05"
fres1_ha$sample <- "Heart Atrial Appendage"
fres1_ha <- merge(fres1_ha, MR1_ha[,c(1,23)], by="exposure", all.x = T)
write.csv(fres1_ha,"./MR/MRresults/eQTL_ha_GWAS1_MRresults_adjusted.csv")

ggplot(fres1_ha, aes(x=effect, y=-log10(pval))) +
  geom_vline(xintercept=1, linetype="dotted") +
  geom_point(data=subset(fres1_ha, !sig)) +
  geom_point(data=subset(fres1_ha, sig), aes(colour=category, size=bon < 0.05)) +
  facet_grid(. ~ outcome, scale="free") +
  geom_label_repel(data=subset(fres1_ha, bon<0.05), aes(label=exposure, fill=category), colour="white", segment.colour="black", point.padding = unit(0.7, "lines"), box.padding = unit(0.7, "lines"), segment.size=0.5, force=2, max.iter=3e3) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.x=element_text(size=20)
  ) +
  scale_colour_brewer(type="qual", palette="Dark2") +
  scale_fill_brewer(type="qual", palette="Dark2") +
  labs(x="OR for CHD per SD increase in Gene Expression", size="Bonferroni",y="P value (-log10)") +
  xlim(c(0.8, 1.1))+
  theme(axis.title.y=element_text(size=18),axis.title.x=element_text(size=18),axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),legend.position="none")
ggsave("./MR/MRresults/volcanoplot_eQTL_ha_CHD1.png", width=8, height=8)

# ha, GWAS2
fres2_ha <- mr_results2_ha
fres2_ha$effect <- exp(fres2_ha$b)
fres2_ha$up_ci <- exp(fres2_ha$b+(1.96*fres2_ha$se))
fres2_ha$lo_ci <- exp(fres2_ha$b-(1.96*fres2_ha$se))
fres2_ha<-fres2_ha[!is.na(fres2_ha$effect),]
fres2_ha$bon <- p.adjust(fres2_ha$pval, method="bonferroni")
fres2_ha$sig <- fres2_ha$bon < 0.05
fres2_ha$category <-"<0.05"
fres2_ha$sample <- "Heart Atrial Appendage"
fres2_ha <- merge(fres2_ha, MR2_ha[,c(1,23)], by="exposure", all.x = T)
write.csv(fres2_ha,"./MR/MRresults/eQTL_ha_GWAS2_MRresults_adjusted.csv")

ggplot(fres2_ha, aes(x=effect, y=-log10(pval))) +
  geom_vline(xintercept=1, linetype="dotted") +
  geom_point(data=subset(fres2_ha, !sig)) +
  geom_point(data=subset(fres2_ha, sig), aes(colour=category, size=bon < 0.05)) +
  facet_grid(. ~ outcome, scale="free") +
  geom_label_repel(data=subset(fres2_ha, bon<0.05), aes(label=exposure, fill=category), colour="white", segment.colour="black", point.padding = unit(0.7, "lines"), box.padding = unit(0.7, "lines"), segment.size=0.5, force=2, max.iter=3e3) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.x=element_text(size=20)
  ) +
  scale_colour_brewer(type="qual", palette="Dark2") +
  scale_fill_brewer(type="qual", palette="Dark2") +
  labs(x="OR for CHD per SD increase in Gene Expression", size="Bonferroni",y="P value (-log10)") +
  xlim(c(0.8, 1.1))+
  theme(axis.title.y=element_text(size=18),axis.title.x=element_text(size=18),axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),legend.position="none")
ggsave("./MR/MRresults/volcanoplot_eQTL_ha_CHD2.png", width=8, height=8)

# hlv, GWAS1
fres1_hlv <- mr_results1_hlv
fres1_hlv$effect <- exp(fres1_hlv$b)
fres1_hlv$up_ci <- exp(fres1_hlv$b+(1.96*fres1_hlv$se))
fres1_hlv$lo_ci <- exp(fres1_hlv$b-(1.96*fres1_hlv$se))
fres1_hlv<-fres1_hlv[!is.na(fres1_hlv$effect),]
fres1_hlv$bon <- p.adjust(fres1_hlv$pval, method="bonferroni")
fres1_hlv$sig <- fres1_hlv$bon < 0.05
fres1_hlv$category <-"<0.05"
fres1_hlv$sample <- "Heart Left Ventricle"
fres1_hlv <- merge(fres1_hlv, MR1_hlv[,c(1,23)], by="exposure", all.x = T)
write.csv(fres1_hlv,"./MR/MRresults/eQTL_hlv_GWAS1_MRresults_adjusted.csv")

ggplot(fres1_hlv, aes(x=effect, y=-log10(pval))) +
  geom_vline(xintercept=1, linetype="dotted") +
  geom_point(data=subset(fres1_hlv, !sig)) +
  geom_point(data=subset(fres1_hlv, sig), aes(colour=category, size=bon < 0.05)) +
  facet_grid(. ~ outcome, scale="free") +
  geom_label_repel(data=subset(fres1_hlv, bon<0.05), aes(label=exposure, fill=category), colour="white", segment.colour="black", point.padding = unit(0.7, "lines"), box.padding = unit(0.7, "lines"), segment.size=0.5, force=2, max.iter=3e3) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.x=element_text(size=20)
  ) +
  scale_colour_brewer(type="qual", palette="Dark2") +
  scale_fill_brewer(type="qual", palette="Dark2") +
  labs(x="OR for CHD per SD increase in Gene Expression", size="Bonferroni",y="P value (-log10)") +
  xlim(c(0.8, 1.1))+
  theme(axis.title.y=element_text(size=18),axis.title.x=element_text(size=18),axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),legend.position="none")
ggsave("./MR/MRresults/volcanoplot_eQTL_hlv_CHD1.png", width=8, height=8)

# hlv, GWAS2
fres2_hlv <- mr_results2_hlv
fres2_hlv$effect <- exp(fres2_hlv$b)
fres2_hlv$up_ci <- exp(fres2_hlv$b+(1.96*fres2_hlv$se))
fres2_hlv$lo_ci <- exp(fres2_hlv$b-(1.96*fres2_hlv$se))
fres2_hlv<-fres2_hlv[!is.na(fres2_hlv$effect),]
fres2_hlv$bon <- p.adjust(fres2_hlv$pval, method="bonferroni")
fres2_hlv$sig <- fres2_hlv$bon < 0.05
fres2_hlv$category <-"<0.05"
fres2_hlv$sample <- "Heart Left Ventricle"
fres2_hlv <- merge(fres2_hlv, MR2_hlv[,c(1,23)], by="exposure", all.x = T)
write.csv(fres2_hlv,"./MR/MRresults/eQTL_hlv_GWAS2_MRresults_adjusted.csv")

ggplot(fres2_hlv, aes(x=effect, y=-log10(pval))) +
  geom_vline(xintercept=1, linetype="dotted") +
  geom_point(data=subset(fres2_hlv, !sig)) +
  geom_point(data=subset(fres2_hlv, sig), aes(colour=category, size=bon < 0.05)) +
  facet_grid(. ~ outcome, scale="free") +
  geom_label_repel(data=subset(fres2_hlv, bon<0.05), aes(label=exposure, fill=category), colour="white", segment.colour="black", point.padding = unit(0.7, "lines"), box.padding = unit(0.7, "lines"), segment.size=0.5, force=2, max.iter=3e3) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.x=element_text(size=20)
  ) +
  scale_colour_brewer(type="qual", palette="Dark2") +
  scale_fill_brewer(type="qual", palette="Dark2") +
  labs(x="OR for CHD per SD increase in Gene Expression", size="Bonferroni",y="P value (-log10)") +
  xlim(c(0.8, 1.1))+
  theme(axis.title.y=element_text(size=18),axis.title.x=element_text(size=18),axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),legend.position="none")
ggsave("./MR/MRresults/volcanoplot_eQTL_hlv_CHD2.png", width=8, height=8)

# table for manuscript
# there is no bonferroni significant MR-eQTL, so we extract the nominal significant findings
sig1 <- rbind(fres1_aa[fres1_aa$pval <= 0.05,],fres1_ac[fres1_ac$pval <= 0.05,],fres1_at[fres1_at$pval <= 0.05,],fres1_ha[fres1_ha$pval <= 0.05,],fres1_hlv[fres1_hlv$pval <= 0.05,])
sig2 <- rbind(fres2_aa[fres2_aa$pval <= 0.05,],fres2_ac[fres2_ac$pval <= 0.05,],fres2_at[fres2_at$pval <= 0.05,],fres2_ha[fres2_ha$pval <= 0.05,],fres2_hlv[fres2_hlv$pval <= 0.05,])

GWAS1 <- read.table("./MR/GWAS1_eQTL_MR_format.txt", header = T)
GWAS2 <- read.table("./MR/GWAS2_eQTL_MR_format.txt", header = T)

results1 <- rbind(fres1_aa, fres1_ac, fres1_at, fres1_ha, fres1_hlv)
results2 <- rbind(fres2_aa, fres2_ac, fres2_at, fres2_ha, fres2_hlv)

eQTL_MR <- read.table("./MR/eQTLs_MR.txt", header = T)
eQTL_MR <- eQTL_MR[,c(2:7,9,8,10,11,16)]
colnames(eQTL_MR)[c(1,6:10)] <- c("exposure","SNP_pos","effect_allele","other_allele","MAF_eQTL","p.eQTL")

sig1 <- sig1[,c(1,4:13,16:17)]
colnames(sig1)[5:11] <- paste(colnames(sig1)[5:11], ".GWAS1", sep = '')
sig1 <- merge(sig1, results2[,c(1,7:13,16)], by=c("exposure","sample"), all.x = T)
colnames(sig1)[14:20] <- paste(colnames(sig1)[14:20], ".GWAS2", sep = '')
colnames(sig1)[c(8,16)] <- c("p.MR.GWAS1","p.MR.GWAS2")
sig1 <- merge(sig1, GWAS1[,c(1,4)], by = "SNP")
colnames(sig1)[21] <- "p.GWAS1"
sig1 <- merge(sig1, GWAS2[,c(1,4)], by = "SNP")
colnames(sig1)[22] <- "p.GWAS2"
sig1 <- merge(sig1[,2:22], eQTL_MR, by=c("exposure","sample"), all.x = T)
sig1 <- sig1[order(sig1$sample),c(1,22:29,2:21,30)]
write.csv(sig1,"./MR/MRresults/eQTL_GWAS1_MRresults_sig.csv")

sig2 <- sig2[,c(1,4:13,16:17)]
colnames(sig2)[5:11] <- paste(colnames(sig2)[5:11], ".GWAS2", sep = '')
sig2 <- merge(sig2, results1[,c(1,7:13,16)], by=c("exposure","sample"), all.x = T)
colnames(sig2)[14:20] <- paste(colnames(sig2)[14:20], ".GWAS1", sep = '')
colnames(sig2)[c(8,16)] <- c("p.MR.GWAS2","p.MR.GWAS1")
sig2 <- merge(sig2, GWAS1[,c(1,4)], by = "SNP")
colnames(sig2)[21] <- "p.GWAS1"
sig2 <- merge(sig2, GWAS2[,c(1,4)], by = "SNP")
colnames(sig2)[22] <- "p.GWAS2"
sig2 <- merge(sig2[,2:22], eQTL_MR, by=c("exposure","sample"), all.x = T)
sig2 <- sig2[order(sig2$sample),c(1,22:29,2:5,13:19,6:12,20:21,30)]
write.csv(sig2,"./MR/MRresults/eQTL_GWAS2_MRresults_sig.csv")

# cross check with MR-mQTL
mQTL_sig1 <- read.csv("./MR/MRresults/mQTL/mQTL_GWAS1_MRresults_sig2.csv")
colnames(mQTL_sig1)[2:31] <- paste(colnames(mQTL_sig1)[2:31], ".mQTL", sep = '')
sig1_overlap <- sig1[sig1$SNP %in% mQTL_sig1$SNP,]
#1
colnames(sig1_overlap)[14:27] <- paste(colnames(sig1_overlap)[14:27], ".eQTL", sep = '')
sig1_overlap <- merge(sig1_overlap, mQTL_sig1, by="SNP", all.x = T)
write.csv(sig1_overlap,"./MR/MRresults/eQTL_mQTL_MRresults_overlap.csv")

mQTL_sig2 <- read.csv("./MR/MRresults/mQTL/mQTL_GWAS2_MRresults_sig2.csv")
sig2_overlap <- sig2[sig2$SNP %in% mQTL_sig2$SNP,]
#0
