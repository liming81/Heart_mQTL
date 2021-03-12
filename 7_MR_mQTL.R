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

# Table S1 for SNP and CpG
tbs1 <- read.csv("./MR/Supp_Table1.csv")

# get beta for MR
CpG <- as.character(unique(tbs1$CpG_site))
bval <- NULL
for (i in 1:22){
  chr <- paste("chr",i, sep = '')
  cpg <- as.character(unique(tbs1$CpG_site[tbs1$chr_CpG == chr]))
  beta <- read.table(paste("./Epigenetic/Beta_CHR",i,".txt", sep = ''), header = T)
  beta_cpg <- beta[beta$site %in% cpg,]
  bval <- rbind(bval, beta_cpg)
}
colnames(bval) <- str_lowercase(colnames(bval))
colnames(bval) <- gsub("x","",colnames(bval), ignore.case =T)
bval <-setNames(data.frame(t(bval[,-1])), bval[,1])
bval$obs <- rownames(bval)
length(tbs1$CpG_site[tbs1$CpG_site %in% colnames(bval)])
# 24188
write.table(bval,file= "./MR/Beta_MR.txt",row.names=F,col.names=T,quote=F,sep='\t')

# get genotypes for MR
SNP <- as.character(unique(tbs1$SNP))
bed <- "./GWAS_new/GWAS_Heart_miss_sub_0.05_miss_snp_0.05_sub83_20181001.bed"
fam <- "./GWAS_new/GWAS_Heart_miss_sub_0.05_miss_snp_0.05_sub83_20181001.fam"
bim <- "./GWAS_new/GWAS_Heart_miss_sub_0.05_miss_snp_0.05_sub83_20181001.bim"
chr_snp <- read.bed(bed, bim, fam)
colnames(chr_snp$snp) <- chr_snp$snp.info$ID
rownames(chr_snp$snp) <- str_lowercase(chr_snp$ind.info$FamID)
var <- intersect(chr_snp$snp.info$ID, SNP)
geno <- as.data.frame(subset(chr_snp$snp,select=var))
geno$obs <- rownames(geno)
length(tbs1$SNP[tbs1$SNP %in% colnames(geno)])
# 24188
write.table(geno,file= "./MR/geno_MR.txt",row.names=F,col.names=T,quote=T,sep='\t')

################################################################
# Start here
################################################################
# Table S1 for SNP and CpG
tbs1 <- read.csv("./MR/Supp_Table1.csv")

# genotypes for MR
geno <- read_excel("./MR/geno_MR.xlsx",sheet = 1)
ID <- geno$obs
rownames(geno) <- geno$obs
geno <- geno[,1:13901]

# beta for MR
bval <- read.table("./MR/Beta_MR.txt", header = T)
bval <- bval[match(ID, bval$obs),]
rownames(bval) <- bval$obs
bval <- bval[,1:1676]

# covariates
covariates <- read.table("./covariates beta value.txt", header=T)
covariates <- covariates[match(ID, covariates$ID),]
rownames(covariates) <- covariates$ID
covariates <- covariates[,c(3,5:14)]

# summary statistics for linear regression
tbs1 <- tbs1[,1:6]
tbs1$beta <- 0
tbs1$se <- 0
tbs1$p_reg <- 0
tbs1$r2 <- 0
tbs1$r2.adj <- 0
for (i in 1:nrow(tbs1)){
  cpg <- as.character(tbs1$CpG_site[i])
  snp <- as.character(tbs1$SNP[i])
  temp <- cbind(subset(bval,select=cpg), subset(geno,select=snp), covariates)
  formula <- as.formula(paste(cpg, "~ .", sep = ''))
  mod <- lm(formula, data=temp)
  tbs1$beta[i] <- summary(mod)$coefficients[2,1]
  tbs1$se[i] <- summary(mod)$coefficients[2,2]
  tbs1$p_reg[i] <- summary(mod)$coefficients[2,4]
  tbs1$r2[i] <- summary(mod)$r.squared
  tbs1$r2.adj[i] <- summary(mod)$adj.r.squared
}
tbs1$p_BH <- p.adjust(tbs1$p_reg, method = "BH")
mQTL_summary <- tbs1[which(abs(tbs1$beta) > 0.1 & (tbs1$r2 > 0.5) & (tbs1$p_BH < 0.05)),]
# 22301
write.table(mQTL_summary, file= "./MR/mQTLs_summary_statistics.txt",row.names=F,col.names=T,quote=T,sep='\t')

# summary statistics for GWAS1
gwas <- read.table("./colocalization/GWAS_summary_statistics.txt", header = T)
gwas$CHR <- paste("chr",gwas$CHR, sep = '')
GWAS1 <- gwas[gwas$study == "GWAS1",]

# summary statistics for GWAS2
GWAS2 <- gwas[gwas$study == "GWAS2",]

# overlapped SNPs
mQTL_summary <- read.table("./MR/mQTLs_summary_statistics.txt", header = T)
snp_mQTL <- as.character(unique(mQTL_summary$SNP))
# 12958
snp_overlap <- snp_mQTL[snp_mQTL %in% GWAS1$SNP & (snp_mQTL %in% GWAS2$SNP)]
# 8677
write.table(as.matrix(snp_overlap,ncol=1), file= "./MR/plink/snp_overlap.txt",row.names=F,col.names=F,quote=T,sep='\t')
mQTL <- mQTL_summary[mQTL_summary$SNP %in% snp_overlap,]
# 13678
write.table(mQTL, file= "./MR/mQTLs_overlap.txt",row.names=F,col.names=T,quote=T,sep='\t')
GWAS1 <- GWAS1[GWAS1$SNP %in% snp_overlap,]
write.table(GWAS1, file= "./MR/GWAS1_overlap.txt",row.names=F,col.names=T,quote=T,sep='\t')
GWAS2 <- GWAS2[GWAS2$SNP %in% snp_overlap,]
write.table(GWAS2, file= "./MR/GWAS2_overlap.txt",row.names=F,col.names=T,quote=T,sep='\t')

# get GWAS bfile with overlapped snps
# plink
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile GWAS_Heart_miss_sub_0.05_miss_snp_0.05_sub83_20181001')
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile GWAS_Heart_miss_sub_0.05_miss_snp_0.05_sub83_20181001 --make-bed --extract snp_overlap.txt --allow-no-sex --out GWAS_overlap')

# double check
bed <- "./MR/plink/GWAS_overlap.bed"
fam <- "./MR/plink/GWAS_overlap.fam"
bim <- "./MR/plink/GWAS_overlap.bim"
chr_snp <- read.bed(bed, bim, fam)
# 83 obs * 8677 snps

# make haploview format
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile GWAS_overlap --recodeHV --chr 1 --allow-no-sex --out mQTL_haploview_1')
# 834
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile GWAS_overlap --recodeHV --chr 2 --allow-no-sex --out mQTL_haploview_2')
# 607
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile GWAS_overlap --recodeHV --chr 3 --allow-no-sex --out mQTL_haploview_3')
# 514
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile GWAS_overlap --recodeHV --chr 4 --allow-no-sex --out mQTL_haploview_4')
# 333
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile GWAS_overlap --recodeHV --chr 5 --allow-no-sex --out mQTL_haploview_5')
# 441
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile GWAS_overlap --recodeHV --chr 6 --allow-no-sex --out mQTL_haploview_6')
# 2021
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile GWAS_overlap --recodeHV --chr 7 --allow-no-sex --out mQTL_haploview_7')
# 441
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile GWAS_overlap --recodeHV --chr 8 --allow-no-sex --out mQTL_haploview_8')
# 446
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile GWAS_overlap --recodeHV --chr 9 --allow-no-sex --out mQTL_haploview_9')
# 109
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile GWAS_overlap --recodeHV --chr 10 --allow-no-sex --out mQTL_haploview_10')
# 420
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile GWAS_overlap --recodeHV --chr 11 --allow-no-sex --out mQTL_haploview_11')
# 319
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile GWAS_overlap --recodeHV --chr 12 --allow-no-sex --out mQTL_haploview_12')
# 361
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile GWAS_overlap --recodeHV --chr 13 --allow-no-sex --out mQTL_haploview_13')
# 259
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile GWAS_overlap --recodeHV --chr 14 --allow-no-sex --out mQTL_haploview_14')
# 164
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile GWAS_overlap --recodeHV --chr 15 --allow-no-sex --out mQTL_haploview_15')
# 181
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile GWAS_overlap --recodeHV --chr 16 --allow-no-sex --out mQTL_haploview_16')
# 274
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile GWAS_overlap --recodeHV --chr 17 --allow-no-sex --out mQTL_haploview_17')
# 385
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile GWAS_overlap --recodeHV --chr 18 --allow-no-sex --out mQTL_haploview_18')
# 82
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile GWAS_overlap --recodeHV --chr 19 --allow-no-sex --out mQTL_haploview_19')
# 259
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile GWAS_overlap --recodeHV --chr 20 --allow-no-sex --out mQTL_haploview_20')
# 128
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile GWAS_overlap --recodeHV --chr 21 --allow-no-sex --out mQTL_haploview_21')
# 48
system('/Users/chenlyu/Downloads/plink_mac/plink --bfile GWAS_overlap --recodeHV --chr 22 --allow-no-sex --out mQTL_haploview_22')
# 51

# check LD r^2 to make sure independent SNPs for each CpG
# get tag SNP for each CpG

# run on terminal
# cd "Documents/LC/Horton/HaploView"
# java -jar Haploview.jar
# tagger; pairwise tagging; r^2 threshold = 0.8; LOD threshold = 3

# get tag SNP list
# 3125 for pairwise tagging, 2821 for aggressive tagging
tagSNP <- NULL
for (i in 1:22){
  list <- read.table(paste("./MR/tagSNP/tagSNP_chr", i, sep = ''))
  tagSNP <- rbind(tagSNP, list)
}
write.table(tagSNP, file= "./MR/tagSNP/tagSNPlist.txt",row.names=F,col.names=F,quote=T,sep='\t')

# get MR datasets
tagSNP<- read.table("./MR/tagSNP/tagSNPlist.txt",header = F)
tagSNP <- as.character(tagSNP[,1])

GWAS1 <- read.table("./MR/GWAS1_overlap.txt", header = T)
GWAS1 <- GWAS1[GWAS1$SNP %in% tagSNP,]
write.table(GWAS1, file= "./MR/GWAS1_MR.txt",row.names=F,col.names=T,quote=T,sep='\t')

GWAS2 <- read.table("./MR/GWAS2_overlap.txt", header = T)
GWAS2 <- GWAS2[GWAS2$SNP %in% tagSNP,]
write.table(GWAS2, file= "./MR/GWAS2_MR.txt",row.names=F,col.names=T,quote=T,sep='\t')

mQTL <- read.table("./MR/mQTLs_overlap.txt", header = T)
mQTL_MR <- mQTL[mQTL$SNP %in% tagSNP,]
# 2821 SNPs, 1475 CpG sites
write.table(mQTL_MR, file= "./MR/mQTLs_MR.txt",row.names=F,col.names=T,quote=T,sep='\t')

# add maf
mQTL_MR <- read.table("./MR/mQTLs_MR.txt", header = T)

bed <- "./MR/plink/GWAS_overlap.bed"
fam <- "./MR/plink/GWAS_overlap.fam"
bim <- "./MR/plink/GWAS_overlap.bim"
chr_snp <- read.bed(bed, bim, fam)
colnames(chr_snp$snp) <- chr_snp$snp.info$ID
rownames(chr_snp$snp) <- str_lowercase(chr_snp$ind.info$FamID)

var <- intersect(chr_snp$snp.info$ID, tagSNP)
geno <- as.data.frame(subset(chr_snp$snp,select=var))
MAF <- as.data.frame(colMeans(geno,na.rm = T)/2)
MAF$SNP <- rownames(MAF)
colnames(MAF)[1] <- "eaf"
mQTL_MR <- merge(mQTL_MR, MAF, by="SNP")
write.table(mQTL_MR, file= "./MR/mQTLs_MR.txt",row.names=F,col.names=T,quote=T,sep='\t')

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

# data formatting
# load MR datasets
mQTL_MR <- read.table("./MR/mQTLs_MR.txt", header = T)
GWAS1_MR <- read.table("./MR/GWAS1_MR.txt", header = T)
GWAS2_MR <- read.table("./MR/GWAS2_MR.txt", header = T)

# exposure datasets
mQTL_MR <- merge(mQTL_MR, GWAS1_MR[,c(1,4,5)], by = "SNP")
colnames(mQTL_MR)[c(2,9,14,15)] <- c("exposure","pval","effect_allele","other_allele")
colnames(mQTL_MR)[c(3:15)] <- paste(colnames(mQTL_MR)[c(3:15)] , ".exposure", sep = '')
mQTL_MR$unit.exposure <- "beta value"
mQTL_MR <- mQTL_MR[,c(1:2,7:9,14:15,13,16)]
write.table(mQTL_MR, file= "./MR/mQTLs_MR_format.txt",row.names=F,col.names=T,quote=T,sep='\t')

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
write.table(GWAS1_MR, file= "./MR/GWAS1_MR_format.txt",row.names=F,col.names=T,quote=T,sep='\t')

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
write.table(GWAS2_MR, file= "./MR/GWAS2_MR_format.txt",row.names=F,col.names=T,quote=T,sep='\t')

# MR
# load MR datasets
mQTL_MR <- read_exposure_data("./MR/mQTLs_MR_format.txt", sep='\t', phenotype_col = "exposure",
                              snp_col = "SNP",
                              beta_col = "beta.exposure",
                              se_col = "se.exposure",
                              eaf_col = "eaf.exposure",
                              effect_allele_col = "effect_allele.exposure",
                              other_allele_col = "other_allele.exposure",
                              pval_col = "pval.exposure",
                              units_col = "units.exposure")
GWAS1_MR <- read_outcome_data("./MR/GWAS1_MR_format.txt", sep='\t', phenotype_col = "outcome",
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
GWAS2_MR <- read_outcome_data("./MR/GWAS2_MR_format.txt", sep='\t', phenotype_col = "outcome",
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

MR1 <- harmonise_data(mQTL_MR, GWAS1_MR, action = 1)
MR2 <- harmonise_data(mQTL_MR, GWAS2_MR, action = 1)

mr_results1 <- mr(MR1, method_list=c("mr_ivw","mr_egger_regression"))
mr_results2 <- mr(MR2, method_list=c("mr_ivw","mr_egger_regression"))

results1 <- cbind.data.frame(mr_results1$exposure,mr_results1$outcome,mr_results1$nsnp,mr_results1$method,mr_results1$b,mr_results1$se,mr_results1$pval)
results2 <- cbind.data.frame(mr_results2$exposure,mr_results2$outcome,mr_results2$nsnp,mr_results2$method,mr_results2$b,mr_results2$se,mr_results2$pval)

write.csv(results1,"./MR/MRresults/mQTL_GWAS1_MRresults.csv")
write.csv(results2,"./MR/MRresults/mQTL_GWAS2_MRresults.csv")

# sensitivity analysis
mr_het1 <- mr_heterogeneity(MR1)
nrow(mr_het1[mr_het1$Q_pval <= 0.01,])
# 6
mr_het2 <- mr_heterogeneity(MR2)
nrow(mr_het2[mr_het2$Q_pval <= 0.01,])
# 7
write.csv(mr_het1,"./MR/MRresults/mQTL_GWAS1_MRheterogeneity.csv")
write.csv(mr_het2,"./MR/MRresults/mQTL_GWAS2_MRheterogeneity.csv")

mr_pt1 <- mr_pleiotropy_test(MR1)
nrow(mr_pt1[mr_pt1$pval <= 0.01,])
# 0
mr_pt2 <-mr_pleiotropy_test(MR2)
nrow(mr_pt2[mr_pt2$pval <= 0.01,])
# 0
write.csv(mr_pt1,"./MR/MRresults/mQTL_GWAS1_MRpleiotropy.csv")
write.csv(mr_pt2,"./MR/MRresults/mQTL_GWAS2_MRpleiotropy.csv")

res_single1 <- mr_singlesnp(MR1)
res_single2 <- mr_singlesnp(MR2)

write.csv(res_single1,"./MR/MRresults/mQTL_GWAS1_MRsinglesnp.csv")
write.csv(res_single2,"./MR/MRresults/mQTL_GWAS2_MRsinglesnp.csv")

# We can also create a volcano plot for multiple MR results 
# First, we can add an Bonferroni correction onto the results (to impose a strict correction for multiple testing)
fres1 <- mr_results1
fres1 <- fres1[(fres1$method=="Inverse variance weighted" | fres1$method=="MR Egger"),]
fres1$effect <- exp(fres1$b)
fres1$up_ci <- exp(fres1$b+(1.96*fres1$se))
fres1$lo_ci <- exp(fres1$b-(1.96*fres1$se))
fres1$sample_size <- 1316
fres1<-fres1[!is.na(fres1$effect),]
fres1$bon <- p.adjust(fres1$pval, method="bonferroni")
fres1$sig <- fres1$bon < 0.05
fres1$category <-"<0.05"
write.csv(fres1,"./MR/MRresults/mQTL_GWAS1_MRresults_adjusted.csv")

ggplot(fres1, aes(x=effect, y=-log10(pval))) +
  geom_vline(xintercept=1, linetype="dotted") +
  geom_point(data=subset(fres1, !sig)) +
  geom_point(data=subset(fres1, sig), aes(colour=category, size=bon < 0.05)) +
  facet_grid(. ~ outcome, scale="free") +
  geom_label_repel(data=subset(fres1, bon<0.05), aes(label=exposure, fill=category), colour="white", segment.colour="black", point.padding = unit(0.7, "lines"), box.padding = unit(0.7, "lines"), segment.size=0.5, force=2, max.iter=3e3) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.x=element_text(size=20)
  ) +
  scale_colour_brewer(type="qual", palette="Dark2") +
  scale_fill_brewer(type="qual", palette="Dark2") +
  labs(x="OR for lung cancer per SD increase in CpG methylation", size="Bonferroni",y="P value (-log10)") +
  xlim(c(0.8, 1.1))+
  theme(axis.title.y=element_text(size=18),axis.title.x=element_text(size=18),axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),legend.position="none")
ggsave("./MR/MRresults/volcanoplot_cpg_CHD1.png", width=8, height=8)

# GWAS2
fres2 <- mr_results2
fres2 <- fres2[(fres2$method=="Inverse variance weighted" | fres2$method=="MR Egger"),]
fres2$effect <- exp(fres2$b)
fres2$up_ci <- exp(fres2$b+(1.96*fres2$se))
fres2$lo_ci <- exp(fres2$b-(1.96*fres2$se))
fres2$sample_size <- 1275
fres2<-fres2[!is.na(fres2$effect),]
fres2$bon <- p.adjust(fres2$pval, method="bonferroni")
fres2$sig <- fres2$bon < 0.05
fres2$category <-"<0.05"
write.csv(fres2,"./MR/MRresults/mQTL_GWAS2_MRresults_adjusted.csv")

ggplot(fres2, aes(x=effect, y=-log10(pval))) +
  geom_vline(xintercept=1, linetype="dotted") +
  geom_point(data=subset(fres2, !sig)) +
  geom_point(data=subset(fres2, sig), aes(colour=category, size=bon < 0.05)) +
  facet_grid(. ~ outcome, scale="free") +
  geom_label_repel(data=subset(fres2, bon<0.05), aes(label=exposure, fill=category), colour="white", segment.colour="black", point.padding = unit(0.7, "lines"), box.padding = unit(0.7, "lines"), segment.size=0.5, force=2, max.iter=3e3) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.x=element_text(size=20)
  ) +
  scale_colour_brewer(type="qual", palette="Dark2") +
  scale_fill_brewer(type="qual", palette="Dark2") +
  labs(x="OR for lung cancer per SD increase in CpG methylation", size="Bonferroni",y="P value (-log10)") +
  xlim(c(0.8, 1.1))+
  theme(axis.title.y=element_text(size=18),axis.title.x=element_text(size=18),axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),legend.position="none")
ggsave("./MR/MRresults/volcanoplot_cpg_CHD2.png", width=8, height=8)

# GWAS1
# Here we have used one of the CpG sites with the greatest number of SNPs (cg22933800) as an example.
dat1 <- MR1[MR1$exposure=="cg22933800",]
mr_ex1 <- mr(dat1, method_list=c("mr_ivw","mr_egger_regression"))
mr_ex1
single_ex1 <- mr_singlesnp(dat1)
single_ex1

png("./MR/MRresults/cg22933800_CHD1_scatter.png")
mr_scatter_plot(mr_ex1, dat1)
dev.off()

# Generate a forest plot of each of the SNP effects, which are then meta-analysed using the IVW and MR-Egger methods
png("./MR/MRresults/cg22933800_CHD1_forest.png")
mr_forest_plot(single_ex1)
dev.off()

# Generate a funnel plot to check asymmetry
png("./MR/MRresults/cg22933800_CHD1_funnel.png")
mr_funnel_plot(single_ex1)
dev.off()

# Run a leave-one-out analysis and generate a plot to test whether any one SNP is driving any pleiotropy or asymmetry in the estimates
res_loo1 <- mr_leaveoneout(dat1)
png("./MR/MRresults/cg22933800_CHD1_loo.png")
mr_leaveoneout_plot(res_loo1)
dev.off()

# GWAS2
dat2 <- MR2[MR2$exposure=="cg22933800",]
mr_ex2 <- mr(dat2, method_list=c("mr_ivw","mr_egger_regression"))
mr_ex2
single_ex2 <- mr_singlesnp(dat2)
single_ex2

png("./MR/MRresults/cg22933800_CHD2_scatter.png")
mr_scatter_plot(mr_ex2, dat2)
dev.off()

# Generate a forest plot of each of the SNP effects, which are then meta-analysed using the IVW and MR-Egger methods
png("./MR/MRresults/cg22933800_CHD2_forest.png")
mr_forest_plot(single_ex2)
dev.off()

# Generate a funnel plot to check asymmetry
png("./MR/MRresults/cg22933800_CHD2_funnel.png")
mr_funnel_plot(single_ex2)
dev.off()

# Run a leave-one-out analysis and generate a plot to test whether any one SNP is driving any pleiotropy or asymmetry in the estimates
res_loo2 <- mr_leaveoneout(dat2)
png("./MR/MRresults/cg22933800_CHD2_loo.png")
mr_leaveoneout_plot(res_loo2)
dev.off()

# table for manuscript
results1 <- read.csv("./MR/MRresults/mQTL_GWAS1_MRresults_adjusted.csv")
results2 <- read.csv("./MR/MRresults/mQTL_GWAS2_MRresults_adjusted.csv")

mr_het1 <- read.csv("./MR/MRresults/mQTL_GWAS1_MRheterogeneity.csv")
mr_het2 <- read.csv("./MR/MRresults/mQTL_GWAS2_MRheterogeneity.csv")

mr_pt1 <- read.csv("./MR/MRresults/mQTL_GWAS1_MRpleiotropy.csv")
mr_pt2 <- read.csv("./MR/MRresults/mQTL_GWAS2_MRpleiotropy.csv")

res_single1 <- read.csv("./MR/MRresults/mQTL_GWAS1_MRsinglesnp.csv")
res_single2 <- read.csv("./MR/MRresults/mQTL_GWAS2_MRsinglesnp.csv")

tbs1 <- read.csv("./MR/Supp_Table1.csv")
GWAS1 <- read.table("./MR/GWAS1_MR_format.txt", header = T)
GWAS2 <- read.table("./MR/GWAS2_MR_format.txt", header = T)
refGene <- read.table("./refGene.merged.bed", header = F)

sig1 <- results1[results1$sig == TRUE,4:17]
sig1  <- merge(sig1, mr_het1[,c(5,6,9)], by = c("exposure","method"))
sig1  <- merge(sig1, mr_pt1[,c(5,8)], by = "exposure")
sig1 <- sig1[,c(1:10,12,15:16)]
sig1 <- merge(sig1, results2[,c(5,6,15)], by = c("exposure","method"))
colnames(sig1)[c(5,7,11:14)] <- c("beta","p.MR","p.MR.bon","p.heterogeneity","p.pleiotropy","p.MR2.bon")
sig1 <- merge(sig1, res_single1[,c(2,7)], by = "exposure")
sig1 <- sig1[! sig1$SNP %in% c("All - MR Egger", "All - Inverse variance weighted"),]
sig1 <- merge(sig1, GWAS1[,c(1,4)], by = "SNP")
colnames(sig1)[16] <- "p.GWAS1"
sig1 <- merge(sig1, GWAS2[,c(1,4)], by = "SNP")
colnames(sig1)[c(2,17)] <- c("CpG_site","p.GWAS2")
sig1 <- merge(sig1, tbs1[,1:6], by = c("SNP","CpG_site"))
sig1 <- sig1[order(sig1$chr_SNP),c(1,20:21,2,18:19,3:17)]
write.csv(sig1,"./MR/MRresults/mQTL_GWAS1_MRresults_sig.csv")

sig1 <- read.csv("./MR/MRresults/mQTL_GWAS1_MRresults_sig.csv")
colnames(sig1)[12:20] <- paste(colnames(sig1)[12:20], ".GWAS1", sep = '')
colnames(sig1)[5] <- "exposure"
sig1 <- merge(sig1, results2[,c(5,6,8:13)], by = c("exposure","method"))
colnames(sig1)[24:29] <- paste(colnames(sig1)[24:29], ".GWAS2", sep = '')
sig1 <- sig1[,c(1:20,22:23,24:29,21)]
sig1  <- merge(sig1, mr_het2[,c(5,6,9)], by = c("exposure","method"))
sig1  <- merge(sig1, mr_pt2[,c(5,8)], by = "exposure")
colnames(sig1)[c(1,23,29:31)] <- c("CpG_site","beta.GWAS2","p.MR.bon.GWAS2","p.heterogeneity.GWAS2","p.pleiotropy.GWAS2")
sig1 <- sig1[,c(3:6,1,7:9,2,10:31)]
write.csv(sig1,"./MR/MRresults/mQTL_GWAS1_MRresults_sig2.csv")

sig2 <- results2[results2$sig == TRUE,4:17]
sig2  <- merge(sig2, mr_het2[,c(5,6,9)], by = c("exposure","method"))
sig2  <- merge(sig2, mr_pt2[,c(5,8)], by = "exposure")
sig2 <- sig2[,c(1:10,12,15:16)]
sig2 <- merge(sig2, results1[,c(5,6,15)], by = c("exposure","method"))
colnames(sig2)[c(5,7,11:14)] <- c("beta","p.MR","p.MR.bon","p.heterogeneity","p.pleiotropy","p.MR1.bon")
sig2 <- merge(sig2, res_single2[,c(2,7)], by = "exposure")
sig2 <- sig2[! sig2$SNP %in% c("All - MR Egger", "All - Inverse variance weighted"),]
sig2 <- merge(sig2, GWAS1[,c(1,4)], by = "SNP")
colnames(sig2)[16] <- "p.GWAS1"
sig2 <- merge(sig2, GWAS2[,c(1,4)], by = "SNP")
colnames(sig2)[c(2,17)] <- c("CpG_site","p.GWAS2")
sig2 <- merge(sig2, tbs1[,1:6], by = c("SNP","CpG_site"))
sig2 <- sig2[order(sig2$CpG_site,sig2$chr_SNP),c(1,20:21,2,18:19,3:17)]
write.csv(sig2,"./MR/MRresults/mQTL_GWAS2_MRresults_sig.csv")

sig2 <- read.csv("./MR/MRresults/mQTL_GWAS2_MRresults_sig.csv")
colnames(sig2)[12:20] <- paste(colnames(sig2)[12:20], ".GWAS2", sep = '')
colnames(sig2)[5] <- "exposure"
sig2 <- merge(sig2, results1[,c(5,6,8:13)], by = c("exposure","method"))
colnames(sig2)[24:29] <- paste(colnames(sig2)[24:29], ".GWAS1", sep = '')
sig2 <- sig2[,c(1:20,22:23,24:29,21)]
sig2  <- merge(sig2, mr_het1[,c(5,6,9)], by = c("exposure","method"))
sig2  <- merge(sig2, mr_pt1[,c(5,8)], by = "exposure")
colnames(sig2)[c(1,23,29:31)] <- c("CpG_site","beta.GWAS1","p.MR.bon.GWAS1","p.heterogeneity.GWAS1","p.pleiotropy.GWAS1")
sig2 <- sig2[,c(3:6,1,7:9,2,10:31)]
write.csv(sig2,"./MR/MRresults/mQTL_GWAS2_MRresults_sig2.csv")

# single SNP result
sig1 <- read.csv("./MR/MRresults/mQTL_GWAS1_MRresults_sig2.csv")
sig2 <- read.csv("./MR/MRresults/mQTL_GWAS2_MRresults_sig2.csv")

res_single1$effect <- exp(res_single1$b)
res_single1$up_ci <- exp(res_single1$b+(1.96*res_single1$se))
res_single1$lo_ci <- exp(res_single1$b-(1.96*res_single1$se))
res_single1 <-res_single1[!is.na(res_single1$effect),]
res_single1$bon <- p.adjust(res_single1$p, method="bonferroni")

res_single2$effect <- exp(res_single2$b)
res_single2$up_ci <- exp(res_single2$b+(1.96*res_single2$se))
res_single2$lo_ci <- exp(res_single2$b-(1.96*res_single2$se))
res_single2 <-res_single2[!is.na(res_single2$effect),]
res_single2$bon <- p.adjust(res_single2$p, method="bonferroni")

colnames(sig1)[5] <- "exposure"
single1 <- merge(sig1[,c(1:11,21:22)],res_single1[,c(2,7:13)],by = c("exposure","SNP"))
colnames(single1)[14:19] <- c("beta.GWAS1","se.GWAS1","p.MR.GWAS1","effect.GWAS1","up_ci.GWAS1","lo_ci.GWAS1")
single1 <- merge(single1,res_single2[,c(2,7:13)],by = c("exposure","SNP"), all.x = T)
colnames(single1)[c(1,20:25)] <- c("CpG_site","beta.GWAS2","se.GWAS2","p.MR.GWAS2","effect.GWAS2","up_ci.GWAS2","lo_ci.GWAS2")
single1 <- single1[,c(2:5,1,6:8,9:25)]
write.csv(single1,"./MR/MRresults/mQTL_GWAS1_MRsinglesnp_sig.csv")

colnames(sig2)[5] <- "exposure"
single2 <- merge(sig2[,c(1:11,21:22)],res_single2[,c(2,7:13)],by = c("exposure","SNP"))
colnames(single2)[14:19] <- c("beta.GWAS2","se.GWAS2","p.MR.GWAS2","effect.GWAS2","up_ci.GWAS2","lo_ci.GWAS2")
single2 <- merge(single2,res_single1[,c(2,7:13)],by = c("exposure","SNP"))
colnames(single2)[c(1,20:25)] <- c("CpG_site","beta.GWAS1","se.GWAS1","p.MR.GWAS1","effect.GWAS1","up_ci.GWAS1","lo_ci.GWAS1")
single2 <- single2[,c(2:5,1,6:8,9:25)]
write.csv(single2,"./MR/MRresults/mQTL_GWAS2_MRsinglesnp_sig.csv")

# supplemental figure 2
cpglist <- c(as.character(unique(sig1$CpG_site)), as.character(unique(sig2$CpG_site)))
pdf("./MR/MRresults/Supp_Figure2.pdf")
for (i in 1:5){
  cpg <- cpglist[i]
  dat1 <- MR1[MR1$exposure==cpg,]
  single_ex1 <- mr_singlesnp(dat1, all_method = "mr_ivw")
  print(mr_forest_plot(single_ex1))
}

for (i in 6:12){
  cpg <- cpglist[i]
  dat2 <- MR2[MR2$exposure==cpg,]
  single_ex2 <- mr_singlesnp(dat2, all_method = "mr_ivw")
  print(mr_forest_plot(single_ex2))
}

dev.off()

