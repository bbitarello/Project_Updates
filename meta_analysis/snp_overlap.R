#!/usr/bin/env Rscript
library(data.table)
library(tictoc)
library(dplyr)
tic('Run script')
tic('Declare variables')
path="to_delete/"
file1="QC_HEIGHT_BBJ_hg38.txt.gz"
file2="1000g_all_EAS_hg38_QC.bim"
file3="HRS_AFR_hg38.bim"
file4=gsub(".bim", ".afreq",file2)
toc()
tic('read in all files')
files<-c(paste0(path, c(file1,file2,file3)), paste0('tmp-dir/', file4))
snps1<-fread(files[1], header=T)$SNP
snps2<-fread(files[2], header=F)[,V2]
snps3<-fread(files[3], header=F)[,V2]
snps4<-fread(files[4], header=T)
toc()
tic('check first overlap')
snps1_2<-intersect(snps1,snps2)
toc()
tic('check total overlap')
snps1_2_3<-intersect(snps1_2,snps3)
toc()
cat(paste0('overlap for ss and LD panel is: ', length(snps1_2), '\n'))
cat(paste0('overlap for ss and LD panel and test data is: ', length(snps1_2_3), '\n'))
tic('Frequency of overlapping varaints')
snps1<-fread(files[1], header=T)
snps1[, overlap:=ifelse(SNP %in% snps1_2_3, TRUE,FALSE)]
snps1<-snps1[overlap==TRUE]
snps4<-snps4[ID %in% snps1_2_3][,MAF:=ifelse(ALT_FREQS>0.5, 1-ALT_FREQS, ALT_FREQS)]
setnames(snps4, 'ID', 'SNP')
final<-merge(snps1, snps4, by='SNP')
final[,DiffMAF:=abs(MAF.x-MAF.y)]

nrow(final[DiffMAF<0.1 & MAF.x>= 0.001])
toc()
toc()
