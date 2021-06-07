#!/usr/bin/env Rscript

#####Rscript 


library(dplyr)
library(tidyverse)
library(VariantAnnotation)

args = commandArgs(trailingOnly=TRUE)
erds<- args[1]

output<- args[2]
vcf = readVcf(erds, "hg38")
df1=as.data.frame(rowRanges(vcf), row.names = NULL) #col 1,2,10 : chrm start filter
df1_subset= df1[c(1,2,10)]
df2=as.data.frame(vcf@info) #col 2,3,4 type, len, end
df2_subset= df2[c(1,3:4)]
cols <- colnames(read.table(pipe(paste0('grep -v "##" ',erds ,' | grep "#"| sed s/#//')), header = TRUE))
cols <- sapply(cols, function(x) gsub("(mrkdp\\.)|(\\.bam)", "", x))
erds_cn <- read.table(erds, stringsAsFactors = FALSE)


file_finale= cbind(df1_subset,df2_subset)


col_order <- c("seqnames", "start", "END",
               "SVTYPE", "SVLEN")
my_data2 <- file_finale[, col_order]
my_data2[my_data2$SVTYPE=="DEL",]$SVLEN= as.numeric(my_data2[my_data2$SVTYPE=="DEL",]$SVLEN)*(-1) 
my_data2$geno= "."


write.table(my_data2,file=output, sep="\t",quote= F, row.names = F, col.names = F)
