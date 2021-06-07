#!/usr/bin/env Rscript

#####Rscript 
###use: Rscript NA12878.diploidSV.vcf NA12878.diploidSV.bed


library(dplyr)
library(tidyverse)
library(VariantAnnotation)

args = commandArgs(trailingOnly=TRUE)
manta<- args[1]
output<- args[2]
vcf = readVcf(manta, "hg38")

df1=as.data.frame(rowRanges(vcf), row.names = NULL) #col 1,2,10 : chrm start filter
df1_subset= df1[c(1,2,10)]
df2=as.data.frame(vcf@info) #col 2,3,4 type, len, end
df2_subset= df2[c(2:4)]
geno_df= as.data.frame((geno(vcf)$GT),row.names = NULL)
file_finale= cbind(df1_subset,df2_subset,geno_df)
row.names(file_finale)=NULL
pass= file_finale%>%
  dplyr::filter(FILTER=="PASS")

col_order <- c("seqnames", "start", "END",
                "SVTYPE", "SVLEN", "NA12878_30x")
my_data2 <- pass[, col_order]
my_data2[my_data2$SVTYPE=="DEL",]$SVLEN= as.numeric(my_data2[my_data2$SVTYPE=="DEL",]$SVLEN)*(-1) 
my_data_final= my_data2[my_data2$SVTYPE!="BND",]


write.table(my_data_final,file=output, sep="\t",quote= F, row.names = F, col.names = F)
