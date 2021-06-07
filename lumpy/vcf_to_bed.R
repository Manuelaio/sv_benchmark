#!/usr/bin/env Rscript

#####Rscript 


library(dplyr)
library(tidyverse)
library(VariantAnnotation)

args = commandArgs(trailingOnly=TRUE)
lumpy<- args[1]

output<- args[2]
vcf = readVcf(lumpy, "hg38")

df1=as.data.frame(rowRanges(vcf), row.names = NULL) #col 1,2,10 : chrm start filter
df1_subset= df1[c(1,2)]
df2=as.data.frame(vcf@info) #col 2,3,4 type, len, end
df2_subset= df2[c(1:3)]
geno_df= as.data.frame((geno(vcf)$GT),row.names = NULL)
file_finale= cbind(df1_subset,df2_subset,geno_df)
row.names(file_finale)=NULL
file_finale[file_finale$SVTYPE=="DEL",]$SVLEN= as.numeric(file_finale[file_finale$SVTYPE=="DEL",]$SVLEN)*(-1) 

#pass= file_finale%>%
 # dplyr::filter(FILTER=="PASS")
#pass$SVLEN= pass$END - pass$start

col_order <- c("seqnames", "start", "END",
               "SVTYPE", "SVLEN", "NA12878_30x")
my_data2 <- file_finale[, col_order]

my_data_final= my_data2[my_data2$SVTYPE!="BND",]


write.table(my_data_final,file=output, sep="\t",quote= F, row.names = F, col.names = F)
