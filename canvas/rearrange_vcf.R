#!/usr/bin/env Rscript

#####Rscript 


library(dplyr)
library(tidyverse)
library(VariantAnnotation)

args = commandArgs(trailingOnly=TRUE)
canvas<- args[1]
output<- args[2]
vcf = readVcf(canvas, "hg38")

df1=as.data.frame(rowRanges(vcf), row.names = NULL) #col 1,2,10 : chrm start filter
df1_subset= df1[c(1,2,8)]
df1_subset$ALT= gsub("<", "", df1_subset$ALT)
df1_subset$ALT= gsub(">", "", df1_subset$ALT)

df2=as.data.frame(vcf@info) #col 2,3,4 type, len, end
df2_subset= df2[c(1, 4)]
geno_df= as.data.frame((geno(vcf)$GT),row.names = NULL)
file_finale= cbind(df1_subset,df2_subset,geno_df)
row.names(file_finale)=NULL
#file_finale[file_finale$SVTYPE=="DEL",]$SVLEN= as.numeric(file_finale[file_finale$SVTYPE=="DEL",]$SVLEN)*(-1) 

#pass= file_finale%>%
# dplyr::filter(FILTER=="PASS")
#pass$SVLEN= pass$END - pass$start

col_order <- c("seqnames", "start", "END",
               "ALT", "REFLEN", "NA12878_30x")
my_data2 <- file_finale[, col_order]

my_data_final= my_data2[my_data2$ALT!="BND",]


my_data_gr=makeGRangesFromDataFrame(my_data_final, ignore.strand = T, keep.extra.columns = T)

seg_dup= read.table("Segmental_dups_hg38_frt_srt", stringsAsFactors = F, col.names = c("chrom", "start","end"))
seg_dup_gr= makeGRangesFromDataFrame(seg_dup, ignore.strand = T, keep.extra.columns = T)

##Intersect goldset and segmental duplications 

tp= findOverlaps(query= my_data_gr, subject = seg_dup_gr, type="any")
intersect=data.frame(my_data_gr[queryHits(tp),],seg_dup_gr[queryHits(tp),])
ann_seg_dup= intersect[!duplicated(intersect[1:4]),]
ann_seg_dup$segdup= "SEG_DUP"
SG_info= ann_seg_dup[c(1,2,3,6,7,8,14)]

my_data_final_SG= dplyr::left_join(my_data2,SG_info)
my_data_final_SG$RR<-"no_annotated"
ann=my_data_final_SG[-c(7)]


write.table(ann,file="NA12878.rearrange.bed", sep="\t",quote= F, row.names = F, col.names = F)
