#!/usr/bin/env Rscript

#####Rscript 


library(dplyr)
library(tidyverse)
library(VariantAnnotation)

args = commandArgs(trailingOnly=TRUE)
erds<- args[1]

output<- args[2]

SG<- args[3]

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



my_data2$SVLEN=ifelse(my_data2$SVLEN=="integer(0)", "0", my_data2$SVLEN)
my_data2$SVLEN=as.numeric(unlist(my_data2$SVLEN))

my_data_gr=makeGRangesFromDataFrame(my_data2, ignore.strand = T, keep.extra.columns = T)

seg_dup= read.table(SG, stringsAsFactors = F, col.names = c("chrom", "start","end"))
seg_dup_gr= makeGRangesFromDataFrame(seg_dup, ignore.strand = T, keep.extra.columns = T)

##Intersect goldset and segmental duplications 

tp= findOverlaps(query= my_data_gr, subject = seg_dup_gr, type="any")
intersect=data.frame(my_data_gr[queryHits(tp),],seg_dup_gr[queryHits(tp),])
if(nrow(intersect)!=0){
ann_seg_dup= intersect[!duplicated(intersect[1:4]),]
ann_seg_dup$segdup= "SEG_DUP"
SG_info= ann_seg_dup[c(1,2,3,6,7,8,14)]
my_data_final_SG= dplyr::left_join(my_data_final,SG_info)
}else{my_data_final$segdup="."
my_data_final_SG=my_data_final
my_data_final_SG$SG="."
}

my_data_final_SG$RR<-"no_annotated"
ann=my_data_final_SG[-c(7)]

write.table(ann,file=output, sep="\t",quote= F, row.names = F, col.names = F)