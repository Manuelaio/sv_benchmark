#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
caller=args[1]
bed_ann=args[2]
output= args[3]
SG= args[4]

library(dplyr)
library(tidyverse)
library(VariantAnnotation)
vcf= caller 
tmp_vcf<-readLines(vcf)
tmp_vcf<-tmp_vcf[-(grep("#CHROM",tmp_vcf)+1):-(length(tmp_vcf))]

cols <- colnames(read.table(pipe(paste0('grep -v "##" ',vcf ,' | grep "#"| sed s/#//')), header = TRUE))
cols <- sapply(cols, function(x) gsub("(mrkdp\\.)|(\\.bam)", "", x))
svaba_uniq1 <- read.table(vcf, stringsAsFactors = FALSE)
#elimino col in pi? che vengono scritte da svaba ? gi? segnalato nelle issue di github 
svaba_uniq= svaba_uniq1[-c(10,11,12)]
colnames(svaba_uniq)= cols
rm(svaba_uniq1)
#substring(svaba_uniq, regexpr("SPAN=", svaba_uniq) + 1
file=read.table(bed_ann,stringsAsFactors = F )
#file=read.table(bed_ann, stringsAsFactors = F)
colnames(file)[2]='POS'
###rearrange vcf 

info_field=data.frame(do.call('rbind', strsplit(as.character(file$V8),';',fixed=TRUE)))
file_df= cbind(file, info_field)
file_df$X2= gsub("SVLEN=", "", file_df$X2)
file_df$X4= gsub("GT=", "", file_df$X4)
file_df$X2=as.numeric(file_df$X2)

file_df$END= file_df$POS + file_df$X2
bed= file_df[c(1:2,13,3,10,12)]
colnames(bed)=c("chrom", "start", "end", "type", "len", "geno")

my_data_gr=makeGRangesFromDataFrame(bed, ignore.strand = T, keep.extra.columns = T)

seg_dup= read.table(SG, stringsAsFactors = F, col.names = c("chrom", "start","end"))
seg_dup_gr= makeGRangesFromDataFrame(seg_dup, ignore.strand = T, keep.extra.columns = T)

##Intersect goldset and segmental duplications 

tp= findOverlaps(query= my_data_gr, subject = seg_dup_gr, type="any")
intersect=data.frame(my_data_gr[queryHits(tp),],seg_dup_gr[queryHits(tp),])
if(nrow(intersect)!=0){
  ann_seg_dup= intersect[!duplicated(intersect[1:4]),]
  ann_seg_dup$segdup= "SEG_DUP"
  SG_info= ann_seg_dup[c(1,2,3,6,7,8,14)]
  my_data_final_SG= dplyr::left_join(bed,SG_info)
}else{bed$segdup="."
my_data_final_SG=bed
my_data_final_SG$SG="."
}



my_data_final_SG$RR<-"no_annotated"
ann=my_data_final_SG[-c(7)]


write.table(ann, file=output, 
            quote = F, row.names = F, col.names = F, sep="\t")

# ann_sv_svaba.vcf


library(dplyr)

df=inner_join(file, svaba_uniq, by='POS')

a= df[!duplicated(df[1:7]),]
df1= a[c(9:17,3)]
#df.f= df1[!duplicated(df1),]
colnames(df1)[10]="sv_type"
colnames(df1)[1]="#CHROM"
df1$INFO<- paste(df1$INFO,df1$sv_type, sep=";SVTYPE=")
#df.f$END= df.f$ALT[gsub(".*:","", df.f$ALT),]
#header=tmp_vcf[-c(length(tmp_vcf))]
#write.table(header, file="/work/emanuela.iovino/intersect_SV/new_analysis/svaba/trio//header_vcf.txt", 
#           quote = F, row.names = F, col.names = T)


write.table(df1, file="/work/emanuela.iovino/intersect_SV/new_analysis/svaba/NA12878/svaba_for_surv.vcf", 
            quote = F, row.names = F, col.names = T, sep="\t")

