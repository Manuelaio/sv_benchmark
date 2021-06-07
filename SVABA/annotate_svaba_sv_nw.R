#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
caller=args[1]
bed_ann=args[2]


vcf= caller 
tmp_vcf<-readLines(vcf)
tmp_vcf<-tmp_vcf[-(grep("#CHROM",tmp_vcf)+1):-(length(tmp_vcf))]

cols <- colnames(read.table(pipe(paste0('grep -v "##" ',vcf ,' | grep "#"| sed s/#//')), header = TRUE))
cols <- sapply(cols, function(x) gsub("(mrkdp\\.)|(\\.bam)", "", x))
svaba_uniq1 <- read.table(vcf, stringsAsFactors = FALSE)
#elimino col in piÃ¹ che vengono scritte da svaba ? giÃ  segnalato nelle issue di github 
svaba_uniq= svaba_uniq1[-c(10,11,12)]
colnames(svaba_uniq)= cols
rm(svaba_uniq1)
#substring(svaba_uniq, regexpr("SPAN=", svaba_uniq) + 1
#file=read.table("/work/emanuela.iovino/intersect_SV/new_analysis/svaba/new_ann_sv.bed")
# ann_sv_svaba.vcf
file=read.table(bed_ann)
colnames(file)[2]='POS'
library(dplyr)

df=inner_join(svaba_uniq,file, by='POS')


df1= df[c(1:12, 14)]
df.f= df1[!duplicated(df1),]
colnames(df.f)[13]="sv_type"
df.f$INFO<- paste(df.f$INFO,df.f$sv_type, sep=";SVTYPE=")
#df.f$END= df.f$ALT[gsub(".*:","", df.f$ALT),]
header=tmp_vcf[-c(length(tmp_vcf))]
write.table(header, file="header_vcf.txt", 
            quote = F, row.names = F, col.names = T)


write.table(df.f, file="svaba_for_surv.vcf", 
            quote = F, row.names = F, col.names = T, sep="\t")

