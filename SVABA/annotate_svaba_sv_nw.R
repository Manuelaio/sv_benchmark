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
#elimino columns in svaba not necessary
               
svaba_uniq= svaba_uniq1[-c(10,11,12)]
colnames(svaba_uniq)= cols
rm(svaba_uniq1)


# ann_sv_svaba.vcf
file=read.table(bed_ann, stringsAsFactors = F)
colnames(file)[2]='POS'
library(dplyr)

df=inner_join(file, svaba_uniq, by='POS')

a= df[!duplicated(df[1:7]),]
df1= a[c(9:17,3)]
#df.f= df1[!duplicated(df1),]
colnames(df1)[10]="sv_type"
colnames(df1)[1]="#CHROM"
df1$INFO<- paste(df1$INFO,df1$sv_type, sep=";SVTYPE=")
#df.f$END= df.f$ALT[gsub(".*:","", df.f$ALT),]



write.table(df1, file="/work/emanuela.iovino/intersect_SV/new_analysis/svaba/NA12878/svaba_for_surv.vcf", 
            quote = F, row.names = F, col.names = T, sep="\t")

