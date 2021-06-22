#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
caller=args[1]
bed_ann=args[2]
output= args[3]

vcf= caller 
tmp_vcf<-readLines(vcf)
tmp_vcf<-tmp_vcf[-(grep("#CHROM",tmp_vcf)+1):-(length(tmp_vcf))]

cols <- colnames(read.table(pipe(paste0('grep -v "##" ',vcf ,' | grep "#"| sed s/#//')), header = TRUE))
cols <- sapply(cols, function(x) gsub("(mrkdp\\.)|(\\.bam)", "", x))
svaba_uniq1 <- read.table(vcf, stringsAsFactors = FALSE)
#elimino col in più che vengono scritte da svaba ? già segnalato nelle issue di github 
svaba_uniq= svaba_uniq1[-c(10,11,12)]
colnames(svaba_uniq)= cols
rm(svaba_uniq1)
#substring(svaba_uniq, regexpr("SPAN=", svaba_uniq) + 1
file=read.table("ann_svaba_sv.vcf",stringsAsFactors = F )
file=read.table(bed_ann, stringsAsFactors = F)
colnames(file)[2]='POS'
###rearrange vcf 

info_field=data.frame(do.call('rbind', strsplit(as.character(file$V8),';',fixed=TRUE)))
file_df= cbind(file, info_field)
file_df$X2= gsub("SVLEN=", "", file_df$X2)
file_df$X4= gsub("GT=", "", file_df$X4)
file_df$X2=as.numeric(file_df$X2)

file_df$END= file_df$POS + file_df$X2
bed= file_df[c(1:2,13,3,10,12)]

write.table(bed, file=output, 
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


write.table(df1, file="svaba_for_surv.vcf", 
            quote = F, row.names = F, col.names = T, sep="\t")

