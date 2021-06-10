library(dplyr)
merge_path= "merge.vcf"
cols <- colnames(read.table(pipe(paste0('grep -v "##" ',merge_path ,' | grep "#"| sed s/#//')), header = TRUE))
cols <- sapply(cols, function(x) gsub("(mrkdp\\.)|(\\.bam)", "", x))
df <- read.table(merge_path, stringsAsFactors = FALSE)
colnames(df)=cols
vcf=read.vcfR(merge_path)
#file=df[grepl("0/1",df$NA12878_30x) | ("0/1",df$NA12878_30x_1) ,]

het=with(df, df[ grepl( '0/1', NA12878_30x) | grepl( '0/1', NA12878_30x_1) | grepl( '0/1', NA12878_30x_2) | 
                   grepl( '0/1', NA12878_30x_3) | grepl('0/1', NA12878_30x_4),])



df$ann=with(df, ifelse(grepl( '0/1', df$NA12878_30x) | grepl( '0/1', df$NA12878_30x_1) | grepl( '0/1', df$NA12878_30x_2) | 
                   grepl( '0/1', df$NA12878_30x_3) | grepl('0/1', df$NA12878_30x_4) | grepl( '1/1', df$NA12878_30x) | grepl( '1/1', df$NA12878_30x_1) | grepl( '1/1', df$NA12878_30x_2) | 
                     grepl( '1/1', df$NA12878_30x_3) | grepl('1/1', df$NA12878_30x_4), "probabile", "incerto" ))
df$FILTER<- paste(df$FILTER,df$ann, sep=";Qual=")

colnames(df)[1]="#CHROM"
df_r= df[-c(15)]
write.table(df_r, file="merge.ann.vcf", quote=F, sep="\t",row.names = F, col.names = T)

