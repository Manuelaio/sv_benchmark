library(dplyr)
library(GenomicRanges)
library(IRanges)

######____________________________________________________________________________________________________________#####
######required ucsc file of Seg_Dup(Segmental_dups_hg38_frt_srt) and rep_mask regions (RepeatMasker_hg38.frt.bed) ###### 
######____________________________________________________________________________________________________________#####


goldset_tmp= read.table("./NA12878_dgv_long_read.merged.tmp.bed", stringsAsFactors = F)
goldset_tmp$V5=sapply(strsplit(goldset_tmp$V5,","), `[`, 1)

goldset_tmp$V6=sapply(strsplit(as.character(goldset_tmp$V6), ","), function(x) paste(unique(x), collapse=","))


goldset_tmp$V4=goldset_tmp$V3-goldset_tmp$V2
goldset= goldset_tmp[c(1,2,3,5,6,4)]
colnames(goldset)= c("chrom","start","end","type", "source","len")

goldset_grange= makeGRangesFromDataFrame(goldset, ignore.strand = T, keep.extra.columns = T)

seg_dup= read.table("Segmental_dups_hg38_frt_srt", stringsAsFactors = F, col.names = c("chrom", "start","end"))
seg_dup_gr= makeGRangesFromDataFrame(seg_dup, ignore.strand = T, keep.extra.columns = T)

##Intersect goldset and segmental duplications 

tp= findOverlaps(query= goldset_grange, subject = seg_dup_gr, type="any")
intersect=data.frame(goldset_grange[queryHits(tp),],seg_dup_gr[queryHits(tp),])
ann_seg_dup= intersect[!duplicated(intersect[1:4]),]
ann_seg_dup$segdup= "SEG_DUP"
SG_info= ann_seg_dup[c(1,2,3,6,7,8,14)]


golset_ann1= dplyr::left_join(goldset,SG_info)

## Intersect gold set and repeted regions 

rep_mask= read.table("RepeatMasker_hg38.frt.bed", stringsAsFactors = F, header = T)
colnames(rep_mask)=c("chrom","start","end","repFamily")
rep_mask_gr= makeGRangesFromDataFrame(rep_mask, ignore.strand = T, keep.extra.columns = T)

tp2= findOverlaps(query= goldset_grange, subject = rep_mask_gr, type="any")
intersect_RR=data.frame(goldset_grange[queryHits(tp2),],rep_mask_gr[queryHits(tp2),])
ann_RR= intersect_RR[!duplicated(intersect_RR[1:4]),]
RR_info= ann_RR[c(1,2,3,6,7,8,14)]


golset_ann2= dplyr::left_join(golset_ann1,RR_info)

golset_ann= golset_ann2[-c(7)]


write.table(golset_ann, file="NA12878_dgv_long_read.merged.annotated.bed", quote = F, sep="\t", col.names = F, row.names = F)

