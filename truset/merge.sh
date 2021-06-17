#!/bin/bash
truset=$1
awk 'BEGIN{OFS="\t"}{if($4=="INS")print $0}' $truset | sort -k 1,1V -k2,2n | mergeBed -c 4,4,5 -o count,collapse,collapse > Trueset_Ins_merged.bed
awk 'BEGIN{OFS="\t"}{if($4=="DEL")print $0}' $truset | sort -k 1,1V -k2,2n | mergeBed -c 4,4,5 -o count,collapse,collapse > Trueset_Del_merged.bed
awk 'BEGIN{OFS="\t"}{if($4=="DUP")print $0}' $truset | sort -k 1,1V -k2,2n | mergeBed -c 4,4,5 -o count,collapse,collapse > Trueset_Dup_merged.bed
awk 'BEGIN{OFS="\t"}{if($4=="INV")print $0}' $truset | sort -k 1,1V -k2,2n | mergeBed -c 4,4,5 -o count,collapse,collapse > Trueset_Dup_merged.bed
awk 'BEGIN{OFS="\t"}{if($4=="INV")print $0}' $truset | sort -k 1,1V -k2,2n | mergeBed -c 4,4,5 -o count,collapse,collapse > Trueset_Inv_merged.bed
awk 'BEGIN{OFS="\t"}{if($4=="DUP")print $0}' $truset | sort -k 1,1V -k2,2n | mergeBed -c 4,4,5 -o count,collapse,collapse > Trueset_Dup_merged.bed

cat Trueset*> NA12878_dgv_long_read.merged.tmp.bed

rm Trueset*
 
Rscript annotate.R
