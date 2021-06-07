## Delly evaluation 

1) run delly caller with default parameter 

``` {r}
delly call -g GRCh38_full_analysis_set_plus_decoy_hla.fa -o  NA12878//tmp.delly.bcf  -x human.hg38.excl.tsv NA12878.alt_bwamem_GRCh38DH.20150706.CEU.illumina_platinum_ped.bam
bcftools view tmp.delly.bcf > NA12878.delly.vcf
``` 

2) convert vcf file in bed for intersect

``` {r}
Rscript vcf_to_bed_delly.R NA12878.delly.vcf NA12878.delly.bed
``` 
3) remove complex regions 

