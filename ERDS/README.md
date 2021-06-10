## ERDS 

1) run ERDS with defoult paramenter 

``` {r}
perl erds_pipeline.pl -b NA12878.alt_bwamem_GRCh38DH.20150706.CEU.illumina_platinum_ped.bam -v NA12878.vcf.gz -r GRCh38_full_analysis_set_plus_decoy_hla.fa -o erds/

```


2) convert vcf file in bed file 

``` {r}
Rscript vcf_to_bed.R NA12878_30x.erds.vcf NA12878_30x.erds.bed
``` 

3) remove complex regions stored in balck list 

``` {r}
intersectBed -a  NA12878_30x.erds.bed -b GRCh38_unified_blacklist.bed -v > NA12878_30x.erds.nbl.bed 

``` 

4) interesct with truset with a reciprocal overlapp 

``` {r}
intersectBed -a  NA12878_30x.erds.nbl.bed  -b NA12878_dgv_long_read.bed -f 0.50 -r -wao > NA12878_30x.erds.nbl.DGV.bed 

``` 


6) calculate precision and recall with Rscript [../benchmarck/precision_recall.R](https://github.com/Manuelaio/sv_benchmark/blob/main/benchmark/precision_recall.R)
