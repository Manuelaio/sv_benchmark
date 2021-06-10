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

3) remove complex regions stored in balck list 

``` {r}
intersectBed -a   NA12878.delly.bed -b GRCh38_unified_blacklist.bed -v > NA12878.delly.nbl.bed 

``` 

4) interesct with truset with a reciprocal overlapp 

``` {r}
intersectBed -a  NA12878.delly.nbl.bed -b NA12878_dgv_long_read.bed -f 0.50 -r -wao > NA12878.delly.nbl.DGV.bed 

``` 

5) rigenotype DEL and DUP with sv2

``` {r}
sv2 -i NA12878.alt_bwamem_GRCh38DH.20150706.CEU.illumina_platinum_ped.bam -b NA12878.delly.nbl.DGV.bed   -snv NA12878.snv.vcf.gz -p NA12878.ped -o ri-genotyping -g hg38 -tmp-dir 

```

6) calculate precision and recall with Rscript [../benchmarck/precision_recall.R](https://github.com/Manuelaio/sv_benchmark/blob/main/benchmark/precision_recall.R)

