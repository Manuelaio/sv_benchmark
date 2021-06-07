## Manta 

1) run Manta with defoult paramenter 

``` {r}
configManta.py --bam NA12878.alt_bwamem_GRCh38DH.20150706.CEU.illumina_platinum_ped.bam --referenceFasta GRCh38_full_analysis_set_plus_decoy_hla.fa --runDir Manta
python /MantaWorkflow/runWorkflow.py
``` 

2) convert vcf file in bed file 

``` {r}
Rscript vcf_to_bed.R NA12878.diploidSV.vcf NA12878.diploidSV.bed
``` 

3) remove complex regions stored in balck list 

``` {r}
intersectBed -a  NA12878.diploidSV.bed -b GRCh38_unified_blacklist.bed -v > NA12878.diploidSV.nbl.bed 

``` 

4) interesct with truset with a reciprocal overlapp 

``` {r}
intersectBed -a  NA12878.diploidSV.nbl.bed  -b NA12878_dgv_long_read.bed -f 0.50 -r -wao > NA12878.diploidSV.nbl.DGV.bed 

``` 

5) rigenotype DEL and DUP with sv2

``` {r}
sv2 -i NA12878.alt_bwamem_GRCh38DH.20150706.CEU.illumina_platinum_ped.bam -b NA12878.diploidSV.nbl.bed  -snv NA12878.snv.vcf.gz -p NA12878.ped -o ri-genotyping -g hg38 -tmp-dir 
```

6) calculate precision and recall 






