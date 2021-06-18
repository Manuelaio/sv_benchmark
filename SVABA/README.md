1) RUN svaba with bash script [run_svaba.sh](https://github.com/Manuelaio/sv_benchmark/blob/main/SVABA/run_svaba.sh) which includes Svaba command, a perl script shared by [Kasugi](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1720-5) and a Rscript to obtain a annotated vcf. 


2) Remove black list regions and intersect witih Gold set 

``` {r}

intersectBed -a NA12878.rearrange.bed -b GRCh38_unified_blacklist.bed -v > NA12878.rearrange.nbl.bed

intersectBed -a NA12878.rearrange.nbl.bed -b NA12878_dgv_long_read.bed -f 0.50 -r -wao >  NA12878.rearrange.nbl.DGV.bed

``` 

3) rigenotype SVs 

``` {r}
sv2 -i NA12878.alt_bwamem_GRCh38DH.20150706.CEU.illumina_platinum_ped.bam -b  NA12878.rearrange.nbl.bed  -snv NA12878.snv.vcf.gz -p NA12878.ped -o ri-genotyping -g hg38 -tmp-dir 

``` 

4) calculate precision and recall with Rscript [../benchmarck/precision_recall.R](https://github.com/Manuelaio/sv_benchmark/blob/main/benchmark/precision_recall.R)

