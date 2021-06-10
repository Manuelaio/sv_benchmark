1) RUN svaba with bash script [run_svaba.sh](https://github.com/Manuelaio/sv_benchmark/blob/main/SVABA/run_svaba.sh) 

2) Run Perl script

3) Run Rscript to rearrange vcf and obtain bed file for intersection

4) Remove black list regions and intersect witih Gold set 

``` {r}

intersectBed -a NA12878.rearrange.bed -b GRCh38_unified_blacklist.bed -v > NA12878.rearrange.nbl.bed

intersectBed -a NA12878.rearrange.nbl.bed -b NA12878_dgv_long_read.bed -f 0.50 -r -wao >  NA12878.rearrange.nbl.DGV.bed

``` 

5) 6) calculate precision and recall with Rscript [../benchmarck/precision_recall.R](https://github.com/Manuelaio/sv_benchmark/blob/main/benchmark/precision_recall.R)

