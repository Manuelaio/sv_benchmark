1) RUN svaba with bash scritp 

2) Run Perl script

3) Run Rscript to rearrange vcf and obtain bed file for intersection

4) Remove black list regions and intersect witih Gold set 

``` {r}

intersectBed -a NA12878.rearrange.bed -b GRCh38_unified_blacklist.bed -v > NA12878.rearrange.nbl.bed

intersectBed -a NA12878.rearrange.nbl.bed -b NA12878_dgv_long_read.bed -f 0.50 -r -wao >  NA12878.rearrange.nbl.DGV.bed

``` 

