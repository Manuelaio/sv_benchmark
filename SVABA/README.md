1) RUN svaba with bash scritp 

2) Run Perl script

3) Run Rscript to rearrange vcf and obtain bed file for intersection

4) Remove black list regions 


intersectBed -a NA12878.rearrange.bed -b GRCh38_unified_blacklist.bed -v > NA12878.rearrange.nbl.bed
