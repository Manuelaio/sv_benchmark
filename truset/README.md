**Truset**

A reference dataset corresponding to NA12878 was obtained from a previous published work of [Kosugi](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1720-5), et al. 2019. 
A reference calset includes:
 
 1. [DGV](http://dgv.tcag.ca/dgv/app/downloads?ref=GRCh37/hg19) variants (DEL,DUP, INS, INV)
 2. [PacBio](https://www.nature.com/articles/nmeth.3454) SVs identifed with long read (DEL, INS)
 3. High confidence NA12878 SV set of the [svclassify](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-2366-2) study (DEL, INS)
 4. Verified nonredundant INV from long read study stored in [InvFEST](http://invfestdb.uab.cat/#:~:text=The%20InvFEST%20database%20stores%20and,the%20resolution%20of%20each%20study.) database (INV)
 
 
 
For each type genomic positions was merged and annotated with segmental duplication and repeted regions using bash script [merge.sh](https://github.com/Manuelaio/sv_benchmark/blob/main/truset/merge.sh)

```r
bash merge.sh NA12878_dgv_long_read.bed

```
The outfile was intersect with ucsc black list with [bedtools](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html) with following command

```r
intersectBed -a NA12878_dgv_long_read.merged.annotated.bed -b GRCh38_unified_blacklist.bed -v > NA12878_dgv_long_read.merged.annotated.nbl.bed

```

The difference between raw file and merged file is resumed in following table 



|    SVS type    |   No merge   |   Merged     |
| :------------- | :----------: | -----------: |
|  DEL           | 9042         | 7292         |
|  DUP           | 2512         | 1755         |
|  INS           | 13525        | 12894        |
|  INV           | 258          |  173         |



A schematic flow is shown in following image 


![alt text](https://github.com/Manuelaio/sv_benchmark/blob/main/truset/refSV.PNG)

