Calculate precision and recall of caller with Rscript

``` {r}
Rscript precision_recall.R manta.vcf NA12878_dgv_long_read.bed output_summary.txt

``` 

The outputs of each caller was merged in a file both of [all SVs](https://github.com/Manuelaio/sv_benchmark/blob/main/benchmark/all_result.txt) and in a file of [cleaned SVs](https://github.com/Manuelaio/sv_benchmark/blob/main/benchmark/all_no_SD.txt) from segmental duplication and using `make_plot.R ` the precion-recall results could be plotted in stored in pdf file. 

The resulting SV calls of the 5 detection tools were merged into one  set using SURVIVOR after SV2 genotyping 

``` {r}
SURVIVOR merge merge.list 100 1 0 0 0 0 merge.vcf

``` 

VCF file is annotated in FILTER column using [vcf_annotation.R](https://github.com/Manuelaio/sv_benchmark/blob/main/benchmark/vcf_annotation.R) wiht additional information:

*  *likely true* when the variants is genotyped as 0/1 or 1/1


*  *uncertain* when the variants is not genotyped, "./." or is genotyped as 0/0




