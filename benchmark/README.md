Calculate precision and recall of caller with Rscript

``` {r}
Rscript precision_recall.R ../Manta/NA12878.diploidSV.nbl.DGV.bed  NA12878_dgv_long_read.bed output_summary_manta.txt
Rscript precision_recall.R ../delly/NA12878.delly.nbl.DGV.bed  NA12878_dgv_long_read.bed output_summary_delly.txt
Rscript precision_recall.R ../lumpy/NA12878.delly.nbl.DGV.bed  NA12878_dgv_long_read.bed output_summary_lumpy.txt
Rscript precision_recall.R ../ERDS/NA12878_30x.erds.nbl.DGV.bed  NA12878_dgv_long_read.bed output_summary_erds.txt
Rscript precision_recall.R ../SVABA/NA12878.rearrange.nbl.DGV.bed  NA12878_dgv_long_read.bed output_summary_svaba.txt
Rscript precision_recall.R ../canvas/NA12878.rearrange.nbl.DGV.bed  NA12878_dgv_long_read.bed output_summary_canvas.txt

``` 

The outputs of each caller was merged in a file both of [all SVs](https://github.com/Manuelaio/sv_benchmark/blob/main/benchmark/all_result.txt) and in a file of [cleaned SVs](https://github.com/Manuelaio/sv_benchmark/blob/main/benchmark/all_no_SD.txt) from segmental duplication and using `make_plot.R ` the precion-recall results could be plotted in stored in pdf file. 

The resulting SV calls of the 5 detection tools were merged into one  set using SURVIVOR after SV2 genotyping 

``` {r}
SURVIVOR merge merge.list 100 1 0 0 0 0 merge.vcf

``` 

VCF file is annotated in FILTER column using [vcf_annotation.R](https://github.com/Manuelaio/sv_benchmark/blob/main/benchmark/vcf_annotation.R) wiht additional information:

*  *likely true* when the variants is genotyped as 0/1 or 1/1


*  *uncertain* when the variants is not genotyped, "./." or is genotyped as 0/0




