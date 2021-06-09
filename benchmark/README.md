Calculate precision and recall of caller with Rscript

``` {r}
Rscript precision_recall.R manta.vcf NA12878_dgv_long_read.bed output_summary.txt

``` 

The outputs of each caller was merged in a single file `all.txt` and using `make_plot.R ` the precion-recall results could be plotted in stored in pdf file. 
