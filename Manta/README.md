## Manta 

1) run Manta with defoult paramenter 

``` {r}
configManta.py --bam NA12878.alt_bwamem_GRCh38DH.20150706.CEU.illumina_platinum_ped.bam --referenceFasta GRCh38_full_analysis_set_plus_decoy_hla.fa --runDir Manta
python /MantaWorkflow/runWorkflow.py
``` 

2) convert vcf file in bed file 

``` {r}
Rscript NA12878.diploidSV.vcf NA12878.diploidSV.bed
``` 
