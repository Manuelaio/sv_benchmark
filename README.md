This is a repository for SVs benchmark and SVs calling from WGS data
Six of most representative tools for SVs detection was included in our analysis **Manta, Delly, Svaba, ERDS, canvas and Lumpy**. Each caller was run with default parameters as report in their main page and extended commands is reported in each folder of this repository. 

The input file of caller is [NA12878.alt_bwamem_GRCh38DH.20150706.CEU.illumina_platinum_ped.cram](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/illumina_platinum_pedigree/data/CEU/NA12878/alignment/NA12878.alt_bwamem_GRCh38DH.20150706.CEU.illumina_platinum_ped.cram) of [1000 genomes](https://www.internationalgenome.org/).


SVs calling (green), processing (orange) and benchmarking (blu) was resumed in following flow.
![alt text](https://github.com/Manuelaio/sv_benchmark/blob/main/flow_sv.jpg)

The required input files are:
  1. Reference file 
  2. BAM file 
  3. SNVs vcf 
  
The calling step produces vcf output for each caller while the processing steps rearrange, annotate, re-genotype and merge calls for prioritization and for benchmark. 
