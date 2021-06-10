

samtools view -h NA12878.alt_bwamem_GRCh38DH.20150706.CEU.illumina_platinum_ped.bam | extractSplitReads_BwaMem -i stdin| samtools view -Sb - > lumpy.tmp.sr.unsort.bam 

samtools sort lumpy.tmp.sr.unsort.bam > lumpy.tmp.sr.sorted.bam 
 
samtools view -b -F 1294 NA12878.alt_bwamem_GRCh38DH.20150706.CEU.illumina_platinum_ped.bam > lumpy.tmp.d.unsort.bam  

samtools sort lumpy.tmp.d.unsort.bam > lumpy.tmp.d.sorted.bam 


lumpyexpress -B NA12878.alt_bwamem_GRCh38DH.20150706.CEU.illumina_platinum_ped.bam -S lumpy.tmp.sr.sorted.bam  -D lumpy.tmp.d.unsort.bam -o output.vcf 

svtyper -B NA12878.alt_bwamem_GRCh38DH.20150706.CEU.illumina_platinum_ped.bam -S lumpy.tmp.sr.sorted.bam -i output.vcf > NA12878.lumpy.vcf
