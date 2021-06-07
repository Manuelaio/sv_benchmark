svaba run -t NA12878.alt_bwamem_GRCh38DH.20150706.CEU.illumina_platinum_ped.bam -p 8 -I -a germline_run -G GRCh38_full_analysis_set_plus_decoy_hla.fa

##annotation BND
perl convert_SvABA_vcf.pl germline_run.svaba.sv.vcf > ann_svaba_sv.vcf 


Rscript annotate_svaba_sv_nw.R /germline_run.svaba.sv.vcf ann_svaba_sv.vcf

cat heder.txt ann_svaba_sv.vcf > ann_sv_svaba_header.vcf
