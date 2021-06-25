#!/usr/bin/env Rscript

####write bash script for SV pipeline
####ex: Rscript sv_pipeline.R -i /Users/emanuelaiovino/Desktop/bam/12302_S18.bam -g hg38 -o /Users/emanuelaiovino/Desktop/


library(optparse)

option_list = list(
    make_option(c("-i", "--input"), type = "character", default = NULL, help = "Input BAM file"),
    make_option(c("-g", "--genome"), type = "character", default = NULL, help = "genome build (hg19 or hg38)"),
    make_option(c("-v", "--snv"), type = "character", default = NULL, help = "SNV vcf "),
    make_option(c("-p", "--ped"), type = "character", default = NULL, help = "ped file"),
    make_option(c("-o", "--output"), type = "character", default = "no_id", help = "Output folder")
)

parseobj = OptionParser(option_list=option_list)
opt = parse_args(parseobj)

if (is.null(opt$input))
    stop(print_help(parseobj))

if (!tolower(opt$genome) %in% c("hg38","hg19"))
    stop("Must specify --genome as hg19 or hg38")



assembly= opt$genome
input=opt$input
snv=opt$snv
dr=opt$output
ped=opt$ped

if(assembly=="hg38"){
   ref="/archive/ngsbo/db/trioCEU_1KGP_resources/GRCh38_full_analysis_set_plus_decoy_hla.fa"}else{ref="/archive/ngsbo/db/hg19/ucsc.hg19.fasta"}

if(assembly=="hg38"){
    ex_site="/archive/ngsbo/db/regions/exclude_site_delly/human.hg38.excl.tsv "}else{ex_site="/archive/ngsbo/db/regions/exclude_site_delly/human.hg19.excl.tsv "}


folder<- paste0(dr,"/delly/")

dir.create(folder)
tmp_d<- paste0(folder,"/tmp/")

dir.create(tmp_d)
srt=paste0 ("#!/bin/bash \n")
srtD=paste0("### DELLY COMMAND \n")
delly=paste0("delly call -g ", ref, " -o  ", tmp_d, "/tmp.delly.bcf ", " -x ", ex_site, " ", input)
srtD1=paste0("###  in ngsra env \n")
dellyGT= paste0("delly call -g ", ref, " -v ", tmp_d, "/tmp.delly.bcf ", " -x ", ex_site, " -o ", tmp_d, "/tmp.delly.gene.bcf ", input )
c=paste0("bcftools view ", tmp_d, "/tmp.delly.gene.bcf > ", folder, "output.delly.vcf")

com_d= c(srt,srtD,delly,srtD1,dellyGT,c)

##1° start
writeLines(com_d)
f.d <- paste0(folder,"/", Sys.Date(), "_sv_delly.sh")
file.create(f.d)
f <- file(f.d, open="w") # or open="a" if appending
write.table(com_d, file = f, sep = " ", append=FALSE, row.names = FALSE, col.names=FALSE, quote=FALSE)
close(f)
##end

##2.start
folder_sv2<- paste0(folder,"/sv2")
dir.create(folder_sv2)
folder_sv2_tmp<- paste0(folder_sv2,"/tmp")
dir.create(folder_sv2_tmp)
bed_d= paste0(folder,"/file.nbl.bed")
sv2.dl=paste0("sv2 -i ", input, " -b ", bed_d, " -snv ", snv, " -p ", ped," -O ",folder_sv2 ,  " -o genotyping -g ", assembly, " -tmp-dir " , folder_sv2_tmp )
writeLines(sv2.dl)
f.sv2 <- paste0(folder,"/", Sys.Date(), "_sv2_sv_delly.sh")
file.create(f.sv2)
dlsv2 <- file(f.sv2, open="w") # or open="a" if appending
write.table(sv2.dl, file = dlsv2, sep = " ", append=FALSE, row.names = FALSE, col.names=FALSE, quote=FALSE)
close(dlsv2)
##end

file.copy("/work/emanuela.iovino/intersect_SV/benchmark/script/vcf_rearrange_delly.R", folder)
file.copy("/work/emanuela.iovino/intersect_SV/benchmark/script/intersect.sh", folder)

#_______________________________________________________________________________#

folder_m<- paste0(dr,"/manta")
dir.create(folder_m)

tmp_l<- paste0(folder_m,"/tmp/")
dir.create(tmp_l)

#1.start
srtM=paste0("### Manta COMMAND - in svtools env \n")
manta=paste0("configManta.py --bam ", input, " --referenceFasta ", ref,  " --runDir ", tmp_l)
manta_py= paste0("python ", tmp_l, "/runWorkflow.py")
com_manta= c(srtM,manta,manta_py)
writeLines(com_manta)
f.m <- paste0(folder_m,"/", Sys.Date(), "_sv_manta.sh")
file.create(f.m)
f.manta <- file(f.m, open="w") # or open="a" if appending
write.table(com_manta, file = f.manta, sep = " ", append=FALSE, row.names = FALSE, col.names=FALSE, quote=FALSE)
close(f.manta)
#end

#2.start
folder_sv2_m<- paste0(folder_m,"/sv2")
dir.create(folder_sv2_m)
folder_sv2_tmp_m<- paste0(folder_sv2_m,"/tmp")
dir.create(folder_sv2_tmp_m)
bed_m= paste0(f.manta,"/file.nbl.bed")

sv2.mt=paste0("sv2 -i ", input, " -b ", bed_m, " -snv ", snv, " -p ", ped," -O ",folder_sv2_m ,  " -o genotyping -g ", assembly, " -tmp-dir " , folder_sv2_tmp_m )
writeLines(sv2.mt)
f.sv2_m <- paste0(folder_m,"/", Sys.Date(), "_sv2_sv_manta.sh")
file.create(f.sv2_m)
Mtsv2 <- file(f.sv2_m, open="w") # or open="a" if appending
write.table(sv2.mt, file = Mtsv2, sep = " ", append=FALSE, row.names = FALSE, col.names=FALSE, quote=FALSE)
close(Mtsv2)
#end 

file.copy("/work/emanuela.iovino/intersect_SV/benchmark/script/vcf_rerrange_manta.R", folder_m)
file.copy("/work/emanuela.iovino/intersect_SV/benchmark/script/intersect.sh", folder_m)
#_________________________________________________________________________________________________________#

folder_s<- paste0(dr,"/svaba")
dir.create(folder_s)

sv_run<- paste0(folder_s,"/svaba_run")
dir.create(sv_run)

#1°start
strS=paste0("#### SVABA COMMAND - in svaba env \n")
svaba=paste0("svaba run -t ", input, " -p 8 -I -a ", sv_run, "/germline_run -G ", ref, "\n" )
ann_svaba=paste0("perl /work/emanuela.iovino/svaba_results/ri_download_NA12878/convert_SvABA_vcf.pl ", 
                 sv_run, "/germline_run.svaba.sv.vcf > ", sv_run,"/ann_svaba_sv.vcf \n")
#Rscript= paste0("Rscript /work/emanuela.iovino/svaba_results/ri_download_NA12878/annotate_svaba_sv_nw.R ",
 #               folder_s,"/germline_run.svaba.sv.vcf ",folder_s,"/ann_svaba_sv.vcf")
#merge=paste0("cat ",folder_s ,"/heder.txt ", folder_s,"/ann_svaba_sv.vcf > ", folder_s, "/ann_sv_svaba_header.vcf")
com_svaba= c(strS,svaba,ann_svaba)
writeLines(com_svaba)
f.svaba <- paste0(folder_s,"/", Sys.Date(), "_sv_svaba.sh")
file.create(f.svaba)
f_svaba <- file(f.svaba, open="w") # or open="a" if appending
write.table(com_svaba, file = f_svaba, sep = " ", append=FALSE, row.names = FALSE, col.names=FALSE, quote=FALSE)
close(f_svaba)
#end

#2.start
folder_sv2_s<- paste0(folder_s,"/sv2")
dir.create(folder_sv2_s)
folder_sv2_tmp_s<- paste0(folder_sv2_s,"/tmp")
dir.create(folder_sv2_tmp_s)
bed_s= paste0(folder_s,"/file.nbl.bed")

sv2.sv=paste0("sv2 -i ", input, " -b ", bed_s, " -snv ", snv, " -p ", ped," -O ",folder_sv2_s ,  " -o genotyping -g ", assembly, " -tmp-dir " , folder_sv2_tmp_s )
writeLines(sv2.sv)
f.sv2_s <- paste0(folder_s,"/", Sys.Date(), "_sv2_sv_svaba.sh")
file.create(f.sv2_s)
SVsv2 <- file(f.sv2_s, open="w") # or open="a" if appending
write.table(sv2.sv, file = SVsv2, sep = " ", append=FALSE, row.names = FALSE, col.names=FALSE, quote=FALSE)
close(SVsv2)
#end 

file.copy("/work/emanuela.iovino/intersect_SV/benchmark/script/vcf_rearrange_svaba.R", folder_s)
file.copy("/work/emanuela.iovino/intersect_SV/benchmark/script/intersect.sh", folder_s)

#_______________________________________________________________________________________________#
folder_l<- paste0(dr,"/lumpy")

dir.create(folder_l)
tmp_l=  paste0(folder_l,"/tmp/")
dir.create(tmp_l)

srtlumpy=paste0("### Lumpy COMMAND - in svtools env \n")
lumpy_1s=paste0("samtools view -h ", input, " | extractSplitReads_BwaMem -i stdin| samtools view -Sb - > ", tmp_l, "/lumpy.tmp.sr.unsort.bam \n")
lumpy_2s= paste0("samtools sort ", tmp_l, "/lumpy.tmp.sr.unsort.bam > ", tmp_l,"/lumpy.tmp.sr.sorted.bam \n ")
lumpy_3s= paste0("samtools view -b -F 1294 ",  input, " > ", tmp_l, "/lumpy.tmp.d.unsort.bam  \n")
lumpy_4= paste0("samtools sort ", tmp_l,"/lumpy.tmp.d.unsort.bam > ", tmp_l, "/lumpy.tmp.d.sorted.bam \n")
lumpy5=paste0("samtools view ", input, " | tail -n+100000 | 
              /shared/conda/miniconda3/pkgs/lumpy-sv-0.3.0-py27hfbaaabd_6/share/lumpy-sv-0.3.0-6/scripts/pairend_distro.py 
              -r 100 -X 4 -N 10000 -o sample.histo \n")
lumpy6= paste0("lumpyexpress -B ", input, " -S ",tmp_l, "/lumpy.tmp.sr.sorted.bam ", " -D ", tmp_l, "/lumpy.tmp.d.unsort.bam",
               " -o ", tmp_l, "/output.vcf \n")
lumpy7= paste0("svtyper -B ", input, " -S ", tmp_l, "/lumpy.tmp.sr.sorted.bam -i ",tmp_l, "/output.vcf > ", folder_l, "/lumpy.vcf" )
com_lumpy= c(srtlumpy, lumpy_1s, lumpy_2s, lumpy_3s, lumpy_4, lumpy5,lumpy6,lumpy7)

#1.
writeLines(com_lumpy)
f.lum <- paste0(folder_l,"/", Sys.Date(), "_sv_lumpy.sh")
file.create(f.lum)
fL <- file(f.lum, open="w") # or open="a" if appending
write.table(com_lumpy, file = fL, sep = " ", append=FALSE, row.names = FALSE, col.names=FALSE, quote=FALSE)
close(fL)
#end

#2.start
folder_sv2_l<- paste0(folder_l,"/sv2")
dir.create(folder_sv2_l)
folder_sv2_tmp_l<- paste0(folder_sv2_l,"/tmp")
dir.create(folder_sv2_tmp_l)
bed_l= paste0(folder_l,"/file.nbl.bed")

sv2.lp=paste0("sv2 -i ", input, " -b ", bed_l, " -snv ", snv, " -p ", ped," -O ",folder_sv2_l ,  " -o genotyping -g ", assembly, " -tmp-dir " , folder_sv2_tmp_l )
writeLines(sv2.lp)
f.sv2_l <- paste0(folder_l,"/", Sys.Date(), "_sv2_sv_lumpy.sh")
file.create(f.sv2_l)
LPsv2 <- file(f.sv2_l, open="w") # or open="a" if appending
write.table(sv2.lp, file = LPsv2, sep = " ", append=FALSE, row.names = FALSE, col.names=FALSE, quote=FALSE)
close(LPsv2)
#end 

file.copy("/work/emanuela.iovino/intersect_SV/benchmark/script/vcf_rearrange_lumpy.R", folder_l)
file.copy("/work/emanuela.iovino/intersect_SV/benchmark/script/intersect.sh", folder_l)

#____________________________________________________________________________________________#

srtERDS=paste0("### ERDS COMMAND - in erds env \n")

folder_e<- paste0(dr,"/erds")

dir.create(folder_e)

erds_1s=paste0("perl /shared/conda/miniconda3/envs/erds/share/erds-1.1-1/erds_pipeline.pl -b ", input, " -v " ,snv, " -r ",  ref, " -o " , folder_e, "/")
#1.
writeLines(erds_1s)
f.e <- paste0(folder_e,"/", Sys.Date(), "_sv_ERDS.sh")
file.create(f.e)
fIL <- file(f.e, open="w") # or open="a" if appending
write.table(erds_1s, file = f.e, sep = " ", append=FALSE, row.names = FALSE, col.names=FALSE, quote=FALSE)
close(fIL)
#end

#2.
folder_sv2_e<- paste0(folder_e,"/sv2")
dir.create(folder_sv2_e)
folder_sv2_tmp_e<- paste0(folder_sv2_e,"/tmp")
dir.create(folder_sv2_tmp_e)
bed_e= paste0(folder_e,"/file.nbl.bed")

sv2.er=paste0("sv2 -i ", input, " -b ", bed_e, " -snv ", snv, " -p ", ped," -O ",folder_sv2_e ,  " -o genotyping -g ", assembly, " -tmp-dir " , folder_sv2_tmp_e )
writeLines(sv2.er)
f.sv2_e <- paste0(folder_e,"/", Sys.Date(), "_sv2_sv_erds.sh")
file.create(f.sv2_e)
ERsv2 <- file(f.sv2_e, open="w") # or open="a" if appending
write.table(sv2.er, file = ERsv2, sep = " ", append=FALSE, row.names = FALSE, col.names=FALSE, quote=FALSE)
close(ERsv2)
#end

file.copy("/work/emanuela.iovino/intersect_SV/benchmark/script/vcf_rearrange_erds.R", folder_e)
file.copy("/work/emanuela.iovino/intersect_SV/benchmark/script/intersect.sh", folder_e)





