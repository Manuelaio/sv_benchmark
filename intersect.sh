#!/bin/bash 

input=$1
output=$2

intersectBed -a $input -b /work/emanuela.iovino/intersect_SV/new_analysis/bad_regions/GRCh38_unified_blacklist.bed -v > file.nbl.bed 
intersectBed -a file.nbl.bed -b /work/emanuela.iovino/intersect_SV/new_analysis/DGV_paiper/Trueset/NA12878_dgv_long_read.merged.annotated.nbl.bed -f 0.50 -r -wao > $output
