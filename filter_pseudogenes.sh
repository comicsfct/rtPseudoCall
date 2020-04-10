#!/bin/bash

#stop in case of error on one of the steps
set -e

annot=$1  #gtf with pseudogenes you want to know if are inside dog or not
input=$2  #output of Dogfinder Get_DoGs
output_prefix=$3  #prefix to name output files
pseudos_removed_from_annot=$4  #file with the pseudogene lines removed from gtf, you want to know if pseudogenes are inside pseudogenes on the opposite strand

if [ -z "$4" ]
  then
    pseudos_removed_from_annot=$annot
fi

#capture pseudogenes inside dogs extension or just downstream - max 30 bases distance
bedtools window -l 30 -r 30 -sw -a $input -b $annot > ${output_prefix}_RTs_ovlap_pseudos.bed

#get only overlaps between dog and pseudogene on the same strand
wc -l ${output_prefix}_RTs_ovlap_pseudos.bed
awk '{if ($6 == $10) print}' ${output_prefix}_RTs_ovlap_pseudos.bed > ${output_prefix}_RTs_same_strand_ovlap_pseudos.bed


#get pseudogenes which are between two dogs in unstranded libraries and thus it is difficult to understand if the reads in the pseudogene region are from the same strand or from the opposite
#echo ${output_prefix}
awk -F $'\t' 'BEGIN {OFS=FS} {if ($6!=$10) print}'   ${output_prefix}_RTs_ovlap_pseudos.bed | while read line; do echo $line | awk -F 'gene_name' '{print $2}' | awk '{print $1}'  | grep -f -  ${output_prefix}_RTs_same_strand_ovlap_pseudos.bed; done > ${output_prefix}_between2RTs.gtf && \
echo ${output_prefix}

#only pseudogenes
awk -F$'\t' 'BEGIN {OFS=FS} {print $7,$8,$9,$10,$11}'  ${output_prefix}_RTs_same_strand_ovlap_pseudos.bed > ${output_prefix}_pseudogenes_same_strand_ovlap.bed
#get pseudogenes which overlap with pseudogenes on the opposite strand
bedtools intersect -wa -wb -loj -a ${output_prefix}_pseudogenes_same_strand_ovlap.bed -b ${pseudos_removed_from_annot}  | awk -F $'\t' 'BEGIN {OFS=FS} {if ($4!=$9) print}' > ${output_prefix}_inside_opposite_strand_pseudogenes.gtf


#echo ${output_prefix}
#filter from the file with the overlap between pseudogenes and dogs of the same strand the pseudogenes present between 2 dogs or that overlap with pseudogenes of opposite strand
awk -F 'gene_name' '{print $2}' ${output_prefix}_inside_opposite_strand_pseudogenes.gtf | awk '{print $1}' | grep -v -f - ${output_prefix}_pseudogenes_same_strand_ovlap.bed | grep -x -v -f ${output_prefix}_between2RTs.gtf  - > ${output_prefix}_samestrand_ovlap_pseudos_filtered.bed
awk -F 'gene_name' '{print $2}' ${output_prefix}_between2RTs.gtf | awk '{print $1}' | grep -v -f - ${output_prefix}_samestrand_ovlap_pseudos_filtered.bed > ${output_prefix}_samestrand_ovlap_pseudos_filtered.temp && mv \
${output_prefix}_samestrand_ovlap_pseudos_filtered.temp ${output_prefix}_samestrand_ovlap_pseudos_filtered.bed 
