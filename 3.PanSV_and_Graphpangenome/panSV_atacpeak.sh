#!/bin/bash
for sv_file in breed_sv/*_sv_CI.bed; do
    breed=$(basename "$sv_file" | cut -d'_' -f1)
    mkdir -p "result/${breed}"
    breed_peak="peak/${breed}_MergePeak.bed"
    bedtools intersect -a "$sv_file" -b "$breed_peak" -wo > "result/${breed}/${breed}_sv_atac_result.txt"
    awk -v breed="$breed" '{print breed "\t" $(NF-1)}' "result/${breed}/${breed}_sv_atac_result.txt" | sort -u > "result/${breed}/${breed}_sv_atac_list.txt"
done