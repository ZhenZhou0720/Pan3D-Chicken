#!/bin/bash
for sv_file in breed_sv/*_sv_CI.bed; do
    breed=$(basename $sv_file | cut -d'_' -f1)
    mkdir -p "result/${breed}"
    bedtools slop -i "$sv_file" -g chrom_sizes.txt -b 5000 > "result/${breed}/${breed}_sv_extended.bed"
    bedtools intersect -a "result/${breed}/${breed}_sv_extended.bed" -b "gene/${breed}_gene.bed" -wo > "result/${breed}/${breed}_sv_gene_result.txt"
    awk -v breed="$breed" '{print breed "\t" $(NF-1)}' "result/${breed}/${breed}_sv_gene_result.txt" | sort -u > "result/${breed}/${breed}_sv_gene_list.txt"
done

