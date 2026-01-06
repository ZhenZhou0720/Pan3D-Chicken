#!/bin/bash
gene_info="$1"
loop_anchor_1="$2"
loop_anchor_2="$3"

bedtools sort -i ""${loop_anchor_1}"" > loop_anchor_1_sorted.bedpe
bedtools sort -i ""${loop_anchor_2}"" > loop_anchor_2_sorted.bedpe
bedtools sort -i ""${gene_info}"" > gene_info_sorted.bed

bedtools intersect -a loop_anchor_1_sorted.bedpe -b gene_info_sorted.bed -loj > panloop_anchor_1_gene.xls
bedtools intersect -a loop_anchor_2_sorted.bedpe -b gene_info_sorted.bed -loj > panloop_anchor_2_gene.xls

rm loop_anchor_1_sorted.bedpe loop_anchor_2_sorted.bedpe gene_info_sorted.bed