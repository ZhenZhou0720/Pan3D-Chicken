#!/bin/bash
SV_FILE="$1"
loop_anchor_1="$2"
loop_anchor_2="$3"

awk '
    BEGIN {
        FS = "\t"
        OFS = "\t"
    }
    NF >= 3 {
        if ($1 == "39") $1 = "W";
        if ($1 == "40") $1 = "Z";
        $1 = "chr" $1;
        print $0
    }
' "${SV_FILE}" > sv_chr.bed

bedtools sort -i "${loop_anchor_1}" > loop_anchor_1_sorted.bedpe
bedtools sort -i "${loop_anchor_2}" > loop_anchor_2_sorted.bedpe
bedtools sort -i sv_chr.bed > sv_sorted.bed

bedtools intersect -a sv_sorted.bed -b loop_anchor_1_sorted.bedpe -wa -wb > 863sv_loop_anchor1.xls
bedtools intersect -a sv_sorted.bed -b loop_anchor_2_sorted.bedpe -wa -wb > 863sv_loop_anchor2.xls

rm loop_anchor_1_sorted.bedpe loop_anchor_2_sorted.bedpe sv_sorted.bed  