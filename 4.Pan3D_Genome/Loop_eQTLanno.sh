#!/bin/bash
QTL_info="$1"
loop_anchor_1="$2"
loop_anchor_2="$3"

bedtools sort -i ""${loop_anchor_1}"" > loop_anchor_1_sorted.bedpe
bedtools sort -i ""${loop_anchor_2}"" > loop_anchor_2_sorted.bedpe
bedtools sort -i ""${QTL_info}"" > QTL_info_sorted.bed

bedtools intersect -a loop_anchor_1_sorted.bedpe -b QTL_info_sorted.bed -loj > panloop_anchor_1_QTL.xls
bedtools intersect -a loop_anchor_2_sorted.bedpe -b QTL_info_sorted.bed -loj > panloop_anchor_2_QTL.xls
rm loop_anchor_1_sorted.bedpe loop_anchor_2_sorted.bedpe QTL_info_sorted.bed