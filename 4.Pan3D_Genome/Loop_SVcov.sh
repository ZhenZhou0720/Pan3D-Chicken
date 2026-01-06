#!/bin/bash
set -e
mkdir -p result

Loop_path="$1"
SV_path="$2"

samples=$(ls "${Loop_path}" | grep 'overlap_50_New_panloopAnchor1.bedpe' | sed 's/\.overlap_50_New_panloopAnchor1\.bedpe//' | sort | uniq)

for sample in ${samples}; do
    loop_anchor1_bedpe="${Loop_path}/${sample}.overlap_50_New_panloopAnchor1.bedpe"
    loop_anchor2_bedpe="${Loop_path}/${sample}.overlap_50_New_panloopAnchor2.bedpe"
    sv_bed="${SV_path}/${sample}_sv_CI.bed"
    bedtools sort -i "${sv_bed}" > "tmp_${sample}_sv.sorted.bed"
    bedtools sort -i "${loop_anchor1_bedpe}" > "tmp_${sample}_loop_anchor1.sorted.bed"
    bedtools sort -i "${loop_anchor2_bedpe}" > "tmp_${sample}_loop_anchor2.sorted.bed"
    bedtools coverage -a "tmp_${sample}_loop_anchor1.sorted.bed" -b "tmp_${sample}_sv.sorted.bed" > "result/${sample}_sv_loopAnchor1_coverage.txt"
    bedtools coverage -a "tmp_${sample}_loop_anchor2.sorted.bed" -b "tmp_${sample}_sv.sorted.bed" > "result/${sample}_sv_loopAnchor2_coverage.txt"
    rm "tmp_${sample}_sv.sorted.bed" "tmp_${sample}_loop_anchor1.sorted.bed" "tmp_${sample}_loop_anchor2.sorted.bed"
done

