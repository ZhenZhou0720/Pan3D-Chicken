#!/bin/bash
set -e
mkdir -p result

TAD_path="$1"
SV_path="$2"
MERGE_TAD_BED="$3"
MERGE_SV_BED="$4"

samples=$(find "${TAD_path}" -maxdepth 1 -type f -name "*_Pan3D_PanTAD.Boundary.bed" | \
          sed -E 's/^.*\/([^\/]+)_Pan3D_PanTAD\.Boundary\.bed$/\1/' | sort -u)

for sample in ${samples}; do
    TAD_bed="${TAD_path}/${sample}_Pan3D_PanTAD.Boundary.bed"
    sv_bed="${SV_path}/${sample}_sv_CI.bed"
    bedtools sort -i "${sv_bed}" > "tmp_${sample}_sv.sorted.bed"
    bedtools sort -i "${TAD_bed}" > "tmp_${sample}_TAD.sorted.bed"
    bedtools coverage -a "tmp_${sample}_TAD.sorted.bed" -b "tmp_${sample}_sv.sorted.bed" > "result/${sample}_sv_TAD_coverage.txt"
    rm "tmp_${sample}_sv.sorted.bed" "tmp_${sample}_TAD.sorted.bed"
done

bedtools sort -i "${MERGE_TAD_BED}" > Pan3D_PanTAD.Boundary.sored.bed
bedtools sort -i "${MERGE_SV_BED}" > sv.sorted.bed
bedtools coverage -a Pan3D_PanTAD.Boundary.sored.bed -b sv.sorted.bed > result/Pan3D_PanTAD.Boundary.SVCoverage.txt
rm Pan3D_PanTAD.Boundary.sored.bed sv.sorted.bed