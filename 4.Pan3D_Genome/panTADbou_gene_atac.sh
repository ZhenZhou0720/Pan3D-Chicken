#!/bin/bash
mkdir -p result
bedtools intersect -a Pan3D_PanTAD.Boundary.bed -b gene_loc_0based_anno.bed -loj > result/All_Pan3D_TADBoundary_Gene.txt

mkdir -p result/ATAC/npz_file 
mkdir -p result/ATAC/matrix_file 
 
samples=($(ls ./ATAC_bw/*_treat_pileup.bw | sed -E 's|./ATAC_bw/(.*)_treat_pileup\.bw|\1|g' | sort))

bw_files=()
for sample in "${samples[@]}"; do
    bw_files+=("./ATAC_bw/${sample}_treat_pileup.bw")
done

multiBigwigSummary BED-file \
    --BED Pan3D_PanTAD.Boundary.bed \
    --bwfiles "${bw_files[@]}" \
    --labels "${samples[@]}" \
    --outFileName ./result/ATAC/npz_file/Pan3D_PanTAD.Boundary.ATAC.npz \
    --outRawCounts ./result/ATAC/matrix_file/Pan3D_PanTAD.Boundary.ATAC.txt \
    --numberOfProcessors 12