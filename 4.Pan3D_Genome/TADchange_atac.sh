#!/bin/bash
mkdir -p result/ 
 
samples=($(ls ./ATAC_bw/*_treat_pileup.bw | sed -E 's|./ATAC_bw/(.*)_treat_pileup\.bw|\1|g' | sort))

bw_files=()
for sample in "${samples[@]}"; do
    bw_files+=("./ATAC_bw/${sample}_treat_pileup.bw")
done

for EVENT in "Stable" "Neo" "Fusion"; do

    multiBigwigSummary BED-file \
        --BED ./change_merge/${EVENT}_merge_change_domain.bed \
        --bwfiles "${bw_files[@]}" \
        --labels "${samples[@]}" \
        --outFileName ./result/${EVENT}_domain.ATAC.npz \
        --outRawCounts ./result/${EVENT}_domain.ATAC.txt \
        --numberOfProcessors 12
        
done


