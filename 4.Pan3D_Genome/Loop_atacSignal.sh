#!/bin/bash
mkdir -p result

ATAC_bw="$1"
loop_anchor_1="$2"
loop_anchor_2="$3"

samples=($(ls ./ATAC_bw/*_treat_pileup.bw | sed -E 's|./ATAC_bw/(.*)_treat_pileup\.bw|\1|g' | sort))
bw_files=()
for sample in "${samples[@]}"; do
    bw_files+=("./ATAC_bw/${sample}_treat_pileup.bw")
done

for FILE in $loop_anchor_1 $loop_anchor_2; do
    base_name=$(basename "${FILE}" .bedpe)
    prefix="${base_name#overlap_50_New_}" 

    multiBigwigSummary BED-file \
        --BED ${FILE} \
        --bwfiles "${bw_files[@]}" \
        --labels "${samples[@]}" \
        --outFileName ./result/${prefix}.ATAC.npz \
        --outRawCounts ./result/${prefix}.ATAC.txt \
        --numberOfProcessors 36
        
done
