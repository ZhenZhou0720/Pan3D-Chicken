#!/bin/bash
mkdir -p result

ATAC_bw="$1"
CRE_file="$2"

samples=($(ls ./ATAC_bw/*_treat_pileup.bw | sed -E 's|./ATAC_bw/(.*)_treat_pileup\.bw|\1|g' | sort))
bw_files=()
for sample in "${samples[@]}"; do
    bw_files+=("./ATAC_bw/${sample}_treat_pileup.bw")
done

for FILE in $CRE_file; do
    prefix=$(basename "${FILE}" .bed)
    
    multiBigwigSummary BED-file \
        --BED ${FILE} \
        --bwfiles "${bw_files[@]}" \
        --labels "${samples[@]}" \
        --outFileName ./result/${prefix}.ATAC.npz \
        --outRawCounts ./result/${prefix}.ATAC.txt \
        --numberOfProcessors 36
        
done
