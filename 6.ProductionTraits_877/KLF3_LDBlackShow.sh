#!/bin/bash
phenotypes=("LW" "EW" "HEW" "DW" "LMW")
for pheno in "${phenotypes[@]}"
do
    LDBlockShow -InVCF ../merge.vcf.gz \
     -OutPut KLF3_${pheno} \
     -InGWAS ../gwas_file/${pheno}.sv.snp.manha_combined_data.txt \
     -TopSite 4:69216244 \
     -Cutline 4.211254 \
     -InGFF ../reference/GGswu.gene.gff3 \
     -Region 4:69206244:69226244 \
     -OutPdf \
     -SeleVar 1
done