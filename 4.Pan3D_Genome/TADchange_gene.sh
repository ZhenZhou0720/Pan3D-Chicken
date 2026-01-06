#!/bin/bash
mkdir -p ./result

bedtools sort -i ./change_merge/Fusion_merge_change_domain.bed > fusion.sorted.bed
bedtools sort -i ./change_merge/Neo_merge_change_domain.bed > neo.sorted.bed
bedtools sort -i ./change_merge/Stable_merge_change_domain.bed > stable.sorted.bed
bedtools sort -i ./gene_loc_0based_anno.bed > gene.sorted.bed

bedtools intersect -a fusion.sorted.bed -b gene.sorted.bed -loj > ./result/Fusion_change_domain_related_gene.bed
bedtools intersect -a neo.sorted.bed -b gene.sorted.bed -loj > ./result/Neo_change_domain_related_gene.bed
bedtools intersect -a stable.sorted.bed -b gene.sorted.bed -loj > ./result/Stable_change_domain_related_gene.bed