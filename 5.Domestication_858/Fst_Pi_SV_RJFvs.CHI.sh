#!/bin/bash
# Fst
vcftools --vcf sv.vcf --weir-fst-pop RJF.txt --weir-fst-pop Chicken.txt --out RJF_CHI
# pi
vcftools --vcf sv.vcf --site-pi --keep RJF.txt --out RJF
vcftools --vcf sv.vcf --site-pi --keep Chicken.txt --out CHI