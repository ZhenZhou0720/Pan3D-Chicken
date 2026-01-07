#!/bin/bash
bedtools slop -i $1 -g chrom_sizes.txt -b 5000 > slop.bed
bedtools intersect -a slop.bed -b $2 -loj > $3
rm slop.bed