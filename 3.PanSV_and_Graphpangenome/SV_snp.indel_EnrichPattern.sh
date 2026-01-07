#!/bin/bash
# SNP
plotHeatmap \
 -m Snp.Peak_matrix.gz \
 --colorMap viridis \
 --heatmapHeight 16 \
 --heatmapWidth 8 \
 --samplesLabel SNP \
 --xAxisLabel SV \
 --yAxisLabel SNP \
 -o SNP.pdf

# InDel
plotHeatmap \
 -m Indel.Peak_matrix.gz \
 --colorMap viridis \
 --heatmapHeight 16 \
 --heatmapWidth 8 \
 --samplesLabel InDel \
 --xAxisLabel SV \
 --yAxisLabel InDel \
 -o InDel.pdf