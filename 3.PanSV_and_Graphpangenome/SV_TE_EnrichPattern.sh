#!/bin/bash
# TE
plotHeatmap \
 -m TE.Peak_matrix.gz \
 --colorMap plasma \
 --heatmapHeight 16 \
 --heatmapWidth 8 \
 --samplesLabel TE \
 --xAxisLabel SV \
 --yAxisLabel TE \
 -o TE.pdf

# LINE
plotHeatmap \
 -m LINE.Peak_matrix.gz \
 --colorMap plasma \
 --heatmapHeight 16 \
 --heatmapWidth 8 \
 --samplesLabel LINE \
 --xAxisLabel SV \
 --yAxisLabel LINE \
 -o LINE.pdf

# SINE
plotHeatmap \
 -m SINE.Peak_matrix.gz \
 --colorMap plasma \
 --heatmapHeight 16 \
 --heatmapWidth 8 \
 --samplesLabel SINE \
 --xAxisLabel SV \
 --yAxisLabel SINE \
 -o SINE.pdf

# LTR
plotHeatmap \
 -m LTR.Peak_matrix.gz \
 --colorMap plasma \
 --heatmapHeight 16 \
 --heatmapWidth 8 \
 --samplesLabel LTR \
 --xAxisLabel SV \
 --yAxisLabel LTR \
 -o LTR.pdf

# DNA 
plotHeatmap \
 -m DNA.Peak_matrix.gz \
 --colorMap plasma \
 --heatmapHeight 16 \
 --heatmapWidth 8 \
 --samplesLabel DNA \
 --xAxisLabel SV \
 --yAxisLabel DNA \
 -o DNA.pdf