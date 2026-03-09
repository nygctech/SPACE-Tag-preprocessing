#!/usr/bin/env bash
source activate spatial_multiome

ref_peaks=$1

# Loop through the remaining arguments (starting from the second argument)
for arg in "${@:2}"; do


  bamCoverage --bam data/0901-1-2_S7_Lall_R1_001.sorted.bam -o data/0901-1-2_S7_Lall_R1_001.bw --normalizeUsing RPGC --effectiveGenomeSize 2150570000 --ignoreForNormalization chrX
done

computeMatrix reference-point --referencePoint center -b 2000 -a 2000 -R results/macs2_callpeak_dedup/2_S7_Lall_R1_001_peaks.narrowPeak -S data/0901-1-2_S7_Lall_R1_001.bw data/SRR12607042_4.fastq.bw data/0919-neb-2_S2_Lall_R1_001.bw data/230509-v8-1_S1_R2_001.trimmed.bw data/SRR16298194_2.trimmed.bw data/231123-12_S9_R2_001.trimmed.bw data/ENCFF713EIC.bigWig --skipZeros -o processed_data/matrix_all_dedup.gz -p 4

plotHeatmap -m processed_data/matrix_all_dedup.gz -out results/coverage_new_dedup.pdf --colorMap Reds --missingDataColor 1

multiBigwigSummary bins -b data/0901-1-2_S7_Lall_R1_001.bw data/SRR12607042_4.fastq.bw data/0919-neb-2_S2_Lall_R1_001.bw data/230509-v8-1_S1_R2_001.trimmed.bw data/SRR16298194_2.trimmed.bw data/231123-12_S9_R2_001.trimmed.bw data/ENCFF713EIC.bigWig -o processed_data/signal_dedup.npz -p 4

plotCorrelation -in processed_data/signal_dedup.npz -c spearman -p heatmap -o results/corr_heatmap_dedup.pdf