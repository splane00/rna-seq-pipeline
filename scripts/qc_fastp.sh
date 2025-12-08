#!/usr/bin/env bash
set -euo pipefail

SAMPLES="metadata/samples.tsv"
TRIM_DIR="data/trimmed"
QC_DIR="data/qc"

ADAPTER_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

mkdir -p "$TRIM_DIR" "$QC_DIR"

# sample   condition   fastq1   fastq2
tail -n +2 "$SAMPLES" | while IFS=$'\t' read -r SAMPLE CONDITION FASTQ1 FASTQ2; do

    echo "Processing $SAMPLE"
    echo "  R1: $FASTQ1"
    echo "  R2: $FASTQ2"

    fastp \
        -i "$FASTQ1" \
        -I "$FASTQ2" \
        -o "$TRIM_DIR/${SAMPLE}_R1.trimmed.fastq.gz" \
        -O "$TRIM_DIR/${SAMPLE}_R2.trimmed.fastq.gz" \
        --adapter_sequence "$ADAPTER_R1" \
        --adapter_sequence_r2 "$ADAPTER_R2" \
        --thread 4 \
        --json "$QC_DIR/${SAMPLE}_fastp.json" \
        --html "$QC_DIR/${SAMPLE}_fastp.html"

done
