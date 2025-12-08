#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------------------------
# MultiQC Aggregation Script
# Aggregates fastp and Salmon QC reports into a single HTML.
# ---------------------------------------------------------

INPUT_DIR="results/qc"
OUTPUT_DIR="results/qc"

echo "Running MultiQC..."
multiqc "$INPUT_DIR" -o "$OUTPUT_DIR"

echo "MultiQC report generated at: $OUTPUT_DIR/multiqc_report.html"
