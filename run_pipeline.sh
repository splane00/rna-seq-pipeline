#!/usr/bin/env bash
set -euo pipefail

#############################################
# RNA-seq Pipeline Master Runner
# Description:
#   Orchestrates the full workflow:
#   1. QC + trimming (fastp)
#   2. MultiQC aggregation
#   3. Salmon index
#   4. Salmon quantification
#   5. DESeq2 input prep (optional placeholder)
#############################################

# Colors for pretty output
GREEN="\033[0;32m"
YELLOW="\033[1;33m"
NC="\033[0m"

echo -e "${GREEN}=== RNA-seq Pipeline Started ===${NC}"

###########################
# Function: check tool
###########################
check_tool() {
    if ! command -v "$1" &> /dev/null; then
        echo -e "${YELLOW}[ERROR] Missing required tool: $1${NC}"
        echo "Install it and rerun the pipeline."
        exit 1
    fi
}

###########################
# Check dependencies
###########################
echo -e "${GREEN}Checking required tools...${NC}"
for tool in fastp multiqc salmon; do
    check_tool "$tool"
done

###########################
# Step 1: QC + trimming
###########################
echo -e "${GREEN}Step 1/4: Running fastp QC + trimming...${NC}"
bash scripts/qc_fastp.sh

###########################
# Step 2: MultiQC aggregation
###########################
echo -e "${GREEN}Step 2/4: Running MultiQC...${NC}"
multiqc data/qc -o data/qc

###########################
# Step 3: Salmon index
###########################
REFERENCE="reference/transcripts.fa"

if [ ! -d reference/salmon_index ]; then
    echo "Salmon index not found — building it..."
    salmon index -t reference/transcripts.fa -i reference/salmon_index --threads 8 --gencode
else
    echo "Salmon index already exists — skipping."
fi

echo -e "${GREEN}Step 3/4: Building Salmon index...${NC}"
bash scripts/salmon_index.sh

###########################
# Step 4: Salmon quantification
###########################
echo -e "${GREEN}Step 4/4: Running Salmon quantification...${NC}"
bash scripts/salmon_quant.sh

###########################
# Optional DESeq2 placeholder
###########################
if [ -f "scripts/deseq2_prep.R" ]; then
    echo -e "${GREEN}Running optional DESeq2 preprocessing...${NC}"
    Rscript scripts/deseq2_prep.R
fi

echo -e "${GREEN}=== Pipeline complete! ===${NC}"
echo "All results available in: results/"
