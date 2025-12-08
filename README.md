# Comprehensive RNA-seq Pipeline with fastp, Salmon, and DESeq2 (SRP075484)
Samantha Lane
M.S. Bioinformatics student, Johns Hopkins University

## Dataset Overview

GSE81698 – “Epigenetic targeting of immune checkpoint PD-L1 by BET bromodomain inhibition”
Identifiers	SRA: SRP075484
BioProject: PRJNA322328
GEO: GSE81698
Study Type	Transcriptome Analysis

Abstract	Epigenetic regulators have emerged as exciting targets for cancer therapy. Additionally, restoration of antitumor immunity by blocking the PD-L1 signaling using antibodies has proven to be beneficial in cancer therapy. Here we show that BET bromodomain inhibition suppresses PD-L1 expression and restores antitumor immunity in ovarian cancer. CD274 (encoding PD-L1) is a direct target of BRD4-mediated gene transcription. In mouse models, treatment with the BET inhibitor JQ1 significantly reduced PD-L1 expression on tumor cells and tumor-associated dendritic cells and macrophages, which correlated with an increase in the activity of antitumor cytotoxic T cells. Together, these data demonstrate an epigenetic approach to block PD-L1 signaling to restore antitumor immunity. Given the fact that BET inhibitors have been proven safe with manageable reversible toxicity in clinical trials, our findings indicate that pharmacological BET inhibitors represent a novel treatment strategy for targeting PD-L1 expression. Overall design: RNA-seq for JQ1 treated and shBRD4 knockdown cells with controls

🧬 RNA-seq Differential Expression Pipeline

Dataset: SRP075484 | Status: Raw → QC → Trimming → Salmon Quant Complete

This repository contains a modular and reproducible RNA-seq differential expression pipeline built using public data from SRP075484 (BioProject PRJNA322328). The project includes four treatment-condition RNA-seq samples (no explicit control). The pipeline uses modern, lightweight tools (fastp, Salmon, MultiQC) and follows current best practices for quantification-based RNA-seq workflows.



## Repository Structure (So Far)
rna-seq-pipeline/
│
├── data/
│   ├── raw/                     # raw FASTQs from SRA
│   ├── trimmed/                 # fastp-trimmed FASTQs
│   ├── quant/                   # Salmon quantification outputs
│   └── reference/               # GENCODE transcriptome + Salmon index
│
└── README.md

## Tools Used
Step	Tool	Purpose
QC	FastQC	Per-sample quality inspection
QC summary	MultiQC	Aggregated QC reporting
Trimming	fastp	Adapter trimming & quality filtering
Quantification	Salmon (v1.10.3)	Alignment-free transcript quantification

The Salmon installation lives in a dedicated conda environment:
conda activate salmon_env


## Pipeline Components

The four paired-end samples from SRP075484 were downloaded using:
fasterq-dump --split-files SRR15074527
fasterq-dump --split-files SRR15074528
fasterq-dump --split-files SRR15074529
fasterq-dump --split-files SRR15074530
Outputs stored in data/raw/

### Quality Control (FastQC + MultiQC)
- FastQC run on all raw FASTQs: fastqc *.fastq --outdir ../qc/raw_fastqc
- MultiQC summary (pending MultiQC install): multiqc ../qc/raw_fastqc -o ../qc/raw_multiqc

### Trimming (fastp)
Each sample was adapter- and quality-trimmed using fastp:
fastp \
  -i SRR15074527_1.fastq \
  -I SRR15074527_2.fastq \
  -o ../trimmed/SRR15074527_1.trimmed.fastq \
  -O ../trimmed/SRR15074527_2.trimmed.fastq \
  --detect_adapter_for_pe \
  --html ../trimmed/SRR15074527.fastp.html \
  --json ../trimmed/SRR15074527.fastp.json
All trimmed reads stored in: data/trimmed/
MultiQC-ready fastp reports also generated.

### Reference Preparation (GENCODE v45)
Downloaded:
gencode.v45.transcripts.fa.gz (used for Salmon index)
gencode.v45.annotation.gtf (used later for tximport + DESeq2)

Stored in: data/reference/

### Salmon Indexing
Index built using:
salmon index \
  -t gencode.v45.transcripts.fa.gz \
  -i salmon_index \
  -k 31
Index stored in:data/reference/salmon_index/

### Salmon Quantification (Alignment-Free)
Quantified all 4 paired-end samples:
for SRR in SRR15074527 SRR15074528 SRR15074529 SRR15074530
do
  salmon quant \
    -i /Users/samilane/Documents/VSCode/Research/rna-seq-pipeline/data/reference/salmon_index \
    -l A \
    -1 ${SRR}_1.trimmed.fastq \
    -2 ${SRR}_2.trimmed.fastq \
    -p 8 \
    -o /Users/samilane/Documents/VSCode/Research/rna-seq-pipeline/data/quant/${SRR}
done

Each sample now has a directory containing:
quant.sf
cmd_info.json
lib_format_counts.json
aux_info/

These quant.sf files will be used for DESeq2.

Next Steps:
1. Install MultiQC in working env for trimmed QC summary
2. Create metadata table (samples.tsv)
3. Import Salmon quantifications into tximport
4. Build DESeq2 dataset
5. Perform:
- differential expression
- PCA
- sample distance heatmap
- volcano plots
- ranked gene lists
- per-condition comparisons

### Notes
All steps are run on macOS using conda environments and Homebrew-installed tools.
Project is designed to be modular, so components can be swapped out (e.g., STAR alignment instead of Salmon).
This pipeline will eventually include full reproducibility and possibly a Snakefile or Nextflow script.

## How to Run


## Requirements


## Results Preview
