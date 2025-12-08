# Comprehensive RNA-seq Pipeline with fastp, Salmon, and DESeq2 (SRP075484)
Samantha Lane
M.S. Bioinformatics student, Johns Hopkins University

## Dataset Overview
GSE81698 – “Epigenetic targeting of immune checkpoint PD-L1 by BET bromodomain inhibition”
Identifiers	SRA: SRP075484
BioProject: PRJNA322328
GEO: GSE81698
Study Type	Transcriptome Analysis
This repository contains a modular and reproducible RNA-seq differential expression pipeline built using public data from SRP075484 (BioProject PRJNA322328). The project includes four treatment-condition RNA-seq samples (no explicit control). The pipeline uses modern, lightweight tools (fastp, Salmon, MultiQC) and follows current best practices for quantification-based RNA-seq workflows.

#### [Reference Page](https://github.com/splane00/rna-seq-pipeline/sources)

## Repository Structure

```text
rna-seq-pipeline/
├── data/
│   ├── raw/
│   ├── trimmed/
│   └── reference/
├── scripts/
│   ├── qc_fastp.sh
│   ├── salmon_index.sh
│   └── salmon_quant.sh
├── results/
│   ├── qc/
│   └── quant/
├── config/
│   └── samples.tsv
├── .gitignore
├── README.md
└── environment.yml
```

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
- FastQC run on all raw FASTQs:  
  fastqc *.fastq --outdir ../qc/raw_fastqc  
- MultiQC summary (pending MultiQC install):  
  multiqc ../qc/raw_fastqc -o ../qc/raw_multiqc  

### Trimming (fastp)  
Each sample was adapter- and quality-trimmed using fastp:  
``` text
fastp \  
  -i SRR15074527_1.fastq \
  -I SRR15074527_2.fastq \
  -o ../trimmed/SRR15074527_1.trimmed.fastq \
  -O ../trimmed/SRR15074527_2.trimmed.fastq \
  --detect_adapter_for_pe \
  --html ../trimmed/SRR15074527.fastp.html \
  --json ../trimmed/SRR15074527.fastp.json
```
All trimmed reads stored in: data/trimmed/
MultiQC-ready fastp reports also generated.

### Reference Preparation (GENCODE v45)
Downloaded:
gencode.v45.transcripts.fa.gz (used for Salmon index)
gencode.v45.annotation.gtf (used later for tximport + DESeq2)

Stored in: data/reference/

### Salmon Indexing
Index built using:
``` text
salmon index \
  -t gencode.v45.transcripts.fa.gz \
  -i salmon_index \
  -k 31
```
Index stored in:data/reference/salmon_index/

### Salmon Quantification (Alignment-Free)
Quantified all 4 paired-end samples:
```text
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
```

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

---

## **DESeq2 Quickstart**

A full DESeq2 walkthrough is available in `docs/DESeq2.md`. The repository also includes a runnable script and an RMarkdown report to perform differential expression, PCA, heatmaps, and volcano plots from Salmon quantifications.

Quick start (create environment and run):

```bash
# create environment from the provided template (adjust channels/versions if needed)
conda env create -f environment.yml
conda activate rna-seq-pipeline

# inside R (if needed) install Bioconductor packages
R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager'); BiocManager::install(c('tximport','DESeq2','apeglm','org.Hs.eg.db'))"

# run non-interactive analysis (outputs -> results/deseq2/ and results/plots/)
Rscript analysis/RNAseq_DESeq2.R
```

Notes:

- Edit `metadata/samples.tsv` before running to match your samples and conditions.
- The script expects Salmon outputs at `data/quant/<SRR>/quant.sf`.
- For gene-level import add a `tx2gene` mapping and set `txOut = FALSE` in `tximport()`.
