# DESeq2 Analysis

This document describes how to run the DESeq2 analysis included in this repository.

Files included

- `metadata/samples.tsv` — sample table (sample, srr, condition, replicate).
- `analysis/RNAseq_DESeq2.Rmd` — interactive RMarkdown walkthrough.
- `analysis/RNAseq_DESeq2.R` — non-interactive script (Rscript) that saves results to `results/deseq2/` and plots to `results/plots/`.

Quick start (conda)

```bash
conda env create -f environment.yml
conda activate rna-seq-pipeline
```

Install R packages (inside R)

```r
install.packages(c("tidyverse", "pheatmap"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("tximport", "tximportData", "DESeq2", "apeglm", "AnnotationDbi", "org.Hs.eg.db", "clusterProfiler"))
```

Run non-interactively

```bash
Rscript analysis/RNAseq_DESeq2.R
```

Notes

- Ensure Salmon `quant.sf` files exist at `data/quant/<SRR>/quant.sf` before running.
- Edit `metadata/samples.tsv` to reflect the correct sample names and conditions.
- For gene-level counts, provide a `tx2gene` table and import with `tximport(..., txOut = FALSE)`.

Outputs

- `results/deseq2/` — CSVs of DE results.
- `results/plots/` — PCA, sample distances, volcano.

See `analysis/RNAseq_DESeq2.Rmd` for interactive plotting, annotation, and further steps (GO enrichment, clusterProfiler examples).
