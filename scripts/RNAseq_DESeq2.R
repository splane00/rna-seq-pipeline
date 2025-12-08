# RNAseq_DESeq2.R
# Minimal, runnable DESeq2 pipeline script using tximport + Salmon outputs
# Run with: Rscript scripts/RNAseq_DESeq2.R

# If packages are missing, uncomment and run the install lines below inside R.
# install.packages(c("tidyverse", "pheatmap"))
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install(c("tximport", "tximportData", "DESeq2", "apeglm", "AnnotationDbi", "org.Hs.eg.db"))

library(tximport)
library(readr)
library(DESeq2)
library(ggplot2)
library(pheatmap)

# Paths
samples_file <- "metadata/samples.tsv"
# Salmon quant results (directories named by sample IDs)
quant_dir <- "results/salmon"
results_dir <- "results/deseq2"
plots_dir <- file.path("results/plots")

# Ensure output directories exist
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# Read sample metadata
samples <- read_tsv(samples_file, show_col_types = FALSE)

# Build file paths to quant.sf files
files <- file.path(quant_dir, samples$sample, "quant.sf")
names(files) <- samples$sample

# Check existence
missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0) {
  stop("Missing quant.sf files:\n", paste(missing_files, collapse = "\n"))
}

# Import with tximport (txOut=TRUE assumes transcript-level; use tx2gene for gene-level)
txi <- tximport(files, type = "salmon", txOut = TRUE)

# Create DESeq2 dataset
coldata <- as.data.frame(samples)
rownames(coldata) <- coldata$sample
coldata$condition <- factor(coldata$condition)

# DESeq2 needs replicates per condition for contrasts
if (any(table(coldata$condition) < 2)) {
  stop("Each condition needs at least 2 replicates for DESeq2. Update metadata or design.")
}

dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)

# Pre-filtering: remove rows with zero or near-zero counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq
dds <- DESeq(dds)

# Transform for visualization
rld <- rlog(dds, blind = FALSE)

# PCA
pdf(file.path(plots_dir, "pca.pdf"), width = 6, height = 5)
print(plotPCA(rld, intgroup = "condition"))
dev.off()

# Sample distance heatmap
sampleDists <- dist(t(assay(rld)))
pdf(file.path(plots_dir, "sample_distance_heatmap.pdf"), width = 6, height = 6)
m <- as.matrix(sampleDists)
rownames(m) <- colnames(m) <- colnames(assay(rld))
pheatmap(m, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists)
dev.off()

# Example pairwise comparison: treatment1 vs treatment2 (modify as needed)
if ("treatment1" %in% coldata$condition && "treatment2" %in% coldata$condition) {
  res_t1_t2 <- results(dds, contrast = c("condition", "treatment1", "treatment2"))
  res_t1_t2 <- lfcShrink(dds, contrast = c("condition", "treatment1", "treatment2"), type = "apeglm")
  res_df <- as.data.frame(res_t1_t2)
  write.csv(res_df, file.path(results_dir, "treatment1_vs_treatment2.csv"))

  # Basic volcano plot
  res_df$gene <- rownames(res_df)
  res_df <- res_df[!is.na(res_df$padj), ]
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(alpha = 0.6) +
    theme_minimal() +
    xlab("log2 fold change") +
    ylab("-log10 adjusted p-value")
  ggsave(filename = file.path(plots_dir, "volcano_t1_vs_t2.png"), plot = p, width = 6, height = 5)
}

message("DESeq2 script finished. Results in: ", results_dir, " and plots in: ", plots_dir)
