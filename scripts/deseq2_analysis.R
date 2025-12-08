#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
})

# ------------------------------------------------------------------
# 1) Load metadata
# -------------------------------------------------------------------
samples <- read_tsv("config/samples.tsv", show_col_types = FALSE)
rownames(samples) <- samples$sample

# Build simple sample_table now (needed for plots & saving later)
sample_table <- samples %>% 
  select(sample, condition) %>% 
  mutate(condition = factor(condition))

# -------------------------------------------------------------------
# 2) Locate quant.sf files
# -------------------------------------------------------------------
files <- file.path("results/salmon", samples$sample, "quant.sf")
names(files) <- samples$sample

# -------------------------------------------------------------------
# 3) Load tx2gene mapping
# -------------------------------------------------------------------
tx2gene <- read_tsv("data/reference/tx2gene.tsv",
                    col_names = FALSE, show_col_types = FALSE)
colnames(tx2gene) <- c("transcript", "gene")

# -------------------------------------------------------------------
# 4) Import quantification with tximport
# -------------------------------------------------------------------
txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE,
  ignoreAfterBar = TRUE
)

# -------------------------------------------------------------------
# 5) Build DESeq2 object with NO design (~1)
# -------------------------------------------------------------------
dds <- DESeqDataSetFromTximport(
  txi,
  colData = samples,
  design = ~ 1
)

# -------------------------------------------------------------------
# 6) Normalization
# -------------------------------------------------------------------
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

# Output directory
outdir <- "results/deseq2"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

write_tsv(as.data.frame(norm_counts),
          file.path(outdir, "normalized_counts.tsv"))

# Save size factors
sf <- sizeFactors(dds)
write_tsv(
  tibble(sample = names(sf), size_factor = sf),
  file.path(outdir, "size_factors.tsv")
)

# -------------------------------------------------------------------
# 7) Variance Stabilizing Transform
# -------------------------------------------------------------------
vsd <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "mean")
vst_mat <- assay(vsd)

write_tsv(
  as.data.frame(vst_mat) %>% 
    rownames_to_column("gene"),
  file.path(outdir, "vst_counts.tsv")
)

# -------------------------------------------------------------------
# 8) PCA plot
# -------------------------------------------------------------------
pca_df <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_df, "percentVar"))

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.7, size = 3) +
  theme_bw() +
  xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ylab(paste0("PC2 (", percentVar[2], "%)"))

ggsave(file.path(outdir, "PCA_plot.png"), p_pca,
       width = 6, height = 5, dpi = 300)

# -------------------------------------------------------------------
# 9) Sample distance heatmap
# -------------------------------------------------------------------
dists <- dist(t(vst_mat))
mat <- as.matrix(dists)

pheatmap(mat,
         clustering_distance_rows = dists,
         clustering_distance_cols = dists,
         main = "Sample-to-Sample Distance (VST)",
         filename = file.path(outdir, "sample_distance_heatmap.png"),
         width = 6, height = 5)

# -------------------------------------------------------------------
# 10) Top 500 most variable genes heatmap
# -------------------------------------------------------------------
top500 <- head(order(matrixStats::rowVars(vst_mat), decreasing = TRUE), 500)

png(file.path(outdir, "heatmap_top500.png"), width = 900, height = 1200)
pheatmap(vst_mat[top500, ], show_rownames = FALSE)
dev.off()

# -------------------------------------------------------------------
# 11) Save sample table
# -------------------------------------------------------------------
write_tsv(sample_table, file.path(outdir, "sample_table_used.tsv"))

# Done
message("DESeq2 normalization, VST, PCA, and QC outputs saved.")
