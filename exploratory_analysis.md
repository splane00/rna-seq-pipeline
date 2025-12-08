# Exploratory Data Analysis

To assess data quality, global transcriptomic structure, and the biological consistency of the four experimental conditions (DMSO, JQ1, BRD4.WT, shBRD4), we performed a series of exploratory analyses on Salmon-quantified RNA-seq data, including variance-stabilizing transformation (VST), principal component analysis (PCA), hierarchical clustering, and sample-to-sample distance visualization.

## Figures
![Top 500 Heatmap](figures/heatmap_top500.png)

![PCA Plot](figures/PCA_plot.png)

![Sample Distance Heatmap](figures/sample_distance_heatmap.png)

## Analysis
### 1. Data Quality and Transformation

Salmon quantifications were imported and summarized to the gene level, followed by DESeq2 variance stabilizing transformation (VST).
VST produces homoscedastic data suitable for visualization and removes the strong mean-variance dependence typical of raw RNA-seq counts.

No samples exhibited aberrant library sizes, extreme outlier gene distributions, or technical artifacts that would indicate low-quality sequencing or preprocessing issues.

### 2. Heatmap of the Top 500 Most Variable Genes

Hierarchical clustering of the top 500 high-variance genes revealed clear separation among all four samples, with two dominant expression signatures:

A transcriptionally active signature (BRD4.WT and, to a lesser extent, DMSO)

A transcriptionally suppressed signature (JQ1 and shBRD4)

The clustering tree showed JQ1 and shBRD4 grouping together, consistent with their shared role in reducing BRD4 function—JQ1 through BET inhibition and shBRD4 through genetic knockdown.
Conversely, BRD4.WT formed a distinct branch reflecting widespread activation of BRD4-dependent transcriptional programs.

The heatmap confirms strong biological signal, minimal noise, and no evidence of swapped or mislabeled samples.


### 3. Principal Component Analysis (PCA)

PCA revealed striking separation among the four conditions, with PC1 (46% of variance) capturing a BRD4 activity gradient:

BRD4.WT positioned at the extreme positive end (maximal BRD4-driven transcription)

shBRD4 at the extreme negative end (strongest transcriptional suppression)

JQ1 between DMSO and shBRD4, reflecting partial BET inhibition

DMSO near the center, representing the untreated baseline

PC2 (34% of variance) captured additional regulatory differences between chemical (JQ1) and genetic (shBRD4) manipulation of BRD4.

Overall, the PCA results demonstrate that the dataset is biologically coherent, with clear, interpretable separation driven by BRD4 perturbation.

### 4. Sample-to-Sample Distance Matrix

The Euclidean distance heatmap of VST-transformed counts further supported these findings:

Closest pair: JQ1 and shBRD4

Second closest: DMSO and BRD4.WT

Most distant: BRD4.WT and shBRD4

These distances reflect the expected transcriptional continuum:

BRD4.WT → DMSO → JQ1 → shBRD4

This ordering aligns with the known roles of BRD4:

Overexpression enhances enhancer activation and transcriptional output

Knockdown suppresses transcription broadly

JQ1 induces intermediate suppression by preventing BRD4’s chromatin binding

The coherence between PCA and distance clustering indicates high-quality data and confirms that each sample reflects its intended biological state.

### 5. Summary of Exploratory Insights

Across all exploratory analyses, the dataset displays:

Strong biological signal despite the absence of replicates

Clear, interpretable separation of sample conditions

Expected relationships based on BRD4 biology

No major QC issues, batch effects, or sample outliers

These results validate the dataset’s suitability for downstream analyses such as gene-ranking, exploratory differential expression, pathway enrichment, and BRD4-dependent gene program characterization.