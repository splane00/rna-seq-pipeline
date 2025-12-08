flowchart TD
    A[FASTQ Files] --> B[QC: FastQC]
    B --> C[Trim Reads]
    C --> D[Align with STAR]
    D --> E[Count Reads with featureCounts]
    E --> F[DESeq2 Differential Expression]
    F --> G[Plots: Volcano, PCA, Heatmap]
