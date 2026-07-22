# ======================================================================
# diagnose_batch_effects.R
# Plots a global PCA combining all 3 datasets to visualize batch variance
# ======================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
})

# Use portable working directory logic
if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd("d:/Zika_wetlab")
}

cat("\n═══ Batch Effect Diagnostics (PCA of all 18 samples) ═══\n")

# ---- 1. Load Count Matrices ----
cat("Loading count matrices from datasets...\n")

load_counts <- function(file, prefix, cols, is_gse265922 = FALSE) {
  if (!file.exists(file)) stop(paste("File missing:", file))
  df <- read.csv(file, row.names = 1, check.names = FALSE)
  
  if (is_gse265922) {
    # GSE265922 has metadata columns first
    df <- df[!is.na(df$ENTREZID) & df$ENTREZID != "", ]
    df <- df[!duplicated(df$ENTREZID), ]
    rownames(df) <- df$ENTREZID
    # The count columns have sample names, we know it's 6 columns
    df <- df[, grep("ZIKV|Control", colnames(df))]
  }
  
  colnames(df) <- paste0(prefix, "_", colnames(df))
  return(df)
}

c1 <- load_counts("GSE146423/results/tables/11_Raw_Counts.csv", "G146423", 1:6)
c2 <- load_counts("GSE233049/results/tables/11_Raw_Counts.csv", "G233049", 1:6)
c3 <- load_counts("GSE265922/results/tables/11_Raw_Counts_Merged.csv", "G265922", 1:6, is_gse265922 = TRUE)

# Intersect genes
common_genes <- intersect(intersect(rownames(c1), rownames(c2)), rownames(c3))
cat(sprintf("Found %d common genes across all 3 datasets.\n", length(common_genes)))

mat <- cbind(
  c1[common_genes, ],
  c2[common_genes, ],
  c3[common_genes, ]
)

# ---- 2. Build metadata ----
conditions <- c(
  rep(c("Control", "ZIKV"), each = 3), # GSE146423
  rep(c("Control", "ZIKV"), each = 3), # GSE233049
  rep("ZIKV", 3), rep("Control", 3)    # GSE265922
)

coldata <- data.frame(
  Dataset = factor(rep(c("GSE146423", "GSE233049", "GSE265922"), each = 6)),
  Condition = factor(conditions, levels = c("Control", "ZIKV")),
  row.names = colnames(mat)
)

# ---- 3. Run DESeq2 and PCA ----
cat("Running joint DESeq2 and VST transformation...\n")
dds <- DESeqDataSetFromMatrix(mat, coldata, ~ Dataset + Condition)
# Pre-filter
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
vsd <- vst(dds, blind = FALSE)

pcaData <- plotPCA(vsd, intgroup = c("Condition", "Dataset"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p <- ggplot(pcaData, aes(PC1, PC2, color = Dataset, shape = Condition)) +
  geom_point(size = 5, alpha = 0.8) +
  scale_shape_manual(values = c(16, 17)) +
  theme_minimal() +
  labs(
    title = "Global PCA: Cross-Dataset Batch Effect Diagnostic",
    subtitle = "Visualizes variance explained by Dataset vs Biological Condition",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance")
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(face = "bold")
  )

dir.create("meta-analysis/plots", recursive = TRUE, showWarnings = FALSE)
ggsave("meta-analysis/plots/Global_PCA_Batch_Diagnostics.png", p, width = 8, height = 6, dpi = 600, bg = "white")

cat("✓ PCA plot saved to meta-analysis/plots/Global_PCA_Batch_Diagnostics.png\n")
cat("Diagnostic Complete.\n")
