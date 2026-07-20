# ======================================================================
# GSE265922: ZIKV-infected A549 — DEG & Enrichment
# GRCh38.p13 NCBI | DESeq2 + apeglm | 600 DPI PNG | 6-sample design
# Input: 6 individual STAR ReadsPerGene.out.tab.gz files
#   Z1/Z2/Z3 = ZIKV infected (n=3)
#   M1/M2/M6 = Mock/Control (n=3)
# ======================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(EnhancedVolcano)
  library(apeglm)
})

setwd("d:/Zika_wetlab/GSE265922")
for (d in c("results/tables", "plots/QC", "plots/DEG", "plots/Enrichment", "plots/GSEA"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- Palettes ----
pal <- list(
  ctrl  = "#2166AC",
  zikv  = "#B2182B",
  up    = "#D73027",
  down  = "#4575B4",
  ns    = "grey80",
  hm    = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  blues = colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
)

# ---- ggplot2 theme ----
theme_pub <- theme_minimal() +
  theme(
    plot.title       = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(size = 13, hjust = 0.5, color = "grey40"),
    axis.title       = element_text(size = 14, face = "bold"),
    axis.text        = element_text(size = 12, color = "black"),
    axis.line        = element_line(color = "black", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    legend.title     = element_text(size = 12, face = "bold"),
    legend.text      = element_text(size = 11),
    legend.position  = "bottom"
  )

# ---- Saving helpers ----
save_png <- function(p, path, file, w, h) {
  ggsave(file.path(path, paste0(file, ".png")), p,
         width = w, height = h, dpi = 600, bg = "white")
}

save_hm <- function(obj, path, file, w, h) {
  png(file.path(path, paste0(file, ".png")), w, h, "in", res = 600)
  grid::grid.newpage()
  grid::grid.draw(obj$gtable)
  dev.off()
}

enrich_plot <- function(obj, title, top_n, path, file, w, h) {
  if (is.null(obj) || nrow(as.data.frame(obj)) == 0) return(invisible())
  p <- dotplot(obj, showCategory = min(top_n, nrow(obj)), title = title) +
    theme_pub +
    theme(axis.text.y = element_text(size = 11), legend.position = "right")
  save_png(p, path, file, w, h)
}

# Safe nrow helper
n_safe <- function(x) if (is.null(x)) 0L else nrow(as.data.frame(x))

# ======================================================================
# 1. LOAD INDIVIDUAL STAR ReadsPerGene.out.tab.gz FILES
# ======================================================================

# ---- Define files and conditions ----
# Z samples = ZIKV infected, M samples = Mock (Control)
star_files <- c(
  "GSM8231986_trimmed_Z1_ReadsPerGene.out.tab.gz",
  "GSM8231987_trimmed_Z2_ReadsPerGene.out.tab.gz",
  "GSM8231988_trimmed_Z3_ReadsPerGene.out.tab.gz",
  "GSM8231989_trimmed_M1_ReadsPerGene.out.tab.gz",
  "GSM8231990_trimmed_M2_ReadsPerGene.out.tab.gz",
  "GSM8231991_trimmed_M6_ReadsPerGene.out.tab.gz"
)

# Sample names — first 3 = ZIKV, last 3 = Control
sample_names <- c("ZIKV_1", "ZIKV_2", "ZIKV_3", "Control_1", "Control_2", "Control_3")

n_zikv <- 3
n_ctrl <- 3

# Check files exist
missing_files <- star_files[!file.exists(star_files)]
if (length(missing_files) > 0) {
  stop("Missing files:\n", paste("  -", missing_files, collapse = "\n"))
}

cat("Reading STAR ReadsPerGene.out.tab.gz files...\n")

# ---- Read each file and extract unstranded counts (column 2) ----
# STAR output format (4 columns, tab-separated):
#   Column 1 = Gene ID
#   Column 2 = Counts unstranded
#   Column 3 = Counts 1st strand
#   Column 4 = Counts 2nd strand
# First 4 lines are summary stats (N_unmapped, N_multimapping, N_noFeature, N_ambiguous)

count_list <- list()
gene_ids   <- NULL

for (i in seq_along(star_files)) {
  cat(sprintf("  [%d/%d] %s", i, length(star_files), star_files[i]))
  
  # Read gzipped STAR output, skip first 4 summary lines
  tmp <- read.table(
    gzfile(star_files[i]),
    header       = FALSE,
    row.names    = NULL,
    sep          = "\t",
    skip         = 4,
    stringsAsFactors = FALSE,
    check.names  = FALSE,
    comment.char = ""
  )
  
  # Verify gene ID consistency across files
  if (i == 1) {
    gene_ids <- tmp[[1]]
    # Show a few gene IDs so we know what format they are
    cat(sprintf("\n  First 5 gene IDs: %s\n",
                paste(head(gene_ids, 5), collapse = ", ")))
  } else if (!identical(gene_ids, tmp[[1]])) {
    warning(sprintf("Gene IDs in %s differ from first file!", star_files[i]))
  }
  
  # Extract unstranded counts (column 2) as integer
  counts_vec <- as.integer(tmp[[2]])
  count_list[[sample_names[i]]] <- counts_vec
  
  cat(sprintf(" → %d genes, %.1fM total reads\n",
              length(counts_vec),
              sum(counts_vec, na.rm = TRUE) / 1e6))
}

# ---- Build count matrix ----
cnt <- as.data.frame(count_list)
rownames(cnt) <- gene_ids

# Remove STAR summary rows (N_unmapped, N_multimapping, etc.)
cnt <- cnt[!grepl("^N_", rownames(cnt)), , drop = FALSE]

cnt <- as.matrix(cnt)
storage.mode(cnt) <- "integer"

cat(sprintf("\nCount matrix assembled: %d genes x %d samples\n", nrow(cnt), ncol(cnt)))
cat("Columns:", paste(colnames(cnt), collapse = ", "), "\n")

# Remove genes with zero counts across all samples
cnt <- cnt[rowSums(cnt) > 0, ]

cat(sprintf("After zero-count filtering: %d genes\n", nrow(cnt)))
cat(sprintf("Library sizes (M): %s\n",
            paste(round(colSums(cnt) / 1e6, 1), collapse = ", ")))

# ======================================================================
# 2. GENE ID MAPPING — Detect ID type & convert to SYMBOL / ENSEMBL
# ======================================================================
cat("\n═══ Gene ID Mapping ═══\n")

# Detect ID type by looking at actual gene IDs
sample_ids <- head(rownames(cnt), 20)
cat("Sample gene IDs:", paste(sample_ids[1:5], collapse = ", "), "\n")

# Try to determine if IDs are Ensembl, Entrez, or something else
is_ensembl <- any(grepl("^ENSG", sample_ids))
is_entrez  <- all(grepl("^[0-9]+$", sample_ids))

if (is_ensembl) {
  cat("→ Detected ENSEMBL gene IDs\n")
  id_type <- "ENSEMBL"
} else if (is_entrez) {
  cat("→ Detected ENTREZ gene IDs\n")
  id_type <- "ENTREZID"
} else {
  # Try to guess — strip version numbers from Ensembl IDs
  cat("→ Gene IDs don't match standard patterns. Trying with version-stripped Ensembl...\n")
  # Many STAR outputs use Ensembl with version: ENSG00000123456.11
  # Let's strip versions and try both
  if (any(grepl("^ENSG", sample_ids, ignore.case = TRUE))) {
    id_type <- "ENSEMBL"
    cat("→ Detected ENSEMBL-like IDs (will strip versions)\n")
  } else {
    id_type <- "SYMBOL"  # fallback
    cat("→ Could not detect ID type. Trying SYMBOL...\n")
  }
}

# ---- Perform mapping based on detected type ----
gene_map <- NULL

if (id_type == "ENSEMBL") {
  # Strip version numbers: ENSG00000123456.11 → ENSG00000123456
  clean_ids <- gsub("\\.[0-9]+$", "", rownames(cnt))
  names(clean_ids) <- rownames(cnt)  # preserve original names
  
  gene_map <- tryCatch(
    bitr(
      unique(clean_ids),
      fromType = "ENSEMBL",
      toType   = c("ENTREZID", "SYMBOL"),
      OrgDb    = org.Hs.eg.db
    ),
    error = function(e) {
      cat("ERROR in ENSEMBL mapping:", e$message, "\n")
      return(NULL)
    }
  )
  
  if (!is.null(gene_map)) {
    # Map back to original IDs (with versions)
    gene_map$ENSEMBL_ORIG <- gene_map$ENSEMBL
    # Create lookup: clean_id → original_id
    id_lookup <- setNames(rownames(cnt), clean_ids)
    gene_map$ENSEMBL_ORIG <- id_lookup[gene_map$ENSEMBL]
  }
  
} else if (id_type == "ENTREZID") {
  gene_map <- tryCatch(
    bitr(
      unique(rownames(cnt)),
      fromType = "ENTREZID",
      toType   = c("SYMBOL", "ENSEMBL"),
      OrgDb    = org.Hs.eg.db
    ),
    error = function(e) {
      cat("ERROR in ENTREZ mapping:", e$message, "\n")
      return(NULL)
    }
  )
  
} else {
  # Try SYMBOL
  gene_map <- tryCatch(
    bitr(
      unique(rownames(cnt)),
      fromType = "SYMBOL",
      toType   = c("ENTREZID", "ENSEMBL"),
      OrgDb    = org.Hs.eg.db
    ),
    error = function(e) {
      cat("ERROR in SYMBOL mapping:", e$message, "\n")
      return(NULL)
    }
  )
}

# ---- Handle mapping results ----
if (is.null(gene_map) || nrow(gene_map) == 0) {
  cat("\n⚠️  WARNING: Gene ID mapping failed! Proceeding with original IDs.\n")
  cat("   DEG table will use raw gene IDs. Enrichment will be skipped.\n")
  
  # Create a minimal gene_map with original IDs
  gene_map <- data.frame(
    ENTREZID = rownames(cnt),
    SYMBOL   = rownames(cnt),
    ENSEMBL  = rownames(cnt),
    stringsAsFactors = FALSE
  )
  
  # Set key column for DESeq2
  key_col <- "ENTREZID"
  use_mapped <- FALSE
} else {
  cat(sprintf("✓ Mapped %d genes successfully\n", nrow(gene_map)))
  
  # Determine which ID column matches our count matrix rownames
  if (id_type == "ENSEMBL") {
    # We stripped versions; filter count matrix using ENSEMBL_ORIG
    keep_genes <- intersect(rownames(cnt), gene_map$ENSEMBL_ORIG)
    key_col <- "ENSEMBL_ORIG"
  } else if (id_type == "ENTREZID") {
    keep_genes <- intersect(rownames(cnt), gene_map$ENTREZID)
    key_col <- "ENTREZID"
  } else {
    keep_genes <- intersect(rownames(cnt), gene_map$SYMBOL)
    key_col <- "SYMBOL"
  }
  
  # Filter count matrix to mapped genes
  cnt <- cnt[keep_genes, ]
  gene_map <- gene_map[match(keep_genes, gene_map[[key_col]]), ]
  cat(sprintf("After gene ID mapping: %d genes\n", nrow(cnt)))
  use_mapped <- TRUE
}

# ======================================================================
# 3. DESeq2
# ======================================================================
col_data <- data.frame(
  condition = factor(
    c(rep("ZIKV", n_zikv), rep("Control", n_ctrl)),
    levels = c("Control", "ZIKV")
  ),
  row.names = colnames(cnt)
)

cat("\nSample metadata:\n")
print(col_data)

dds <- DESeqDataSetFromMatrix(
  countData = cnt,
  colData   = col_data,
  design    = ~ condition
)

# Pre-filter: keep genes with ≥ 10 counts in at least 3 samples (smallest group size)
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
cat(sprintf("After prefilter (≥10 counts in ≥3 samples): %d genes\n", nrow(dds)))

# Run DESeq2
cat("\nRunning DESeq2...\n")
dds <- DESeq(dds)

# LFC shrinkage with apeglm
cat("Applying apeglm shrinkage...\n")
res <- lfcShrink(dds, coef = "condition_ZIKV_vs_Control", type = "apeglm")

# ======================================================================
# 4. BUILD DEG TABLE
# ======================================================================

# Decide which gene ID column to use for joining
# gene_map has: ENTREZID, SYMBOL, ENSEMBL (and possibly ENSEMBL_ORIG)
# The rownames of counts/dds match the original IDs

if (use_mapped) {
  # Create a join-able ID column in gene_map to match res rownames
  if (id_type == "ENSEMBL") {
    gene_map$JOIN_ID <- gene_map$ENSEMBL_ORIG
  } else if (id_type == "ENTREZID") {
    gene_map$JOIN_ID <- gene_map$ENTREZID
  } else {
    gene_map$JOIN_ID <- gene_map$SYMBOL
  }
} else {
  gene_map$JOIN_ID <- gene_map$ENTREZID
}

# Build the DEG table
deg <- as.data.frame(res) %>%
  rownames_to_column("GeneID") %>%
  left_join(
    gene_map %>% dplyr::select(JOIN_ID, SYMBOL, ENSEMBL, ENTREZID),
    by = c("GeneID" = "JOIN_ID")
  ) %>%
  mutate(
    sig = padj < 0.05 & abs(log2FoldChange) > 1,
    dir = case_when(
      sig & log2FoldChange > 1  ~ "Up",
      sig & log2FoldChange < -1 ~ "Down",
      TRUE                       ~ "NS"
    ),
    nlp = {
      x <- -log10(padj)
      ifelse(is.infinite(x), max(x[is.finite(x)]) + 10, x)
    }
  ) %>%
  arrange(padj)

# Add mean normalized counts per group
deg$Mean_Ctrl <- rowMeans(counts(dds, normalized = TRUE)[deg$GeneID, 1:n_ctrl])
deg$Mean_ZIKV <- rowMeans(counts(dds, normalized = TRUE)[deg$GeneID, (n_ctrl + 1):(n_ctrl + n_zikv)])

n_up   <- sum(deg$dir == "Up")
n_down <- sum(deg$dir == "Down")

cat(sprintf(
  "\n═══ DEG RESULTS ═══\n  Total: %d | Up: %d | Down: %d\n\n",
  n_up + n_down, n_up, n_down
))

# Use SYMBOL if available, otherwise GeneID
deg$Label <- ifelse(is.na(deg$SYMBOL) | deg$SYMBOL == "", deg$GeneID, deg$SYMBOL)

cat("Top 15 Upregulated:\n")
deg %>%
  filter(dir == "Up") %>%
  head(15) %>%
  dplyr::select(Label, log2FoldChange, padj, Mean_Ctrl, Mean_ZIKV) %>%
  print()

cat("\nTop 15 Downregulated:\n")
deg %>%
  filter(dir == "Down") %>%
  head(15) %>%
  dplyr::select(Label, log2FoldChange, padj, Mean_Ctrl, Mean_ZIKV) %>%
  print()

# ======================================================================
# 5. QC PLOTS
# ======================================================================
vsd <- vst(dds, blind = FALSE)

# 5a. PCA (with 95% confidence ellipses)
pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))

p_pca <- ggplot(pca, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 6, alpha = 0.9) +
  stat_ellipse(aes(fill = condition), alpha = 0.15, level = 0.95,
               geom = "polygon", show.legend = FALSE) +
  geom_text_repel(
    size        = 5,
    fontface    = "bold",
    show.legend = FALSE,
    box.padding = 0.8
  ) +
  scale_color_manual(values = c(Control = pal$ctrl, ZIKV = pal$zikv)) +
  scale_fill_manual(values = c(Control = pal$ctrl, ZIKV = pal$zikv)) +
  labs(
    title    = "PCA — GSE265922: ZIKV vs Control A549",
    subtitle = sprintf("%d expressed genes | n=3 per group", nrow(dds)),
    x        = paste0("PC1: ", pv[1], "%"),
    y        = paste0("PC2: ", pv[2], "%")
  ) +
  coord_fixed(1) +
  theme_pub

save_png(p_pca, "plots/QC", "01_PCA", 9, 8)

# 5b. Sample-to-sample distances
dists <- dist(t(assay(vsd)))
dm    <- as.matrix(dists)
dimnames(dm) <- list(
  paste0(vsd$condition, " — ", colnames(vsd)),
  paste0(vsd$condition, " — ", colnames(vsd))
)

hm_dist <- pheatmap(
  dm,
  clustering_distance_rows    = dists,
  clustering_distance_cols    = dists,
  annotation_col = data.frame(
    Condition = vsd$condition,
    row.names = colnames(vsd)
  ),
  annotation_colors = list(
    Condition = c(Control = pal$ctrl, ZIKV = pal$zikv)
  ),
  main             = "Sample-to-Sample Distances (VST) — GSE265922",
  color            = pal$blues,
  display_numbers  = TRUE,
  number_format    = "%.0f",
  fontsize         = 13,
  fontsize_number  = 11,
  silent           = TRUE
)

save_hm(hm_dist, "plots/QC", "02_Sample_Distances", 9, 8)

# 5c. Dispersion estimates
png(file.path("plots/QC", "03_Dispersion.png"), 8, 6, "in", res = 600)
plotDispEsts(dds, main = "DESeq2: Dispersion Estimates — GSE265922", cex = 0.6)
dev.off()

# 5d. Expression density
dens_df <- as.data.frame(counts(dds, normalized = TRUE)) %>%
  rownames_to_column("GeneID") %>%
  pivot_longer(
    cols      = -GeneID,
    names_to  = "Sample",
    values_to = "Count"
  ) %>%
  mutate(
    Group     = ifelse(grepl("Control", Sample), "Control", "ZIKV"),
    log2Count = log2(Count + 1)
  )

p_dens <- ggplot(dens_df, aes(log2Count, fill = Group, color = Group)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  scale_fill_manual(values = c(Control = pal$ctrl, ZIKV = pal$zikv)) +
  scale_color_manual(values = c(Control = pal$ctrl, ZIKV = pal$zikv)) +
  labs(
    title    = "Expression Density — GSE265922",
    x        = expression(Log[2] * "(Counts+1)"),
    subtitle = "log2(DESeq2-normalized counts + 1) | n=3 per group",
    y        = "Density"
  ) +
  theme_pub

save_png(p_dens, "plots/QC", "04_Expression_Density", 8, 6)

# 5e. Library sizes
lib <- data.frame(
  Sample  = colnames(cnt),
  Reads_M = colSums(cnt) / 1e6,
  Group   = factor(c(rep("ZIKV", n_zikv), rep("Control", n_ctrl)))
)

p_lib <- ggplot(lib, aes(Sample, Reads_M, fill = Group)) +
  geom_bar(stat = "identity", alpha = 0.85, width = 0.7) +
  geom_text(
    aes(label = sprintf("%.1fM", Reads_M)),
    vjust = -0.5,
    size  = 4
  ) +
  scale_fill_manual(values = c(Control = pal$ctrl, ZIKV = pal$zikv)) +
  labs(
    title    = "Library Sizes — GSE265922 (3 Control vs 3 ZIKV)",
    x        = "",
    y        = "Million Reads"
  ) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

save_png(p_lib, "plots/QC", "05_Library_Sizes", 8, 6)
cat("QC: ✓ 01-05\n")

# ======================================================================
# 6. DEG PLOTS
# ======================================================================
top30 <- deg %>%
  filter(sig) %>%
  slice_min(padj, n = 30)

p_volcano <- ggplot(deg, aes(log2FoldChange, nlp, color = dir, alpha = dir)) +
  geom_point(size = 1.8) +
  scale_color_manual(
    values = c(Up = pal$up, Down = pal$down, NS = pal$ns),
    name   = NULL
  ) +
  scale_alpha_manual(
    values = c(Up = 0.8, Down = 0.8, NS = 0.25),
    guide  = "none"
  ) +
  geom_hline(yintercept = -log10(0.05), lty = "longdash", color = "grey50") +
  geom_vline(xintercept = c(-1, 1), lty = "longdash", color = "grey50") +
  geom_text_repel(
    data        = top30,
    aes(label   = Label),
    size        = 3.2,
    fontface    = "italic",
    max.overlaps = 30,
    box.padding  = 0.6,
    segment.size = 0.3
  ) +
  labs(
    title    = "ZIKV vs Control — A549 Cells (GSE265922)",
    subtitle = sprintf(
      "%d Up · %d Down · %d DEGs | n=3 per group",
      n_up, n_down, n_up + n_down
    ),
    x = expression(Log[2] ~ "Fold Change"),
    y = expression(-Log[10] ~ "P-adj")
  ) +
  theme_pub

save_png(p_volcano, "plots/DEG", "06_Volcano", 10, 8.5)

# EnhancedVolcano
top_labels <- top30$Label[1:min(25, nrow(top30))]
png(file.path("plots/DEG", "07_EnhancedVolcano.png"), 11, 9, "in", res = 600)
EnhancedVolcano(
  deg,
  lab          = deg$Label,
  x            = "log2FoldChange",
  y            = "padj",
  title        = "ZIKV vs Control (A549) — DESeq2 | GSE265922",
  subtitle     = sprintf(
    "%d up · %d down · %d total | n=3 per group",
    n_up, n_down, n_up + n_down
  ),
  caption      = sprintf("padj < 0.05, |log2FC| > 1 | %d genes", nrow(deg)),
  pCutoff      = 0.05,
  FCcutoff     = 1,
  pointSize    = 2,
  labSize      = 3.5,
  colAlpha     = 0.6,
  legendPosition    = "bottom",
  drawConnectors    = TRUE,
  max.overlaps      = 25,
  selectLab         = top_labels,
  border            = "full",
  borderWidth       = 0.8
)
dev.off()

# MA plot
png(file.path("plots/DEG", "08_MA_Plot.png"), 9, 7, "in", res = 600)
DESeq2::plotMA(res, main = "MA Plot — GSE265922", ylim = c(-6, 6), alpha = 0.05)
dev.off()

# P-value histogram
p_pval <- ggplot(deg, aes(pvalue)) +
  geom_histogram(bins = 50, fill = "grey60", color = "grey30", alpha = 0.7, na.rm = TRUE) +
  geom_vline(xintercept = 0.05, lty = "dashed", color = pal$up, linewidth = 1) +
  labs(
    title    = "P-value Distribution — GSE265922",
    subtitle = sprintf("%d genes tested | n=3 per group", nrow(deg)),
    x        = "Raw P-value",
    y        = "Frequency"
  ) +
  theme_pub

save_png(p_pval, "plots/DEG", "09_Pvalue_Histogram", 8, 6)
cat("DEG: ✓ 06-09\n")

# ======================================================================
# 7. HEATMAP — Top 50 DEGs
# ======================================================================
n_deg  <- sum(deg$sig)
if (n_deg > 0) {
  n_hm   <- min(50, n_deg)
  top50  <- deg %>%
    filter(sig) %>%
    slice_min(padj, n = n_hm)
  
  hm_mat <- assay(vsd)[top50$GeneID, ]
  lbls   <- top50$Label
  lbls[is.na(lbls) | lbls == ""] <- top50$GeneID[is.na(lbls) | lbls == ""]
  dup <- duplicated(lbls)
  if (any(dup)) {
    for (i in which(dup)) {
      lbls[i] <- paste0(lbls[i], "_", i)
    }
  }
  rownames(hm_mat) <- lbls
  
  hm_top <- pheatmap(
    hm_mat,
    scale                  = "row",
    annotation_col = data.frame(
      Condition = vsd$condition,
      row.names = colnames(vsd)
    ),
    annotation_colors = list(
      Condition = c(Control = pal$ctrl, ZIKV = pal$zikv)
    ),
    main                   = sprintf(
      "Top %d DEGs — ZIKV vs Control A549 (GSE265922)", nrow(hm_mat)
    ),
    show_rownames          = TRUE,
    fontsize_row           = 6,
    fontsize_col           = 12,
    clustering_method      = "ward.D2",
    clustering_distance_rows   = "correlation",
    clustering_distance_cols   = "correlation",
    color                  = pal$hm,
    border_color           = NA,
    silent                 = TRUE
  )
  
  save_hm(hm_top, "plots/DEG", "10_Heatmap_Top50", 10, 13)
  cat("✓ 10_Heatmap.png\n")
} else {
  cat("⚠️  No significant DEGs found — skipping heatmap\n")
}

# ======================================================================
# 8. ENRICHMENT (GO + KEGG)
# ======================================================================
cat("\n═══ Enrichment ═══\n")

# Check if we have ENTREZ IDs for enrichment
if (use_mapped && "ENTREZID" %in% colnames(gene_map) && any(!is.na(gene_map$ENTREZID))) {
  
  bg  <- unique(deg$GeneID)
  up  <- unique(deg$GeneID[deg$dir == "Up"])
  dn  <- unique(deg$GeneID[deg$dir == "Down"])
  all <- unique(deg$GeneID[deg$sig])
  
  # Map GeneIDs to Entrez IDs
  id_map <- setNames(gene_map$ENTREZID, gene_map$JOIN_ID)
  bg_entrez  <- na.omit(unique(id_map[bg]))
  up_entrez  <- na.omit(unique(id_map[up]))
  dn_entrez  <- na.omit(unique(id_map[dn]))
  all_entrez <- na.omit(unique(id_map[all]))
  
  cat(sprintf("Background: %d | Up: %d | Down: %d | All sig: %d (Entrez IDs)\n",
              length(bg_entrez), length(up_entrez), length(dn_entrez), length(all_entrez)))
  
  # Helper functions
  run_go <- function(g, o) {
    if (length(g) < 5) return(NULL)
    enrichGO(
      gene          = g,
      universe      = bg_entrez,
      OrgDb         = org.Hs.eg.db,
      ont           = o,
      keyType       = "ENTREZID",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      readable      = TRUE
    )
  }
  
  run_kegg <- function(g) {
    if (length(g) < 5) return(NULL)
    enrichKEGG(
      gene          = g,
      universe      = bg_entrez,
      organism      = "hsa",
      keyType       = "kegg",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2
    )
  }
  
  # Run enrichments
  go_bp_up  <- tryCatch(run_go(up_entrez, "BP"), error = function(e) { cat("GO BP Up error:", e$message, "\n"); NULL })
  go_bp_dn  <- tryCatch(run_go(dn_entrez, "BP"), error = function(e) { cat("GO BP Down error:", e$message, "\n"); NULL })
  go_mf_up  <- tryCatch(run_go(up_entrez, "MF"), error = function(e) { cat("GO MF Up error:", e$message, "\n"); NULL })
  go_mf_dn  <- tryCatch(run_go(dn_entrez, "MF"), error = function(e) { cat("GO MF Down error:", e$message, "\n"); NULL })
  go_cc_up  <- tryCatch(run_go(up_entrez, "CC"), error = function(e) { cat("GO CC Up error:", e$message, "\n"); NULL })
  go_cc_dn  <- tryCatch(run_go(dn_entrez, "CC"), error = function(e) { cat("GO CC Down error:", e$message, "\n"); NULL })
  kegg_up   <- tryCatch(run_kegg(up_entrez),     error = function(e) { cat("KEGG Up error:", e$message, "\n"); NULL })
  kegg_dn   <- tryCatch(run_kegg(dn_entrez),     error = function(e) { cat("KEGG Down error:", e$message, "\n"); NULL })
  kegg_all  <- tryCatch(run_kegg(all_entrez),    error = function(e) { cat("KEGG All error:", e$message, "\n"); NULL })
  
  cat(sprintf(
    "GO BP: ↑%d ↓%d | GO MF: ↑%d ↓%d | GO CC: ↑%d ↓%d | KEGG: ↑%d ↓%d All:%d\n",
    n_safe(go_bp_up), n_safe(go_bp_dn),
    n_safe(go_mf_up), n_safe(go_mf_dn),
    n_safe(go_cc_up), n_safe(go_cc_dn),
    n_safe(kegg_up),  n_safe(kegg_dn),
    n_safe(kegg_all)
  ))
  
  # Dot plots
  enrich_plot(go_bp_up, "GO BP — Upregulated (GSE265922)", 20,
              "plots/Enrichment", "11_GO_BP_Up", 13, 9)
  enrich_plot(go_mf_up, "GO MF — Upregulated (GSE265922)", 15,
              "plots/Enrichment", "12_GO_MF_Up", 12, 7.5)
  enrich_plot(kegg_up,  "KEGG — Upregulated (GSE265922)", 15,
              "plots/Enrichment", "13_KEGG_Up", 12, 7.5)
  enrich_plot(kegg_all, "KEGG — All DEGs (GSE265922)", 15,
              "plots/Enrichment", "14_KEGG_All", 12, 7.5)
  
  # GO BP Up vs Down comparison
  if (!is.null(go_bp_up) && !is.null(go_bp_dn) &&
      nrow(as.data.frame(go_bp_up)) > 0 &&
      nrow(as.data.frame(go_bp_dn)) > 0) {
    comp <- bind_rows(
      as.data.frame(go_bp_up) %>% mutate(Reg = "Upregulated"),
      as.data.frame(go_bp_dn) %>% mutate(Reg = "Downregulated")
    ) %>%
      group_by(Reg) %>%
      slice_min(p.adjust, n = 15) %>%
      ungroup()
    
    p_comp <- ggplot(
      comp,
      aes(
        x    = reorder(str_wrap(Description, 55), -log10(p.adjust)),
        y    = -log10(p.adjust),
        fill = Reg
      )
    ) +
      geom_bar(
        stat     = "identity",
        position = position_dodge(0.8),
        width    = 0.7,
        alpha    = 0.9
      ) +
      coord_flip() +
      scale_fill_manual(
        values = c(Upregulated = pal$up, Downregulated = pal$down)
      ) +
      labs(
        title    = "GO BP: Up vs Down (GSE265922)",
        subtitle = "n=3 per group",
        x        = "",
        y        = expression(-Log[10] * "(P-adj)")
      ) +
      theme_pub
    
    save_png(p_comp, "plots/Enrichment", "15_GO_BP_Up_vs_Down", 15, 8)
  }
  
  # GO BP Enrichment Network
  if (!is.null(go_bp_up) && nrow(as.data.frame(go_bp_up)) > 5) {
    sim_go <- enrichplot::pairwise_termsim(go_bp_up)
    p_net  <- enrichplot::emapplot(
      sim_go,
      showCategory = min(30, nrow(as.data.frame(go_bp_up))),
      layout       = "nicely"
    ) +
      labs(title    = "GO BP Enrichment Map — GSE265922",
           subtitle = "n=3 per group") +
      theme_pub +
      theme(legend.position = "right")
    
    save_png(p_net, "plots/Enrichment", "16_GO_BP_Network", 14, 12)
  }
  cat("Enrichment: ✓ 11-16\n")
  
} else {
  cat("⚠️  Skipping enrichment — no valid Entrez IDs available\n")
  go_bp_up <- go_bp_dn <- go_mf_up <- go_mf_dn <- go_cc_up <- go_cc_dn <- NULL
  kegg_up <- kegg_dn <- kegg_all <- NULL
}

# ======================================================================
# 9. GSEA
# ======================================================================
cat("═══ GSEA ═══\n")

if (use_mapped && exists("bg_entrez") && length(bg_entrez) > 100) {
  
  # Build ranked gene list using Entrez IDs
  rnk <- deg %>%
    filter(!is.na(ENTREZID)) %>%
    arrange(desc(log2FoldChange)) %>%
    distinct(ENTREZID, .keep_all = TRUE) %>%
    pull(log2FoldChange, ENTREZID) %>%
    sort(decreasing = TRUE)
  
  cat(sprintf("GSEA ranked list: %d genes\n", length(rnk)))
  
  gsea <- tryCatch(
    gseGO(
      geneList      = rnk,
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",
      keyType       = "ENTREZID",
      minGSSize     = 10,
      maxGSSize     = 500,
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH",
      verbose       = FALSE
    ),
    error = function(e) {
      cat("GSEA error:", e$message, "\n")
      return(NULL)
    }
  )
  
  cat(sprintf("GSEA GO BP: %d gene sets\n", n_safe(gsea)))
  
  if (!is.null(gsea) && nrow(as.data.frame(gsea)) > 0) {
    
    # Running-score plot
    png(file.path("plots/GSEA", "17_GSEA_RunningScore.png"), 16, 12, "in", res = 600)
    ids <- head(which(gsea@result$p.adjust < 0.05), 4)
    if (length(ids) > 0) {
      print(
        enrichplot::gseaplot2(
          gsea, ids,
          title        = "GSEA — GSE265922 (n=3 per group)",
          pvalue_table = FALSE,
          subplots     = 1:2
        )
      )
    }
    dev.off()
    
    if (nrow(as.data.frame(gsea)) >= 5) {
      # Ridge plot
      p_ridge <- enrichplot::ridgeplot(gsea, showCategory = 20, fill = "p.adjust") +
        labs(
          title    = "GSEA GO BP — GSE265922",
          subtitle = "n=3 per group",
          x        = "Running Enrichment Score"
        ) +
        theme_minimal() +
        theme(
          plot.title  = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.title  = element_text(size = 14, face = "bold"),
          axis.text   = element_text(size = 10, color = "black")
        )
      
      ggsave(
        file.path("plots/GSEA", "18_GSEA_RidgePlot.png"),
        p_ridge,
        width = 14, height = 10, dpi = 600, bg = "white"
      )
      
      # Dot plot
      p_dot <- enrichplot::dotplot(
        gsea,
        showCategory = 20,
        title        = "GSEA — GSE265922 (n=3 per group)"
      ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text  = element_text(size = 11, color = "black")
        )
      
      ggsave(
        file.path("plots/GSEA", "19_GSEA_Dotplot.png"),
        p_dot,
        width = 16, height = 10, dpi = 600, bg = "white"
      )
      
      # ── Enrichment Map (network) ────────────────────────────────────
      if (nrow(as.data.frame(gsea)) >= 10) {
        tryCatch({
          gsea_pair <- enrichplot::pairwise_termsim(gsea, showCategory = min(50, nrow(as.data.frame(gsea))))
          p_emap <- enrichplot::emapplot(gsea_pair, showCategory = min(30, nrow(as.data.frame(gsea))),
                             layout = "nicely", node_label = "category") +
            scale_color_gradient2(low = pal$gsea_down, mid = "white",
                                  high = pal$gsea_up, midpoint = 0, name = "NES") +
            labs(title = "GSEA GO BP — Enrichment Map") +
            theme_minimal() + theme(
              plot.title       = element_text(size = 16, face = "bold", hjust = 0.5),
              legend.position  = "right",
              legend.text      = element_text(size = 9),
              legend.title     = element_text(size = 10)
            )
          ggsave(file.path("plots/GSEA", "20_GSEA_EnrichmentMap.png"), p_emap, width = 16, height = 12, dpi = 600, bg = "white")
          cat("✓ 20_GSEA_EnrichmentMap.png\n")
        }, error = function(e) cat(sprintf("⊙ 20_GSEA_EnrichmentMap FAILED: %s\n", e$message)))
      }
      
      # ── Gene-Concept Network ────────────────────────────────────────
      if (nrow(as.data.frame(gsea)) >= 3) {
        tryCatch({
          p_cnet <- enrichplot::cnetplot(gsea, showCategory = min(5, nrow(as.data.frame(gsea))),
                             foldChange = rnk,
                             node_label = "all") +
            scale_color_gradient2(low = pal$gsea_down, mid = "grey90",
                                  high = pal$gsea_up, midpoint = 0, name = "log2FC") +
            labs(title = "GSEA GO BP — Gene-Concept Network",
                 subtitle = "Top 5 gene sets | Nodes coloured by fold-change") +
            theme_minimal() + theme(
              plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
              legend.position = "right"
            )
          ggsave(file.path("plots/GSEA", "21_GSEA_CnetPlot.png"), p_cnet, width = 14, height = 10, dpi = 600, bg = "white")
          cat("✓ 21_GSEA_CnetPlot.png\n")
        }, error = function(e) cat(sprintf("⊙ 21_GSEA_CnetPlot FAILED: %s\n", e$message)))
      }
      
      # ── Leading-edge Heatplot ───────────────────────────────────────
      if (nrow(as.data.frame(gsea)) >= 5) {
        tryCatch({
          p_heat <- enrichplot::heatplot(gsea, showCategory = min(10, nrow(as.data.frame(gsea))),
                             foldChange = rnk) +
            scale_fill_gradient2(low = pal$gsea_down, mid = "white",
                                 high = pal$gsea_up, midpoint = 0, name = "log2FC") +
            labs(title = "GSEA GO BP — Leading-Edge Heatmap",
                 subtitle = "Top 10 gene sets | log2FC of core enrichment genes") +
            theme_minimal() + theme(
              plot.title  = element_text(size = 16, face = "bold", hjust = 0.5),
              axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
              axis.text.y = element_text(size = 9),
              legend.position = "right"
            )
          ggsave(file.path("plots/GSEA", "22_GSEA_Heatplot.png"), p_heat, width = 16, height = 10, dpi = 600, bg = "white")
          cat("✓ 22_GSEA_Heatplot.png\n")
        }, error = function(e) cat(sprintf("⊙ 22_GSEA_Heatplot FAILED: %s\n", e$message)))
      }
      
    }
    cat("GSEA: ✓ 17-22\n")
  } else {
    cat("GSEA: No significant gene sets found\n")
  }
  
} else {
  cat("⚠️  Skipping GSEA — no valid Entrez IDs\n")
  gsea <- NULL
}

# ======================================================================
# 10. SAVE TABLES
# ======================================================================
nc <- as.data.frame(counts(dds, normalized = TRUE)) %>%
  rownames_to_column("GeneID") %>%
  left_join(
    gene_map %>% dplyr::select(JOIN_ID, SYMBOL, ENSEMBL, ENTREZID),
    by = c("GeneID" = "JOIN_ID")
  ) %>%
  dplyr::select(GeneID, SYMBOL, ENSEMBL, ENTREZID, everything()) %>%
  arrange(GeneID)

write.csv(deg,                      "results/tables/01_DEG_Complete.csv",       row.names = FALSE)
write.csv(filter(deg, sig),         "results/tables/02_DEG_Significant.csv",     row.names = FALSE)
write.csv(filter(deg, dir == "Up"), "results/tables/03_DEG_Upregulated.csv",    row.names = FALSE)
write.csv(filter(deg, dir == "Down"),"results/tables/04_DEG_Downregulated.csv", row.names = FALSE)
write.csv(nc,                       "results/tables/05_Normalized_Counts.csv",   row.names = FALSE)

if (!is.null(go_bp_up) && nrow(as.data.frame(go_bp_up)) > 0)
  write.csv(as.data.frame(go_bp_up), "results/tables/06_GO_BP_Up.csv", row.names = FALSE)
if (!is.null(go_mf_up) && nrow(as.data.frame(go_mf_up)) > 0)
  write.csv(as.data.frame(go_mf_up), "results/tables/07_GO_MF_Up.csv", row.names = FALSE)
if (!is.null(kegg_up) && nrow(as.data.frame(kegg_up)) > 0)
  write.csv(as.data.frame(kegg_up),  "results/tables/08_KEGG_Up.csv",   row.names = FALSE)
if (!is.null(kegg_all) && nrow(as.data.frame(kegg_all)) > 0)
  write.csv(as.data.frame(kegg_all), "results/tables/09_KEGG_All.csv",  row.names = FALSE)
if (!is.null(gsea) && nrow(as.data.frame(gsea)) > 0)
  write.csv(as.data.frame(gsea),     "results/tables/10_GSEA_GO_BP.csv", row.names = FALSE)

# Save raw count matrix
raw_cnt_out <- as.data.frame(cnt) %>%
  rownames_to_column("GeneID") %>%
  left_join(
    gene_map %>% dplyr::select(JOIN_ID, SYMBOL, ENSEMBL, ENTREZID),
    by = c("GeneID" = "JOIN_ID")
  ) %>%
  dplyr::select(GeneID, SYMBOL, ENSEMBL, ENTREZID, everything()) %>%
  arrange(GeneID)
write.csv(raw_cnt_out, "results/tables/11_Raw_Counts_Merged.csv", row.names = FALSE)

cat("✓ Tables saved\n")

# ======================================================================
# 11. SUMMARY
# ======================================================================
cat(sprintf(
  "
%s
  GSE265922 — COMPLETE | DESeq2 + apeglm | ZIKV vs Control A549
  Design: n=%d Control vs n=%d ZIKV (3 biological replicates each)
  Input: %d STAR ReadsPerGene.out.tab.gz files merged
  Gene ID type: %s
%s
  DEGs   : %d (↑ %d  |  ↓ %d)   padj < 0.05, |log2FC| > 1
  GO BP  : ↑ %d  ↓ %d    GO MF : ↑ %d  ↓ %d    GO CC : ↑ %d  ↓ %d
  KEGG   : ↑ %d  ↓ %d  All:%d    GSEA  : %d gene sets
  19 PNG @ 600 DPI  |  11 CSV tables
%s
",
  strrep("=", 66), n_ctrl, n_zikv, length(star_files),
  ifelse(exists("id_type"), id_type, "unknown"),
  strrep("=", 66),
  n_up + n_down, n_up, n_down,
  n_safe(go_bp_up), n_safe(go_bp_dn),
  n_safe(go_mf_up), n_safe(go_mf_dn),
  n_safe(go_cc_up), n_safe(go_cc_dn),
  n_safe(kegg_up),  n_safe(kegg_dn),
  n_safe(kegg_all), n_safe(gsea),
  strrep("=", 66)
))

