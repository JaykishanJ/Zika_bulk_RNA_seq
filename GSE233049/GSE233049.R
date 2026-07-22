# ======================================================================
# GSE233049: ZIKV-infected A549 — Comprehensive DEG & Enrichment
# GRCh38.p13 NCBI | DESeq2 + apeglm | 600 DPI PNG | Publication-ready
# FULLY FIXED: GSEA plots, emapplot, infinite values, 1:many mapping
# ======================================================================
suppressPackageStartupMessages({
  library(tidyverse); library(DESeq2); library(ggrepel); library(pheatmap)
  library(RColorBrewer); library(clusterProfiler); library(org.Hs.eg.db)
  library(enrichplot); library(EnhancedVolcano); library(apeglm)
})
setwd("d:/Zika_wetlab/GSE233049")
for (d in c("results/tables", "plots/QC", "plots/DEG", "plots/Enrichment", "plots/GSEA"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ── Palette & Theme ───────────────────────────────────────────────────
pal <- list(
  ctrl      = "#2166AC",
  zikv      = "#B2182B",
  up        = "#D73027",
  down      = "#4575B4",
  ns        = "grey80",
  hm        = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  blues     = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
  gsea_div  = c("#4575B4","#91BFDB","#E0F3F8","#FEE090","#FC8D59","#D73027"),
  gsea_up   = "#D73027",
  gsea_down = "#4575B4"
)

theme_pub <- theme_minimal() + theme(
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

save_png <- function(p, path, file, w, h) {
  ggsave(file.path(path, paste0(file, ".png")), p,
         width = w, height = h, dpi = 600, limitsize = FALSE)
}

save_hm <- function(obj, path, file, w, h) {
  png(file.path(path, paste0(file, ".png")), w, h, "in", res = 600)
  grid::grid.newpage(); grid::grid.draw(obj$gtable); dev.off()
}

enrich_plot <- function(obj, title, top_n, path, file, w, h) {
  if (is.null(obj) || nrow(as.data.frame(obj)) == 0) return(invisible())
  p <- dotplot(obj, showCategory = min(top_n, nrow(obj)), title = title) +
    theme_pub + theme(axis.text.y = element_text(size = 11), legend.position = "right")
  save_png(p, path, file, w, h)
}

# ── 1. Load & Prepare ─────────────────────────────────────────────────
cnt <- read.table("GSE233049_raw_counts_GRCh38.p13_NCBI.tsv.gz",
                  header = TRUE, row.names = 1, sep = "\t",
                  stringsAsFactors = FALSE, check.names = FALSE)
# FIX FOR PROBLEM 1: The original code used `cnt <- cnt[, 1:6]`, which incorrectly 
# selected WT_Mock (1:3) and KO_Mock (4:6) samples, leading to an analysis with no ZIKV infection.
# The correct columns for WT_ZIKV are 7:9, so we now select `c(1:3, 7:9)`.
cnt <- cnt[, c(1:3, 7:9)]
colnames(cnt) <- c("Control_1","Control_2","Control_3","ZIKV_1","ZIKV_2","ZIKV_3")
cnt <- as.matrix(cnt); storage.mode(cnt) <- "integer"
stopifnot(all(cnt == floor(cnt), na.rm = TRUE))
cnt <- cnt[rowSums(cnt) > 0, ]

cat(sprintf("Loaded: %d genes x %d samples\n", nrow(cnt), ncol(cnt)))
cat(sprintf("Library sizes (M): %s\n",
            paste(round(colSums(cnt)/1e6, 1), collapse=", ")))

# ── 2. Gene ID Mapping (FIXED: deduplicate 1:many mappings) ────────────
map <- bitr(rownames(cnt), fromType = "ENTREZID",
            toType = c("SYMBOL","ENSEMBL"), OrgDb = org.Hs.eg.db) %>%
  distinct(ENTREZID, .keep_all = TRUE)          # ← FIX: keep only first match per ENTREZID

g <- intersect(rownames(cnt), map$ENTREZID)
cnt <- cnt[g, ]
map <- map[match(g, map$ENTREZID), ]
cat(sprintf("After symbol mapping: %d genes\n", nrow(cnt)))

# ── 3. DESeq2 ─────────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(cnt,
                              data.frame(condition = factor(rep(c("Control","ZIKV"), each = 3),
                                                            levels = c("Control","ZIKV")),
                                         row.names = colnames(cnt)), ~ condition)
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
cat(sprintf("After prefilter: %d genes\n", nrow(dds)))

dds <- DESeq(dds)
res <- lfcShrink(dds, coef = "condition_ZIKV_vs_Control", type = "apeglm")

# ── Build DEG table with ROBUST infinite-padj handling ─────────────────
deg <- as.data.frame(res) %>%
  rownames_to_column("GeneID") %>%
  left_join(rename(map, GeneID = ENTREZID), by = "GeneID") %>%
  mutate(
    sig = padj < 0.05 & abs(log2FoldChange) > 1,
    dir = case_when(
      sig & log2FoldChange > 1  ~ "Up",
      sig & log2FoldChange < -1 ~ "Down",
      TRUE                       ~ "NS"
    ),
    # FIX: cap infinite -log10(padj) at a sensible value
    nlp = {
      x <- -log10(padj)
      finite_max <- max(x[is.finite(x)], na.rm = TRUE)
      ifelse(is.infinite(x), finite_max * 1.02, x)
    }
  ) %>%
  arrange(padj)

deg$Mean_Ctrl <- rowMeans(counts(dds, normalized = TRUE)[deg$GeneID, 1:3])
deg$Mean_ZIKV <- rowMeans(counts(dds, normalized = TRUE)[deg$GeneID, 4:6])

n_up   <- sum(deg$dir == "Up")
n_down <- sum(deg$dir == "Down")

cat(sprintf("\n═══ DEGs: %d (↑%d ↓%d) ═══\n", n_up + n_down, n_up, n_down))
cat("\nTop 15 Upregulated:\n")
deg %>% filter(dir == "Up") %>% head(15) %>%
  dplyr::select(SYMBOL, log2FoldChange, padj, Mean_Ctrl, Mean_ZIKV) %>% print()
cat("\nAll Downregulated:\n")
deg %>% filter(dir == "Down") %>%
  dplyr::select(SYMBOL, log2FoldChange, padj, Mean_Ctrl, Mean_ZIKV) %>% print()

# ── 4. QC Plots ───────────────────────────────────────────────────────
vsd <- vst(dds, blind = FALSE)

# 4a. PCA
pca <- plotPCA(vsd, "condition", returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))
save_png(
  ggplot(pca, aes(PC1, PC2, color = condition, label = name)) +
    geom_point(size = 6, alpha = 0.9) +
    geom_text_repel(size = 5, fontface = "bold", show.legend = FALSE, box.padding = 0.8) +
    scale_color_manual(values = c(Control = pal$ctrl, ZIKV = pal$zikv)) +
    labs(title = "Principal Component Analysis",
         subtitle = sprintf("ZIKV vs Control A549 (%d genes)", nrow(dds)),
         x = paste0("PC1: ", pv[1], "%"), y = paste0("PC2: ", pv[2], "%")) +
    coord_fixed(1) + theme_pub, "plots/QC", "01_PCA", 9, 8)
cat("✓ 01_PCA.png\n")

# 4b. Sample Distances
dists <- dist(t(assay(vsd))); dm <- as.matrix(dists)
dimnames(dm) <- list(
  paste0(vsd$condition," — ",colnames(vsd)),
  paste0(vsd$condition," — ",colnames(vsd))
)
ann_col <- data.frame(Condition = vsd$condition, row.names = colnames(vsd))
save_hm(
  pheatmap(dm, clustering_distance_rows = dists,
           clustering_distance_cols = dists, annotation_col = ann_col,
           annotation_colors = list(Condition = c(Control = pal$ctrl, ZIKV = pal$zikv)),
           main = "Sample-to-Sample Distances (VST)", color = pal$blues,
           display_numbers = TRUE, number_format = "%.0f",
           fontsize = 13, fontsize_number = 11, silent = TRUE),
  "plots/QC", "02_Sample_Distances", 9, 8)
cat("✓ 02_Sample_Distances.png\n")

# 4c. Dispersion
png(file.path("plots/QC","03_Dispersion.png"), 8, 6, "in", res = 600)
plotDispEsts(dds, main = "DESeq2: Dispersion Estimates", cex = 0.6); dev.off()
cat("✓ 03_Dispersion.png\n")

# 4d. Expression Density
save_png(
  as.data.frame(counts(dds, normalized = TRUE)) %>% rownames_to_column("GeneID") %>%
    pivot_longer(-GeneID, names_to = "Sample", values_to = "Count") %>%
    mutate(Group = ifelse(grepl("Control", Sample), "Control", "ZIKV"),
           log2Count = log2(Count + 1)) %>%
    ggplot(aes(log2Count, fill = Group, color = Group)) +
    geom_density(alpha = 0.3, linewidth = 0.8) +
    scale_fill_manual(values = c(Control = pal$ctrl, ZIKV = pal$zikv)) +
    scale_color_manual(values = c(Control = pal$ctrl, ZIKV = pal$zikv)) +
    labs(title = "Expression Density — Normalized Counts",
         subtitle = "log2(DESeq2 normalized counts + 1)",
         x = expression(Log[2] * "(Norm. Counts + 1)"), y = "Density") + theme_pub,
  "plots/QC", "04_Expression_Density", 8, 6)
cat("✓ 04_Expression_Density.png\n")

# 4e. Library Sizes
lib <- data.frame(Sample = colnames(cnt), Reads_M = colSums(cnt)/1e6,
                  Group = factor(rep(c("Control","ZIKV"), each = 3)))
save_png(
  ggplot(lib, aes(Sample, Reads_M, fill = Group)) +
    geom_bar(stat = "identity", alpha = 0.85, width = 0.7) +
    geom_text(aes(label = sprintf("%.1fM", Reads_M)), vjust = -0.5, size = 4) +
    scale_fill_manual(values = c(Control = pal$ctrl, ZIKV = pal$zikv)) +
    labs(title = "Sequencing Library Sizes", x = "", y = "Million Reads") +
    theme_pub + theme(axis.text.x = element_text(angle = 30, hjust = 1)),
  "plots/QC", "05_Library_Sizes", 8, 6)
cat("✓ 05_Library_Sizes.png\n")

# ── 5. DEG Plots ──────────────────────────────────────────────────────
top30 <- deg %>% filter(sig) %>% slice_min(padj, n = 30)

# 5a. Volcano
save_png(
  ggplot(deg, aes(log2FoldChange, nlp, color = dir, alpha = dir)) +
    geom_point(size = 1.8) +
    scale_color_manual(values = c(Up = pal$up, Down = pal$down, NS = pal$ns), name = NULL) +
    scale_alpha_manual(values = c(Up = 0.8, Down = 0.8, NS = 0.25), guide = "none") +
    geom_hline(yintercept = -log10(0.05), lty = "longdash", color = "grey50") +
    geom_vline(xintercept = c(-1, 1), lty = "longdash", color = "grey50") +
    geom_text_repel(data = top30, aes(label = SYMBOL), size = 3.2, fontface = "italic",
                    max.overlaps = 30, box.padding = 0.6, segment.size = 0.3) +
    labs(title = "ZIKV vs Control — A549 Cells",
         subtitle = sprintf("%d Up  ·  %d Down  ·  %d DEGs", n_up, n_down, n_up + n_down),
         x = expression(Log[2]~"Fold Change"),
         y = expression(-Log[10]~"P-adj")) + theme_pub,
  "plots/DEG", "06_Volcano", 10, 8.5)
cat("✓ 06_Volcano.png\n")

# 5b. EnhancedVolcano
p_ev <- EnhancedVolcano(deg, lab = deg$SYMBOL, x = "log2FoldChange", y = "padj",
                title = "ZIKV vs Control (A549) — DESeq2 | GRCh38.p13",
                subtitle = sprintf("%d up  ·  %d down  ·  %d total", n_up, n_down, n_up + n_down),
                caption = sprintf("padj < 0.05, |log2FC| > 1 | %d genes", nrow(deg)),
                pCutoff = 0.05, FCcutoff = 1, pointSize = 2, labSize = 3.5, colAlpha = 0.6,
                legendPosition = "bottom", drawConnectors = TRUE, max.overlaps = 25,
                selectLab = top30$SYMBOL[1:min(25, nrow(top30))],
                border = "full", borderWidth = 0.8)
save_png(p_ev, "plots/DEG", "07_EnhancedVolcano", 11, 9)
cat("✓ 07_EnhancedVolcano.png\n")

# 5c. MA Plot
png(file.path("plots/DEG","08_MA_Plot.png"), 9, 7, "in", res = 600)
DESeq2::plotMA(res, main = "MA Plot — apeglm-shrunken LFC", ylim = c(-6, 6), alpha = 0.05)
dev.off()
cat("✓ 08_MA_Plot.png\n")

# 5d. P-value Histogram
save_png(
  ggplot(deg, aes(pvalue)) +
    geom_histogram(bins = 50, fill = "grey60", color = "grey30", alpha = 0.7, na.rm = TRUE) +
    geom_vline(xintercept = 0.05, lty = "dashed", color = pal$up, linewidth = 1) +
    labs(title = "P-value Distribution", subtitle = sprintf("%d genes tested", nrow(deg)),
         x = "Raw P-value", y = "Frequency") + theme_pub,
  "plots/DEG", "09_Pvalue_Histogram", 8, 6)
cat("✓ 09_Pvalue_Histogram.png\n")

# ── 6. Heatmap ────────────────────────────────────────────────────────
top50 <- deg %>% filter(sig) %>% slice_min(padj, n = 50)
hm_mat <- assay(vsd)[top50$GeneID, ]
lbls <- top50$SYMBOL
lbls[is.na(lbls) | lbls == ""] <- top50$GeneID[is.na(lbls) | lbls == ""]
dup <- duplicated(lbls)
if (any(dup)) for (i in which(dup)) lbls[i] <- paste0(lbls[i], "_", i)
rownames(hm_mat) <- lbls
ann <- data.frame(Condition = vsd$condition, row.names = colnames(vsd))

save_hm(
  pheatmap(hm_mat, scale = "row", annotation_col = ann,
           annotation_colors = list(Condition = c(Control = pal$ctrl, ZIKV = pal$zikv)),
           main = "Top 50 DEGs — ZIKV vs Control A549", show_rownames = TRUE,
           fontsize_row = 6, fontsize_col = 12, clustering_method = "ward.D2",
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           color = pal$hm, border_color = NA, silent = TRUE),
  "plots/DEG", "10_Heatmap_Top50", 10, 13)
cat("✓ 10_Heatmap_Top50.png\n")

# ── 7. Enrichment ─────────────────────────────────────────────────────
cat("\n═══ Running Enrichment ═══\n")

bg  <- unique(deg$GeneID)
up  <- unique(deg$GeneID[deg$dir == "Up"])
dn  <- unique(deg$GeneID[deg$dir == "Down"])
all <- unique(deg$GeneID[deg$sig])

run_go <- function(g, o) {
  enrichGO(g, universe = bg, OrgDb = org.Hs.eg.db, ont = o,
           keyType = "ENTREZID", pAdjustMethod = "BH",
           pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
}
run_kegg <- function(g) {
  enrichKEGG(g, universe = bg, organism = "hsa",
             keyType = "ncbi-geneid", pAdjustMethod = "BH",
             pvalueCutoff = 0.05, qvalueCutoff = 0.2)
}

go_bp_up <- run_go(up, "BP")
go_bp_dn <- run_go(dn, "BP")   # may be empty — that's OK
go_mf_up <- run_go(up, "MF")
go_mf_dn <- run_go(dn, "MF")
go_cc_up <- run_go(up, "CC")
go_cc_dn <- run_go(dn, "CC")
kegg_up  <- run_kegg(up)
kegg_dn  <- run_kegg(dn)
kegg_all <- run_kegg(all)

cat(sprintf("GO BP: ↑%d ↓%d | GO MF: ↑%d ↓%d | GO CC: ↑%d ↓%d | KEGG: ↑%d ↓%d All:%d\n",
            nrow(go_bp_up), nrow(go_bp_dn), nrow(go_mf_up), nrow(go_mf_dn),
            nrow(go_cc_up), nrow(go_cc_dn), nrow(kegg_up), nrow(kegg_dn), nrow(kegg_all)))

# ── 8. Enrichment Plots ───────────────────────────────────────────────
enrich_plot(go_bp_up, "GO BP — Upregulated", 20, "plots/Enrichment", "11_GO_BP_Up", 13, 9)
enrich_plot(go_mf_up, "GO MF — Upregulated", 15, "plots/Enrichment", "12_GO_MF_Up", 12, 7.5)
enrich_plot(kegg_up,  "KEGG — Upregulated",  15, "plots/Enrichment", "13_KEGG_Up", 12, 7.5)
enrich_plot(kegg_all, "KEGG — All DEGs",      15, "plots/Enrichment", "14_KEGG_All", 12, 7.5)
cat("✓ 11-14 Enrichment dotplots\n")

# 8e. GO BP Up vs Down comparison — only if BOTH have terms
if (nrow(as.data.frame(go_bp_up)) > 0 && nrow(as.data.frame(go_bp_dn)) > 0) {
  comp <- bind_rows(
    as.data.frame(go_bp_up) %>% mutate(Reg = "Upregulated"),
    as.data.frame(go_bp_dn) %>% mutate(Reg = "Downregulated")
  ) %>%
    group_by(Reg) %>% slice_min(p.adjust, n = 15) %>% ungroup()
  save_png(
    ggplot(comp, aes(
      reorder(str_wrap(Description, 55), -log10(p.adjust)),
      -log10(p.adjust), fill = Reg)) +
      geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.9) +
      coord_flip() +
      scale_fill_manual(values = c(Upregulated = pal$up, Downregulated = pal$down)) +
      labs(title = "GO BP: Up vs Downregulated", x = "",
           y = expression(-Log[10]*"(P-adj)")) + theme_pub,
    "plots/Enrichment", "15_GO_BP_Up_vs_Down", 15, 8)
  cat("✓ 15_GO_BP_Up_vs_Down.png\n")
} else {
  cat("⊙ 15_GO_BP_Up_vs_Down.png SKIPPED — no downregulated GO BP terms\n")
}

# 8f. GO BP Network — FIXED: tryCatch around emapplot
if (nrow(as.data.frame(go_bp_up)) > 5) {
  res_emap <- tryCatch({
    ego_pair <- pairwise_termsim(go_bp_up)
    # Pre-filter to avoid overloaded network
    n_show <- min(30, nrow(as.data.frame(go_bp_up)))
    p <- emapplot(ego_pair, showCategory = n_show,
                  layout = "nicely", node_label = "category") +
      labs(title = "GO BP Enrichment Map — Upregulated") +
      theme_pub + theme(
        legend.position  = "right",
        legend.text      = element_text(size = 9),
        legend.title     = element_text(size = 10)
      )
    save_png(p, "plots/Enrichment", "16_GO_BP_Network", 14, 12)
    TRUE
  }, error = function(e) {
    cat(sprintf("⊙ 16_GO_BP_Network.png FAILED: %s\n", e$message))
    return(FALSE)
  })
  if (isTRUE(res_emap)) cat("✓ 16_GO_BP_Network.png\n")
} else {
  cat("⊙ 16_GO_BP_Network.png SKIPPED — ≤5 GO BP terms\n")
}

# ── 9. GSEA — FULLY REWRITTEN ─────────────────────────────────────────
cat("\n═══ Running GSEA ═══\n")

set.seed(42)   # reproducibility

# Build ranked gene list: sorted by log2FoldChange (most up → most down)
rnk <- deg %>%
  mutate(
    p_safe = ifelse(pvalue == 0 | is.na(pvalue), 1e-300, pvalue),
    rank_metric = sign(log2FoldChange) * -log10(p_safe)
  ) %>%
  filter(!is.na(rank_metric)) %>%
  distinct(GeneID, .keep_all = TRUE) %>%
  arrange(desc(rank_metric)) %>%
  pull(rank_metric, GeneID) %>%
  sort(decreasing = TRUE)

cat(sprintf("Ranked list: %d genes (range: %.2f to %.2f)\n",
            length(rnk), min(rnk), max(rnk)))

gsea <- gseGO(
  rnk,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  keyType       = "ENTREZID",
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  seed          = TRUE,
  verbose       = FALSE
)

n_gsea <- nrow(as.data.frame(gsea))
cat(sprintf("GSEA GO BP: %d significant gene sets\n", n_gsea))

# ── 9a. Running Score ──────────────────────────────────────────────────
if (n_gsea > 0) {
  
  # Pick top 4 pathways by |NES| (strongest activation + strongest suppression)
  gsea_df <- as.data.frame(gsea) %>% mutate(absNES = abs(NES)) %>% arrange(desc(absNES))
  top_ids <- gsea_df$ID[1:min(4, n_gsea)]
  
  png(file.path("plots/GSEA","17_GSEA_RunningScore.png"), 14, 10, "in", res = 600)
  # FIXED: gseaplot2 with explicit IDs (not row indices from which())
  print(gseaplot2(gsea, geneSetID = top_ids, base_size = 14,
            title    = "GSEA GO BP — Top 4 Pathways by |NES|",
            color    = c(pal$gsea_up, pal$gsea_down, "#FB8C00", "#7B1FA2"),
            pvalue_table = FALSE))
  dev.off()
  cat(sprintf("✓ 17_GSEA_RunningScore.png (%d pathways)\n", length(top_ids)))
  
  # ── 9b. Ridge Plot ──────────────────────────────────────────────────
  if (n_gsea >= 5) {
    tryCatch({
      # FIXED: ridgeplot with proper p.adjust colour mapping
      p_ridge <- ridgeplot(gsea, showCategory = min(20, n_gsea)) +
        scale_fill_gradientn(
          colours = pal$gsea_div,
          trans   = "reverse",          # low p.adj = vibrant colour
          name    = "Adjusted\nP-value"
        ) +
        labs(
          title    = "GSEA GO BP — ZIKV vs Control A549",
          subtitle = sprintf("Ridge density of %d significant gene sets", n_gsea)
        ) +
        theme_pub + theme(
          legend.position = "right",
          axis.text.y     = element_text(size = 10)
        )
      save_png(p_ridge, "plots/GSEA", "18_GSEA_RidgePlot", 13, 9)
      cat("✓ 18_GSEA_RidgePlot.png\n")
    }, error = function(e) cat(sprintf("⊙ 18_GSEA_RidgePlot FAILED: %s\n", e$message)))
  }
  
  # ── 9c. Dotplot ─────────────────────────────────────────────────────
  if (n_gsea >= 3) {
    tryCatch({
      p_dot <- dotplot(gsea, showCategory = min(20, n_gsea),
                       title = "GSEA GO BP — Significant Gene Sets",
                       split = ".sign") +
        facet_grid(. ~ .sign, scales = "free_y", space = "free_y") +
        scale_color_gradientn(
          colours = pal$gsea_div,
          trans   = "reverse",
          name    = "Adjusted\nP-value"
        ) +
        theme_pub + theme(
          strip.background = element_rect(fill = "grey95", color = "grey60"),
          strip.text       = element_text(size = 12, face = "bold"),
          axis.text.y      = element_text(size = 10)
        )
      save_png(p_dot, "plots/GSEA", "19_GSEA_Dotplot", 14, 9)
      cat("✓ 19_GSEA_Dotplot.png\n")
    }, error = function(e) cat(sprintf("⊙ 19_GSEA_Dotplot FAILED: %s\n", e$message)))
  }
  
  # ── 9d. Enrichment Map (network) ────────────────────────────────────
  if (n_gsea >= 10) {
    tryCatch({
      gsea_pair <- pairwise_termsim(gsea, showCategory = min(50, n_gsea))
      p_emap <- emapplot(gsea_pair, showCategory = min(30, n_gsea),
                         layout = "nicely", node_label = "category") +
        scale_color_gradient2(low = pal$gsea_down, mid = "white",
                              high = pal$gsea_up, midpoint = 0, name = "NES") +
        labs(title = "GSEA GO BP — Enrichment Map") +
        theme_pub + theme(
          legend.position  = "right",
          legend.text      = element_text(size = 9),
          legend.title     = element_text(size = 10)
        )
      save_png(p_emap, "plots/GSEA", "20_GSEA_EnrichmentMap", 16, 12)
      cat("✓ 20_GSEA_EnrichmentMap.png\n")
    }, error = function(e) cat(sprintf("⊙ 20_GSEA_EnrichmentMap FAILED: %s\n", e$message)))
  }
  
  # ── 9e. Gene-Concept Network ────────────────────────────────────────
  if (n_gsea >= 3) {
    tryCatch({
      p_cnet <- cnetplot(gsea, showCategory = min(5, n_gsea),
                         foldChange = rnk,
                         node_label = "all") +
        scale_color_gradient2(low = pal$gsea_down, mid = "grey90",
                              high = pal$gsea_up, midpoint = 0, name = "log2FC") +
        labs(title = "GSEA GO BP — Gene-Concept Network",
             subtitle = "Top 5 gene sets | Nodes coloured by fold-change") +
        theme_pub + theme(legend.position = "right")
      save_png(p_cnet, "plots/GSEA", "21_GSEA_CnetPlot", 14, 10)
      cat("✓ 21_GSEA_CnetPlot.png\n")
    }, error = function(e) cat(sprintf("⊙ 21_GSEA_CnetPlot FAILED: %s\n", e$message)))
  }
  
  # ── 9f. Leading-edge Heatplot ───────────────────────────────────────
  if (n_gsea >= 5) {
    tryCatch({
      p_heat <- heatplot(gsea, showCategory = min(10, n_gsea),
                         foldChange = rnk) +
        scale_fill_gradient2(low = pal$gsea_down, mid = "white",
                             high = pal$gsea_up, midpoint = 0, name = "log2FC") +
        labs(title = "GSEA GO BP — Leading-Edge Heatmap",
             subtitle = "Top 10 gene sets | log2FC of core enrichment genes") +
        theme_pub + theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 9),
          legend.position = "right"
        )
      save_png(p_heat, "plots/GSEA", "22_GSEA_Heatplot", 16, 10)
      cat("✓ 22_GSEA_Heatplot.png\n")
    }, error = function(e) cat(sprintf("⊙ 22_GSEA_Heatplot FAILED: %s\n", e$message)))
  }
  
} else {
  cat("⊙ GSEA yielded 0 significant gene sets — all GSEA plots skipped.\n")
  cat("  Hint: try relaxing pvalueCutoff to 0.1 or setting it to Inf,\n")
  cat("  then filtering by enrichmentScore in post-processing.\n")
}

# ── 10. Save Tables ───────────────────────────────────────────────────
nc <- as.data.frame(counts(dds, normalized = TRUE)) %>%
  rownames_to_column("GeneID") %>%
  left_join(rename(map, GeneID = ENTREZID), by = "GeneID") %>%
  dplyr::select(GeneID, SYMBOL, ENSEMBL, everything()) %>%
  arrange(GeneID)

write.csv(deg,                        "results/tables/01_DEG_Complete.csv",      row.names = FALSE)
write.csv(filter(deg, sig),           "results/tables/02_DEG_Significant.csv",    row.names = FALSE)
write.csv(filter(deg, dir == "Up"),   "results/tables/03_DEG_Upregulated.csv",    row.names = FALSE)
write.csv(filter(deg, dir == "Down"), "results/tables/04_DEG_Downregulated.csv",  row.names = FALSE)
write.csv(nc,                         "results/tables/05_Normalized_Counts.csv",  row.names = FALSE)
if (nrow(as.data.frame(go_bp_up)) > 0) write.csv(as.data.frame(go_bp_up), "results/tables/06_GO_BP_Up.csv",   row.names = FALSE)
if (nrow(as.data.frame(go_mf_up)) > 0) write.csv(as.data.frame(go_mf_up), "results/tables/07_GO_MF_Up.csv",   row.names = FALSE)
if (nrow(as.data.frame(kegg_up))  > 0) write.csv(as.data.frame(kegg_up),  "results/tables/08_KEGG_Up.csv",     row.names = FALSE)
if (nrow(as.data.frame(kegg_all)) > 0) write.csv(as.data.frame(kegg_all), "results/tables/09_KEGG_All.csv",    row.names = FALSE)
if (n_gsea > 0)                     write.csv(as.data.frame(gsea),        "results/tables/10_GSEA_GO_BP.csv",  row.names = FALSE)
write.csv(as.data.frame(cnt),                                         "results/tables/11_Raw_Counts.csv")
cat("✓ All tables saved\n")

# ── 11. Summary ───────────────────────────────────────────────────────
cat(sprintf("
%s
  GSE233049 — COMPLETE  |  GRCh38.p13  |  DESeq2 + apeglm
%s
  DEGs   : %d (↑ %d  |  ↓ %d)   padj < 0.05, |log2FC| > 1
  GO BP  : ↑ %d  ↓ %d    GO MF : ↑ %d  ↓ %d    GO CC : ↑ %d  ↓ %d
  KEGG   : ↑ %d  ↓ %d  All:%d    GSEA  : %d gene sets
  Plots  : up to 22 PNG @ 600 DPI   |   Tables : 11 CSV
%s
",
            strrep("=", 66), strrep("=", 66),
            n_up + n_down, n_up, n_down,
            nrow(as.data.frame(go_bp_up)), nrow(as.data.frame(go_bp_dn)),
            nrow(as.data.frame(go_mf_up)), nrow(as.data.frame(go_mf_dn)),
            nrow(as.data.frame(go_cc_up)), nrow(as.data.frame(go_cc_dn)),
            nrow(as.data.frame(kegg_up)), nrow(as.data.frame(kegg_dn)),
            nrow(as.data.frame(kegg_all)), n_gsea,
            strrep("=", 66)))

