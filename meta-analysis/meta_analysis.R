# ======================================================================
# ZIKV Meta-Analysis & In Silico Validation
# Combines DEGs from GSE146423, GSE233049, GSE265922
# ======================================================================

# ---- 1. Install & Load Required Packages ----
packages <- c("tidyverse", "pheatmap", "RColorBrewer", "ggrepel", "ggVennDiagram", "enrichR")
bioc_packages <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "STRINGdb", "igraph")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = "https://cloud.r-project.org")
}
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggVennDiagram)
  library(enrichR)
  library(STRINGdb)
  library(igraph)
})

if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
for (d in c("results", "plots")) dir.create(d, recursive = TRUE, showWarnings = FALSE)

theme_pub <- theme_minimal() + theme(
  plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
  axis.title = element_text(size = 14, face = "bold"),
  axis.text = element_text(size = 12, color = "black")
)

# ---- 2. Load Datasets ----
cat("\n═══ Loading Individual Analyses ═══\n")
paths <- list(
  GSE146423 = "../GSE146423/results/tables/01_DEG_Complete.csv",
  GSE233049 = "../GSE233049/results/tables/01_DEG_Complete.csv",
  GSE265922 = "../GSE265922/results/tables/01_DEG_Complete.csv"
)

for (p in names(paths)) {
  if (!file.exists(paths[[p]])) stop(paste("File missing:", paths[[p]], "- Make sure individual analyses have completed."))
}

deg_list <- lapply(paths, read.csv)
up_genes <- list()
dn_genes <- list()
all_lfc <- list()

for (name in names(deg_list)) {
  df <- deg_list[[name]]
  # Fix for labels
  if (!"SYMBOL" %in% colnames(df) && "Label" %in% colnames(df)) {
    df$SYMBOL <- df$Label
  }
  df <- df[!is.na(df$SYMBOL) & df$SYMBOL != "", ]
  
  up_genes[[name]] <- df$SYMBOL[df$padj < 0.05 & df$log2FoldChange > 1 & !is.na(df$padj)]
  dn_genes[[name]] <- df$SYMBOL[df$padj < 0.05 & df$log2FoldChange < -1 & !is.na(df$padj)]
  
  lfc_map <- setNames(df$log2FoldChange, df$SYMBOL)
  all_lfc[[name]] <- lfc_map
}

# ---- 3. Meta-Analysis Intersection ----
common_up <- Reduce(intersect, up_genes)
common_dn <- Reduce(intersect, dn_genes)

# FIX FOR PROBLEM 4 & 2: Added explicit calculation and export of 2-of-3 dataset overlaps 
# so that the main "198/355" headline numbers are reproducible within the main pipeline 
# instead of relying on the external scratch script (temp_query.R).
# Calculate 2-of-3 overlap
all_up_counts <- table(unlist(up_genes))
up_2of3 <- names(all_up_counts[all_up_counts >= 2])
all_dn_counts <- table(unlist(dn_genes))
dn_2of3 <- names(all_dn_counts[all_dn_counts >= 2])

cat(sprintf("Upregulated across all 3 datasets: %d genes\n", length(common_up)))
cat(sprintf("Downregulated across all 3 datasets: %d genes\n", length(common_dn)))
cat(sprintf("Upregulated in at least 2 of 3 datasets: %d genes\n", length(up_2of3)))
cat(sprintf("Downregulated in at least 2 of 3 datasets: %d genes\n", length(dn_2of3)))

write.csv(data.frame(Gene = common_up), "results/Common_Upregulated_3of3.csv", row.names = FALSE)
write.csv(data.frame(Gene = common_dn), "results/Common_Downregulated_3of3.csv", row.names = FALSE)
write.csv(data.frame(Gene = up_2of3), "results/Common_Upregulated_2of3.csv", row.names = FALSE)
write.csv(data.frame(Gene = dn_2of3), "results/Common_Downregulated_2of3.csv", row.names = FALSE)

# Venn Diagrams
venn_up <- ggVennDiagram(up_genes, category.names = names(up_genes)) + 
  scale_fill_gradient(low = "white", high = "#D73027") + 
  labs(title = "Upregulated DEGs Intersection")
ggsave("plots/01_Venn_Up.png", venn_up, width = 8, height = 6, dpi = 600, bg="white")

venn_dn <- ggVennDiagram(dn_genes, category.names = names(dn_genes)) + 
  scale_fill_gradient(low = "white", high = "#4575B4") + 
  labs(title = "Downregulated DEGs Intersection")
ggsave("plots/02_Venn_Down.png", venn_dn, width = 8, height = 6, dpi = 600, bg="white")

# ---- 4. Heatmap of Common Signature ----
all_common <- c(common_up, common_dn)
if(length(all_common) > 1) {
  hm_mat <- matrix(NA, nrow = length(all_common), ncol = 3, dimnames = list(all_common, names(all_lfc)))
  for (name in names(all_lfc)) hm_mat[, name] <- all_lfc[[name]][all_common]
  hm_mat[is.na(hm_mat)] <- 0
  
  png("plots/03_Common_Signature_Heatmap.png", width = 8, height = max(6, length(all_common)*0.15), units = "in", res = 600)
  pheatmap(hm_mat, 
           color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
           cluster_cols = FALSE, clustering_method = "ward.D2",
           main = "Log2 Fold Change of Core ZIKV Signature across datasets",
           fontsize_row = max(4, 12 - (length(all_common) / 20)))
  dev.off()
}

# ---- 5. Functional Enrichment of Common Genes ----
cat("\n═══ Functional Enrichment (Core Signature) ═══\n")
common_entrez <- bitr(all_common, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

if(length(common_entrez) > 0) {
  ego <- enrichGO(gene = common_entrez, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
  if(!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    p_go <- dotplot(ego, showCategory = 15, title = "GO BP Enrichment of Core ZIKV Signature") + theme_pub
    ggsave("plots/04_GO_Enrichment.png", p_go, width = 10, height = 8, dpi = 600, bg="white")
    write.csv(as.data.frame(ego), "results/Common_GO_BP.csv", row.names = FALSE)
  }
}

# ---- 6. In Silico Validation 1: enrichR TF Master Regulators ----
cat("\n═══ In Silico Validation: Master Regulators ═══\n")
setEnrichrSite("Enrichr")
dbs <- c("TRRUST_Transcription_Factors_2019", "ChEA_2016")
enrichr_res <- enrichr(all_common, dbs)

for(db in dbs) {
  res <- enrichr_res[[db]]
  res_sig <- res[res$Adjusted.P.value < 0.05, ]
  write.csv(res, paste0("results/TF_Regulators_", db, ".csv"), row.names = FALSE)
  
  if(nrow(res_sig) > 0) {
    plot_df <- res_sig[1:min(10, nrow(res_sig)), ]
    plot_df$Term <- reorder(plot_df$Term, -log10(plot_df$Adjusted.P.value))
    p_tf <- ggplot(plot_df, aes(x = Term, y = -log10(Adjusted.P.value))) +
      geom_col(fill = "#8073ac") + coord_flip() + theme_pub +
      labs(title = paste("Master Transcription Factors (", db, ")"), x = "Transcription Factor")
    ggsave(paste0("plots/05_TF_", db, ".png"), p_tf, width = 8, height = 6, dpi = 600, bg="white")
  }
}

# ---- 7. In Silico Validation 2: PPI Network (STRINGdb) ----
cat("\n═══ In Silico Validation: PPI Hub Genes ═══\n")
options(timeout = max(300, getOption("timeout")))
tryCatch({
  string_db <- STRINGdb$new(version="11.5", species=9606, score_threshold=400, input_directory="")
  mapped <- string_db$map(data.frame(gene = all_common), "gene", removeUnmappedRows = TRUE)

  if(nrow(mapped) > 0) {
    hits <- mapped$STRING_id
    sub_net <- string_db$get_subnetwork(hits)
    
    deg <- degree(sub_net)
    hub_genes <- sort(deg, decreasing = TRUE)
    
    hub_df <- data.frame(STRING_id = names(hub_genes), Degree = as.numeric(hub_genes))
    hub_df <- merge(hub_df, mapped[, c("gene", "STRING_id")], by = "STRING_id")
    hub_df <- hub_df[order(hub_df$Degree, decreasing = TRUE), ]
    write.csv(hub_df, "results/PPI_Hub_Genes.csv", row.names = FALSE)
    
    cat("Top 5 Hub Genes in ZIKV Signature:\n")
    print(head(hub_df, 5))
    
    if(length(hub_genes) > 1) {
      # Map STRING IDs to gene symbols for labels
      V(sub_net)$label <- mapped$gene[match(V(sub_net)$name, mapped$STRING_id)]
      
      # Highlight top 10% hub genes in red, others in blue
      deg_thresh <- quantile(degree(sub_net), 0.9)
      V(sub_net)$color <- ifelse(degree(sub_net) >= deg_thresh, "#D73027", "#4575B4")
      
      png("plots/06_PPI_Network.png", width = 10, height = 10, units = "in", res = 600)
      set.seed(123)
      plot(sub_net, 
           layout = layout_with_fr,
           vertex.size = 5 + (degree(sub_net) / max(degree(sub_net))) * 8, # Size by degree
           vertex.label = V(sub_net)$label,
           vertex.label.cex = 0.8,
           vertex.label.color = "black",
           vertex.label.dist = 1.2,
           vertex.color = V(sub_net)$color,
           vertex.frame.color = "white",
           edge.width = 0.5,
           edge.color = "gray80",
           main = "Literature-Grade PPI Network of Core ZIKV Signature")
      legend("bottomright", legend=c("Top 10% Hub Genes", "Peripheral Genes"), col=c("#D73027", "#4575B4"), pch=19, bty="n", cex=1.2)
      dev.off()
    }
  }
}, error = function(e) {
  cat("STRINGdb validation skipped due to network timeout or small gene list:\n", e$message, "\n")
})

cat("\n═══ Meta-Analysis Complete ═══\n")
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
