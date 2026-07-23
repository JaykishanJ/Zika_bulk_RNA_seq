# ======================================================================
# miRNA target analysis (R version) — same logic as miRNA_target_analysis.py
# Intersect miRTarBase validated targets of miR-134-5p / 146a-5p / 146b-5p
# with the 66-gene conserved ZIKV core; build annotated tables + 5-gene panel.
#
# INPUTS
#   ../00_Raw_miRNA_Data/hsa_MTI.csv
#   ../meta-analysis/results/Common_Upregulated_3of3.csv
#   ../GSE*/results/tables/01_DEG_Complete.csv
#
# RUN FROM this folder:  Rscript miRNA_target_analysis.R
# (uses data.table::fread for the ~337 MB miRTarBase file)
# ======================================================================
suppressPackageStartupMessages({
  library(data.table)   # fast reader for the large MTI file
  library(tidyverse)
})

proj      <- ".."
mti_path  <- file.path(proj, "00_Raw_miRNA_Data", "hsa_MTI.csv")
core_path <- file.path(proj, "meta-analysis", "results", "Common_Upregulated_3of3.csv")
datasets  <- c("GSE146423", "GSE233049", "GSE265922")

mirs  <- c("hsa-miR-134-5p", "hsa-miR-146a-5p", "hsa-miR-146b-5p")
panel <- c("IRF7", "IFIT3", "RSAD2", "CCL5", "SAMD9L")

# ---- 1. core genes ----
core <- read.csv(core_path)[[1]]
cat(sprintf("Core genes: %d\n", length(core)))

# ---- 2. load miRTarBase, keep our 3 miRNAs ----
mti <- fread(mti_path, colClasses = "character")
setnames(mti, trimws(names(mti)))
mti <- mti[miRNA %in% mirs]
mti[, `Target Gene` := trimws(`Target Gene`)]
mti[, PMID := gsub("\\.0$", "", `References (PMID)`)]
cat(sprintf("miRTarBase rows for the 3 miRNAs: %d\n", nrow(mti)))

# ---- 3. interactions with the core (+ evidence tier) ----
inter <- map_dfr(mirs, function(mir) {
  sub  <- mti[miRNA == mir]
  func <- unique(sub[grepl("Functional MTI", `Support Type`) &
                     !grepl("Non-Functional", `Support Type`)]$`Target Gene`)
  hits <- sort(intersect(core, unique(sub$`Target Gene`)))
  map_dfr(hits, function(g) {
    ev <- sub[`Target Gene` == g]
    tibble(miRNA = mir, core_gene = g,
           evidence  = ifelse(g %in% func, "Functional (strong)", "Non-functional/weak"),
           n_records = nrow(ev),
           experiments = paste(sort(unique(unlist(strsplit(ev$Experiments, "//")))), collapse = "; "),
           PMIDs = paste(sort(unique(ev$PMID[ev$PMID != ""])), collapse = "; "))
  })
})
write.csv(inter, "miRNA_target_interactions.csv", row.names = FALSE)
for (mir in mirs) cat(sprintf("  %s: %d core targets\n", mir, sum(inter$miRNA == mir)))

# ---- 4. per-dataset log2FC ----
lfc <- map(setNames(datasets, datasets), function(ds) {
  d <- read.csv(file.path(proj, ds, "results", "tables", "01_DEG_Complete.csv"))
  if (!"SYMBOL" %in% names(d) || all(is.na(d$SYMBOL))) d$SYMBOL <- d$Label
  setNames(d$log2FoldChange, d$SYMBOL)
})
tgt <- map(setNames(mirs, gsub("hsa-", "", mirs)),
           ~ inter$core_gene[inter$miRNA == .x])

# ---- 5. annotated 66-gene table ----
ann <- tibble(gene = core)
for (ds in datasets) ann[[paste0(ds, "_log2FC")]] <- round(unname(lfc[[ds]][core]), 2)
ann$miR_146a_5p    <- ifelse(core %in% tgt[["miR-146a-5p"]], "yes", "")
ann$miR_146b_5p    <- ifelse(core %in% tgt[["miR-146b-5p"]], "yes", "")
ann$miR_134_5p     <- ifelse(core %in% tgt[["miR-134-5p"]],  "yes", "")
ann$in_5gene_panel <- ifelse(core %in% panel, "yes", "")
write.csv(ann, "66_core_genes_miRNA_annotated.csv", row.names = FALSE)

# ---- 6. 5-gene panel ----
mir_of <- c(IRF7 = "146a", IFIT3 = "146a+146b", RSAD2 = "146a",
            CCL5 = "146a", SAMD9L = "146a+146b")
pan <- tibble(gene = panel, miRNA = mir_of[panel])
for (ds in datasets) pan[[ds]] <- round(unname(lfc[[ds]][panel]), 2)
write.csv(pan, "5_gene_validation_panel.csv", row.names = FALSE)

cat(sprintf("\nSummary: 146a=%d/66, 146b=%d/66, 134=%d/66\n",
            length(tgt[["miR-146a-5p"]]), length(tgt[["miR-146b-5p"]]),
            length(tgt[["miR-134-5p"]])))
