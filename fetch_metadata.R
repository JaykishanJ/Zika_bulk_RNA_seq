# ======================================================================
# fetch_metadata.R
# Fetches experimental metadata for datasets from GEO
# ======================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(tidyverse)
})

if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd("d:/Zika_wetlab")
}

cat("\n═══ Fetching Metadata from NCBI GEO ═══\n")

geo_ids <- c("GSE146423", "GSE233049", "GSE265922")

metadata <- data.frame()

for (id in geo_ids) {
  cat(sprintf("Fetching %s...\n", id))
  tryCatch({
    gse <- getGEO(id, GSEMatrix = TRUE, getGPL = FALSE)
    if (length(gse) > 0) {
      pd <- pData(gse[[1]])
      
      # Try to safely extract columns
      titles <- if ("title" %in% colnames(pd)) pd$title else rownames(pd)
      source <- if ("source_name_ch1" %in% colnames(pd)) pd$source_name_ch1 else "Unknown"
      chars  <- if ("characteristics_ch1" %in% colnames(pd)) pd$characteristics_ch1 else "None"
      
      df <- data.frame(
        Dataset = id,
        Sample = rownames(pd),
        Title = titles,
        Source = source,
        Characteristics = chars
      )
      metadata <- bind_rows(metadata, df)
      cat(sprintf("  ✓ Added %d samples.\n", nrow(df)))
    }
  }, error = function(e) {
    cat("  x Error fetching", id, ":", e$message, "\n")
  })
}

dir.create("meta-analysis/results", recursive = TRUE, showWarnings = FALSE)
write.csv(metadata, "meta-analysis/results/Dataset_Metadata.csv", row.names = FALSE)
cat("✓ Metadata saved to meta-analysis/results/Dataset_Metadata.csv\n")
