paths <- list(
  G14 = "d:/Zika_wetlab/GSE146423/results/tables/01_DEG_Complete.csv",
  G23 = "d:/Zika_wetlab/GSE233049/results/tables/01_DEG_Complete.csv",
  G26 = "d:/Zika_wetlab/GSE265922/results/tables/01_DEG_Complete.csv"
)
degs <- lapply(paths, read.csv)
up_genes <- list()
for(n in names(degs)) {
  df <- degs[[n]]
  if(!"SYMBOL" %in% colnames(df)) df$SYMBOL <- df$Label
  up_genes[[n]] <- df$SYMBOL[df$padj < 0.05 & df$log2FoldChange > 1 & !is.na(df$padj) & !is.na(df$SYMBOL) & df$SYMBOL != ""]
}
all_up <- table(unlist(up_genes))
in_two <- names(all_up[all_up >= 2])
cat("Genes up in at least 2 datasets:", length(in_two), "\n")
classic <- c("IFIT1","IFIT2","IFIT3","ISG15","MX1","MX2","OAS1","OAS2","OAS3","CXCL10","RSAD2","STAT1","IRF7","DDX58","PARP9","HERC5")
found <- intersect(in_two, classic)
cat(paste(found, collapse=", "), "\n")
