# ======================================================================
# Publication-Ready PPI Network Construction & Hub Gene Analysis
# Companion script to meta_analysis.R (ZIKV Meta-Analysis)
#
# INPUTS  : results/Common_Upregulated_3of3.csv
#           results/Common_Downregulated_3of3.csv
#           (falls back to *_2of3.csv if 3of3 lists are too small)
# OUTPUTS : results/ppi/  -> node table, edge table, hub ranking, modules
#           plots/ppi/    -> publication figures (PNG 600 dpi + PDF vector)
#           Cytoscape     -> .graphml, edge/node .tsv (Cytoscape-importable)
# ======================================================================

# ---- 1. Packages ------------------------------------------------------
cran_pkgs <- c("tidyverse", "igraph", "ggraph", "tidygraph", "ggrepel",
               "RColorBrewer", "scales", "httr", "readr")
bioc_pkgs <- c("STRINGdb")

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(ggrepel)
  library(RColorBrewer)
  library(scales)
  library(httr)
})

if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

dir.create("results/ppi", recursive = TRUE, showWarnings = FALSE)
dir.create("plots/ppi",   recursive = TRUE, showWarnings = FALSE)

# ---- 2. USER PARAMETERS ----------------------------------------------
PARAMS <- list(
  species          = 9606,      # 9606 = Homo sapiens
  string_version   = "12.0",    # STRING release ("11.5" also valid)
  score_threshold  = 700,       # 400=medium, 700=high, 900=highest confidence
  min_genes        = 10,        # if 3of3 list < this, fall back to 2of3
  top_hubs_label   = 30,        # how many hub labels to draw on the figure
  hub_top_n        = 20,        # hubs reported in the bar plot / table
  remove_isolated  = TRUE,      # drop nodes with degree 0
  largest_component_only = FALSE,# TRUE = plot only the giant component
  seed             = 42
)
set.seed(PARAMS$seed)
options(timeout = max(600, getOption("timeout")))

# ---- 3. Load gene lists from meta-analysis ---------------------------
cat("\n=== Loading DEG lists ===\n")

read_gene_list <- function(path) {
  if (!file.exists(path)) return(character(0))
  x <- read.csv(path, stringsAsFactors = FALSE)
  unique(as.character(x[[1]]))[!is.na(x[[1]]) & nzchar(x[[1]])]
}

up <- read_gene_list("results/Common_Upregulated_3of3.csv")
dn <- read_gene_list("results/Common_Downregulated_3of3.csv")

if (length(c(up, dn)) < PARAMS$min_genes) {
  cat("3-of-3 core list is small; falling back to 2-of-3 lists.\n")
  up <- read_gene_list("results/Common_Upregulated_2of3.csv")
  dn <- read_gene_list("results/Common_Downregulated_2of3.csv")
}

gene_df <- bind_rows(
  tibble(gene = up, direction = "Up"),
  tibble(gene = dn, direction = "Down")
) %>% distinct(gene, .keep_all = TRUE)

stopifnot(nrow(gene_df) > 2)
cat(sprintf("Input genes: %d (%d up, %d down)\n",
            nrow(gene_df), sum(gene_df$direction == "Up"),
            sum(gene_df$direction == "Down")))

# ---- 4. Retrieve STRING interactions ---------------------------------
# Primary: STRING REST API (returns symbols + all evidence channels directly).
# Fallback: STRINGdb Bioconductor package.
cat("\n=== Querying STRING ===\n")

get_string_api <- function(genes, species, threshold, version) {
  url <- "https://string-db.org/api/tsv/network"
  res <- httr::POST(
    url,
    body = list(
      identifiers    = paste(genes, collapse = "\r"),
      species        = species,
      required_score = threshold,
      caller_identity = "R_PPI_pipeline"
    ),
    encode = "form"
  )
  if (httr::status_code(res) != 200) stop("STRING API returned an error.")
  readr::read_tsv(httr::content(res, "text", encoding = "UTF-8"),
                  show_col_types = FALSE)
}

edges_raw <- tryCatch(
  get_string_api(gene_df$gene, PARAMS$species,
                 PARAMS$score_threshold, PARAMS$string_version),
  error = function(e) {
    cat("API failed (", e$message, ") - using STRINGdb package.\n")
    library(STRINGdb)
    sdb <- STRINGdb$new(version = PARAMS$string_version,
                        species = PARAMS$species,
                        score_threshold = PARAMS$score_threshold,
                        input_directory = "")
    mapped <- sdb$map(as.data.frame(gene_df), "gene", removeUnmappedRows = TRUE)
    net <- sdb$get_interactions(mapped$STRING_id)
    key <- setNames(mapped$gene, mapped$STRING_id)
    tibble(
      preferredName_A = key[net$from],
      preferredName_B = key[net$to],
      score           = net$combined_score / 1000
    )
  }
)

edges <- edges_raw %>%
  transmute(from = preferredName_A,
            to   = preferredName_B,
            combined_score = as.numeric(score)) %>%
  filter(!is.na(from), !is.na(to), from != to) %>%
  # collapse reciprocal duplicates (A-B and B-A)
  mutate(a = pmin(from, to), b = pmax(from, to)) %>%
  group_by(a, b) %>%
  summarise(combined_score = max(combined_score), .groups = "drop") %>%
  rename(from = a, to = b)

if (nrow(edges) == 0) stop("No interactions passed the score threshold. Lower PARAMS$score_threshold.")
cat(sprintf("Interactions retrieved: %d\n", nrow(edges)))

# ---- 5. Build graph & topology ---------------------------------------
g <- graph_from_data_frame(edges, directed = FALSE,
                           vertices = data.frame(name = unique(c(edges$from, edges$to))))
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE,
              edge.attr.comb = list(combined_score = "max"))

if (PARAMS$remove_isolated) g <- delete_vertices(g, which(degree(g) == 0))
if (PARAMS$largest_component_only) {
  comp <- components(g)
  g <- induced_subgraph(g, which(comp$membership == which.max(comp$csize)))
}

# --- MCC (Maximal Clique Centrality, cytoHubba's top-performing method) ---
mcc_score <- function(graph) {
  cl <- max_cliques(graph, min = 2)
  sapply(V(graph)$name, function(v) {
    idx <- which(sapply(cl, function(c) v %in% names(c)))
    if (length(idx) == 0) return(0)
    sum(factorial(pmax(sapply(cl[idx], length) - 1, 0)))
  })
}

V(g)$Direction   <- gene_df$direction[match(V(g)$name, gene_df$gene)]
V(g)$Degree      <- degree(g)
V(g)$Betweenness <- betweenness(g, normalized = TRUE)
V(g)$Closeness   <- closeness(g, normalized = TRUE)
V(g)$Eigenvector <- eigen_centrality(g)$vector
V(g)$MCC         <- mcc_score(g)
V(g)$Module      <- as.character(membership(cluster_louvain(g)))

node_tbl <- tibble(
  Gene        = V(g)$name,
  Direction   = V(g)$Direction,
  Degree      = V(g)$Degree,
  Betweenness = round(V(g)$Betweenness, 4),
  Closeness   = round(V(g)$Closeness, 4),
  Eigenvector = round(V(g)$Eigenvector, 4),
  MCC         = V(g)$MCC,
  Module      = V(g)$Module
) %>%
  mutate(Hub_Rank = rank(-Degree, ties.method = "min")) %>%
  arrange(Hub_Rank)

edge_tbl <- as_data_frame(g, what = "edges") %>%
  rename(Source = from, Target = to, Combined_Score = combined_score) %>%
  mutate(Interaction = "pp")

# ---- 6. Network statistics -------------------------------------------
net_stats <- tibble(
  Metric = c("Nodes", "Edges", "Average degree", "Network density",
             "Clustering coefficient", "Avg path length",
             "Connected components", "Modules (Louvain)", "Modularity"),
  Value  = c(vcount(g), ecount(g),
             round(mean(degree(g)), 2), round(edge_density(g), 4),
             round(transitivity(g, type = "global"), 4),
             round(mean_distance(g), 3),
             components(g)$no,
             length(unique(V(g)$module)),
             round(modularity(cluster_louvain(g)), 4))
)
print(net_stats, n = Inf)

write_csv(node_tbl,  "results/ppi/PPI_Node_Table.csv")
write_csv(edge_tbl,  "results/ppi/PPI_Edge_Table.csv")
write_csv(net_stats, "results/ppi/PPI_Network_Statistics.csv")
write_csv(head(node_tbl, PARAMS$hub_top_n), "results/ppi/PPI_Top_Hub_Genes.csv")

# ---- 7. Cytoscape export ---------------------------------------------
cat("\n=== Exporting for Cytoscape ===\n")

# (a) GraphML - retains all node/edge attributes, one-click import
write_graph(g, "results/ppi/PPI_Network_Cytoscape.graphml", format = "graphml")

# (b) SIF - simple interaction format
write_tsv(edge_tbl %>% select(Source, Interaction, Target),
          "results/ppi/PPI_Network.sif", col_names = FALSE)

# (c) Plain tables for File > Import > Network / Table from File
write_tsv(edge_tbl, "results/ppi/Cytoscape_Edges.tsv")
write_tsv(node_tbl, "results/ppi/Cytoscape_Nodes.tsv")

# ---- 8. Publication figure: main network ------------------------------
cat("\n=== Rendering figures ===\n")

tg <- as_tbl_graph(g) %>%
  mutate(label_flag = rank(-Degree, ties.method = "min") <= PARAMS$top_hubs_label)

set.seed(PARAMS$seed)
p_net <- ggraph(tg, layout = "fr") +
  geom_edge_link(aes(width = combined_score, alpha = combined_score),
                 colour = "grey65", show.legend = TRUE) +
  geom_node_point(aes(size = Degree, fill = Degree, shape = Direction),
                  colour = "white", stroke = 0.35) +
  geom_node_text(aes(label = ifelse(label_flag, name, "")),
                 repel = TRUE, size = 3.1, fontface = "bold",
                 max.overlaps = Inf, segment.size = 0.25,
                 segment.colour = "grey55", box.padding = 0.35) +
  scale_edge_width(range = c(0.15, 1.2), name = "Confidence") +
  scale_edge_alpha(range = c(0.25, 0.75), guide = "none") +
  scale_size_continuous(range = c(2, 11), name = "Degree") +
  scale_fill_gradientn(colours = c("#4575B4", "#91BFDB", "#FEE090", "#FC8D59", "#D73027"),
                       name = "Degree") +
  scale_shape_manual(values = c(Up = 21, Down = 24), name = "Regulation",
                     na.value = 22) +
  guides(fill = guide_colourbar(order = 1),
         size = guide_legend(order = 2, override.aes = list(fill = "grey40")),
         shape = guide_legend(order = 3, override.aes = list(size = 4, fill = "grey40"))) +
  labs(title = "Protein–Protein Interaction Network of the Core ZIKV Signature",
       subtitle = sprintf("STRING v%s | confidence \u2265 %.2f | %d nodes, %d edges",
                          PARAMS$string_version, PARAMS$score_threshold / 1000,
                          vcount(g), ecount(g))) +
  theme_graph(base_family = "sans") +
  theme(plot.title    = element_text(size = 15, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, colour = "grey30", hjust = 0.5),
        legend.position = "right",
        legend.title  = element_text(size = 10, face = "bold"),
        legend.text   = element_text(size = 9))

ggsave("plots/ppi/01_PPI_Network.png", p_net, width = 12, height = 10, dpi = 600, bg = "white")
ggsave("plots/ppi/01_PPI_Network.pdf", p_net, width = 12, height = 10, device = cairo_pdf)

# ---- 9. Module-coloured network --------------------------------------
mod_keep <- names(sort(table(V(g)$module), decreasing = TRUE))
tg_mod <- as_tbl_graph(g) %>%
  mutate(Module = factor(Module, levels = mod_keep),
         label_flag = rank(-Degree, ties.method = "min") <= PARAMS$top_hubs_label)

set.seed(PARAMS$seed)
p_mod <- ggraph(tg_mod, layout = "fr") +
  geom_edge_link(colour = "grey75", alpha = 0.4, width = 0.25) +
  geom_node_point(aes(size = Degree, fill = Module), shape = 21,
                  colour = "white", stroke = 0.35) +
  geom_node_text(aes(label = ifelse(label_flag, name, "")),
                 repel = TRUE, size = 3, fontface = "bold", max.overlaps = Inf) +
  scale_size_continuous(range = c(2, 10)) +
  scale_fill_brewer(palette = "Set2", name = "Module") +
  labs(title = "Functional Modules (Louvain Community Detection)",
       subtitle = sprintf("%d modules | modularity = %.3f",
                          length(unique(V(g)$module)),
                          modularity(cluster_louvain(g)))) +
  theme_graph(base_family = "sans") +
  theme(plot.title    = element_text(size = 15, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, colour = "grey30", hjust = 0.5))

ggsave("plots/ppi/02_PPI_Modules.png", p_mod, width = 12, height = 10, dpi = 600, bg = "white")
ggsave("plots/ppi/02_PPI_Modules.pdf", p_mod, width = 12, height = 10, device = cairo_pdf)

# ---- 10. Hub gene bar plot -------------------------------------------
hub_plot_df <- node_tbl %>%
  slice_head(n = PARAMS$hub_top_n) %>%
  mutate(Gene = fct_reorder(Gene, Degree))

p_hub <- ggplot(hub_plot_df, aes(x = Gene, y = Degree, fill = Degree)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = Degree), hjust = -0.25, size = 3.2, fontface = "bold") +
  coord_flip(clip = "off") +
  scale_fill_gradientn(colours = c("#91BFDB", "#FEE090", "#D73027"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(title = sprintf("Top %d Hub Genes by Node Degree", PARAMS$hub_top_n),
       x = NULL, y = "Node degree") +
  theme_minimal(base_size = 12) +
  theme(plot.title  = element_text(size = 15, face = "bold", hjust = 0.5),
        axis.text.y = element_text(face = "bold", colour = "black"),
        axis.text.x = element_text(colour = "black"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank())

ggsave("plots/ppi/03_Top_Hub_Genes.png", p_hub, width = 7, height = 8, dpi = 600, bg = "white")
ggsave("plots/ppi/03_Top_Hub_Genes.pdf", p_hub, width = 7, height = 8, device = cairo_pdf)

# ---- 11. Hub metric concordance heatmap -------------------------------
hub_long <- node_tbl %>%
  slice_head(n = PARAMS$hub_top_n) %>%
  select(Gene, Degree, Betweenness, Closeness, Eigenvector, MCC) %>%
  pivot_longer(-Gene, names_to = "Metric", values_to = "Value") %>%
  group_by(Metric) %>%
  mutate(Scaled = rescale(Value)) %>%
  ungroup() %>%
  mutate(Gene = factor(Gene, levels = rev(node_tbl$Gene[seq_len(PARAMS$hub_top_n)])),
         Metric = factor(Metric, levels = c("Degree", "MCC", "Betweenness",
                                            "Closeness", "Eigenvector")))

p_metric <- ggplot(hub_long, aes(Metric, Gene, fill = Scaled)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  scale_fill_gradientn(colours = c("#FFFFFF", "#FEE090", "#FC8D59", "#D73027"),
                       name = "Scaled\ncentrality") +
  labs(title = "Centrality Concordance Across Hub Metrics", x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title  = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text.y = element_text(face = "bold", colour = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),
        panel.grid  = element_blank())

ggsave("plots/ppi/04_Hub_Metric_Heatmap.png", p_metric, width = 6.5, height = 8, dpi = 600, bg = "white")
ggsave("plots/ppi/04_Hub_Metric_Heatmap.pdf", p_metric, width = 6.5, height = 8, device = cairo_pdf)

# ---- 12. Degree distribution (scale-free check) -----------------------
dd <- tibble(Degree = degree(g)) %>% count(Degree, name = "Frequency")
tryCatch({
  fit <- fit_power_law(degree(g)[degree(g) > 0])
  p_dd <- ggplot(dd, aes(Degree, Frequency)) +
    geom_point(size = 2.4, colour = "#2166AC", alpha = 0.8) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                colour = "#D73027", linewidth = 0.7, linetype = "dashed") +
    scale_x_log10() + scale_y_log10() +
    labs(title = "Degree Distribution",
         subtitle = sprintf("Power-law fit: alpha = %.2f, KS p = %.3f",
                            fit$alpha, fit$KS.p),
         x = "Node degree (log10)", y = "Frequency (log10)") +
    theme_minimal(base_size = 12) +
    theme(plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, colour = "grey30", hjust = 0.5),
          axis.text     = element_text(colour = "black"))
  ggsave("plots/ppi/05_Degree_Distribution.png", p_dd, width = 6.5, height = 5, dpi = 600, bg = "white")
}, error = function(e) {
  cat("Power law fit failed, skipping degree distribution plot.\n")
})

# ---- 13. Top-hub subnetwork (clean figure for main text) --------------
top_nodes <- node_tbl$Gene[seq_len(min(PARAMS$hub_top_n, nrow(node_tbl)))]
g_sub <- induced_subgraph(g, which(V(g)$name %in% top_nodes))
g_sub <- delete_vertices(g_sub, which(degree(g_sub) == 0))

if (vcount(g_sub) > 2) {
  set.seed(PARAMS$seed)
  p_sub <- ggraph(as_tbl_graph(g_sub), layout = "kk") +
    geom_edge_link(aes(width = combined_score), colour = "grey60", alpha = 0.55) +
    geom_node_point(aes(size = Degree, fill = Degree), shape = 21,
                    colour = "white", stroke = 0.5) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3.6,
                   fontface = "bold", max.overlaps = Inf) +
    scale_edge_width(range = c(0.3, 1.8), name = "Confidence") +
    scale_size_continuous(range = c(5, 16), name = "Degree") +
    scale_fill_gradientn(colours = c("#FEE090", "#FC8D59", "#D73027"), name = "Degree") +
    labs(title = sprintf("Core Hub Subnetwork (top %d genes)", length(top_nodes))) +
    theme_graph(base_family = "sans") +
    theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))

  ggsave("plots/ppi/06_Hub_Subnetwork.png", p_sub, width = 9, height = 8, dpi = 600, bg = "white")
  ggsave("plots/ppi/06_Hub_Subnetwork.pdf", p_sub, width = 9, height = 8, device = cairo_pdf)
}

# ---- 14. Session info -------------------------------------------------
cat("\n=== PPI Analysis Complete ===\n")
cat("Tables : results/ppi/\n")
cat("Figures: plots/ppi/\n")
cat("Cytoscape: results/ppi/PPI_Network_Cytoscape.graphml\n")
writeLines(capture.output(sessionInfo()), "results/ppi/sessionInfo_PPI.txt")
