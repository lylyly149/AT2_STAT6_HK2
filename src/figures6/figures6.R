#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Figures 6: Monocle pseudotime heatmap and STAT6-positive alluvial plot
# - Reproduces the Monocle pseudotime heatmap (3 clusters) from scRNA seq code.R
# - Exports cluster assignments as state_DEGs_pseudotime.csv
# - Builds an alluvial plot linking STAT6 positively regulated genes to
#   their assigned heatmap cluster (or "not variable" if absent)
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle)
  library(Matrix)
  library(pheatmap)
  library(ggplot2)
  library(ggalluvial)
  library(dplyr)
  library(yaml)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

# ------------------------------------------------------------------------------
# Paths and configuration
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
script_dir <- dirname(script_path)
project_root <- normalizePath(file.path(script_dir, "../.."))

config_path <- file.path(project_root, "config", "paths.local.yaml")
if (!file.exists(config_path)) {
  config_path <- file.path(project_root, "config", "paths.yaml")
}
cfg <- yaml::read_yaml(config_path)

set.seed(cfg$parameters$seed %||% 123)

sce_at2_path <- cfg$raw_data$seurat_at2 %||% stop("config raw_data$seurat_at2 not set")
stat6_geneset_path <- cfg$raw_data$stat6_geneset %||% stop("config raw_data$stat6_geneset not set")

fig_dir <- file.path(project_root, "figures", "figures6")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading AT2 Seurat object from:", sce_at2_path, "\n")
sce <- readRDS(sce_at2_path)

if (is.null(sce@meta.data$nomipf)) {
  if ("orig.ident" %in% colnames(sce@meta.data)) {
    sce@meta.data$nomipf <- ifelse(substr(sce@meta.data$orig.ident, 1, 3) == "IPF", "IPF", "NOM")
  } else {
    stop("nomipf not found and orig.ident unavailable to derive IPF/NOM status.")
  }
}

# ------------------------------------------------------------------------------
# Build Monocle CellDataSet (matching scRNA seq code.R logic)
# ------------------------------------------------------------------------------
counts_mat <- sce@assays$RNA@counts
if (!inherits(counts_mat, "dgCMatrix")) {
  counts_mat <- as(counts_mat, "dgCMatrix")
}
pd <- new("AnnotatedDataFrame", data = sce@meta.data)
fd <- new("AnnotatedDataFrame",
          data = data.frame(gene_short_name = rownames(counts_mat),
                            row.names = rownames(counts_mat)))

monocle_cds <- newCellDataSet(counts_mat,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.01,
                              expressionFamily = negbinomial.size())

Mono.cds <- estimateSizeFactors(monocle_cds)
Mono.cds <- estimateDispersions(Mono.cds)
disp_table <- dispersionTable(Mono.cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
Mono.cds <- setOrderingFilter(Mono.cds, unsup_clustering_genes$gene_id)

Mono.cds <- reduceDimension(Mono.cds, max_components = 2, method = "DDRTree")
Mono.cds <- orderCells(Mono.cds)
if ("State" %in% colnames(pData(Mono.cds)) && any(pData(Mono.cds)$State == 6)) {
  Mono.cds <- orderCells(Mono.cds, root_state = 6)
} else {
  cat("State 6 not found; keeping default root state.\n")
}

# ------------------------------------------------------------------------------
# Differential expression by State and pseudotime heatmap
# ------------------------------------------------------------------------------
cat("Running differentialGeneTest by State...\n")
states_de <- differentialGeneTest(Mono.cds, fullModelFormulaStr = "~State")
states_de <- states_de[order(states_de$qval), ]

state_de_genes <- rownames(states_de)
state_de_genes <- state_de_genes[state_de_genes %in% unsup_clustering_genes$gene_id]

if (length(state_de_genes) == 0) {
  stop("No state DE genes passed the mean_expression filter; cannot plot heatmap.")
}

heatmap_file <- file.path(fig_dir, "figure6_pseudotime_heatmap.pdf")
cat("Plotting pseudotime heatmap to:", heatmap_file, "\n")
pdf(heatmap_file, width = 8, height = 10)
heatmap_obj <- plot_pseudotime_heatmap(Mono.cds[state_de_genes, ],
                                       num_clusters = 3,
                                       cores = 1,
                                       show_rownames = FALSE,
                                       return_heatmap = TRUE)
dev.off()

# Extract clusters from heatmap tree
if (!is.null(heatmap_obj$tree_row)) {
  clusters <- cutree(heatmap_obj$tree_row, k = 3)
  clustering <- data.frame(
    Gene = names(clusters),
    Cluster = as.character(clusters),
    stringsAsFactors = FALSE
  )
} else {
  stop("Heatmap object does not contain tree_row for clustering.")
}

cluster_csv <- file.path(fig_dir, "state_DEGs_pseudotime.csv")
write.csv(clustering, cluster_csv, row.names = FALSE)
cat("Saved cluster assignments to:", cluster_csv, "\n")

# ------------------------------------------------------------------------------
# Alluvial plot: STAT6 positively regulated genes -> heatmap cluster
# ------------------------------------------------------------------------------
cat("Building alluvial plot using STAT6 positively regulated genes from:",
    stat6_geneset_path, "\n")
stat6_sets <- read.csv(stat6_geneset_path, header = TRUE, check.names = FALSE,
                       stringsAsFactors = FALSE)
if (!("STAT6 positively regulated genes" %in% colnames(stat6_sets))) {
  stop("STAT6 positively regulated genes column not found in geneset file.")
}
stat6_pos <- stat6_sets[["STAT6 positively regulated genes"]]
stat6_pos <- stat6_pos[!is.na(stat6_pos) & stat6_pos != ""]

cluster_lookup <- setNames(clustering$Cluster, clustering$Gene)
cluster_assign <- ifelse(stat6_pos %in% names(cluster_lookup),
                         cluster_lookup[stat6_pos],
                         "not variable")

alluvial_df <- data.frame(
  STAT6_Positive_Gene = stat6_pos,
  Cluster = factor(cluster_assign, levels = c("1", "2", "3", "not variable")),
  stringsAsFactors = FALSE
)

alluvial_file <- file.path(fig_dir, "figure6_alluvial_STAT6_positive_clusters.pdf")
cat("Plotting alluvial diagram to:", alluvial_file, "\n")

alluvial_plot <- ggplot(alluvial_df,
                        aes(axis1 = "STAT6 positively regulated genes",
                            axis2 = Cluster,
                            y = 1)) +
  geom_alluvium(aes(fill = Cluster), width = 1/12, knot.pos = 0.35, alpha = 0.8) +
  geom_stratum(width = 1/8, fill = "grey80", color = "grey40") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("STAT6 positively regulated genes", "Pseudotime cluster"),
                   expand = c(.1, .1)) +
  labs(fill = "Cluster",
       y = "Genes",
       title = "STAT6 positively regulated genes across pseudotime clusters") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 11),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13))

ggsave(alluvial_file, alluvial_plot, width = 8, height = 6, dpi = 300)
cat("Saved alluvial plot to:", alluvial_file, "\n")

cat("\nDone. Outputs written to:", fig_dir, "\n")
