#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Figure 2D: Differential expression of AT2 cells (IPF vs NOM)
# - Load AT2 subset (sce_at2.rds from config)
# - Ensure IPF/NOM group labels
# - Run FindMarkers with contrast IPF - NOM
# - Save DEG table to CSV
# - Run GO BP enrichment on significant genes (padj < 0.05 & |logFC| > 0.5) with no GO cutoffs
# - Save GO table to CSV
# ------------------------------------------------------------------------------

set.seed(123)

suppressPackageStartupMessages({
  library(Seurat)
  library(yaml)
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

# Project root and config
args <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
script_dir <- dirname(script_path)
project_root <- normalizePath(file.path(script_dir, "../.."))

config_path <- file.path(project_root, "config", "paths.local.yaml")
if (!file.exists(config_path)) {
  config_path <- file.path(project_root, "config", "paths.yaml")
}
cfg <- yaml::read_yaml(config_path)

# Paths
sce_at2_path <- cfg$raw_data$seurat_at2 %||% stop("config raw_data$seurat_at2 not set")
out_dir <- file.path(project_root, "figures", "figures2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_csv <- file.path(out_dir, "figures2d_deg_ipf_vs_nom_at2.csv")
go_csv <- file.path(out_dir, "figures2d_go_bp_ipf_vs_nom_at2.csv")

cat("Loading AT2 Seurat object from:", sce_at2_path, "\n")
sce_at2 <- readRDS(sce_at2_path)

# Handle multi-layer assays if present
if (inherits(sce_at2[["RNA"]], "Assay5")) {
  cat("Joining multi-layer assay into single counts layer for Seurat operations\n")
  sce_at2 <- JoinLayers(sce_at2)
}

# Ensure group labels
if (!"nomipf" %in% colnames(sce_at2@meta.data)) {
  if (!"orig.ident" %in% colnames(sce_at2@meta.data)) {
    stop("Neither nomipf nor orig.ident found in meta.data; cannot assign IPF/NOM labels.")
  }
  sce_at2@meta.data[["nomipf"]] <- ifelse(substr(sce_at2@meta.data[["orig.ident"]], 1, 3) == "IPF", "IPF", "NOM")
}

Idents(sce_at2) <- sce_at2@meta.data$nomipf
groups <- levels(Idents(sce_at2))

if (!all(c("IPF", "NOM") %in% groups)) {
  stop("Both IPF and NOM groups are required for FindMarkers; found: ", paste(groups, collapse = ", "))
}

cat("Running FindMarkers (IPF - NOM)...\n")
deg <- FindMarkers(
  sce_at2,
  ident.1 = "IPF",
  ident.2 = "NOM",
  assay = "RNA",
  test.use = "wilcox"
)

deg <- deg %>%
  mutate(gene = rownames(.)) %>%
  dplyr::select(gene, dplyr::everything())

write.csv(deg, out_csv, row.names = FALSE)
cat("Saved DEG table to:", out_csv, "\n")

# ------------------------------------------------------------------------------
# GO enrichment on significant genes (padj < 0.05 & |logFC| > 0.5)
# ------------------------------------------------------------------------------
sig_genes <- deg %>%
  filter(!is.na(p_val_adj), p_val_adj < 0.05, abs(avg_log2FC) > 0.5) %>%
  pull(gene) %>%
  unique()

cat(sprintf("Significant genes for GO (padj<0.05 & |logFC|>0.5): %d\n", length(sig_genes)))

if (length(sig_genes) == 0) {
  stop("No significant genes; GO enrichment cannot be performed.")
}

all_genes <- deg$gene %>% unique()

universe_entrez <- mapIds(org.Hs.eg.db, keys = all_genes, column = "ENTREZID",
                          keytype = "SYMBOL", multiVals = "first") %>%
  na.omit() %>%
  unique()

sig_entrez <- mapIds(org.Hs.eg.db, keys = sig_genes, column = "ENTREZID",
                     keytype = "SYMBOL", multiVals = "first") %>%
  na.omit() %>%
  unique()

cat(sprintf("Mapped %d/%d significant genes to Entrez IDs\n",
            length(sig_entrez), length(sig_genes)))

if (length(sig_entrez) == 0) {
  stop("No significant genes could be mapped to Entrez IDs; GO enrichment cannot be performed.")
}

ego <- suppressMessages(enrichGO(
  gene = sig_entrez,
  universe = universe_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = TRUE
))

ego_df <- as.data.frame(ego)

if (is.null(ego_df) || nrow(ego_df) == 0) {
  stop("enrichGO returned no results.")
}

write.csv(ego_df, go_csv, row.names = FALSE)
cat("Saved GO BP enrichment to:", go_csv, "\n")
