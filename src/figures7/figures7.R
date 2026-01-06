#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Figures 7: Spatial plots for STAT6/lesion signatures in SIO2_56
# - Feature plots for Col1a1/Col1a2/Col3a1
# - Spatial maps for ctrllesion, lesion_stat6, ctrl_stat6
# Outputs are written to the SIO2_56 data directory.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(ggpubr)
  library(msigdbr)
  library(yaml)
})

set.seed(123)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

# Paths
script_path <- normalizePath(sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]))
script_dir <- dirname(script_path)
project_root <- normalizePath(file.path(script_dir, "../.."))
project_parent <- normalizePath(file.path(project_root, ".."))
config_path <- file.path(project_root, "config", "paths.local.yaml")
if (!file.exists(config_path)) {
  config_path <- file.path(project_root, "config", "paths.yaml")
}
cfg <- yaml::read_yaml(config_path)

data_dir <- cfg$raw_data$spatial_sio2_56 %||% file.path(project_parent, "spatial trans", "dataset", "SIO2_56")
fig_dir <- cfg$figures$figures7 %||% file.path(project_root, "figures", "figures7")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

out_exp <- file.path(fig_dir, "exp col1a2 col1a1 col3a1.pdf")
out_ctrllesion <- file.path(fig_dir, "ctrllesion spatial.pdf")
out_stat6_lesion <- file.path(fig_dir, "stat6 lesion.pdf")
out_stat6_ctrl <- file.path(fig_dir, "stat6 ctrl spatial.pdf")
out_scatter <- file.path(fig_dir, "cor_stat6_emt_glycolysis_scatter.pdf")
out_crossk <- file.path(fig_dir, "crossK_pcf_stat6_emt_glycolysis.pdf")

cat("Loading spatial data from:", data_dir, "\n")
counts <- Read10X_h5(file.path(data_dir, "filtered_feature_bc_matrix.h5"))
sce <- CreateSeuratObject(counts, assay = "Spatial")

# Load tissue positions
pos <- read.csv(file.path(data_dir, "spatial", "tissue_positions_list.csv"), stringsAsFactors = FALSE)
rownames(pos) <- pos$barcode
pos <- pos[colnames(sce), ]  # align
sce@meta.data$in_tissue <- pos$in_tissue
sce@meta.data$array_row <- pos$array_row
sce@meta.data$array_col <- pos$array_col
sce@meta.data$imagerow <- pos$pxl_row_in_fullres
sce@meta.data$imagecol <- pos$pxl_col_in_fullres

# SCTransform and clustering
sce <- SCTransform(sce, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
sce <- RunPCA(sce, assay = "SCT", verbose = FALSE)
sce <- FindNeighbors(sce, reduction = "pca", dims = 1:20, verbose = FALSE)
sce <- FindClusters(sce, resolution = 0.8, verbose = FALSE)
sce <- RunUMAP(sce, reduction = "pca", dims = 1:20)
sce@meta.data$seurat_clusters <- as.character(sce@meta.data$seurat_clusters)

# IFAT2 flag
sce@meta.data$IFAT2 <- ifelse(sce@assays$SCT@data["Sftpc", ] != 0 &
                                sce@assays$SCT@data["Epcam", ] != 0,
                              "AT2", "non-AT2")

# Gene sets: STAT6 positive, EMT, Glycolysis (convert human -> rat)
stat6_geneset_path <- cfg$raw_data$stat6_geneset %||% file.path(project_parent, "scRNAseq", "STAT6 geneset.csv")
stat6_df <- read.csv(stat6_geneset_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
stat6_pos <- stat6_df[["STAT6 positively regulated genes"]]
stat6_pos <- stat6_pos[!is.na(stat6_pos) & stat6_pos != ""]

convert_genes <- function(human_genes) {
  human_genes <- human_genes[!is.na(human_genes) & human_genes != ""]
  if (requireNamespace("homologene", quietly = TRUE)) {
    conv <- homologene::homologene(human_genes, inTax = 9606, outTax = 10116)
    mapping <- setNames(as.character(conv[,"10116"]), conv[,"9606"])
    rat_genes <- mapping[human_genes]
    rat_genes <- rat_genes[!is.na(rat_genes)]
    if (length(rat_genes) > 0) return(unique(rat_genes))
  }
  # Fallback: return human genes if conversion unavailable
  unique(human_genes)
}

m_df <- msigdbr(species = "Homo sapiens")
emt_set <- m_df %>%
  filter(gs_name == "GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_TO_MESENCHYMAL_TRANSITION") %>%
  pull(gene_symbol) %>% unique()
gly_set <- m_df %>%
  filter(gs_name == "HALLMARK_GLYCOLYSIS") %>%
  pull(gene_symbol) %>% unique()

stat6_genes <- convert_genes(stat6_pos)
emt_genes <- convert_genes(emt_set)
gly_genes <- convert_genes(gly_set)

sce <- AddModuleScore(sce,
                      features = list(stat6_genes, emt_genes, gly_genes),
                      name = "ModuleScore",
                      assay = "SCT",
                      nbin = 24,
                      ctrl = 100)
score_cols <- grep("^ModuleScore", colnames(sce@meta.data), value = TRUE)
colnames(sce@meta.data)[match(score_cols, colnames(sce@meta.data))] <- c("STAT6", "EMT", "Glycolysis")

sce@meta.data$STAT6_binary <- ifelse(sce@meta.data$STAT6 > median(sce@meta.data$STAT6, na.rm = TRUE), "High", "Low")
sce@meta.data$EMT_binary <- ifelse(sce@meta.data$EMT > median(sce@meta.data$EMT, na.rm = TRUE), "High", "Low")
sce@meta.data$Gly_binary <- ifelse(sce@meta.data$Glycolysis > median(sce@meta.data$Glycolysis, na.rm = TRUE), "High", "Low")

# ctrl vs lesion clusters
ctrl_clusters <- c("0", "2", "3", "4", "5", "7", "8")
lesion_clusters <- c("1", "6", "9")
sce@meta.data$ctrllesion <- ifelse(sce@meta.data$seurat_clusters %in% ctrl_clusters, "ctrl", "lesion")

sce@meta.data$lesion_stat6 <- ifelse(sce@meta.data$seurat_clusters %in% ctrl_clusters, "ctrl",
                                     ifelse(sce@meta.data$IFAT2 == "AT2" &
                                              sce@meta.data$STAT6_binary == "High", "High", "Low"))
sce@meta.data$ctrl_stat6 <- ifelse(sce@meta.data$seurat_clusters %in% lesion_clusters, "lesion",
                                   ifelse(sce@meta.data$IFAT2 == "AT2" &
                                            sce@meta.data$STAT6_binary == "High", "High", "Low"))

# Coordinates for ggplot
meta <- sce@meta.data

# Feature plots: Col1a1/Col1a2/Col3a1
feat_plot <- function(gene) {
  df <- FetchData(sce, vars = c(gene, "imagerow", "imagecol"), slot = "data", assay = "SCT")
  colnames(df)[1] <- "expr"
  ggplot(df, aes(x = imagerow, y = imagecol, color = expr)) +
    geom_point(size = 2.0, alpha = 0.8) +
    scale_color_gradientn(colors = c("white", "mistyrose", "salmon", "darkred")) +
    coord_fixed() +
    theme_void() +
    ggtitle(gene) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}
p_col1a1 <- feat_plot("Col1a1")
p_col1a2 <- feat_plot("Col1a2")
p_col3a1 <- feat_plot("Col3a1")

ggsave(out_exp, (p_col1a1 | p_col1a2) / p_col3a1, width = 12, height = 12, dpi = 300)
cat("Saved:", out_exp, "\n")

# ctrllesion spatial map
p_ctrllesion <- ggplot(meta, aes(x = imagerow, y = imagecol, color = ctrllesion)) +
  geom_point(size = 2.15, alpha = 0.8) +
  scale_color_manual(values = c("ctrl" = "#66C2A5", "lesion" = "#CD5C5C")) +
  coord_fixed() +
  theme_void() +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 1))
ggsave(out_ctrllesion, p_ctrllesion, width = 8, height = 8, dpi = 300)
cat("Saved:", out_ctrllesion, "\n")

# lesion_stat6 spatial map
meta$lesion_stat6 <- factor(meta$lesion_stat6, levels = c("ctrl", "Low", "High"))
p_lesion_stat6 <- ggplot(meta, aes(x = imagerow, y = imagecol, color = lesion_stat6)) +
  geom_point(size = 2.15, alpha = 0.8) +
  scale_color_manual(values = c("ctrl" = "grey90", "Low" = "purple", "High" = "#CD5C5C")) +
  coord_fixed() +
  theme_void() +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 1))
ggsave(out_stat6_lesion, p_lesion_stat6, width = 8, height = 8, dpi = 300)
cat("Saved:", out_stat6_lesion, "\n")

# ctrl_stat6 spatial map
meta$ctrl_stat6 <- factor(meta$ctrl_stat6, levels = c("lesion", "Low", "High"))
p_ctrl_stat6 <- ggplot(meta, aes(x = imagerow, y = imagecol, color = ctrl_stat6)) +
  geom_point(size = 2.15, alpha = 0.8) +
  scale_color_manual(values = c("lesion" = "purple", "Low" = "grey90", "High" = "#66C2A5")) +
  coord_fixed() +
  theme_void() +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 1))
ggsave(out_stat6_ctrl, p_ctrl_stat6, width = 8, height = 8, dpi = 300)
cat("Saved:", out_stat6_ctrl, "\n")

# ------------------------------------------------------------------------------
# Correlation scatter: STAT6 vs EMT and STAT6 vs Glycolysis (AT2 cells)
# ------------------------------------------------------------------------------
at2_cells <- rownames(meta)[meta$IFAT2 == "AT2"]
if (length(at2_cells) >= 3) {
  expr_df <- FetchData(sce, vars = c("STAT6", "EMT", "Glycolysis"), cells = at2_cells)
  expr_df <- expr_df[complete.cases(expr_df), ]
  if (nrow(expr_df) >= 3) {
    p_scatter_emt <- ggplot(expr_df, aes(x = STAT6, y = EMT)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", se = TRUE, color = "blue") +
      stat_cor(method = "pearson",
               label.x = min(expr_df$STAT6),
               label.y = max(expr_df$EMT)) +
      labs(title = "STAT6 vs EMT (AT2)",
           x = "STAT6 module score",
           y = "EMT module score") +
      theme_bw()
    p_scatter_gly <- ggplot(expr_df, aes(x = STAT6, y = Glycolysis)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", se = TRUE, color = "red") +
      stat_cor(method = "pearson",
               label.x = min(expr_df$STAT6),
               label.y = max(expr_df$Glycolysis)) +
      labs(title = "STAT6 vs Glycolysis (AT2)",
           x = "STAT6 module score",
           y = "Glycolysis module score") +
      theme_bw()
    ggsave(out_scatter, p_scatter_emt / p_scatter_gly, width = 7, height = 8, dpi = 300)
    cat("Saved:", out_scatter, "\n")
  } else {
    cat("Insufficient AT2 cells with non-NA STAT6/EMT/Glycolysis for scatter plots.\n")
  }
} else {
  cat("No AT2 cells found; skip scatter plots.\n")
}

# ------------------------------------------------------------------------------
# Spatstat Cross-K / PCF for STAT6 vs EMT and STAT6 vs Glycolysis (all spots)
# ------------------------------------------------------------------------------
meta$STAT6_binary <- factor(meta$STAT6_binary, levels = c("Low", "High"))
meta$EMT_binary <- factor(meta$EMT_binary, levels = c("Low", "High"))
meta$Gly_binary <- factor(meta$Gly_binary, levels = c("Low", "High"))

if (requireNamespace("spatstat.geom", quietly = TRUE) &&
    requireNamespace("spatstat.core", quietly = TRUE)) {

  win <- spatstat.geom::owin(xrange = range(meta$imagerow, na.rm = TRUE),
                             yrange = range(meta$imagecol, na.rm = TRUE))

  compute_cross <- function(df, mark1, mark2, label1, label2) {
    df_sub <- df[df[[mark1]] == "High" | df[[mark2]] == "High", ]
    if (nrow(df_sub) < 10) return(NULL)
    df_sub$mark <- ifelse(df_sub[[mark1]] == "High", label1, label2)
    df_sub$mark <- factor(df_sub$mark, levels = c(label1, label2))
    spatstat.geom::ppp(x = df_sub$imagerow,
                       y = df_sub$imagecol,
                       window = win,
                       marks = df_sub$mark)
  }

  pp_emt <- compute_cross(meta, "STAT6_binary", "EMT_binary", "STAT6", "EMT")
  pp_gly <- compute_cross(meta, "STAT6_binary", "Gly_binary", "STAT6", "Gly")

  pdf(out_crossk, width = 10, height = 8)
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  if (!is.null(pp_emt)) {
    K_res <- spatstat.core::Kcross(pp_emt, i = "STAT6", j = "EMT", correction = "Ripley")
    pcf_res <- spatstat.core::pcf(pp_emt, i = "STAT6", j = "EMT", correction = "Ripley")
    plot(K_res, main = "Cross K: STAT6 vs EMT", xlab = "r", ylab = "K(r)")
    plot(pcf_res, main = "Cross PCF: STAT6 vs EMT", xlab = "r", ylab = "pcf(r)")
  } else {
    plot.new(); title("Insufficient points for STAT6 vs EMT")
    plot.new(); title("Insufficient points for STAT6 vs EMT")
  }
  if (!is.null(pp_gly)) {
    K_res <- spatstat.core::Kcross(pp_gly, i = "STAT6", j = "Gly", correction = "Ripley")
    pcf_res <- spatstat.core::pcf(pp_gly, i = "STAT6", j = "Gly", correction = "Ripley")
    plot(K_res, main = "Cross K: STAT6 vs Glycolysis", xlab = "r", ylab = "K(r)")
    plot(pcf_res, main = "Cross PCF: STAT6 vs Glycolysis", xlab = "r", ylab = "pcf(r)")
  } else {
    plot.new(); title("Insufficient points for STAT6 vs Glycolysis")
    plot.new(); title("Insufficient points for STAT6 vs Glycolysis")
  }
  dev.off()
  cat("Saved:", out_crossk, "\n")
} else {
  cat("spatstat.geom/core not available; skipping Cross-K/PCF plot. Install spatstat.geom and spatstat.core if needed.\n")
}

cat("\nDone. Outputs written to:", fig_dir, "\n")
