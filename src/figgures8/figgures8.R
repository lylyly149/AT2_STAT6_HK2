#!/usr/bin/env Rscript
# ==============================================================================
# Comment 7 Figure 1: scFEA glycolysis flux analysis in AT2 cells (IPF vs NOM)
# ==============================================================================
#
# This script generates figures for scFEA metabolic flux analysis:
#   - Figure 1A: Volcano plot showing glycolysis module flux changes (Cohen's d)
#   - Figure 1B: Boxplot of glycolysis scores (NOM vs IPF)
#   - Figure 1C: Dot plot of glycolysis flux modules M1-M6
#
# Input: scFEA flux results from AT2 cells
# Runtime environment: conda activate r_4_3 (R 4.3+)
# ==============================================================================

set.seed(42)

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(ggrepel)
  library(Seurat)
  library(optparse)
  library(yaml)
})

# Get script directory and project root
args <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
script_dir <- dirname(script_path)
project_root <- normalizePath(file.path(script_dir, "../.."))

# Load configuration
config <- yaml::read_yaml(file.path(project_root, "config/paths.yaml"))

# Set working directory to project root
setwd(project_root)

# Create output directories
fig_dir <- "figures/figgures8"
processed_dir <- "data/processed"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n====================\n")
cat("Comment 7 Figure 1: AT2 glycolysis flux analysis\n")
cat("====================\n\n")

# ==============================================================================
# Parameters
# ==============================================================================
group_a <- "NOM"  # Reference group
group_b <- "IPF"  # Comparison group
group_column <- "nomipf"
padj_cutoff <- config$parameters$stats$padj_cutoff

# ==============================================================================
# Step 1: Load scFEA flux data and metadata
# ==============================================================================
cat("=== Loading scFEA flux data ===\n")

# Read data (adjust paths based on config)
flux_file <- config$processed_data$scfea_flux_at2
metadata_file <- config$processed_data$scfea_input_at2_metadata
module_info_file <- config$raw_data$scfea_module_info

read_table <- function(path) read_csv(path, show_col_types = FALSE)

flux <- read_table(flux_file)
meta <- read_table(metadata_file)
mod_info <- read_table(module_info_file)

# Prepare cell_id column
first_col <- names(flux)[1]
if (first_col %in% c("cell_id", "Unnamed: 0")) {
  names(flux)[1] <- "cell_id"
} else {
  warning("Treating first column '", first_col, "' as cell_id.")
  names(flux)[1] <- "cell_id"
}

# Filter cells
flux <- flux %>% semi_join(meta, by = "cell_id")
meta <- meta %>%
  filter(cell_id %in% flux$cell_id) %>%
  mutate(group = .data[[group_column]]) %>%
  filter(group %in% c(group_a, group_b))
flux <- flux %>% filter(cell_id %in% meta$cell_id)

cat(sprintf("Retained %d cells for analysis\n", nrow(meta)))

# ==============================================================================
# Step 2: Identify glycolysis/TCA modules
# ==============================================================================
cat("\n=== Processing module information ===\n")

mod_info <- mod_info %>%
  mutate(label_full = paste0(Module_id, ": ", Compound_IN_name, "\u2192", Compound_OUT_name))

modules_in_flux <- intersect(mod_info$Module_id, names(flux))
cat(sprintf("Found %d modules in flux data\n", length(modules_in_flux)))

flux_gly <- flux %>%
  select(cell_id, all_of(modules_in_flux)) %>%
  inner_join(meta %>% select(cell_id, group, orig.ident), by = "cell_id")

# ==============================================================================
# Step 3: Calculate differential flux statistics
# ==============================================================================
cat("\n=== Calculating differential flux ===\n")

calc_stats <- function(df, module) {
  vals <- df[[module]]
  g <- df$group
  a_vals <- vals[g == group_a]
  b_vals <- vals[g == group_b]
  p <- tryCatch(wilcox.test(b_vals, a_vals)$p.value, error = function(e) NA_real_)
  mean_a <- mean(a_vals, na.rm = TRUE)
  mean_b <- mean(b_vals, na.rm = TRUE)
  min_mean <- min(mean_a, mean_b, na.rm = TRUE)
  shift <- if (is.finite(min_mean) && min_mean <= 0) abs(min_mean) + 1e-6 else 1e-6
  logfc <- log2(mean_b + shift) - log2(mean_a + shift)  # IPF - NOM
  pooled_std <- sqrt(((length(a_vals) - 1) * stats::var(a_vals) + (length(b_vals) - 1) * stats::var(b_vals)) /
                       (length(a_vals) + length(b_vals) - 2))
  cohen_d <- (mean_b - mean_a) / pooled_std
  tibble(
    Module_id = module,
    mean_group_a = mean_a,
    mean_group_b = mean_b,
    log2FC = logfc,
    cohen_d = cohen_d,
    pvalue = p
  )
}

stats_tbl <- bind_rows(lapply(modules_in_flux, function(m) calc_stats(flux_gly, m))) %>%
  mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  left_join(mod_info %>% select(Module_id, SM_anno, label_full), by = "Module_id") %>%
  arrange(padj)

# Save statistical results
write_csv(stats_tbl, file.path(processed_dir, "scfea_flux_diff_at2_all.csv"))
write_csv(filter(stats_tbl, !is.na(padj) & padj < padj_cutoff),
          file.path(processed_dir, "scfea_flux_diff_at2_sig.csv"))

# ==============================================================================
# Figure 1A: Volcano plot (Cohen's d)
# ==============================================================================
cat("\n=== Generating Figure 1A: Volcano plot ===\n")

gly_only <- c("M_1", "M_2", "M_3", "M_4", "M_5", "M_6")
label_lookup <- mod_info$label_full
names(label_lookup) <- mod_info$Module_id

volcano_df <- stats_tbl %>%
  mutate(
    padj_plot = ifelse(!is.na(padj) & padj <= 0, .Machine$double.eps, padj),
    neglog10_padj = -log10(padj_plot),
    label_plot = ifelse(Module_id %in% gly_only, label_lookup[Module_id], NA_character_),
    label_flag = !is.na(label_plot)
  )

ylim_top <- volcano_df %>%
  filter(is.finite(neglog10_padj)) %>%
  summarise(top = max(neglog10_padj, na.rm = TRUE)) %>%
  pull(top)
if (!is.finite(ylim_top)) ylim_top <- 10
ylim_top <- ylim_top * 1.1

p_volcano <- ggplot(volcano_df, aes(x = cohen_d, y = neglog10_padj)) +
  geom_point(data = subset(volcano_df, !label_flag), color = "grey70", size = 2.5, alpha = 0.85) +
  geom_point(data = subset(volcano_df, label_flag), color = "#d73027", size = 2.6, alpha = 0.95) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    title = paste0("Glycolysis/TCA modules: ", group_b, " vs ", group_a),
    x = "Cohen's d (IPF - NOM)",
    y = "-log10(padj)"
  ) +
  theme_bw() +
  coord_cartesian(ylim = c(0, ylim_top))

if (any(!is.na(volcano_df$label_plot))) {
  p_volcano <- p_volcano +
    geom_text_repel(
      data = subset(volcano_df, !is.na(label_plot)),
      aes(label = label_plot),
      max.overlaps = Inf,
      min.segment.length = 0,
      box.padding = 0.4,
      segment.color = "grey40",
      size = 3
    )
}

ggsave(file.path(fig_dir, "volcano_glycolysis_modules_cohend.png"), p_volcano, width = 7, height = 5.5, dpi = 300)
ggsave(file.path(fig_dir, "volcano_glycolysis_modules_cohend.pdf"), p_volcano, width = 7, height = 5.5)

# ==============================================================================
# Figure 1B: Glycolysis score boxplot
# ==============================================================================
cat("\n=== Generating Figure 1B: Glycolysis score boxplot ===\n")

# Calculate glycolysis score (sum of M1-M6 flux)
glycolysis_score <- flux_gly %>%
  mutate(glycolysis_score = rowSums(select(., all_of(gly_only)), na.rm = TRUE)) %>%
  select(cell_id, group, glycolysis_score)

p_boxplot <- ggplot(glycolysis_score, aes(x = group, y = glycolysis_score, fill = group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.3) +
  scale_fill_manual(values = c("NOM" = "#5aa4ae", "IPF" = "#ea5514")) +
  labs(
    title = "Glycolysis flux score (M1-M6)",
    x = "Group",
    y = "Glycolysis score"
  ) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "glycolysis_score_boxplot.png"), p_boxplot, width = 5, height = 4, dpi = 300)
ggsave(file.path(fig_dir, "glycolysis_score_boxplot.pdf"), p_boxplot, width = 5, height = 4)

# ==============================================================================
# Figure 1C: Dot plot of glycolysis modules M1-M6
# ==============================================================================
cat("\n=== Generating Figure 1C: Dot plot of M1-M6 ===\n")

if (all(gly_only %in% modules_in_flux)) {
  # Build Seurat object on flux values
  mat_for_seurat <- t(as.matrix(flux_gly %>% select(all_of(modules_in_flux))))
  colnames(mat_for_seurat) <- flux_gly$cell_id
  meta_df <- flux_gly %>% select(cell_id, group) %>% as.data.frame()
  rownames(meta_df) <- meta_df$cell_id
  seu <- CreateSeuratObject(counts = mat_for_seurat, assay = "FLUX")
  DefaultAssay(seu) <- "FLUX"
  seu$group <- meta_df[colnames(seu), "group", drop = TRUE]
  Idents(seu) <- factor(seu$group, levels = c(group_a, group_b))

  module_labels <- mod_info %>%
    filter(Module_id %in% gly_only) %>%
    mutate(label = paste0(Module_id, ": ", Compound_IN_name, "\u2192", Compound_OUT_name)) %>%
    select(Module_id, label)

  features_plot <- gsub("_", "-", gly_only)
  label_vec <- setNames(module_labels$label, features_plot)

  dot <- DotPlot(seu, features = features_plot, cols = c("RdYlBu")) +
    geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
    guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"))) +
    scale_x_discrete(labels = label_vec) +
    labs(title = "Glycolysis flux modules (M1-M6)", x = "Module", y = "Group") +
    theme_bw() +
    theme(
      panel.grid.major = element_line(color = "grey85"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    coord_flip()

  ggsave(file.path(fig_dir, "dotplot_glycolysis_M1_M6.png"), dot, width = 5, height = 8, dpi = 300)
  ggsave(file.path(fig_dir, "dotplot_glycolysis_M1_M6.pdf"), dot, width = 5, height = 8)
}

cat("\n=== Comment 7 Figure 1 analysis complete ===\n")
cat("All figures saved to:", fig_dir, "\n")
