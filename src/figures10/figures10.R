#!/usr/bin/env Rscript
# ==============================================================================
# Comment 7 Figure 3: Metabolomics analysis of STAT6 KO vs WT mice
# ==============================================================================
#
# This script generates figures for targeted metabolomics analysis:
#   - Figure 3A: Heatmap of glycolysis metabolites across treatment groups
#   - Figure 3B: Bubble plots of pathway enrichment for each comparison
#           (KO_vs_WT_BLM, KO_vs_WT_CS, KO_vs_WT_CTRL)
#   - Figure 3C: Bar plots of key metabolites (glucose-6-phosphate, pyruvate)
#
# Input: Targeted metabolomics data from STAT6 KO and WT mice
# Runtime environment: conda activate r_4_3 (R 4.3+)
# ==============================================================================

set.seed(42)

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(pheatmap)
  library(ggplot2)
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
fig_dir <- "figures/figures10"
processed_dir <- "data/processed"
enrichment_dir <- file.path(fig_dir, "enrichment")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(enrichment_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n====================\n")
cat("Comment 7 Figure 3: STAT6 KO metabolomics analysis\n")
cat("====================\n\n")

# ==============================================================================
# Step 1: Load metabolomics data
# ==============================================================================
cat("=== Loading metabolomics data ===\n")

metab_raw_path <- config$raw_data$metabolomics_raw
sample_info_path <- config$raw_data$metabolomics_sample_info

metab_raw <- read_csv(metab_raw_path, locale = locale(encoding = "GBK"), show_col_types = FALSE)
sample_info0 <- read_csv(sample_info_path, show_col_types = FALSE)

# Validate required columns
if (!"Metabolite" %in% names(metab_raw)) {
  stop("Input metabolomics matrix must include a 'Metabolite' column.")
}

analysis_col <- "Sample Analysis Name"
group_col <- "Group Name"
if (!(analysis_col %in% names(sample_info0) && group_col %in% names(sample_info0))) {
  stop("Sample info file must include 'Sample Analysis Name' and 'Group Name'.")
}

# Match samples
sample_cols <- intersect(colnames(metab_raw), sample_info0[[analysis_col]])
if (length(sample_cols) == 0) {
  stop("No overlapping sample columns between metabolomics table and metadata.")
}

sample_info <- sample_info0 %>%
  rename(SampleAnalysisName = all_of(analysis_col), Group = all_of(group_col)) %>%
  filter(SampleAnalysisName %in% sample_cols) %>%
  mutate(sampleID = SampleAnalysisName) %>%
  separate(Group, into = c("Genotype", "Treatment"), sep = "_", remove = FALSE, fill = "right")
sample_info <- sample_info[match(sample_cols, sample_info$sampleID), ]

cat(sprintf("Detected %d samples and %d metabolites\n", nrow(sample_info), nrow(metab_raw)))

# ==============================================================================
# Step 2: Data preprocessing
# ==============================================================================
cat("\n=== Data preprocessing ===\n")

# Handle duplicate metabolite names
metab_ids <- metab_raw$Metabolite
if (anyDuplicated(metab_ids)) {
  metab_ids <- make.unique(metab_ids)
  cat("Non-unique metabolite names detected; applying make.unique()\n")
}
metabolite_base <- metab_raw$Metabolite
names(metabolite_base) <- metab_ids
metab_raw$Metabolite <- metab_ids

# Extract expression matrix
expr_raw <- metab_raw %>% select(all_of(sample_cols)) %>% as.matrix()
rownames(expr_raw) <- metab_ids

# Impute missing values (half-minimum)
expr_imputed <- expr_raw
for (i in seq_len(nrow(expr_imputed))) {
  values <- expr_imputed[i, ]
  if (any(is.na(values))) {
    min_val <- suppressWarnings(min(values, na.rm = TRUE))
    if (is.finite(min_val)) {
      values[is.na(values)] <- min_val * 0.5
    } else {
      values[is.na(values)] <- 0
    }
  }
  expr_imputed[i, ] <- values
}

# Log2 transformation and Z-score normalization
expr_log2 <- log2(expr_imputed + 1)
expr_z <- t(scale(t(expr_log2)))
expr_z[is.na(expr_z)] <- 0

# ==============================================================================
# Figure 3A: Heatmap of glycolysis metabolites
# ==============================================================================
cat("\n=== Generating Figure 3A: Metabolite heatmap ===\n")

# Define glycolysis-related metabolites for labeling
glycolysis_hits <- c(
  "glucose", "glucose 1-phosphate", "glucose 6-phosphate",
  "fructose", "fructose 6-phosphate", "glyceraldehyde",
  "glucosamine 6-phosphate", "galactose 1-phosphate", "ribose 5-phosphate",
  "pyruvic acid", "pyruvate", "lactate", "lactic acid"
)

# Prepare annotation
ann_col <- sample_info %>%
  select(sampleID, Genotype, Treatment) %>%
  column_to_rownames("sampleID")

# Order samples by group
group_order <- c("WT_CTRL", "KO_CTRL", "WT_CS", "KO_CS", "WT_BLM", "KO_BLM")
sample_info <- sample_info %>%
  mutate(Group = paste(Genotype, Treatment, sep = "_")) %>%
  mutate(Group = factor(Group, levels = group_order)) %>%
  arrange(Group, Genotype, Treatment)

expr_z <- expr_z[, sample_info$sampleID, drop = FALSE]
ann_col <- sample_info %>%
  select(sampleID, Genotype, Treatment) %>%
  column_to_rownames("sampleID")

# Highlight glycolysis metabolites
labels_row <- metabolite_base[rownames(expr_z)]
labels_row <- ifelse(tolower(labels_row) %in% glycolysis_hits, labels_row, "")

ann_colors <- list(
  Genotype = c(WT = "#4B9CD3", KO = "#D33F6A"),
  Treatment = c(CTRL = "#A0A0A0", CS = "#F1C40F", BLM = "#1ABC9C")
)

pheatmap(
  expr_z,
  annotation_col = ann_col,
  show_rownames = TRUE,
  labels_row = labels_row,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  clustering_distance_rows = "euclidean",
  annotation_colors = ann_colors,
  border_color = NA,
  annotation_legend = TRUE,
  annotation_names_col = TRUE,
  filename = file.path(fig_dir, "heatmap_STAT6KO_glycolysis.png"),
  width = 9,
  height = 6
)
pheatmap(
  expr_z,
  annotation_col = ann_col,
  show_rownames = TRUE,
  labels_row = labels_row,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  clustering_distance_rows = "euclidean",
  annotation_colors = ann_colors,
  border_color = NA,
  annotation_legend = TRUE,
  annotation_names_col = TRUE,
  filename = file.path(fig_dir, "heatmap_STAT6KO_glycolysis.pdf"),
  width = 9,
  height = 6
)

# ==============================================================================
# Step 3: Differential abundance analysis (limma)
# ==============================================================================
cat("\n=== Performing differential abundance analysis ===\n")

sample_info$Genotype <- factor(sample_info$Genotype, levels = c("WT", "KO"))
sample_info$Treatment <- factor(sample_info$Treatment, levels = c("CTRL", "CS", "BLM"))

design <- model.matrix(~ Genotype * Treatment, data = sample_info)
colnames(design) <- make.names(colnames(design))
fit <- lmFit(expr_log2, design)
fit <- eBayes(fit)

contrast_mat <- makeContrasts(
  KO_CTRL_vs_WT_CTRL = GenotypeKO,
  KO_CS_vs_WT_CS = GenotypeKO + GenotypeKO.TreatmentCS,
  KO_BLM_vs_WT_BLM = GenotypeKO + GenotypeKO.TreatmentBLM,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_mat)
fit2 <- eBayes(fit2)

metab_annot <- metab_raw %>% select(-all_of(sample_cols))

# Extract and save results
save_with_annotation <- function(tbl, annot_tbl, out_path) {
  annotated <- tbl %>%
    rownames_to_column("Metabolite") %>%
    left_join(annot_tbl, by = "Metabolite")
  write_csv(annotated, out_path)
}

res_ctrl <- topTable(fit2, coef = "KO_CTRL_vs_WT_CTRL", number = Inf, adjust.method = "fdr")
res_cs <- topTable(fit2, coef = "KO_CS_vs_WT_CS", number = Inf, adjust.method = "fdr")
res_blm <- topTable(fit2, coef = "KO_BLM_vs_WT_BLM", number = Inf, adjust.method = "fdr")

save_with_annotation(res_ctrl, metab_annot, file.path(processed_dir, "diff_KO_CTRL_vs_WT_CTRL_all_metabolites.csv"))
save_with_annotation(res_cs, metab_annot, file.path(processed_dir, "diff_KO_CS_vs_WT_CS_all_metabolites.csv"))
save_with_annotation(res_blm, metab_annot, file.path(processed_dir, "diff_KO_BLM_vs_WT_BLM_all_metabolites.csv"))

cat("Saved differential abundance results\n")

# ==============================================================================
# Figure 3B: Bubble plots for pathway enrichment
# ==============================================================================
cat("\n=== Generating Figure 3B: Pathway enrichment bubble plots ===\n")

plot_bubble <- function(df, label, png_path, pdf_path, top_n = 10) {
  if (nrow(df) == 0) {
    warning("No rows to plot for ", label, "; skipping.")
    return()
  }
  df <- df[order(df$Impact, decreasing = TRUE), , drop = FALSE]
  df <- head(df, top_n)
  df$Pathway <- factor(df$Pathway, levels = rev(df$Pathway))

  plt <- ggplot(df, aes(x = Impact, y = Pathway)) +
    geom_point(aes(size = Hits, color = NegLog10P)) +
    scale_color_gradient(low = "#6a51a3", high = "#d73027") +
    scale_size(range = c(3, 10)) +
    guides(
      color = guide_colorbar(order = 1),
      size = guide_legend(order = 2)
    ) +
    theme_bw() +
    labs(
      title = label,
      x = "Impact",
      y = "Pathway",
      color = "-log10(p)",
      size = "Hits"
    )

  ggsave(png_path, plot = plt, width = 6, height = 4.5, dpi = 300)
  ggsave(pdf_path, plot = plt, width = 6, height = 4.5)
}

scfea_dir <- config$raw_data$scfea_dir
enrichment_files <- list(
  "KO_vs_WT_BLM"  = file.path(scfea_dir, "results/metabolomics_stat6/enrichment/KEGG_pathway_down_KO_vs_WT_BLM.csv"),
  "KO_vs_WT_CS"   = file.path(scfea_dir, "results/metabolomics_stat6/enrichment/KEGG_pathway_down_KO_vs_WT_CS.csv"),
  "KO_vs_WT_CTRL" = file.path(scfea_dir, "results/metabolomics_stat6/enrichment/KEGG_pathway_down_KO_vs_WT_CTRL.csv")
)

for (label in names(enrichment_files)) {
  path <- enrichment_files[[label]]
  if (file.exists(path)) {
    df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    colnames(df)[1] <- "Pathway"
    df$Raw_p <- as.numeric(df[["Raw p"]])
    df$NegLog10P <- -log10(df$Raw_p)
    df$Impact <- as.numeric(df$Impact)
    df$Hits <- as.numeric(df$Hits)

    plot_bubble(
      df = df,
      label = label,
      png_path = file.path(enrichment_dir, paste0("bubble_", label, ".png")),
      pdf_path = file.path(enrichment_dir, paste0("bubble_", label, ".pdf"))
    )
  } else {
    cat(sprintf("  Enrichment file not found: %s\n", path))
  }
}

# ==============================================================================
# Figure 3C: Bar plots of key metabolites
# ==============================================================================
cat("\n=== Generating Figure 3C: Key metabolite bar plots ===\n")

plot_metabolite_bar <- function(metabolite_name, expr_df, sample_info_df, output_prefix) {
  # Find metabolite row (case-insensitive partial match)
  row_idx <- grep(metabolite_name, rownames(expr_df), ignore.case = TRUE)
  if (length(row_idx) == 0) {
    cat(sprintf("  Metabolite '%s' not found, skipping\n", metabolite_name))
    return()
  }

  plot_data <- data.frame(
    sampleID = colnames(expr_df),
    expression = as.numeric(expr_df[row_idx[1], ]),
    stringsAsFactors = FALSE
  ) %>%
    left_join(sample_info_df, by = "sampleID") %>%
    mutate(Group = paste(Genotype, Treatment, sep = "_"))

  plot_data$Group <- factor(plot_data$Group, levels = group_order)

  # Summary statistics
  summary_df <- plot_data %>%
    group_by(Group) %>%
    summarise(
      mean = mean(expression, na.rm = TRUE),
      sd = sd(expression, na.rm = TRUE),
      n = n(),
      se = sd / sqrt(n),
      .groups = "drop"
    )

  p <- ggplot(summary_df, aes(x = Group, y = mean, fill = Group)) +
    geom_bar(stat = "identity", alpha = 0.8, width = 0.7) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  width = 0.2, color = "black") +
    geom_jitter(data = plot_data, aes(x = Group, y = expression),
                width = 0.15, size = 2, alpha = 0.6, inherit.aes = FALSE) +
    scale_fill_manual(values = c(
      "WT_CTRL" = "#4B9CD3", "KO_CTRL" = "#D33F6A",
      "WT_CS" = "#4B9CD3", "KO_CS" = "#D33F6A",
      "WT_BLM" = "#4B9CD3", "KO_BLM" = "#D33F6A"
    )) +
    labs(
      title = rownames(expr_df)[row_idx[1]],
      x = "Group",
      y = "log2(abundance + 1)"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  ggsave(file.path(fig_dir, paste0(output_prefix, "_bar_log2.png")), p, width = 6, height = 4, dpi = 300)
  ggsave(file.path(fig_dir, paste0(output_prefix, "_bar_log2.pdf")), p, width = 6, height = 4)
}

# Plot key metabolites
plot_metabolite_bar("glucose 6-phosphate", expr_log2, sample_info, "glucose6p")
plot_metabolite_bar("pyruvic", expr_log2, sample_info, "pyruvic_acid")

cat("\n=== Comment 7 Figure 3 analysis complete ===\n")
cat("All figures saved to:", fig_dir, "\n")
