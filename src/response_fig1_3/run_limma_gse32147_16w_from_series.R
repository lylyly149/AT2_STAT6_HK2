#!/usr/bin/env Rscript
# ==============================================================================
# GSE32147 bulk RNA-seq (silica model, 16 week): limma from GEO series matrix
# ==============================================================================
# - Uses local GEO series matrix and platform annotation
# - Filters to 16-week samples, contrasts 15mg silica vs control
# - Outputs limmaresult to data/processed
# ==============================================================================

set.seed(123)

suppressPackageStartupMessages({
  library(yaml)
  library(limma)
  library(GEOquery)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
})

# ------------------------------------------------------------------------------#
# Paths and setup
# ------------------------------------------------------------------------------#
args <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
script_dir <- dirname(script_path)
project_root <- normalizePath(file.path(script_dir, "../.."))

config <- yaml::read_yaml(file.path(project_root, "config/paths.yaml"))
setwd(project_root)

series_path <- config$raw_data$gse32147_raw$series_matrix
platform_path <- config$raw_data$gse32147_raw$platform
output_csv <- config$processed_data$gse32147_limmaresult
output_pdf <- file.path(dirname(output_csv), "volcano_gse32147_16w_from_series.pdf")

dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)

cat("\n====================\n")
cat("GSE32147 limma (16 week) from series matrix\n")
cat("====================\n\n")

cat("Series matrix: ", series_path, "\n", sep = "")
cat("Platform:      ", platform_path, "\n", sep = "")
cat("Output CSV:    ", output_csv, "\n\n", sep = "")
cat("Volcano PDF:   ", output_pdf, "\n\n", sep = "")

# ------------------------------------------------------------------------------#
# Load series and platform from local files
# ------------------------------------------------------------------------------#
gset <- getGEO(filename = series_path, getGPL = FALSE)
exprset <- exprs(gset)
pdata <- pData(gset)

gpl <- getGEO(filename = platform_path)
gpl_anno <- Table(gpl)
gpl_anno$Symbol <- str_split_fixed(gpl_anno$Symbol, pattern = "_", 2)[, 1]

# ------------------------------------------------------------------------------#
# Build expression matrix with gene symbols
# ------------------------------------------------------------------------------#
n_samples <- ncol(exprset)
ex <- cbind(exprset, ID = rownames(exprset))
exprset_new <- merge(ex, gpl_anno, by = "ID", incomparables = TRUE)
exprset_new <- data.frame(exprset_new[, 2:(n_samples + 1)], symbol = exprset_new$Symbol)

matrix <- na.omit(exprset_new)
matrix[, 1:n_samples] <- apply(matrix[, 1:n_samples], 2, as.numeric)

matrix2 <- aggregate(x = matrix[, 1:n_samples],
                     by = list(matrix$symbol),
                     FUN = mean)
matrix <- data.frame(matrix2[, 2:ncol(matrix2)], row.names = matrix2$Group.1)
matrix <- matrix[-1, ]
matrix <- normalizeBetweenArrays(matrix)

# ------------------------------------------------------------------------------#
# Subset 16-week samples and run limma
# ------------------------------------------------------------------------------#
sel <- grepl("16 week", pdata$source_name_ch1)
pdata_16w <- pdata[sel, ]
sel_names <- rownames(pdata_16w)
matrix <- matrix[, sel_names]

group_info <- ifelse(grepl("15mg", pdata_16w$source_name_ch1), "Sio2", "control")
names(group_info) <- sel_names

f <- factor(group_info, levels = c("Sio2", "control"))
design <- model.matrix(~0 + f)
colnames(design) <- c("Sio2", "control")
rownames(design) <- sel_names

contrast.matrix <- makeContrasts(`Sio2-control` = Sio2 - control, levels = design)

fit <- lmFit(matrix, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
nrDEG <- na.omit(topTable(fit2, adjust = "BH", coef = 1, n = Inf))

write.csv(nrDEG, file = output_csv)
cat(sprintf("Wrote results: %s\n", output_csv))
cat(sprintf("Rows: %d | adj.P.Val < 0.01: %d | P.Value < 0.01: %d\n",
            nrow(nrDEG),
            sum(nrDEG$adj.P.Val < 0.01),
            sum(nrDEG$P.Value < 0.01)))

# ------------------------------------------------------------------------------#
# Volcano plot (padj threshold only; up/down by sign of logFC)
# ------------------------------------------------------------------------------#
padj_thr <- 0.01
# Rat symbols for plotting labels
genes_to_label <- c("Ccl7","Arg1","Col16a1","Col17a1","Col15a1","Col16a3")

plot_df <- nrDEG
plot_df$gene <- rownames(plot_df)
plot_df$status <- ifelse(plot_df$adj.P.Val < padj_thr & plot_df$logFC > 0, "Up",
                         ifelse(plot_df$adj.P.Val < padj_thr & plot_df$logFC < 0, "Down", "Not changed"))
plot_df$neg_log10_padj <- -log10(plot_df$adj.P.Val)

cnt <- table(plot_df$status)
for(s in c("Down","Not changed","Up")) if(is.na(cnt[s])) cnt[s] <- 0
labels <- c(
  paste0("Down regulated(", cnt["Down"], ")"),
  paste0("Not changed(",   cnt["Not changed"], ")"),
  paste0("Up regulated(",  cnt["Up"], ")")
)
names(labels) <- c("Down","Not changed","Up")

label_df <- plot_df[tolower(plot_df$gene) %in% tolower(genes_to_label), ]

p <- ggplot(plot_df, aes(x = logFC, y = neg_log10_padj, color = status)) +
  geom_point(size = 2, alpha = 0.9) +
  geom_hline(yintercept = -log10(padj_thr), linetype = 2, color = "grey40") +
  geom_text_repel(data = label_df, aes(label = gene),
                  size = 4, max.overlaps = Inf, color = "black") +
  scale_color_manual(
    values = c("Down" = "#00bfc4", "Not changed" = "grey70", "Up" = "#f8766d"),
    breaks = c("Down","Not changed","Up"),
    labels = labels
  ) +
  labs(x = "log2(Fold Change)\nSio2_vs_control (16w)",
       y = "-log10(adj.P.value)",
       color = NULL) +
  theme_bw(base_size = 14) +
  theme(legend.position = "top")

ggsave(output_pdf, p, width = 6, height = 6)
cat(sprintf("Wrote volcano: %s\n", output_pdf))
cat("\n=== Done ===\n")
