#!/usr/bin/env Rscript
# ==============================================================================
# GSE16846 bulk RNA-seq (BLM model): limma from GEO series matrix
# ==============================================================================
# - Uses local GEO series matrix and platform annotation
# - Replicates the notrend limma workflow used in original analysis
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

series_path <- config$raw_data$gse16846_raw$series_matrix
platform_path <- config$raw_data$gse16846_raw$platform
output_csv <- config$processed_data$gse16846_limmaresult
output_pdf <- file.path(dirname(output_csv), "volcano_gse16846_from_series.pdf")

dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)

cat("\n====================\n")
cat("GSE16846 limma from series matrix\n")
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
gpl_anno$Symbol <- str_split_fixed(gpl_anno$`Gene Symbol`, pattern = "_", 2)[, 1]
gpl_anno <- cbind(gpl_anno, symbol = gpl_anno$Symbol)

# ------------------------------------------------------------------------------#
# Build expression matrix with gene symbols
# ------------------------------------------------------------------------------#
ex <- cbind(exprset, ID = rownames(exprset))
exprset_new <- merge(ex, gpl_anno, by = "ID", incomparables = TRUE)
exprset_new <- data.frame(exprset_new[, 2:7], symbol = exprset_new$Symbol)

matrix <- na.omit(exprset_new)
matrix[, 1:6] <- apply(matrix[, 1:6], 2, as.numeric)

matrix2 <- aggregate(x = matrix[, 1:6],
                     by = list(matrix$symbol),
                     FUN = mean)
matrix <- data.frame(matrix2[, 2:ncol(matrix2)], row.names = matrix2$Group.1)
matrix <- matrix[-1, ]
matrix <- normalizeBetweenArrays(matrix)

# ------------------------------------------------------------------------------#
# Design and limma
# ------------------------------------------------------------------------------#
pdata$treat <- ifelse(pdata$`treatment:ch1` == "Bleomycin", "bleo", "ctrl")

col_names <- colnames(matrix)
group <- data.frame(group = pdata[col_names, "treat"], row.names = col_names)

f <- factor(group[, 1], levels = c("ctrl", "bleo"))
design <- model.matrix(~0 + f)
colnames(design) <- c("ctrl", "bleo")
rownames(design) <- col_names

contrast.matrix <- makeContrasts(bleo - ctrl, levels = design)

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
# Volcano plot (p-value threshold only; up/down by sign of logFC)
# ------------------------------------------------------------------------------#
p_thr <- 0.01
# Mouse symbols for plotting labels
genes_to_label <- c("Acta2","Col6a2","Col18a1","Col1a1","Col1a2","Tgfb3","Col15a1")

plot_df <- nrDEG
plot_df$gene <- rownames(plot_df)
plot_df$status <- ifelse(plot_df$P.Value < p_thr & plot_df$logFC > 0, "Up",
                         ifelse(plot_df$P.Value < p_thr & plot_df$logFC < 0, "Down", "Not changed"))
plot_df$neg_log10_p <- -log10(plot_df$P.Value)

cnt <- table(plot_df$status)
for(s in c("Down","Not changed","Up")) if(is.na(cnt[s])) cnt[s] <- 0
labels <- c(
  paste0("Down regulated(", cnt["Down"], ")"),
  paste0("Not changed(",   cnt["Not changed"], ")"),
  paste0("Up regulated(",  cnt["Up"], ")")
)
names(labels) <- c("Down","Not changed","Up")

label_df <- plot_df[tolower(plot_df$gene) %in% tolower(genes_to_label), ]

p <- ggplot(plot_df, aes(x = logFC, y = neg_log10_p, color = status)) +
  geom_point(size = 2, alpha = 0.9) +
  geom_hline(yintercept = -log10(p_thr), linetype = 2, color = "grey40") +
  geom_text_repel(data = label_df, aes(label = gene),
                  size = 4, max.overlaps = Inf, color = "black") +
  scale_color_manual(
    values = c("Down" = "#00bfc4", "Not changed" = "grey70", "Up" = "#f8766d"),
    breaks = c("Down","Not changed","Up"),
    labels = labels
  ) +
  labs(x = "log2(Fold Change)\nBleomycin_vs_ctrl",
       y = "-log10(P.value)",
       color = NULL) +
  theme_bw(base_size = 14) +
  theme(legend.position = "top")

ggsave(output_pdf, p, width = 6, height = 6)
cat(sprintf("Wrote volcano: %s\n", output_pdf))
cat("\n=== Done ===\n")
