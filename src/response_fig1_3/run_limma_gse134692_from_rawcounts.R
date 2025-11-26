#!/usr/bin/env Rscript
# ==============================================================================
# GSE134692 bulk RNA-seq: Reproduce limma-voom DEG results from raw counts
# ==============================================================================
# - Uses local raw counts, design table, and gene annotation
# - Replicates the filtering/voom/limma steps from TCGAdownload.R
# - Outputs limmaresult.csv-equivalent to data/processed
# ==============================================================================

set.seed(123)

suppressPackageStartupMessages({
  library(yaml)
  library(limma)
  library(edgeR)
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

counts_path <- config$raw_data$gse134692_raw$counts
design_path <- config$raw_data$gse134692_raw$design
anno_path <- config$raw_data$gse134692_raw$annotation
output_csv <- config$processed_data$gse134692_limmaresult
output_pdf <- file.path(dirname(output_csv), "volcano_gse134692_from_rawcounts.pdf")

dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)

cat("\n====================\n")
cat("GSE134692 limma-voom from raw counts\n")
cat("====================\n\n")

cat("Counts file:   ", counts_path, "\n", sep = "")
cat("Design file:   ", design_path, "\n", sep = "")
cat("Annotation:    ", anno_path, "\n", sep = "")
cat("Output CSV:    ", output_csv, "\n\n", sep = "")
cat("Volcano PDF:   ", output_pdf, "\n\n", sep = "")

# ------------------------------------------------------------------------------#
# Load data
# ------------------------------------------------------------------------------#
counts <- read.table(counts_path, header = TRUE, sep = "\t", check.names = FALSE)
rownames(counts) <- counts$ensembl
counts <- counts[, -1]

anno <- read.table(anno_path, header = TRUE, sep = "\t", check.names = FALSE)
pdata <- read.table(design_path, header = TRUE, sep = "\t", check.names = FALSE)

# ------------------------------------------------------------------------------#
# Prepare matrix (replicates TCGAdownload.R logic)
# ------------------------------------------------------------------------------#
same <- intersect(rownames(counts), anno$ensembl)
rownames(anno) <- anno$ensembl
anno <- anno[, -1]

count_1 <- cbind(anno[same, ], counts[same, ])
count_1 <- count_1[, -c(2:9)]

rt <- as.matrix(count_1)
rownames(rt) <- rt[, 1]
COUNT_1 <- rt[, 2:ncol(rt)]
dimnames <- list(rownames(COUNT_1), colnames(COUNT_1))
data <- matrix(as.numeric(as.matrix(COUNT_1)), nrow = nrow(COUNT_1), dimnames = dimnames)
data <- avereps(data)

COUNT <- data[rowMeans(data) > 0, ]
COUNT <- round(COUNT)

# Drop ALI samples
pdata <- pdata[pdata$DiseaseStatus != "ALI", ]
col_names <- pdata$sample_id
COUNT <- COUNT[, col_names]

cat(sprintf("After preprocessing: %d genes x %d samples\n", nrow(COUNT), ncol(COUNT)))
cat(sprintf("Group counts: %s\n\n",
            paste(names(table(pdata$DiseaseStatus)), table(pdata$DiseaseStatus), collapse = "; ")))

# Gene filter (<=10 counts in >=80% of samples are removed)
qualified_genes <- c()
for (g in rownames(COUNT)) {
  qualification <- COUNT[g, ] <= 10
  if (sum(qualification) < 0.8 * length(COUNT)) qualified_genes <- c(qualified_genes, g)
}
COUNT_FILTER <- COUNT[qualified_genes, ]
cat(sprintf("After count filter: %d genes retained\n\n", nrow(COUNT_FILTER)))

# ------------------------------------------------------------------------------#
# limma-voom differential expression
# ------------------------------------------------------------------------------#
group <- factor(pdata$DiseaseStatus, levels = c("Normal", "IPF"))
design <- model.matrix(~0 + group)
colnames(design) <- c("Normal", "IPF")
rownames(design) <- colnames(COUNT_FILTER)

y <- DGEList(counts = COUNT_FILTER, group = group)
keep_gene <- rowSums(cpm(y) > 1) >= 2
y <- y[keep_gene, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

v <- voom(y, design, plot = FALSE, normalize = "quantile")

fit <- lmFit(v, design)
cont.matrix <- makeContrasts(`IPF-Normal` = IPF - Normal, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

nrDEG_limma_voom <- topTable(fit2, coef = "IPF-Normal", adjust = "BH", number = Inf)
nrDEG_limma_voom <- na.omit(nrDEG_limma_voom)

write.csv(nrDEG_limma_voom, file = output_csv)
cat(sprintf("Wrote results: %s\n", output_csv))
cat(sprintf("Rows: %d | adj.P.Val < 0.01: %d | P.Value < 0.01: %d\n",
            nrow(nrDEG_limma_voom),
            sum(nrDEG_limma_voom$adj.P.Val < 0.01),
            sum(nrDEG_limma_voom$P.Value < 0.01)))
# ------------------------------------------------------------------------------#
# Volcano plot (padj threshold only; up/down by sign of logFC)
# ------------------------------------------------------------------------------#
padj_thr <- 0.01
genes_to_label <- c("ACTA2","COL6A2","COL18A1","COL1A1","COL1A2","TGFB3","COL15A1")

plot_df <- nrDEG_limma_voom
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

label_df <- plot_df[plot_df$gene %in% genes_to_label, ]

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
  labs(x = "log2(Fold Change)\nIPF_vs_Ctrl",
       y = "-log10(adj.P.value)",
       color = NULL) +
  theme_bw(base_size = 14) +
  theme(legend.position = "top")

ggsave(output_pdf, p, width = 6, height = 6)
cat(sprintf("Wrote volcano: %s\n", output_pdf))
cat("\n=== Done ===\n")
