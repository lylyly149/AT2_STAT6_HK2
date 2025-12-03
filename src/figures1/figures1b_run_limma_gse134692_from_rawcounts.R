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
  library(dplyr)
  library(tibble)
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Hs.eg.db)
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
output_dir <- file.path("figures", "figures1")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_csv <- file.path(output_dir, "figures1b_limmaresult_gse134692_from_rawcounts.csv")
output_pdf <- file.path(output_dir, "figures1b_volcano_gse134692_from_rawcounts.pdf")

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

# ------------------------------------------------------------------------------
# Enrichment on significant genes + bubble plot (Reactome/GO BP/KEGG merged)
# ------------------------------------------------------------------------------
padj_cut <- 0.05
sig_genes <- rownames(nrDEG_limma_voom)[nrDEG_limma_voom$adj.P.Val < padj_cut]

if (length(sig_genes) == 0) {
  cat("No significant genes (adj.P.Val < 0.05); skipping enrichment.\n")
  } else {
  all_symbol <- rownames(nrDEG_limma_voom)
  universe_entrez <- mapIds(org.Hs.eg.db, keys = all_symbol, column = "ENTREZID",
                            keytype = "SYMBOL", multiVals = "first") |> na.omit() |> unique()
  sig_entrez <- mapIds(org.Hs.eg.db, keys = sig_genes, column = "ENTREZID",
                       keytype = "SYMBOL", multiVals = "first") |> na.omit() |> unique()

  if (length(sig_entrez) == 0) {
    cat("No significant genes could be mapped to Entrez IDs; skipping enrichment.\n")
  } else {
    target_terms <- c(
      "T cell differentiation",
      "regulation of inflammatory response",
      "Extracellular matrix organization",
      "Collagen degradation",
      "cytokine-mediated signaling pathway",
      "Activation of Matrix Metalloproteinases"
    )

    ego <- suppressMessages(enrichGO(
      gene = sig_entrez,
      universe = universe_entrez,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.1,
      qvalueCutoff = 0.1,
      readable = TRUE
    ))

    ekegg <- tryCatch({
      suppressMessages(enrichKEGG(
        gene = sig_entrez,
        universe = universe_entrez,
        organism = "hsa",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.1,
        qvalueCutoff = 0.1
      ))
    }, error = function(e) {
      message("KEGG enrichment failed (skipping): ", conditionMessage(e))
      NULL
    })

    ereact <- suppressMessages(enrichPathway(
      gene = sig_entrez,
      universe = universe_entrez,
      organism = "human",
      pvalueCutoff = 0.1,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.1,
      readable = TRUE
    ))

    collect_res <- function(enrich_obj, db_label) {
      if (is.null(enrich_obj)) return(NULL)
      df <- as.data.frame(enrich_obj)
      if (nrow(df) == 0) return(NULL)
      df$Database <- db_label
      df
    }

    df_all <- bind_rows(
      collect_res(ego, "GO:BP"),
      collect_res(ekegg, "KEGG"),
      collect_res(ereact, "Reactome")
    )

    if (is.null(df_all) || nrow(df_all) == 0) {
      cat("Enrichment returned no terms; skipping bubble plot.\n")
    } else {
      df_targets <- df_all |>
        filter(Description %in% target_terms) |>
        distinct(Description, Database, .keep_all = TRUE)

      if (nrow(df_targets) == 0) {
        cat("Target pathways not found in enrichment results; skipping bubble plot.\n")
      } else {
        # enforce display order from target_terms
        df_targets$Description <- factor(df_targets$Description,
                                         levels = target_terms[match(tolower(target_terms), tolower(target_terms))])
        df_targets <- df_targets |>
          mutate(
            GeneRatioNum = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
            negLog10P = -log10(pvalue)
          ) |>
          arrange(negLog10P)

        bubble_plot <- ggplot(df_targets, aes(x = GeneRatioNum, y = Description)) +
          geom_point(aes(size = Count, color = negLog10P)) +
          scale_color_gradient(low = "#a1d99b", high = "#238b45") +
          scale_size_continuous(range = c(3, 10)) +
          labs(
            x = "Gene Ratio",
            y = NULL,
            color = "-log10(p value)",
            size = "Count",
            title = "Significant genes enrichment (dotplot, selected pathways)"
          ) +
          theme_bw(base_size = 12) +
          theme(
            legend.position = "right",
            plot.title = element_text(face = "bold", hjust = 0)
          )

        enrich_csv <- file.path(output_dir, "figures1e_enrichment_sig_selected_pathways.csv")
        bubble_png <- file.path(output_dir, "figures1e_enrichment_sig_dotplot.png")
        bubble_pdf <- file.path(output_dir, "figures1e_enrichment_sig_dotplot.pdf")

        write.csv(df_targets, enrich_csv, row.names = FALSE)
        ggsave(bubble_png, bubble_plot, width = 7, height = 5, dpi = 300)
        ggsave(bubble_pdf, bubble_plot, width = 7, height = 5)

        cat(sprintf("Wrote enrichment table: %s\n", enrich_csv))
        cat(sprintf("Wrote bubble plots: %s and %s\n", bubble_png, bubble_pdf))
      }
    }
  }
}
cat("\n=== Done ===\n")
