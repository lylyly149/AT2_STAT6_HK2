#!/usr/bin/env Rscript
# ==============================================================================
# Response Figure 1,3: Venn diagrams of DEGs from multiple pulmonary fibrosis datasets
# ==============================================================================
#
# This script generates Venn diagrams comparing:
#   - GSE134692 (human IPF) DEGs vs PF target genes
#   - GSE16846 (mouse BLM model) DEGs vs PF target genes
#   - GSE32147 (rat CS model) DEGs vs PF target genes
#   - Three-way intersection of all datasets
#
# Runtime environment: conda activate r_4_3 (R 4.3+)
# ==============================================================================

set.seed(123)

# Load required libraries
suppressPackageStartupMessages({
  library(VennDiagram)
  library(biomaRt)
  library(yaml)
  library(dplyr)
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

# Create output directory
fig_dir <- file.path("figures", "figures1")
processed_dir <- "data/processed"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n====================\n")
cat("Response Figure 1,3: Venn diagrams for PF target gene intersection\n")
cat("====================\n\n")

# ==============================================================================
# Step 1: Load DEG results from bulk RNA-seq analyses
# ==============================================================================
cat("=== Loading DEG results ===\n")

# Load DEG results using paths from config
deg_gse134692 <- read.csv(config$raw_data$bulk_deg$gse134692, row.names=1, check.names=F)
deg_gse16846 <- read.csv(config$raw_data$bulk_deg$gse16846, row.names=1, check.names=F)
deg_gse32147 <- read.csv(config$raw_data$bulk_deg$gse32147, row.names=1, check.names=F)

cat(sprintf("Loaded datasets:\n"))
cat(sprintf("  GSE134692 (IPF): %d genes\n", nrow(deg_gse134692)))
cat(sprintf("  GSE16846 (BLM mouse): %d genes\n", nrow(deg_gse16846)))
cat(sprintf("  GSE32147 (CS rat): %d genes\n", nrow(deg_gse32147)))

# ==============================================================================
# Step 2: Filter significant DEGs
# ==============================================================================
cat("\n=== Filtering significant DEGs ===\n")

# Apply statistical cutoffs
padj_cutoff <- config$parameters$stats$padj_cutoff
pval_cutoff <- config$parameters$stats$pval_cutoff

deg_gse134692_sig <- deg_gse134692[deg_gse134692$adj.P.Val < 0.01, ]
deg_gse16846_sig <- deg_gse16846[deg_gse16846$P.Value < 0.01, ]
deg_gse32147_sig <- deg_gse32147[deg_gse32147$adj.P.Val < 0.01, ]

cat(sprintf("Filtered significant DEGs:\n"))
cat(sprintf("  GSE134692: %d DEGs (adj.P.Val < 0.01)\n", nrow(deg_gse134692_sig)))
cat(sprintf("  GSE16846: %d DEGs (P.Value < 0.01)\n", nrow(deg_gse16846_sig)))
cat(sprintf("  GSE32147: %d DEGs (adj.P.Val < 0.01)\n", nrow(deg_gse32147_sig)))

# ==============================================================================
# Step 3: Convert mouse/rat gene symbols to human orthologs
# ==============================================================================
cat("\n=== Converting species gene symbols to human ===\n")

# Setup biomaRt connections
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org/")
rat <- useMart("ensembl", dataset="rnorvegicus_gene_ensembl", host="https://dec2021.archive.ensembl.org/")

# Function to convert rat genes to human orthologs
convert_rat_to_human <- function(rat_genes) {
  result <- data.frame(Rat_Symbol=rat_genes, Human_Symbol=NA, stringsAsFactors=F)
  for(i in 1:length(rat_genes)) {
    tryCatch({
      conversion <- getLDS(attributes="rgd_symbol", filters="rgd_symbol",
                          values=rat_genes[i], mart=rat,
                          attributesL="hgnc_symbol", martL=human, uniqueRows=T)
      if(nrow(conversion) > 0) result$Human_Symbol[i] <- conversion[1, 2]
    }, error=function(e) {})
    if(i %% 100 == 0) cat(sprintf("  %d/%d genes processed\n", i, length(rat_genes)))
  }
  return(result)
}

# Convert GSE16846 (rat/mouse BLM model) to human
cat("Converting GSE16846 genes (mouse/rat to human)...\n")
gse16846_conversion <- convert_rat_to_human(rownames(deg_gse16846_sig))
deg_gse16846_sig$Human_Symbol <- gse16846_conversion$Human_Symbol
deg_gse16846_sig <- deg_gse16846_sig[!is.na(deg_gse16846_sig$Human_Symbol), ]
cat(sprintf("  Retained %d genes with human orthologs\n", nrow(deg_gse16846_sig)))

# Convert GSE32147 (rat) to human
cat("Converting GSE32147 genes (rat to human)...\n")
gse32147_conversion <- convert_rat_to_human(rownames(deg_gse32147_sig))
deg_gse32147_sig$Human_Symbol <- gse32147_conversion$Human_Symbol
deg_gse32147_sig <- deg_gse32147_sig[!is.na(deg_gse32147_sig$Human_Symbol), ]
cat(sprintf("  Retained %d genes with human orthologs\n", nrow(deg_gse32147_sig)))

# GSE134692 is already human
deg_gse134692_sig$Human_Symbol <- rownames(deg_gse134692_sig)
cat(sprintf("GSE134692: Using %d human genes\n\n", nrow(deg_gse134692_sig)))

# Save conversion results
write.csv(gse16846_conversion, file.path(processed_dir, "gse16846_ortholog_conversion.csv"), row.names=FALSE)
write.csv(gse32147_conversion, file.path(processed_dir, "gse32147_ortholog_conversion.csv"), row.names=FALSE)

# ==============================================================================
# Step 4: Load PF target genes
# ==============================================================================
cat("=== Loading PF target genes ===\n")

pf_target <- read.csv(config$raw_data$pf_target_genes, header=FALSE, stringsAsFactors=F)
pf_target_genes <- pf_target$V1
cat(sprintf("PF target genes: %d\n\n", length(pf_target_genes)))

# ==============================================================================
# Step 5: Create gene lists (normalized to uppercase)
# ==============================================================================
genes_134692 <- toupper(deg_gse134692_sig$Human_Symbol)
genes_16846 <- toupper(deg_gse16846_sig$Human_Symbol)
genes_32147 <- toupper(deg_gse32147_sig$Human_Symbol)
genes_pf_target <- toupper(pf_target_genes)

cat("Gene counts per dataset (after conversion):\n")
cat(sprintf("  GSE134692 (IPF): %d\n", length(genes_134692)))
cat(sprintf("  GSE16846 (BLM): %d\n", length(genes_16846)))
cat(sprintf("  GSE32147 (CS): %d\n", length(genes_32147)))
cat(sprintf("  PF Target: %d\n\n", length(genes_pf_target)))

# ==============================================================================
# Step 6: Calculate intersections
# ==============================================================================
cat("=== Calculating intersections ===\n")

intersect_134692_pf <- intersect(genes_134692, genes_pf_target)
intersect_16846_pf <- intersect(genes_16846, genes_pf_target)
intersect_32147_pf <- intersect(genes_32147, genes_pf_target)
intersect_all_three <- Reduce(intersect, list(genes_134692, genes_16846, genes_32147))

cat("Intersections:\n")
cat(sprintf("  GSE134692 & PF target: %d genes\n", length(intersect_134692_pf)))
cat(sprintf("  GSE16846 & PF target: %d genes\n", length(intersect_16846_pf)))
cat(sprintf("  GSE32147 & PF target: %d genes\n", length(intersect_32147_pf)))
cat(sprintf("  All three datasets: %d genes\n\n", length(intersect_all_three)))

# Save intersection gene lists
write.csv(data.frame(Gene=intersect_134692_pf), file.path(processed_dir, "intersection_GSE134692_PF_target.csv"), row.names=FALSE)
write.csv(data.frame(Gene=intersect_16846_pf), file.path(processed_dir, "intersection_GSE16846_PF_target.csv"), row.names=FALSE)
write.csv(data.frame(Gene=intersect_32147_pf), file.path(processed_dir, "intersection_GSE32147_PF_target.csv"), row.names=FALSE)
write.csv(data.frame(Gene=intersect_all_three), file.path(processed_dir, "intersection_all_three_datasets.csv"), row.names=FALSE)

# ==============================================================================
# Venn Diagram 1: GSE134692 vs PF target
# ==============================================================================
cat("=== Generating Venn Diagram 1: GSE134692 vs PF target ===\n")

pdf(file.path(fig_dir, "figures1a_venn_diagram_GSE134692_vs_PF_target.pdf"), width=8, height=8)
venn.plot1 <- venn.diagram(
  x = list(
    GSE134692 = genes_134692,
    "PF Target" = genes_pf_target
  ),
  filename = NULL,
  category.names = c("GSE134692\n(IPF)", "PF Target\nGenes"),
  fill = c("#FF6B6B", "#FFD93D"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.fontface = "bold",
  fontfamily = "sans",
  cat.fontfamily = "sans",
  margin = 0.1
)
grid.draw(venn.plot1)
dev.off()
cat("  Saved figures1a_venn_diagram_GSE134692_vs_PF_target.pdf\n")

# ==============================================================================
# Venn Diagram 2: GSE16846 vs PF target
# ==============================================================================
cat("=== Generating Venn Diagram 2: GSE16846 vs PF target ===\n")

pdf(file.path(fig_dir, "figures1a_venn_diagram_GSE16846_vs_PF_target.pdf"), width=8, height=8)
venn.plot2 <- venn.diagram(
  x = list(
    GSE16846 = genes_16846,
    "PF Target" = genes_pf_target
  ),
  filename = NULL,
  category.names = c("GSE16846\n(BLM)", "PF Target\nGenes"),
  fill = c("#4ECDC4", "#FFD93D"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.fontface = "bold",
  fontfamily = "sans",
  cat.fontfamily = "sans",
  margin = 0.1
)
grid.draw(venn.plot2)
dev.off()
cat("  Saved figures1a_venn_diagram_GSE16846_vs_PF_target.pdf\n")

# ==============================================================================
# Venn Diagram 3: GSE32147 vs PF target
# ==============================================================================
cat("=== Generating Venn Diagram 3: GSE32147 vs PF target ===\n")

pdf(file.path(fig_dir, "figures1a_venn_diagram_GSE32147_vs_PF_target.pdf"), width=8, height=8)
venn.plot3 <- venn.diagram(
  x = list(
    GSE32147 = genes_32147,
    "PF Target" = genes_pf_target
  ),
  filename = NULL,
  category.names = c("GSE32147\n(CS)", "PF Target\nGenes"),
  fill = c("#95E1D3", "#FFD93D"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.fontface = "bold",
  fontfamily = "sans",
  cat.fontfamily = "sans",
  margin = 0.1
)
grid.draw(venn.plot3)
dev.off()
cat("  Saved figures1a_venn_diagram_GSE32147_vs_PF_target.pdf\n")

# ==============================================================================
# Venn Diagram 4: Three-way intersection of all datasets
# ==============================================================================
cat("=== Generating Venn Diagram 4: Three datasets intersection ===\n")

pdf(file.path(fig_dir, "figures1a_venn_diagram_three_datasets.pdf"), width=8, height=8)
venn.plot4 <- venn.diagram(
  x = list(
    GSE134692 = genes_134692,
    GSE16846 = genes_16846,
    GSE32147 = genes_32147
  ),
  filename = NULL,
  category.names = c("GSE134692\n(IPF)", "GSE16846\n(BLM)", "GSE32147\n(CS)"),
  fill = c("#FF6B6B", "#4ECDC4", "#95E1D3"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.fontface = "bold",
  fontfamily = "sans",
  cat.fontfamily = "sans",
  margin = 0.1
)
grid.draw(venn.plot4)
dev.off()
cat("  Saved figures1a_venn_diagram_three_datasets.pdf\n")

cat("\n=== Response Figure 1,3 analysis complete ===\n")
cat("All figures saved to:", fig_dir, "\n")
