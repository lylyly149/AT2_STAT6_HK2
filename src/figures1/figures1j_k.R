# ------------------------------------------------------------------------------
# figures1j_k.R
# Generate PF-target DEG intersections and heatmaps for GSE134692, GSE32147 (16w),
# and GSE16846, mirroring final_analysis.R but using figure1-specific limma inputs
# and writing outputs with the prefix "figures1jk_" under figures/figures1.
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrastr)
  library(GEOquery)
  library(limma)
  library(pheatmap)
  library(biomaRt)
  library(dplyr)
  library(stringr)
  library(ggrepel)
  library(readr)
  library(TFTF)
  library(GOSemSim)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(dynamicTreeCut)
  library(WGCNA)
  library(yaml)
})

# ---------------------------------------------------------------------------
# Paths and prefixes
# ---------------------------------------------------------------------------
base_dir <- file.path("D:", "researchproject", "STAT6 and EMT",
                      "AT2_STAT6_HK2", "figures", "figures1")

# Fallback for WSL/Linux path
if (!dir.exists(base_dir) && dir.exists("/mnt/d/researchproject/STAT6 and EMT/AT2_STAT6_HK2/figures/figures1")) {
  base_dir <- "/mnt/d/researchproject/STAT6 and EMT/AT2_STAT6_HK2/figures/figures1"
}

if (!dir.exists(base_dir)) {
  stop("figures/figures1 directory not found. Please check the path.")
}

# Locate project root and config
project_root <- normalizePath(file.path(base_dir, "..", ".."), winslash = "/", mustWork = TRUE)
config_path <- file.path(project_root, "config", "paths.local.yaml")
if (!file.exists(config_path)) {
  config_path <- file.path(project_root, "config", "paths.yaml")
}
cfg <- yaml::read_yaml(config_path)

output_prefix <- "figures1jk_"

deg_gse134692_path <- file.path(base_dir, "figures1b_limmaresult_gse134692_from_rawcounts.csv")
deg_gse32147_path  <- file.path(base_dir, "figures1c_limmaresult_gse32147_16w_from_series.csv")
deg_gse16846_path  <- file.path(base_dir, "figures1d_limmaresult_gse16846_from_series.csv")

pf_hs_path <- "/mnt/d/researchproject/STAT6 and EMT/pharmacy/PF target Hs.csv"
pf_rn_path <- "/mnt/d/researchproject/STAT6 and EMT/pharmacy/PFtarget Rn.csv"

setwd(base_dir)

cat("\n====================\n")
cat("Step 1: Filter DEGs and intersect with PF targets\n")
cat("====================\n\n")

# Read limma results (precomputed for figures1)
deg_gse134692 <- read.csv(deg_gse134692_path, row.names = 1, check.names = FALSE)
deg_gse16846  <- read.csv(deg_gse16846_path,  row.names = 1, check.names = FALSE)
deg_gse32147  <- read.csv(deg_gse32147_path,  row.names = 1, check.names = FALSE)

# Filter DEGs (same thresholds as final_analysis.R)
deg_gse134692_sig <- deg_gse134692[deg_gse134692$adj.P.Val < 0.01, ]
deg_gse16846_sig  <- deg_gse16846[deg_gse16846$P.Value   < 0.01, ]
deg_gse32147_sig  <- deg_gse32147[deg_gse32147$adj.P.Val < 0.01, ]

# Convert to uppercase gene symbols
deg_gse134692_sig$Gene <- toupper(rownames(deg_gse134692_sig))
deg_gse16846_sig$Gene  <- toupper(rownames(deg_gse16846_sig))
deg_gse32147_sig$Gene  <- toupper(rownames(deg_gse32147_sig))

# PF target lists (unchanged paths)
pf_target_hs <- read.csv(pf_hs_path, header = TRUE, stringsAsFactors = FALSE)
pf_target_rn <- read.csv(pf_rn_path, header = TRUE, stringsAsFactors = FALSE)
pf_target_hs_upper <- toupper(pf_target_hs[, 1])
pf_target_rn_upper <- toupper(pf_target_rn[, 1])

# Intersect with PF targets
intersect_gse134692 <- deg_gse134692_sig[deg_gse134692_sig$Gene %in% pf_target_hs_upper, ]
intersect_gse16846  <- deg_gse16846_sig[deg_gse16846_sig$Gene   %in% pf_target_rn_upper, ]
intersect_gse32147  <- deg_gse32147_sig[deg_gse32147_sig$Gene   %in% pf_target_rn_upper, ]

cat(paste("GSE134692:", nrow(intersect_gse134692), "PF-related DEGs\n"))
cat(paste("GSE16846:",  nrow(intersect_gse16846),  "PF-related DEGs\n"))
cat(paste("GSE32147:",  nrow(intersect_gse32147),  "PF-related DEGs\n"))

cat("\n====================\n")
cat("Step 2: Generate heatmaps\n")
cat("====================\n\n")

# ===== GSE134692 Heatmap =====
cat("Processing GSE134692...\n")

count <- read.table("/mnt/d/researchproject/STAT6 and EMT/pharmacy/12.17 revise/another IPF dataset/GSE134692_raw_counts.txt/GSE134692_raw_counts.txt",
                    header = TRUE, sep = "\t", check.names = FALSE)
rownames(count) <- count$ensembl
count <- count[, -1]

anno <- read.table("/mnt/d/researchproject/STAT6 and EMT/pharmacy/12.17 revise/another IPF dataset/GSE134692_gene_annotation.txt/GSE134692_gene_annotation.txt",
                   header = TRUE, sep = "", check.names = FALSE)

same <- intersect(rownames(count), anno$ensembl)
rownames(anno) <- anno$ensembl
anno <- anno[, -1]
count_1 <- cbind(anno[same, ], count[same, ])
count_1 <- count_1[, -c(2:9)]

rt <- as.matrix(count_1)
rownames(rt) <- rt[, 1]
COUNT_1 <- rt[, 2:ncol(rt)]
dimnames <- list(rownames(COUNT_1), colnames(COUNT_1))
data <- matrix(as.numeric(as.matrix(COUNT_1)), nrow = nrow(COUNT_1), dimnames = dimnames)
data <- avereps(data)
matrix_134692 <- data[rowMeans(data) > 0, ]
matrix_134692 <- normalizeBetweenArrays(matrix_134692)

heatmap_info <- read.csv("/mnt/d/researchproject/STAT6 and EMT/pharmacy/12.17 revise/another IPF dataset/exp for heatmap.csv",
                         header = FALSE, row.names = 1, check.names = FALSE)
group_info_134692 <- as.character(heatmap_info[1, ])
sample_names_134692 <- as.character(heatmap_info[2, ])

matrix_134692 <- matrix_134692[, sample_names_134692]
names(group_info_134692) <- sample_names_134692

genes_heatmap_134692 <- intersect(rownames(matrix_134692), rownames(intersect_gse134692))
heat_134692 <- matrix_134692[genes_heatmap_134692, ]

annotation_col_134692 <- data.frame(group = group_info_134692, row.names = sample_names_134692)
choose_matrix_134692 <- t(scale(t(heat_134692)))
bk <- seq(-1.4, 1.4, 0.01)

pdf(file.path(base_dir, paste0(output_prefix, "heatmap_gse134692_PF_targets.pdf")), width = 10, height = 12)
pheatmap(choose_matrix_134692, cluster_cols = TRUE, scale = "row",
         annotation_col = annotation_col_134692, show_colnames = FALSE, show_rownames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(length(bk) - 1),
         legend_breaks = seq(-1, 1, 1), breaks = bk, fontsize_row = 8)
dev.off()
cat(paste("  Saved", paste0(output_prefix, "heatmap_gse134692_PF_targets.pdf"),
          "(", length(genes_heatmap_134692), "genes)\n"))

# ===== GSE16846 Heatmap =====
cat("Processing GSE16846...\n")

gse <- getGEO(filename = "/mnt/d/researchproject/STAT6 and EMT/pharmacy/12.17 revise/another BLEO dataset/GSE16846_series_matrix.txt.gz",
              getGPL = FALSE)
gpl <- getGEO(filename = "/mnt/d/researchproject/STAT6 and EMT/pharmacy/12.17 revise/another BLEO dataset/GPL339.soft")
gset <- gse
exprset <- exprs(gset)
gpl_anno <- Table(gpl)
gpl_anno$Symbol <- str_split_fixed(gpl_anno$`Gene Symbol`, pattern = "_", 2)[, 1]

ex <- cbind(exprset, ID = row.names(exprset))
exprset_new <- merge(ex, gpl_anno, by = "ID", incomparables = TRUE)
n_samples <- ncol(exprset)
exprset_new <- data.frame(exprset_new[, 2:(n_samples + 1)], symbol = exprset_new$Symbol)
matrix_16846 <- na.omit(exprset_new)
matrix_16846[, 1:n_samples] <- apply(matrix_16846[, 1:n_samples], 2, FUN = as.numeric)

matrix2 <- aggregate(x = matrix_16846[, 1:(ncol(matrix_16846) - 1)],
                     by = list(matrix_16846$symbol), FUN = mean)
matrix_16846 <- data.frame(matrix2[, 2:ncol(matrix2)], row.names = matrix2$Group.1)
matrix_16846 <- matrix_16846[-1, ]
matrix_16846 <- normalizeBetweenArrays(matrix_16846)

heatmap_info <- read.csv("/mnt/d/researchproject/STAT6 and EMT/pharmacy/12.17 revise/another BLEO dataset/potential target exp for heatmap.csv",
                         header = FALSE, row.names = 1, check.names = FALSE)
group_info_16846 <- as.character(heatmap_info[1, ])
sample_names_16846 <- as.character(heatmap_info[2, ])

matrix_16846 <- matrix_16846[, sample_names_16846]

genes_heatmap_16846 <- intersect(rownames(matrix_16846), rownames(intersect_gse16846))
heat_16846 <- matrix_16846[genes_heatmap_16846, ]

annotation_col_16846 <- data.frame(group = group_info_16846, row.names = sample_names_16846)
choose_matrix_16846 <- t(scale(t(heat_16846)))

pdf(file.path(base_dir, paste0(output_prefix, "heatmap_gse16846_PF_targets.pdf")), width = 10, height = 12)
pheatmap(choose_matrix_16846, cluster_cols = TRUE, scale = "row",
         annotation_col = annotation_col_16846, show_colnames = FALSE, show_rownames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(length(bk) - 1),
         legend_breaks = seq(-1, 1, 1), breaks = bk, fontsize_row = 8)
dev.off()
cat(paste("  Saved", paste0(output_prefix, "heatmap_gse16846_PF_targets.pdf"),
          "(", length(genes_heatmap_16846), "genes)\n"))

# ===== GSE32147 Heatmap (16 week only) =====
cat("Processing GSE32147 (16 week samples)...\n")

gse <- getGEO(filename = "/mnt/d/researchproject/STAT6 and EMT/pharmacy/GSE32147_series_matrix.txt.gz",
              getGPL = FALSE)
gpl <- getGEO(filename = "/mnt/d/researchproject/STAT6 and EMT/pharmacy/GPL6101.soft")
gset <- gse
exprset <- exprs(gset)
pdata_32147 <- pData(gset)
gpl_anno <- Table(gpl)
gpl_anno$Symbol <- str_split_fixed(gpl_anno$Symbol, pattern = "_", 2)[, 1]

ex <- cbind(exprset, ID = row.names(exprset))
exprset_new <- merge(ex, gpl_anno, by = "ID", incomparables = TRUE)
exprset_new <- data.frame(exprset_new[, 2:61], symbol = exprset_new$Symbol)
matrix_32147 <- na.omit(exprset_new)
matrix_32147[, 1:60] <- apply(matrix_32147[, 1:60], 2, FUN = as.numeric)

matrix2 <- aggregate(x = matrix_32147[, 1:60], by = list(matrix_32147$symbol), FUN = mean)
matrix_32147 <- data.frame(matrix2[, 2:ncol(matrix2)], row.names = matrix2$Group.1)
matrix_32147 <- matrix_32147[-1, ]
matrix_32147 <- normalizeBetweenArrays(matrix_32147)

pdata_32147_16w <- pdata_32147[grepl("16 week", pdata_32147$source_name_ch1), ]
matrix_32147 <- matrix_32147[, rownames(pdata_32147_16w)]

group_info_32147 <- ifelse(grepl("15mg", pdata_32147_16w$source_name_ch1), "Silica", "Control")
names(group_info_32147) <- rownames(pdata_32147_16w)

genes_heatmap_32147 <- intersect(rownames(matrix_32147), rownames(intersect_gse32147))
heat_32147 <- matrix_32147[genes_heatmap_32147, ]

annotation_col_32147 <- data.frame(group = group_info_32147, row.names = rownames(pdata_32147_16w))
choose_matrix_32147 <- t(scale(t(heat_32147)))

pdf(file.path(base_dir, paste0(output_prefix, "heatmap_gse32147_PF_targets.pdf")), width = 10, height = 12)
pheatmap(choose_matrix_32147, cluster_cols = TRUE, scale = "row",
         annotation_col = annotation_col_32147, show_colnames = FALSE, show_rownames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(length(bk) - 1),
         legend_breaks = seq(-1, 1, 1), breaks = bk, fontsize_row = 8)
dev.off()
cat(paste("  Saved", paste0(output_prefix, "heatmap_gse32147_PF_targets.pdf"),
          "(", length(genes_heatmap_32147), "genes)\n"))

cat("\n====================\n")
cat("Step 3: Convert rat symbols to human symbols\n")
cat("====================\n\n")

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
rat   <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

convert_rat_to_human <- function(rat_genes) {
  result <- data.frame(Rat_Symbol = rat_genes, Human_Symbol = NA, stringsAsFactors = FALSE)
  for (i in seq_along(rat_genes)) {
    tryCatch({
      conversion <- getLDS(attributes = "rgd_symbol", filters = "rgd_symbol",
                           values = rat_genes[i], mart = rat,
                           attributesL = "hgnc_symbol", martL = human, uniqueRows = TRUE)
      if (nrow(conversion) > 0) result$Human_Symbol[i] <- conversion[1, 2]
    }, error = function(e) {})
    if (i %% 10 == 0) cat(paste0("  ", i, "/", length(rat_genes), " processed\n"))
  }
  return(result)
}

cat("Converting GSE16846 genes...\n")
gse16846_conversion <- convert_rat_to_human(rownames(intersect_gse16846))
intersect_gse16846$Human_Symbol <- gse16846_conversion$Human_Symbol
intersect_gse16846 <- intersect_gse16846[!is.na(intersect_gse16846$Human_Symbol), ]
cat(paste("  Retained", nrow(intersect_gse16846), "genes\n"))

cat("Converting GSE32147 genes...\n")
gse32147_conversion <- convert_rat_to_human(rownames(intersect_gse32147))
intersect_gse32147$Human_Symbol <- gse32147_conversion$Human_Symbol
intersect_gse32147 <- intersect_gse32147[!is.na(intersect_gse32147$Human_Symbol), ]
cat(paste("  Retained", nrow(intersect_gse32147), "genes\n"))

intersect_gse134692$Human_Symbol <- rownames(intersect_gse134692)

cat("\n====================\n")
cat("Step 4: Save results\n")
cat("====================\n\n")

write.csv(intersect_gse134692, file.path(base_dir, paste0(output_prefix, "gse134692_PF_target_DEGs.csv")), row.names = TRUE)
write.csv(intersect_gse16846,  file.path(base_dir, paste0(output_prefix, "gse16846_PF_target_DEGs.csv")),  row.names = TRUE)
write.csv(intersect_gse32147,  file.path(base_dir, paste0(output_prefix, "gse32147_PF_target_DEGs.csv")),  row.names = TRUE)

cat("\n=== Analysis Complete ===\n")
cat(paste("GSE134692:", nrow(intersect_gse134692), "PF-related DEGs\n"))
cat(paste("GSE16846:",  nrow(intersect_gse16846),  "PF-related DEGs\n"))
cat(paste("GSE32147:",  nrow(intersect_gse32147),  "PF-related DEGs\n"))
cat("\nFiles generated:\n")
cat(paste0("  - ", output_prefix, "gse134692_PF_target_DEGs.csv\n"))
cat(paste0("  - ", output_prefix, "gse16846_PF_target_DEGs.csv\n"))
cat(paste0("  - ", output_prefix, "gse32147_PF_target_DEGs.csv\n"))
cat(paste0("  - ", output_prefix, "heatmap_gse134692_PF_targets.pdf\n"))
cat(paste0("  - ", output_prefix, "heatmap_gse16846_PF_targets.pdf\n"))
cat(paste0("  - ", output_prefix, "heatmap_gse32147_PF_targets.pdf\n"))

# ==============================================================================
# Step 5: Create union dataset (in-memory, no CSV) — equivalent to create_union_dataset.R
# ==============================================================================
cat("\n====================\n")
cat("Step 5: Create union dataset (no CSV output)\n")
cat("====================\n\n")

# Add Source column and select overlapping columns
add_source_and_subset <- function(df, source_label) {
  df$Source <- source_label
  cols_to_keep <- c("Human_Symbol", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "Source")
  df_subset <- df[, intersect(cols_to_keep, colnames(df))]
  return(df_subset)
}

deg_134692_subset <- add_source_and_subset(intersect_gse134692, "IPF")
deg_16846_subset  <- add_source_and_subset(intersect_gse16846,  "BLM")
deg_32147_subset  <- add_source_and_subset(intersect_gse32147,  "CS")

union_all <- rbind(
  deg_134692_subset,
  deg_16846_subset,
  deg_32147_subset
)

# Drop rows with missing Human_Symbol
union_all <- union_all[!is.na(union_all$Human_Symbol), ]

# Add Gene_ID and order (mirroring create_union_dataset.R)
union_all$Gene_ID <- paste0(union_all$Human_Symbol, "_", union_all$Source)
cols_order <- c("Gene_ID", "Human_Symbol", setdiff(colnames(union_all), c("Gene_ID", "Human_Symbol")))
union_all <- union_all[, cols_order]
union_all <- union_all[order(union_all$Human_Symbol, union_all$Source), ]

cat("Union created with all records\n")
cat("Total records: ", nrow(union_all), "\n", sep = "")
cat("\nBreakdown by source:\n")
cat("IPF: ", sum(union_all$Source == "IPF"), "\n", sep = "")
cat("BLM: ", sum(union_all$Source == "BLM"), "\n", sep = "")
cat("CS: ", sum(union_all$Source == "CS"), "\n", sep = "")
cat("\nTotal unique genes: ", length(unique(union_all$Human_Symbol)), "\n", sep = "")
gene_counts <- table(union_all$Human_Symbol)
cat("Genes appearing in 1 dataset: ", sum(gene_counts == 1), "\n", sep = "")
cat("Genes appearing in 2 datasets: ", sum(gene_counts == 2), "\n", sep = "")
cat("Genes appearing in 3 datasets: ", sum(gene_counts == 3), "\n", sep = "")

# ==============================================================================
# Step 6: Gene similarity clustering — equivalent to gene_similarity_clustering.R
# ==============================================================================
cat("\n====================\n")
cat("Step 6: Gene similarity clustering\n")
cat("====================\n\n")

unique_genes <- unique(union_all$Human_Symbol)

if (length(unique_genes) == 0) {
  stop("No genes available for clustering.")
}

# Map symbols to Entrez IDs
gene_mapping <- bitr(unique_genes,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)

cat("Successfully mapped genes:", nrow(gene_mapping), "\n")
cat("Unmapped genes:", length(unique_genes) - nrow(gene_mapping), "\n\n")

mapping_out <- file.path(base_dir, paste0(output_prefix, "gene_symbol_to_entrezid_mapping.csv"))
write.csv(gene_mapping, mapping_out, row.names = FALSE)

# GO semantic data
hsGO <- godata("org.Hs.eg.db", ont = "BP", computeIC = TRUE)

entrez_ids <- gene_mapping$ENTREZID

cat("Calculating semantic similarity matrix for", length(entrez_ids), "genes...\n")
sim_matrix <- mgeneSim(entrez_ids, semData = hsGO, measure = "Wang",
                       combine = "BMA", verbose = TRUE)

# Align symbols to matrix
matrix_entrez <- rownames(sim_matrix)
gene_mapping_filtered <- gene_mapping[gene_mapping$ENTREZID %in% matrix_entrez, ]
gene_mapping_filtered <- gene_mapping_filtered[match(matrix_entrez, gene_mapping_filtered$ENTREZID), ]

rownames(sim_matrix) <- gene_mapping_filtered$SYMBOL
colnames(sim_matrix) <- gene_mapping_filtered$SYMBOL

sim_out <- file.path(base_dir, paste0(output_prefix, "gene_similarity_matrix.csv"))
write.csv(sim_matrix, sim_out)
cat("Similarity matrix saved to", sim_out, "\n\n")

# Hierarchical clustering and dynamic tree cut
dist_matrix <- as.dist(1 - sim_matrix)
hc <- hclust(dist_matrix, method = "average")

modules <- cutreeDynamic(dendro = hc, distM = as.matrix(dist_matrix),
                         deepSplit = 1, pamRespectsDendro = FALSE,
                         minClusterSize = 15)

module_assignment <- data.frame(
  Gene = gene_mapping_filtered$SYMBOL,
  EntrezID = gene_mapping_filtered$ENTREZID,
  Module = modules
)
module_assignment <- module_assignment[order(module_assignment$Module), ]

module_out <- file.path(base_dir, paste0(output_prefix, "gene_module_assignment.csv"))
write.csv(module_assignment, module_out, row.names = FALSE)

cat("Number of modules identified:", max(modules), "\n")
cat("Module sizes:\n")
print(table(modules))
cat("\n")

# ---------------------------------------------------------------------------
# Convert module format (equivalent to convert_module_format.R) with prefix
# ---------------------------------------------------------------------------
n_modules <- max(modules)
module_colors <- labels2colors(1:n_modules)

color_map <- data.frame(
  Module = 1:n_modules,
  module_color = module_colors
)

module_data_formatted <- merge(module_assignment, color_map, by = "Module", all.x = TRUE)
module_data_formatted$module_name <- paste0("module_", module_data_formatted$Module)

gene_module <- data.frame(
  gene_id = module_data_formatted$Gene,
  module_name = module_data_formatted$module_name,
  module_color = module_data_formatted$module_color
)

gene_module <- gene_module[order(gene_module$module_name, gene_module$gene_id), ]

gene_module_out <- file.path(base_dir, paste0(output_prefix, "gene_module.csv"))
write.csv(gene_module, gene_module_out, row.names = FALSE, quote = FALSE)

cat("Module color mapping:\n")
for (i in 1:n_modules) {
  n_genes <- sum(gene_module$module_name == paste0("module_", i))
  cat(sprintf("  module_%d -> %-10s (%d genes)\n", i, color_map$module_color[i], n_genes))
}

cat("Gene module file saved to ", gene_module_out, "\n\n", sep = "")

# ==============================================================================
# Step 7: Filter edge file (equivalent to fix_edge_file.R) using configured path
# ==============================================================================
cat("====================\n")
cat("Step 7: Filter edge file\n")
cat("====================\n\n")

edge_file_path <- tryCatch(cfg$raw_data$network$edge_file, error = function(e) NULL)
if (is.null(edge_file_path)) {
  edge_file_path <- "/mnt/d/researchproject/STAT6 and EMT/wd/edge_file.csv"
}
if (!file.exists(edge_file_path) && file.exists("D:/researchproject/STAT6 and EMT/wd/edge_file.csv")) {
  edge_file_path <- "D:/researchproject/STAT6 and EMT/wd/edge_file.csv"
}

cat("Reading edge file from: ", edge_file_path, "\n", sep = "")
edge_file <- read.csv(edge_file_path, header = TRUE)

cat("Original edge count:", nrow(edge_file), "\n")
cat("Columns in edge file:", paste(colnames(edge_file), collapse = ", "), "\n\n")

valid_genes <- gene_module$gene_id

edge_file_filtered <- edge_file[
  edge_file$from %in% valid_genes & edge_file$to %in% valid_genes,
]

cat("Filtered edge count:", nrow(edge_file_filtered), "\n")
cat("Removed edges:", nrow(edge_file) - nrow(edge_file_filtered), "\n")

if ("combined_score" %in% colnames(edge_file_filtered)) {
  cat("\nCombined score statistics:\n")
  cat("  Min:", min(edge_file_filtered$combined_score, na.rm = TRUE), "\n")
  cat("  Max:", max(edge_file_filtered$combined_score, na.rm = TRUE), "\n")
  cat("  Mean:", mean(edge_file_filtered$combined_score, na.rm = TRUE), "\n")
} else {
  cat("\nNote: No combined_score column found in edge file\n")
}
cat("\n")

edge_filtered_out <- file.path(base_dir, paste0(output_prefix, "edge_file_filtered.csv"))
edge_backup_out   <- file.path(base_dir, paste0(output_prefix, "edge_file_original_backup.csv"))

write.csv(edge_file_filtered, edge_filtered_out, row.names = FALSE)
write.csv(edge_file, edge_backup_out, row.names = FALSE)

cat("Edge files saved:\n")
cat("  1. ", edge_filtered_out, " - Cleaned edge file\n", sep = "")
cat("  2. ", edge_backup_out,   " - Backup of original\n", sep = "")
cat("\nNow you can use the filtered edge file for downstream network plotting.\n\n")

# ==============================================================================
# Step 8: Basic network layout (equivalent to 3-1_basic_network_plot.R)
# ==============================================================================
cat("====================\n")
cat("Step 8: Basic network layout\n")
cat("====================\n\n")

# ks_large_network function (minimal inline copy from 2-ks_large_network_函数.R)
ks_large_network <- function(edge_df, node_df,
                             cluster_spread = 0.06,
                             pull_toward_center = 0.25,
                             vertex_size = 3,
                             repel_radius = 0.2,
                             repel_strength = 0.05,
                             repel_iter = 50) {
  nodes <- node_df$gene
  module <- as.factor(node_df$module)
  n_nodes <- length(nodes)
  module_levels <- levels(module)
  n_modules <- length(module_levels)
  mod_sizes <- as.numeric(table(module))

  r_mod <- sqrt(runif(n_modules))
  theta_mod <- runif(n_modules, 0, 2 * pi)
  centers <- cbind(r_mod * cos(theta_mod), r_mod * sin(theta_mod))

  for (iter in 1:repel_iter) {
    for (i in 1:(n_modules - 1)) {
      for (j in (i + 1):n_modules) {
        dx <- centers[i, 1] - centers[j, 1]; dy <- centers[i, 2] - centers[j, 2]
        dist <- sqrt(dx^2 + dy^2)
        if (dist < repel_radius) {
          push <- (repel_radius - dist) * repel_strength
          if (dist > 1e-6) { ux <- dx / dist; uy <- dy / dist } else { ux <- runif(1, -1, 1); uy <- runif(1, -1, 1) }
          centers[i, 1] <- centers[i, 1] + push * ux; centers[i, 2] <- centers[i, 2] + push * uy
          centers[j, 1] <- centers[j, 1] - push * ux; centers[j, 2] <- centers[j, 2] - push * uy
        }
      }
    }
    rr <- sqrt(rowSums(centers^2))
    idx <- which(rr > 1)
    if (length(idx) > 0) centers[idx, ] <- centers[idx, ] / rr[idx]
  }

  pool_size <- max(20000, round(10 * n_nodes))
  theta_pool <- runif(pool_size, 0, 2 * pi)
  r_pool <- sqrt(runif(pool_size))
  r_pool <- r_pool + 0.15 * (1 - r_pool)
  r_pool[r_pool > 1] <- 1
  pool_coords <- cbind(r_pool * cos(theta_pool), r_pool * sin(theta_pool))

  dist_mat <- sapply(1:n_modules, function(j) {
    dx <- pool_coords[, 1] - centers[j, 1]; dy <- pool_coords[, 2] - centers[j, 2]; dx^2 + dy^2
  })
  assigned <- max.col(-dist_mat)

  chosen_idx <- integer(0); node_order <- integer(0)
  for (m in 1:n_modules) {
    members <- which(module == module_levels[m]); need <- length(members)
    cand <- which(assigned == m)
    if (length(cand) < need) {
      n_extra <- need - length(cand)
      theta_extra <- runif(n_extra, 0, 2 * pi)
      r_extra <- runif(n_extra, 0, 0.05 + sqrt(sum(centers[m, ]^2)))
      extra_pts <- cbind(centers[m, 1] + r_extra * cos(theta_extra),
                         centers[m, 2] + r_extra * sin(theta_extra))
      rr <- sqrt(extra_pts[, 1]^2 + extra_pts[, 2]^2)
      extra_pts[rr > 1, ] <- extra_pts[rr > 1, ] / rr[rr > 1]
      pool_coords <- rbind(pool_coords, extra_pts)
      cand <- c(cand, (nrow(pool_coords) - n_extra + 1):nrow(pool_coords))
    }
    pick <- sample(cand, need)
    chosen_idx <- c(chosen_idx, pick)
    node_order <- c(node_order, members)
  }

  coords <- pool_coords[chosen_idx, ]
  module_index <- match(module[node_order], module_levels)
  centers_rep <- centers[module_index, , drop = FALSE]
  coords <- coords * (1 - pull_toward_center) + centers_rep * pull_toward_center +
    matrix(rnorm(n_nodes * 2, sd = cluster_spread), ncol = 2)
  rr <- sqrt(coords[, 1]^2 + coords[, 2]^2)
  coords[rr > 1, ] <- coords[rr > 1, ] / rr[rr > 1]

  node_coords <- data.frame(
    gene = nodes[node_order],
    module = module[node_order],
    x = coords[, 1],
    y = coords[, 2]
  )

  g <- igraph::graph_from_data_frame(edge_df, vertices = node_df, directed = FALSE)
  g <- igraph::simplify(g)

  return(list(node_coords = node_coords, g = g))
}

# Prepare node data (harmonize columns)
node_df <- gene_module
colnames(node_df)[colnames(node_df) == "gene_id"] <- "gene"
colnames(node_df)[colnames(node_df) == "module_name"] <- "module"

edge_input <- edge_file_filtered

set.seed(12345)

layout_result <- ks_large_network(edge_input, node_df,
                                  cluster_spread = 0.06,
                                  pull_toward_center = 0.20,
                                  vertex_size = 4,
                                  repel_radius = 0.25,
                                  repel_strength = 0.08,
                                  repel_iter = 50)

out_rdata_path <- file.path(base_dir, paste0(output_prefix, "out.RData"))
node_coords <- layout_result$node_coords
g <- layout_result$g
save(node_coords, g, file = out_rdata_path)

cat("\nOutput files for network layout:\n")
cat("  - ", out_rdata_path, " (node_coords, g)\n", sep = "")

# ==============================================================================
# Step 9: Regulation-colored network data (equivalent to 3-4_regulation_network.R)
# ==============================================================================
cat("====================\n")
cat("Step 9: Regulation-colored network data\n")
cat("====================\n\n")

# Reload to avoid stale env
load(out_rdata_path)

# Edges and coords
edges <- igraph::as_data_frame(g, what = "edges")

# Read original edge file to get combined_score (use filtered to align genes)
edge_file_raw <- edge_input

if ("combined_score" %in% colnames(edge_file_raw)) {
  edges <- merge(edges, edge_file_raw[, c("from", "to", "combined_score")],
                 by = c("from", "to"), all.x = TRUE)
  cat("Merged with edge file to get combined_score\n")
}

node_names <- node_coords$gene
coords_mat <- cbind(node_coords$x, node_coords$y)
rownames(coords_mat) <- node_names

edge.data <- data.frame(
  from.x = coords_mat[match(edges$from, node_names), 1],
  from.y = coords_mat[match(edges$from, node_names), 2],
  to.x   = coords_mat[match(edges$to, node_names), 1],
  to.y   = coords_mat[match(edges$to, node_names), 2]
)

if ("combined_score" %in% colnames(edges)) {
  edge.data$combined_score <- edges$combined_score
  cat("Combined score found. Mapping edge width based on score.\n")
  cat("  Original edges:", nrow(edge.data), "\n")
  cat("  Score range:", min(edge.data$combined_score, na.rm = TRUE), "to",
      max(edge.data$combined_score, na.rm = TRUE), "\n")

  edge.data <- edge.data[edge.data$combined_score >= 0.6, ]
  cat("  After filtering (score >= 0.6):", nrow(edge.data), "edges\n")
} else {
  edge.data$combined_score <- 0.5
  cat("No combined_score found. Using default edge width.\n")
}
cat("\n")

cat("Adding gene source and regulation info...\n\n")

# Use in-memory union_all from Step 5
source_data <- union_all
cat("Loaded source data:", nrow(source_data), "records\n")

deg_data_ranked <- source_data %>%
  group_by(Source) %>%
  mutate(
    padj_percentile = rank(-log10(adj.P.Val), ties.method = "average") / n(),
    logfc_percentile = rank(abs(logFC), ties.method = "average") / n()
  ) %>%
  ungroup()

cat("Percentile calculation complete.\n")
cat("  Datasets:", paste(unique(deg_data_ranked$Source), collapse = ", "), "\n\n")

gene_info <- deg_data_ranked %>%
  group_by(Human_Symbol) %>%
  summarise(
    source_group = paste(sort(unique(Source)), collapse = "+"),
    n_sources = n_distinct(Source),
    padj_percentile_mean = mean(padj_percentile, na.rm = TRUE),
    logfc_percentile_mean = mean(logfc_percentile, na.rm = TRUE),
    logfc_signs = paste(unique(sign(logFC)), collapse = ","),
    all_positive = all(logFC > 0),
    all_negative = all(logFC < 0),
    mixed_direction = n_distinct(sign(logFC[logFC != 0])) > 1,
    .groups = "drop"
  )

gene_info$regulation <- ifelse(gene_info$mixed_direction, "Inconsistent",
                               ifelse(gene_info$all_positive, "Up-regulated",
                                      ifelse(gene_info$all_negative, "Down-regulated", "Inconsistent")))

gene_info$source_norm <- (gene_info$n_sources - 1) / 2
gene_info$combined_score <- (gene_info$source_norm * 0.5 +
                               gene_info$logfc_percentile_mean * 0.4 +
                               gene_info$padj_percentile_mean * 0.1)
gene_info$combined_score_transformed <- gene_info$combined_score^1.5
gene_info$node_size <- 3 + gene_info$combined_score_transformed * 9
gene_info$node_alpha <- 0.6 + gene_info$combined_score_transformed * 0.4

cat("Regulation status distribution:\n")
print(table(gene_info$regulation))
cat("\n")

cat("Score summary:\n")
cat("  Source range:", min(gene_info$n_sources), "to", max(gene_info$n_sources), "\n")
cat("  Node size range:", round(min(gene_info$node_size), 2),
    "to", round(max(gene_info$node_size), 2), "\n")
cat("  Node alpha range:", round(min(gene_info$node_alpha), 2),
    "to", round(max(gene_info$node_alpha), 2), "\n\n")

node_coords$regulation <- gene_info$regulation[match(node_coords$gene, gene_info$Human_Symbol)]
node_coords$node_size <- gene_info$node_size[match(node_coords$gene, gene_info$Human_Symbol)]
node_coords$node_alpha <- gene_info$node_alpha[match(node_coords$gene, gene_info$Human_Symbol)]

unmatched <- sum(is.na(node_coords$regulation))
if (unmatched > 0) {
  median_size <- median(gene_info$node_size, na.rm = TRUE)
  median_alpha <- median(gene_info$node_alpha, na.rm = TRUE)
  cat("WARNING:", unmatched, "genes in network not found in source data\n")
  cat("Assigning median size:", median_size, "and alpha:", median_alpha, "\n\n")
  node_coords$regulation[is.na(node_coords$regulation)] <- "Unknown"
  node_coords$node_size[is.na(node_coords$node_size)] <- median_size
  node_coords$node_alpha[is.na(node_coords$node_alpha)] <- median_alpha
}

cat("Regulation distribution in network:\n")
print(table(node_coords$regulation))
cat("\n")

regulation_out <- file.path(base_dir, paste0(output_prefix, "node_coords_with_regulation.csv"))
write.csv(node_coords, regulation_out, row.names = FALSE)

cat("Saved regulation-annotated nodes to ", regulation_out, "\n\n", sep = "")

# Persist enriched node metadata (size/alpha/regulation) for downstream steps
save(node_coords, g, file = out_rdata_path)

# ==============================================================================
# Step 10: TF prediction (equivalent to network_tf_prediction_batch_v2.R)
# ==============================================================================
cat("====================\n")
cat("Step 10: TF prediction and importance scoring\n")
cat("====================\n\n")

node_data <- read_csv(regulation_out, show_col_types = FALSE)
cat("读取到", nrow(node_data), "个节点\n")
cat("节点分数(node_size)范围:", round(min(node_data$node_size), 2), "到",
    round(max(node_data$node_size), 2), "\n\n")

all_genes <- node_data$gene
gene_scores <- setNames(node_data$node_size, node_data$gene)

cat("准备预测", length(all_genes), "个基因的上游转录因子\n\n")

batch_size <- 50
num_batches <- ceiling(length(all_genes) / batch_size)

cat("将", length(all_genes), "个基因分成", num_batches, "批处理\n")
cat("每批", batch_size, "个基因\n\n")

all_tf_target_pairs <- list()

for (i in 1:num_batches) {
  start_idx <- (i - 1) * batch_size + 1
  end_idx <- min(i * batch_size, length(all_genes))
  batch_genes <- all_genes[start_idx:end_idx]

  cat("批次", i, "/", num_batches, ": 处理基因", start_idx, "-", end_idx,
      "(", length(batch_genes), "个基因)\n")

  tryCatch({
    result <- predict_TF(
      datasets = c("hTFtarget", "KnockTF", "ENCODE", "CHEA", "TRRUST", "GTRD", "ChIP_Atlas"),
      target = batch_genes,
      app = FALSE
    )

    prediction_dbs <- c("hTFtarget", "KnockTF", "ENCODE", "CHEA", "TRRUST", "GTRD", "ChIP_Atlas")

    batch_pairs <- list()
    for (db_name in prediction_dbs) {
      if (!is.null(result[[db_name]]) && is.data.frame(result[[db_name]]) && nrow(result[[db_name]]) > 0) {
        df <- result[[db_name]]
        if ("TF" %in% colnames(df) && "Target" %in% colnames(df)) {
          tf_target <- df %>%
            select(TF, Target) %>%
            mutate(database = db_name, batch = i)
          batch_pairs[[db_name]] <- tf_target
        }
      }
    }

    if (length(batch_pairs) > 0) {
      batch_combined <- bind_rows(batch_pairs)
      all_tf_target_pairs[[i]] <- batch_combined
      cat("  -> 找到", nrow(batch_combined), "条TF-target记录\n")
    } else {
      cat("  -> 未找到TF-target记录\n")
    }

  }, error = function(e) {
    cat("  -> 批次", i, "出错:", e$message, "\n")
  })

  if (i < num_batches) {
    Sys.sleep(2)
  }
}

cat("\n所有批次处理完成！\n\n")

cat("合并所有批次结果...\n")
all_pairs <- bind_rows(all_tf_target_pairs)
cat("总共", nrow(all_pairs), "条原始TF-target记录\n\n")

tf_target_all_out <- file.path(base_dir, paste0(output_prefix, "network_tf_target_all_raw.csv"))
write.csv(all_pairs, tf_target_all_out, row.names = FALSE)

cat("========== 筛选高可信度TF-Target对 ==========\n\n")

tf_target_support <- all_pairs %>%
  group_by(TF, Target) %>%
  summarise(
    db_support_count = n_distinct(database),
    supporting_databases = paste(sort(unique(database)), collapse = ";"),
    .groups = "drop"
  )

cat("去重后共有", nrow(tf_target_support), "个唯一TF-Target对\n")

high_confidence_pairs <- tf_target_support %>%
  filter(db_support_count >= 3)

cat("至少3个数据库支持的TF-Target对:", nrow(high_confidence_pairs), "\n")
cat("数据库支持数分布:\n")
print(table(high_confidence_pairs$db_support_count))
cat("\n")

tf_target_high_out <- file.path(base_dir, paste0(output_prefix, "tf_target_high_confidence.csv"))
write.csv(high_confidence_pairs, tf_target_high_out, row.names = FALSE)

cat("========== 计算转录因子重要性分数 ==========\n\n")

high_confidence_pairs$target_node_score <- gene_scores[high_confidence_pairs$Target]

unmatched <- sum(is.na(high_confidence_pairs$target_node_score))
if (unmatched > 0) {
  cat("警告:", unmatched, "个target基因未找到node_score，将使用中位数填充\n")
  median_score <- median(gene_scores, na.rm = TRUE)
  high_confidence_pairs$target_node_score[is.na(high_confidence_pairs$target_node_score)] <- median_score
}

tf_importance <- high_confidence_pairs %>%
  group_by(TF) %>%
  summarise(
    target_count = n_distinct(Target),
    total_pairs = n(),
    avg_db_support = mean(db_support_count),
    max_db_support = max(db_support_count),
    avg_target_score = mean(target_node_score),
    max_target_score = max(target_node_score),
    tf_importance_score = sum(target_node_score),
    top_targets = paste(head(Target[order(-target_node_score)], 10), collapse = ";"),
    all_databases = paste(sort(unique(unlist(strsplit(supporting_databases, ";")))), collapse = ";"),
    .groups = "drop"
  ) %>%
  arrange(desc(tf_importance_score))

cat("发现", nrow(tf_importance), "个转录因子调控网络基因\n\n")

cat("TF重要性分数统计:\n")
cat("  分数范围:", round(min(tf_importance$tf_importance_score), 2), "到",
    round(max(tf_importance$tf_importance_score), 2), "\n")
cat("  平均靶基因数:", round(mean(tf_importance$target_count), 1), "\n")
cat("  平均数据库支持数:", round(mean(tf_importance$avg_db_support), 2), "\n\n")

tf_importance_out <- file.path(base_dir, paste0(output_prefix, "tf_importance_score_based.csv"))
write.csv(tf_importance, tf_importance_out, row.names = FALSE)

cat("Top 30 最重要的转录因子:\n")
print(head(tf_importance, 30))

cat("\n========== 生成详细TF-Target关系表 ==========\n\n")

# ==============================================================================
# Step 11: TF coverage analysis (equivalent to tf_network_coverage_analysis.R)
# ==============================================================================
cat("====================\n")
cat("Step 11: TF coverage analysis\n")
cat("====================\n\n")

tf_importance <- read_csv(tf_importance_out, show_col_types = FALSE)
network_high_conf <- read_csv(tf_target_high_out, show_col_types = FALSE)

cat("网络分析发现", nrow(tf_importance), "个转录因子\n")
cat("网络中共有", nrow(network_high_conf), "个高可信度TF-Target对（≥3数据库支持）\n\n")

network_tf_targets <- network_high_conf %>%
  group_by(TF) %>%
  summarise(
    network_target_count = n_distinct(Target),
    network_targets = paste(sort(unique(Target)), collapse = ";"),
    .groups = "drop"
  )

cat("网络中有", nrow(network_tf_targets), "个TF有高可信度靶基因\n\n")

tf_list <- unique(network_tf_targets$TF)
cat("将查询", length(tf_list), "个转录因子的全部靶基因\n\n")

all_tf_global_targets <- list()

for (i in seq_along(tf_list)) {
  tf_name <- tf_list[i]
  cat("(", i, "/", length(tf_list), ") 查询", tf_name, "的靶基因...")
  tryCatch({
    result <- predict_target(
      datasets = c("hTFtarget", "KnockTF", "ENCODE", "CHEA", "TRRUST", "GTRD", "ChIP_Atlas"),
      tf = tf_name,
      app = FALSE
    )

    prediction_dbs <- c("hTFtarget", "KnockTF", "ENCODE", "CHEA", "TRRUST", "GTRD", "ChIP_Atlas")
    tf_targets_list <- list()
    for (db_name in prediction_dbs) {
      if (!is.null(result[[db_name]]) && is.data.frame(result[[db_name]]) && nrow(result[[db_name]]) > 0) {
        df <- result[[db_name]]
        if ("Target" %in% colnames(df)) {
          df$database <- db_name
          df$TF <- tf_name
          tf_targets_list[[db_name]] <- df %>% select(TF, Target, database)
        }
      }
    }

    if (length(tf_targets_list) > 0) {
      tf_all_targets <- bind_rows(tf_targets_list)
      tf_target_support <- tf_all_targets %>%
        group_by(TF, Target) %>%
        summarise(
          db_support_count = n_distinct(database),
          supporting_databases = paste(sort(unique(database)), collapse = ";"),
          .groups = "drop"
        ) %>%
        filter(db_support_count >= 3)

      all_tf_global_targets[[tf_name]] <- tf_target_support
      cat(" 找到", nrow(tf_target_support), "个高可信度靶基因\n")
    } else {
      cat(" 无结果\n")
    }
  }, error = function(e) {
    cat(" 出错:", e$message, "\n")
  })

  if (i < length(tf_list)) Sys.sleep(1)
}

cat("\n全部TF查询完成！\n\n")

if (length(all_tf_global_targets) > 0) {
  global_targets_combined <- bind_rows(all_tf_global_targets)

  global_tf_targets <- global_targets_combined %>%
    group_by(TF) %>%
    summarise(
      global_target_count = n_distinct(Target),
      global_targets = paste(sort(unique(Target)), collapse = ";"),
      .groups = "drop"
    )

  cat("全局统计:\n")
  cat("  TF数量:", nrow(global_tf_targets), "\n")
  cat("  靶基因总数:", n_distinct(global_targets_combined$Target), "\n\n")

  global_targets_out <- file.path(base_dir, paste0(output_prefix, "tf_global_targets_high_confidence.csv"))
  global_summary_out <- file.path(base_dir, paste0(output_prefix, "tf_global_target_summary.csv"))
  write.csv(global_targets_combined, global_targets_out, row.names = FALSE)
  write.csv(global_tf_targets, global_summary_out, row.names = FALSE)

  coverage_analysis <- network_tf_targets %>%
    left_join(global_tf_targets, by = "TF") %>%
    mutate(
      coverage_rate = network_target_count / global_target_count,
      coverage_percentage = round(coverage_rate * 100, 2),
      missing_target_count = global_target_count - network_target_count
    ) %>%
    arrange(desc(coverage_rate))

  cat("覆盖率范围:",
      round(min(coverage_analysis$coverage_percentage), 2), "% 到",
      round(max(coverage_analysis$coverage_percentage), 2), "%\n")
  cat("平均覆盖率:", round(mean(coverage_analysis$coverage_percentage), 2), "%\n")
  cat("中位覆盖率:", round(median(coverage_analysis$coverage_percentage), 2), "%\n\n")

  cat("Top 20 覆盖率最高的TF:\n")
  print(coverage_analysis %>%
          select(TF, network_target_count, global_target_count, coverage_percentage) %>%
          head(20))

  cat("\n\nTop 20 覆盖率最低的TF:\n")
  print(coverage_analysis %>%
          select(TF, network_target_count, global_target_count, coverage_percentage) %>%
          tail(20))

  coverage_out <- file.path(base_dir, paste0(output_prefix, "tf_network_coverage_analysis.csv"))
  write.csv(coverage_analysis, coverage_out, row.names = FALSE)

  tf_importance_with_cov <- tf_importance %>%
    left_join(coverage_analysis %>%
                select(TF, network_target_count, global_target_count,
                       coverage_rate, coverage_percentage, missing_target_count),
              by = "TF") %>%
    arrange(desc(tf_importance_score), desc(coverage_rate))

  tf_importance_cov_out <- file.path(base_dir, paste0(output_prefix, "tf_importance_with_coverage.csv"))
  write.csv(tf_importance_with_cov, tf_importance_cov_out, row.names = FALSE)

  cat("整合完成！Top 30 TF（按重要性分数）:\n")
  print(tf_importance_with_cov %>%
          select(TF, tf_importance_score, target_count, network_target_count,
                 global_target_count, coverage_percentage) %>%
          head(30))

  # Visualizations (optional)
  coverage_pdf <- file.path(base_dir, paste0(output_prefix, "tf_coverage_distribution.pdf"))
  coverage_png <- file.path(base_dir, paste0(output_prefix, "tf_coverage_distribution.png"))
  p2_pdf <- file.path(base_dir, paste0(output_prefix, "tf_network_vs_global_targets.pdf"))
  p2_png <- file.path(base_dir, paste0(output_prefix, "tf_network_vs_global_targets.png"))

  p1 <- ggplot(coverage_analysis, aes(x = coverage_percentage)) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7, color = "white") +
    geom_vline(aes(xintercept = median(coverage_percentage)),
               color = "red", linetype = "dashed", size = 1) +
    labs(title = "Distribution of TF Network Coverage Rate",
         subtitle = paste("Red line: median =", round(median(coverage_analysis$coverage_percentage), 2), "%"),
         x = "Coverage Rate (%)",
         y = "Number of TFs") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(coverage_pdf, p1, width = 10, height = 6)
  ggsave(coverage_png, p1, width = 10, height = 6, dpi = 300)

  p2 <- ggplot(coverage_analysis, aes(x = global_target_count, y = network_target_count)) +
    geom_point(aes(color = coverage_percentage, size = coverage_percentage), alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_text(data = coverage_analysis %>% arrange(desc(coverage_percentage)) %>% head(10),
              aes(label = TF), hjust = -0.1, vjust = -0.1, size = 2.5, check_overlap = TRUE) +
    scale_color_gradient(low = "blue", high = "red", name = "Coverage (%)") +
    scale_size_continuous(range = c(2, 8), name = "Coverage (%)") +
    labs(title = "Network vs Global Target Genes for Each TF",
         subtitle = "Dashed line: 100% coverage",
         x = "Global High-Confidence Target Count",
         y = "Network High-Confidence Target Count") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(p2_pdf, p2, width = 12, height = 8)
  ggsave(p2_png, p2, width = 12, height = 8, dpi = 300)

  cat("\nCoverage analysis outputs:\n")
  cat("  - ", global_targets_out, "\n", sep = "")
  cat("  - ", global_summary_out, "\n", sep = "")
  cat("  - ", coverage_out, "\n", sep = "")
  cat("  - ", tf_importance_cov_out, "\n", sep = "")
  cat("  - ", coverage_pdf, "\n", sep = "")
  cat("  - ", coverage_png, "\n", sep = "")
  cat("  - ", p2_pdf, "\n", sep = "")
  cat("  - ", p2_png, "\n", sep = "")

} else {
  cat("No global TF targets found; skipping coverage analysis outputs.\n")
}

# ==============================================================================
# Step 12: Build module_info with terms (manual annotations provided)
# ==============================================================================
cat("\n====================\n")
cat("Step 12: Build module_info with terms\n")
cat("====================\n\n")

# Compute module centers from node_coords
module_info <- node_coords %>%
  group_by(module) %>%
  summarise(
    x_center = mean(x),
    y_center = mean(y),
    .groups = "drop"
  )

module_terms <- data.frame(
  module = paste0("module_", c(1, 10, 11, 12, 13, 14, 15, 16, 2, 3, 4, 5, 6, 7, 8, 9)),
  terms = c(
    "epithelial cell proliferation",
    "renal system development",
    "cardiac muscle tissue morphogenesis",
    "cilium organization",
    "lipid transport",
    "humoral immune response",
    "regulation of protein stability",
    "extracellular matrix disassembly",
    "respiratory system development",
    "positive regulation of cytokine production",
    "negative regulation of blood coagulation",
    "monoatomic anion transport",
    "chemokine-mediated signaling pathway",
    "regulation of tube diameter",
    "reactive oxygen species metabolic process",
    "wound healing"
  ),
  stringsAsFactors = FALSE
)

module_info <- module_info %>%
  left_join(module_terms, by = "module") %>%
  arrange(module)

missing_terms <- sum(is.na(module_info$terms))
if (missing_terms > 0) {
  cat("WARNING: ", missing_terms, " modules missing terms annotations\n", sep = "")
}

module_info_out <- file.path(base_dir, paste0(output_prefix, "module_info.csv"))
write.csv(module_info, module_info_out, row.names = FALSE)

cat("Module info with terms saved to ", module_info_out, "\n\n", sep = "")

# ==============================================================================
# Step 13: STAT6-annotated network plot (equivalent to 3-3_STAT6_annotated_network.R)
# ==============================================================================
cat("====================\n")
cat("Step 13: STAT6-annotated network plot\n")
cat("====================\n\n")

# Load network data
load(out_rdata_path)

# STAT6 targets
stat6_known_targets <- c(
  "IL4R", "FCER2", "HLA-DRA", "IGHE", "CCL11", "CCL17", "CCL22", "CCL26",
  "PDCD1LG2", "SOCS1", "GATA3", "MAF", "ARG1", "RETN", "MRC1", "ALOX15",
  "CHI3L1", "POSTN", "PPARG", "BCL2", "BCL2L1", "CDKN1A", "CIITA"
)

tf_targets <- read.csv(file.path(base_dir, paste0(output_prefix, "tf_global_targets_high_confidence.csv")), header = TRUE)
stat6_db_targets <- tf_targets %>%
  filter(TF == "STAT6") %>%
  pull(Target)

all_stat6_targets <- unique(c(stat6_known_targets, stat6_db_targets))

network_genes <- node_coords$gene
stat6_targets_in_network <- intersect(all_stat6_targets, network_genes)

cat("STAT6 targets in network:", length(stat6_targets_in_network), "\n")
cat("  Genes:", paste(head(stat6_targets_in_network, 20), collapse = ", "))
if (length(stat6_targets_in_network) > 20) cat(", ...")
cat("\n\n")

node_coords$is_stat6_target <- node_coords$gene %in% stat6_targets_in_network

# Edges with combined_score
edges <- igraph::as_data_frame(g, what = "edges")
edge_file_raw <- edge_input
if ("combined_score" %in% colnames(edge_file_raw)) {
  edges <- merge(edges, edge_file_raw[, c("from", "to", "combined_score")],
                 by = c("from", "to"), all.x = TRUE)
}

coords_mat <- cbind(node_coords$x, node_coords$y)
rownames(coords_mat) <- node_coords$gene
edge.data <- data.frame(
  from.x = coords_mat[match(edges$from, node_coords$gene), 1],
  from.y = coords_mat[match(edges$from, node_coords$gene), 2],
  to.x   = coords_mat[match(edges$to, node_coords$gene), 1],
  to.y   = coords_mat[match(edges$to, node_coords$gene), 2]
)

if ("combined_score" %in% colnames(edges)) {
  edge.data$combined_score <- edges$combined_score
  cat("  Edge score range:", min(edge.data$combined_score, na.rm = TRUE), "to",
      max(edge.data$combined_score, na.rm = TRUE), "\n")
  edge.data <- edge.data[edge.data$combined_score >= 0.6, ]
  cat("  Edges after score >= 0.6:", nrow(edge.data), "\n\n")
} else {
  edge.data$combined_score <- 0.5
  cat("No combined_score found. Using default edge width.\n\n")
}

# Module colors
n_modules <- length(unique(node_coords$module))
color.module <- c("#ecb888", "#af88bb", "#a032cb", "#efbed6", "#fc496a",
                  "#b6d37f", "#589336", "#7fd68e", "#52c465", "#3372e0",
                  "#84d7f6", "#5394c3", "#6376b3", "#7f6cd7", "#c4ceff",
                  "#fc9d40", "#5c95e0", "#cd7560", "#ff70e4")
if (n_modules > length(color.module)) {
  color.module <- colorRampPalette(color.module)(n_modules)
}
names(color.module) <- unique(node_coords$module)
color.module["module_14"] <- "#00AA88"
color.module["module_16"] <- "#FFD700"

# Regulation border colors
color.regulation <- c(
  "Up-regulated"   = "#E63946",
  "Down-regulated" = "#4A7BA7",
  "Inconsistent"   = "grey70",
  "Unknown"        = "#CCCCCC"
)

stat6_nodes <- node_coords %>% filter(is_stat6_target)
non_stat6_nodes <- node_coords %>% filter(!is_stat6_target)

cat("STAT6 targets:", nrow(stat6_nodes), "\n")
cat("Non-STAT6 targets:", nrow(non_stat6_nodes), "\n\n")

module_info_use <- read.csv(module_info_out, header = TRUE)

module_info_use <- module_info_use %>%
  mutate(angle = atan2(y_center, x_center)) %>%
  arrange(angle)

min_angle <- 0.2
for (i in 2:nrow(module_info_use)) {
  if ((module_info_use$angle[i] - module_info_use$angle[i - 1]) < min_angle) {
    module_info_use$angle[i] <- module_info_use$angle[i - 1] + min_angle
  }
}

r_text <- 1.15
module_info_use <- module_info_use %>%
  mutate(
    x_text = cos(angle) * r_text,
    y_text = sin(angle) * r_text,
    dist_to_center = sqrt((x_text - x_center)^2 + (y_text - y_center)^2),
    draw_line = dist_to_center > 0.5
  )

stat6_nodes$label_x <- stat6_nodes$x
stat6_nodes$label_y <- stat6_nodes$y

p <- ggplot() +
  geom_segment(aes(x = from.x, y = from.y, xend = to.x, yend = to.y,
                   linewidth = combined_score,
                   alpha = combined_score),
               color = "grey72", data = edge.data) +
  scale_linewidth_continuous(range = c(0, 0.05), guide = "none") +
  scale_alpha_continuous(range = c(0.2, 1.0), guide = "none") +
  geom_point_rast(aes(x = x, y = y, color = module, size = node_size, alpha = node_alpha),
                  data = non_stat6_nodes, show.legend = FALSE) +
  geom_point(aes(x = x, y = y, size = node_size * 1.5, color = regulation),
             data = stat6_nodes, shape = 21, stroke = 2.5, fill = NA,
             show.legend = FALSE) +
  geom_point_rast(aes(x = x, y = y, color = module, size = node_size, alpha = node_alpha),
                  data = stat6_nodes, show.legend = FALSE) +
  geom_text_repel(aes(x = label_x, y = label_y, label = gene),
                  data = stat6_nodes,
                  size = 3.5, fontface = "bold", color = "black",
                  bg.color = "white", bg.r = 0.1,
                  box.padding = 0.5, point.padding = 0.3,
                  segment.color = "grey50", segment.size = 0.3,
                  max.overlaps = 30, force = 2,
                  min.segment.length = 0) +
  scale_size_continuous(range = c(3, 12), guide = "none") +
  scale_alpha_continuous(range = c(0.6, 1.0), guide = "none") +
  geom_segment(data = dplyr::filter(module_info_use, draw_line),
               aes(x = x_center, y = y_center,
                   xend = x_text, yend = y_text),
               color = "grey30") +
  geom_text(data = module_info_use,
            aes(x = x_text,
                y = y_text,
                label = str_wrap(terms, width = 25),
                color = module),
            hjust = 0.5, vjust = 0.5,
            fontface = "bold", show.legend = FALSE, size = 5) +
  scale_colour_manual(values = c(color.module, color.regulation)) +
  theme_void() +
  theme(plot.margin = margin(10, 10, 10, 10)) +
  coord_cartesian(xlim = c(-1.71, 1.71), ylim = c(-1.5, 1.5), clip = "off")

plot_png <- file.path(base_dir, "network_STAT6_annotated.png")
plot_pdf <- file.path(base_dir, "network_STAT6_annotated.pdf")

cat("Saving STAT6-annotated plot...\n")
ggsave(plot_png, p, width = 16, height = 14, dpi = 300, device = "png")
ggsave(plot_pdf, p, width = 16, height = 14, device = "pdf")

stat6_output <- stat6_nodes %>%
  select(gene, module, regulation, node_size, node_alpha, x, y) %>%
  arrange(regulation, gene)

stat6_out_csv <- file.path(base_dir, "STAT6_targets_in_network.csv")
write.csv(stat6_output, file = stat6_out_csv, row.names = FALSE)

cat("\nOutput files:\n")
cat("  1. ", plot_png, " - Network with STAT6 targets highlighted (PNG)\n", sep = "")
cat("  2. ", plot_pdf, " - Network with STAT6 targets highlighted (PDF)\n", sep = "")
cat("  3. ", stat6_out_csv, " - List of STAT6 targets with regulation info\n", sep = "")
cat("\nDone!\n")

# Dendrogram plot
dendro_out <- file.path(base_dir, paste0(output_prefix, "gene_clustering_dendrogram.pdf"))
module_colors <- labels2colors(modules)

pdf(dendro_out, width = 12, height = 8)
plotDendroAndColors(hc, module_colors, "Modules",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene clustering dendrogram and module assignment")
dev.off()
cat("Dendrogram saved to", dendro_out, "\n")

# Heatmap if not too large
if (nrow(sim_matrix) <= 200) {
  heatmap_out <- file.path(base_dir, paste0(output_prefix, "gene_similarity_heatmap.pdf"))
  pdf(heatmap_out, width = 12, height = 12)
  pheatmap(sim_matrix[order(modules), order(modules)],
           cluster_rows = FALSE, cluster_cols = FALSE,
           show_rownames = FALSE, show_colnames = FALSE,
           color = colorRampPalette(c("white", "yellow", "red"))(100),
           main = "Gene similarity matrix (ordered by modules)")
  dev.off()
  cat("Heatmap saved to", heatmap_out, "\n")
} else {
  cat("Matrix too large for heatmap (", nrow(sim_matrix), " genes)\n", sep = "")
}

cat("\nGene similarity clustering outputs:\n")
cat("  - ", mapping_out, "\n", sep = "")
cat("  - ", sim_out, "\n", sep = "")
cat("  - ", module_out, "\n", sep = "")
cat("  - ", dendro_out, "\n", sep = "")
