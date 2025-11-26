#!/usr/bin/env Rscript
# ==============================================================================
# Export Seurat object counts and metadata for scFEA
# ==============================================================================
#
# Usage:
#   Rscript src/scfea/export_scfea_inputs.R --config config/scfea_export_at2.yaml
#   Rscript src/scfea/export_scfea_inputs.R --config config/scfea_export_ipf.yaml
# Requires: conda r_4_3 environment with Seurat.
# ==============================================================================

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

parse_meta_fields <- function(x) {
  if (is.null(x) || length(x) == 0) return(NULL)
  if (is.character(x) && length(x) == 1) return(trimws(unlist(strsplit(x, ","))))
  as.character(x)
}

prepare_feature_rename <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.list(x)) x <- unlist(x)
  if (length(x) == 0) return(NULL)
  x
}

validate_filters <- function(filters) {
  if (is.null(filters)) return(list())
  if (!is.list(filters)) stop("Filters must be a list of column/include/exclude definitions")
  lapply(filters, function(flt) {
    column <- flt$column %||% flt$col
    if (is.null(column) || !nzchar(column)) stop("Each filter must define a 'column'")
    include <- flt$include
    exclude <- flt$exclude
    if (!is.null(include) && !is.vector(include)) include <- as.vector(include)
    if (!is.null(exclude) && !is.vector(exclude)) exclude <- as.vector(exclude)
    list(
      column = as.character(column),
      include = if (is.null(include)) NULL else as.character(include),
      exclude = if (is.null(exclude)) NULL else as.character(exclude)
    )
  })
}

apply_filters_to_object <- function(obj, filters) {
  if (length(filters) == 0) return(obj)
  meta <- obj@meta.data
  keep_cells <- rownames(meta)
  for (flt in filters) {
    column <- flt$column
    if (!column %in% colnames(meta)) stop("Filter column not found in metadata: ", column)
    values <- meta[keep_cells, column, drop = TRUE]
    if (!is.null(flt$include)) keep_cells <- keep_cells[values %in% flt$include]
    values <- meta[keep_cells, column, drop = TRUE]
    if (!is.null(flt$exclude)) keep_cells <- keep_cells[!values %in% flt$exclude]
    if (length(keep_cells) == 0) stop("No cells remain after applying filter on column: ", column)
  }
  message("Cells remaining after filters: ", length(keep_cells))
  subset(obj, cells = keep_cells)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(yaml)
})

parse_args <- function(args) {
  out <- list()
  skip_next <- FALSE
  for (i in seq_along(args)) {
    if (skip_next) {
      skip_next <- FALSE
      next
    }
    token <- args[[i]]
    if (!startsWith(token, "--")) next
    token <- substring(token, 3)
    if (identical(token, "help")) {
      out$help <- TRUE
      next
    }
    if (grepl("=", token, fixed = TRUE)) {
      kv <- strsplit(token, "=", fixed = TRUE)[[1]]
      out[[kv[1]]] <- kv[2]
    } else {
      if (i == length(args)) stop(sprintf("Missing value for argument '--%s'", token))
      out[[token]] <- args[[i + 1]]
      skip_next <- TRUE
    }
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

config <- list()
if (!is.null(args$config)) {
  if (!file.exists(args$config)) stop("Config file not found: ", args$config)
  config <- yaml::read_yaml(args$config)
}
filters <- validate_filters(config$filters)
feature_output_names <- prepare_feature_rename(config$feature_output_names)

if (!is.null(args$help) || (length(args) == 0 && length(config) == 0)) {
  cat("Required arguments (CLI or config):\n",
      "  --seurat_file <path>\n",
      "  --counts_output <path>\n",
      "  --metadata_output <path>\n",
      "Optional arguments:\n",
      "  --config <config file path>\n",
      "  --assay <assay name, default RNA>\n",
      "  --slot <data slot, default counts>\n",
      "  --meta_fields <comma separated list of metadata columns>\n",
      "  --feature_genes <comma separated gene list>\n",
      sep = "")
  quit(status = 0)
}

get_option <- function(name, default = NULL) {
  if (!is.null(args[[name]])) return(args[[name]])
  if (!is.null(config[[name]])) return(config[[name]])
  default
}

seurat_file <- get_option("seurat_file")
counts_output <- get_option("counts_output")
metadata_output <- get_option("metadata_output")
if (is.null(seurat_file) || is.null(counts_output) || is.null(metadata_output)) {
  stop("Missing required arguments: seurat_file, counts_output, metadata_output")
}

assay_name <- get_option("assay", "RNA")
slot_name <- get_option("slot", "counts")

feature_genes <- parse_meta_fields(args$feature_genes) %||%
  parse_meta_fields(config$feature_genes)
feature_assay <- get_option("feature_assay", assay_name)
feature_slot <- get_option("feature_slot", "data")
feature_prefix <- get_option("feature_prefix", NULL)

meta_fields <- parse_meta_fields(args$meta_fields) %||%
  parse_meta_fields(config$meta_fields) %||%
  c("orig.ident", "nomipf", "celltype", "seurat_clusters",
    "nCount_RNA", "nFeature_RNA", "percent.mt")
if (length(filters) > 0) {
  filter_columns <- vapply(filters, function(flt) flt$column, character(1))
  meta_fields <- unique(c(meta_fields, filter_columns))
}

message("Loading Seurat object from: ", seurat_file)
obj <- readRDS(seurat_file)

if (length(filters) > 0) {
  message("Applying metadata filters...")
  obj <- apply_filters_to_object(obj, filters)
}

counts <- GetAssayData(obj, assay = assay_name, slot = slot_name)
dir.create(dirname(counts_output), recursive = TRUE, showWarnings = FALSE)
message("Writing counts to: ", counts_output)
write.csv(as.matrix(counts), counts_output)

meta <- obj@meta.data
meta$cell_id <- rownames(meta)
available_fields <- intersect(meta_fields, colnames(meta))
missing_fields <- setdiff(meta_fields, available_fields)
if (length(missing_fields) > 0) {
  warning("Missing metadata fields: ", paste(missing_fields, collapse = ", "))
}
meta_df <- meta[, c("cell_id", available_fields), drop = FALSE]

if (!is.null(feature_genes) && length(feature_genes) > 0) {
  if (!feature_assay %in% names(obj@assays)) {
    stop("Assay not found for feature extraction: ", feature_assay)
  }
  message("Fetching feature values for: ", paste(feature_genes, collapse = ", "))
  assay_features <- rownames(obj[[feature_assay]])
  missing_features <- setdiff(feature_genes, assay_features)
  if (length(missing_features) > 0) {
    warning("Features not found in assay '", feature_assay, "': ",
            paste(missing_features, collapse = ", "))
  }
  available_features <- intersect(feature_genes, assay_features)
  if (length(available_features) > 0) {
    feature_data <- FetchData(
      obj,
      vars = available_features,
      assay = feature_assay,
      slot = feature_slot,
      cells = colnames(obj)
    )
    feature_data$cell_id <- rownames(feature_data)
    rename_lookup <- feature_output_names %||% character(0)
    if (length(rename_lookup) > 0) {
      for (nm in names(rename_lookup)) {
        idx <- which(colnames(feature_data) == nm)
        if (length(idx) > 0) colnames(feature_data)[idx] <- rename_lookup[[nm]]
      }
    }
    if (!is.null(feature_prefix) && nzchar(feature_prefix)) {
      gene_idx <- which(colnames(feature_data) != "cell_id")
      colnames(feature_data)[gene_idx] <- paste0(feature_prefix, colnames(feature_data)[gene_idx])
    }
    gene_cols <- setdiff(colnames(feature_data), "cell_id")
    if (length(intersect(gene_cols, colnames(meta_df))) > 0) {
      stop("Feature columns overlap with existing metadata columns. Adjust feature_prefix or feature_output_names.")
    }
    ordered_idx <- match(meta_df$cell_id, feature_data$cell_id)
    feature_values <- feature_data[ordered_idx, gene_cols, drop = FALSE]
    meta_df <- cbind(meta_df, feature_values)
  } else {
    warning("None of the requested features were available in assay '", feature_assay, "'.")
  }
}
dir.create(dirname(metadata_output), recursive = TRUE, showWarnings = FALSE)
message("Writing metadata to: ", metadata_output)
write.csv(meta_df, metadata_output, row.names = FALSE)

message("Extraction complete.")
