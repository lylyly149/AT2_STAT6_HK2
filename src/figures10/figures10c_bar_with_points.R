#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(svglite)
  library(ragg)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- args[grep("--file=", args)]
if (length(script_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("--file=", "", script_arg[1])))
  project_root <- normalizePath(file.path(script_dir, "../.."))
} else {
  project_root <- getwd()
}
setwd(project_root)

fig_dir <- file.path("figures", "figures10")
processed_dir <- file.path("data", "processed")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

metab_raw_path <- file.path("data", "config_data", "metab", "metab_raw.csv")
sample_info_path <- file.path("data", "config_data", "metab", "metab_sample_info.csv")

metab_raw <- read_csv(metab_raw_path, locale = locale(encoding = "GBK"), show_col_types = FALSE)
sample_info <- read_csv(sample_info_path, show_col_types = FALSE) %>%
  rename(
    sampleID = all_of("Sample Analysis Name"),
    Group = all_of("Group Name")
  ) %>%
  separate(Group, into = c("Genotype", "Treatment"), sep = "_", remove = FALSE)

treatment_order <- c("CTRL", "CS", "BLM")
genotype_order <- c("WT", "KO")
sample_info <- sample_info %>%
  mutate(
    Treatment = factor(Treatment, levels = treatment_order),
    Genotype = factor(Genotype, levels = genotype_order)
  )

metabolites <- tibble(
  metabolite_query = c("glucose 6-phosphate", "pyruvic"),
  metabolite_label = c("Glucose 6-phosphate", "Pyruvic acid"),
  output_id = c("glucose6p", "pyruvic_acid")
)

extract_metabolite <- function(query, label, output_id) {
  idx <- grep(query, metab_raw$Metabolite, ignore.case = TRUE)
  if (length(idx) == 0) {
    stop("Metabolite not found: ", query)
  }

  row <- metab_raw[idx[1], , drop = FALSE]
  raw_values <- as.numeric(row[, sample_info$sampleID])
  tibble(
    sampleID = sample_info$sampleID,
    metabolite = label,
    output_id = output_id,
    raw_intensity = raw_values,
    log2_intensity = log2(raw_values + 1)
  ) %>%
    left_join(sample_info, by = "sampleID")
}

plot_data <- bind_rows(Map(
  extract_metabolite,
  metabolites$metabolite_query,
  metabolites$metabolite_label,
  metabolites$output_id
))

write_csv(plot_data, file.path(processed_dir, "s10c_bar_with_points_source.csv"))

summary_data <- plot_data %>%
  group_by(metabolite, output_id, Treatment, Genotype) %>%
  summarise(
    mean = mean(log2_intensity, na.rm = TRUE),
    sd = sd(log2_intensity, na.rm = TRUE),
    n = dplyr::n(),
    sem = sd / sqrt(n),
    .groups = "drop"
  )

write_csv(summary_data, file.path(processed_dir, "s10c_bar_with_points_summary.csv"))

format_p <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 1e-4, "<0.0001", sprintf("%.4f", p)))
}

p_value_data <- plot_data %>%
  group_by(metabolite, output_id, Treatment) %>%
  summarise(
    comparison = "WT_vs_KO",
    p_value = t.test(
      log2_intensity[Genotype == "WT"],
      log2_intensity[Genotype == "KO"]
    )$p.value,
    p_label = format_p(p_value),
    .groups = "drop"
  )

write_csv(p_value_data, file.path(processed_dir, "s10c_bar_with_points_pvalues.csv"))

theme_s10c <- function(base_size = 5.8) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.3, colour = "black"),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 0.6),
      legend.key.size = grid::unit(2.4, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      plot.title = element_text(size = base_size + 0.3, face = "bold", hjust = 0),
      plot.margin = margin(3, 4, 3, 4),
      panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey88")
    )
}

add_p_bar <- function(plot, x_start, x_end, y, label, text_y_offset = 0.10) {
  plot +
    annotate("segment", x = x_start, xend = x_end, y = y, yend = y, linewidth = 0.25) +
    annotate("text", x = (x_start + x_end) / 2, y = y + text_y_offset, label = label, size = 1.75)
}

make_panel <- function(metabolite_name, y_limit, p_bars) {
  df <- filter(plot_data, metabolite == metabolite_name)
  sm <- filter(summary_data, metabolite == metabolite_name)

  dodge <- position_dodge(width = 0.62)
  point_pos <- position_jitterdodge(jitter.width = 0.10, jitter.height = 0, dodge.width = 0.62)

  p <- ggplot(sm, aes(x = Treatment, y = mean, fill = Genotype)) +
    geom_col(position = dodge, width = 0.52, colour = "grey35", linewidth = 0.25, alpha = 0.95) +
    geom_errorbar(
      aes(ymin = mean - sem, ymax = mean + sem),
      position = dodge,
      width = 0.14,
      linewidth = 0.35,
      colour = "black"
    ) +
    geom_point(
      data = df,
      aes(x = Treatment, y = log2_intensity, group = Genotype),
      position = point_pos,
      inherit.aes = FALSE,
      shape = 21,
      size = 1.05,
      stroke = 0.25,
      colour = "black",
      fill = "white"
    ) +
    scale_fill_manual(values = c(WT = "#A6A6A6", KO = "#4D9BD3")) +
    scale_y_continuous(
      limits = c(0, y_limit),
      expand = expansion(mult = c(0.015, 0.02))
    ) +
    labs(
      title = metabolite_name,
      x = NULL,
      y = "log2(Intensity + 1)"
    ) +
    coord_cartesian(clip = "off") +
    theme_s10c()

  for (i in seq_len(nrow(p_bars))) {
    p <- add_p_bar(p, p_bars$x_start[i], p_bars$x_end[i], p_bars$y[i], p_bars$label[i])
  }

  p
}

# Significance labels are omitted so they can be added manually during layout.
g6p_bars <- tibble(
  x_start = numeric(),
  x_end = numeric(),
  y = numeric(),
  label = character()
)

pyruvic_bars <- tibble(
  x_start = numeric(),
  x_end = numeric(),
  y = numeric(),
  label = character()
)

p_g6p <- make_panel("Glucose 6-phosphate", y_limit = 5.95, p_bars = g6p_bars)
p_pyruvic <- make_panel("Pyruvic acid", y_limit = 11.45, p_bars = pyruvic_bars)

combined <- ggarrange(
  p_g6p,
  p_pyruvic,
  ncol = 1,
  nrow = 2,
  align = "v",
  heights = c(1, 1.08)
)

save_figure <- function(plot, prefix, width_mm, height_mm, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  agg_png(paste0(prefix, ".png"), width = width_in, height = height_in, units = "in", res = dpi)
  print(plot)
  dev.off()

  ggsave(paste0(prefix, ".pdf"), plot, width = width_in, height = height_in, device = cairo_pdf)

  svglite(paste0(prefix, ".svg"), width = width_in, height = height_in)
  print(plot)
  dev.off()

  agg_tiff(paste0(prefix, ".tiff"), width = width_in, height = height_in, units = "in", res = dpi, compression = "lzw")
  print(plot)
  dev.off()
}

save_figure(combined, file.path(fig_dir, "S10c_bar_with_points"), width_mm = 88, height_mm = 112)
save_figure(p_g6p, file.path(fig_dir, "S10c_glucose6p_bar_with_points"), width_mm = 88, height_mm = 54)
save_figure(p_pyruvic, file.path(fig_dir, "S10c_pyruvic_acid_bar_with_points"), width_mm = 88, height_mm = 58)

save_figure(combined, file.path(fig_dir, "S10c_bar_with_points_taller"), width_mm = 72, height_mm = 124)
save_figure(p_g6p, file.path(fig_dir, "S10c_glucose6p_bar_with_points_taller"), width_mm = 72, height_mm = 64)
save_figure(p_pyruvic, file.path(fig_dir, "S10c_pyruvic_acid_bar_with_points_taller"), width_mm = 72, height_mm = 68)

save_figure(combined, file.path(fig_dir, "S10c_bar_with_points_taller_nosig"), width_mm = 72, height_mm = 124)
save_figure(p_g6p, file.path(fig_dir, "S10c_glucose6p_bar_with_points_taller_nosig"), width_mm = 72, height_mm = 64)
save_figure(p_pyruvic, file.path(fig_dir, "S10c_pyruvic_acid_bar_with_points_taller_nosig"), width_mm = 72, height_mm = 68)

cat("S10c bar plots with replicate points saved to: ", fig_dir, "\n", sep = "")
