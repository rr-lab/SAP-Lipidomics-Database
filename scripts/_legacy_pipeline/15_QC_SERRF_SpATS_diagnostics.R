#!/usr/bin/env Rscript

# Updated S1 QC figure from the rebuilt pre-SERRF, SERRF, and SpATS pipeline.
# Panels G/H use the same horizontal display convention: the longer field axis
# is always plotted horizontally, avoiding a compressed vertical LIN map.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

project_root <- "/Users/nirwantandukar/Documents/Github/SoLD_paper"
new_data_dir <- Sys.getenv("QC_DATA_DIR", unset = "/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database/data/new_data")
out_file <- Sys.getenv("QC_FIG_OUT", unset = "/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database/fig/new_figures/SuppFig_QC_SERRF_SpATS_diagnostics.png")
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

# Shared project styling, with compact adjustments for the eight-panel layout.
plot_theme <- theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, colour = "black"),
    axis.line = element_line(colour = "black", linewidth = .5),
    panel.grid = element_blank(), legend.position = "top", legend.title = element_blank(),
    legend.text = element_text(size = 12), plot.margin = margin(15, 15, 15, 15)
  )

theme_qc <- plot_theme + theme(
  plot.title = element_text(size = 13, hjust = 0),
  plot.subtitle = element_text(size = 8.5, colour = "grey30"),
  axis.title = element_text(size = 12, face = "bold"),
  axis.text = element_text(size = 10, colour = "black"),
  legend.position = "bottom", legend.text = element_text(size = 10),
  plot.margin = margin(5, 7, 5, 7)
)

classify_type <- function(x) {
  case_when(
    str_detect(x, regex("QC", ignore_case = TRUE)) ~ "QC",
    str_detect(x, regex("ISTD", ignore_case = TRUE)) ~ "ISTD",
    str_detect(x, regex("CHECK", ignore_case = TRUE)) ~ "Check",
    str_detect(x, regex("BLANK", ignore_case = TRUE)) ~ "Blank",
    TRUE ~ "Sample"
  )
}

create_run_order_plot <- function(raw_file) {
  x <- read_csv(raw_file, show_col_types = FALSE, name_repair = "minimal")
  run_cols <- names(x)[str_detect(names(x), "Run[0-9]+")]
  m <- as.matrix(x[, run_cols, drop = FALSE])
  storage.mode(m) <- "numeric"
  dat <- tibble(
    Injection = run_cols,
    Run = as.integer(str_extract(run_cols, "(?<=Run)[0-9]+")),
    Type = classify_type(run_cols),
    TIC = colSums(m, na.rm = TRUE)
  ) %>% filter(Type %in% c("Sample", "Check", "QC", "Blank", "ISTD")) %>% arrange(Run)
  ggplot(dat, aes(Run, TIC, colour = Type)) +
    geom_line(data = filter(dat, Type == "Sample"), aes(group = 1), linewidth = 0.25, alpha = 0.35) +
    geom_point(size = 1.1, alpha = 0.8) +
    scale_colour_manual(values = c(Sample = "#440154", Check = "#FDE725", QC = "#E64B35", Blank = "#777777", ISTD = "#00A087")) +
    scale_y_continuous(labels = scientific) +
    labs(x = "Run order", y = "TIC") + theme_qc
}

read_serrf_input <- function(path) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  injections <- as.character(x[3, -c(1, 2)])
  types <- as.character(x[1, -c(1, 2)])
  feature_ids <- as.character(x[-c(1, 2, 3), 2])
  mat <- as.matrix(x[-c(1, 2, 3), -c(1, 2), drop = FALSE])
  storage.mode(mat) <- "numeric"
  colnames(mat) <- injections
  rownames(mat) <- feature_ids
  list(mat = mat, meta = tibble(Injection = injections, Type = case_when(types == "qc" ~ "QC", types == "validate" ~ "ISTD", TRUE ~ classify_type(injections))))
}

read_serrf_normalized <- function(path) {
  x <- read_csv(path, show_col_types = FALSE, name_repair = "minimal")
  ids <- as.character(x[[1]])
  mat <- as.matrix(x[, -1, drop = FALSE]); storage.mode(mat) <- "numeric"; rownames(mat) <- ids
  list(mat = mat, meta = tibble(Injection = colnames(mat), Type = classify_type(colnames(mat))))
}

create_rsd_plot <- function(path) {
  x <- read_csv(path, show_col_types = FALSE)
  d <- tibble(Stage = factor(c("Before SERRF", "After SERRF"), levels = c("Before SERRF", "After SERRF")),
              RSD = c(median(x$QC_none, na.rm = TRUE) * 100, median(x$QC_SERRF, na.rm = TRUE) * 100),
              Pct = c(mean(x$QC_none < .30, na.rm = TRUE) * 100, mean(x$QC_SERRF < .30, na.rm = TRUE) * 100))
  ggplot(d, aes(Stage, RSD, fill = Stage)) +
    geom_col(width = .62, colour = "black", linewidth = .25) +
    geom_text(aes(label = sprintf("%.1f%%", RSD)), vjust = -.35, size = 3.1, fontface = "bold") +
    scale_fill_manual(values = c("Before SERRF" = "#9E9E9E", "After SERRF" = "#009E73")) +
    labs(x = NULL, y = "Median QC-RSD (%)") +
    expand_limits(y = max(d$RSD) * 1.15) + theme_qc + theme(legend.position = "none")
}

create_pca_plot <- function(input_path, normalized_path) {
  raw <- read_serrf_input(input_path); serrf <- read_serrf_normalized(normalized_path)
  ids <- intersect(rownames(raw$mat), rownames(serrf$mat))
  injections <- intersect(raw$meta$Injection, serrf$meta$Injection)
  raw_meta <- raw$meta %>% filter(Injection %in% injections, Type %in% c("Sample", "Check", "QC", "ISTD"))
  injections <- raw_meta$Injection
  raw_x <- t(log10(pmax(raw$mat[ids, injections, drop = FALSE], 0) + 1))
  serrf_x <- t(log10(pmax(serrf$mat[ids, injections, drop = FALSE], 0) + 1))
  combo <- rbind(raw_x, serrf_x)
  combo[!is.finite(combo)] <- NA_real_
  med <- apply(combo, 2, median, na.rm = TRUE); med[!is.finite(med)] <- 0
  for (j in seq_len(ncol(combo))) combo[is.na(combo[, j]), j] <- med[j]
  keep <- apply(combo, 2, sd) > 0
  pca <- prcomp(combo[, keep, drop = FALSE], center = TRUE, scale. = TRUE)
  ve <- pca$sdev^2 / sum(pca$sdev^2)
  n <- length(injections)
  d <- tibble(Stage = rep(c("Raw", "SERRF"), each = n), Type = rep(raw_meta$Type, 2), PC1 = pca$x[, 1], PC2 = pca$x[, 2])
  ggplot(d, aes(PC1, PC2, colour = Stage, shape = Type)) +
    geom_point(size = 1.25, alpha = .75) +
    scale_colour_manual(values = c(Raw = "#7A7A7A", SERRF = "#009E73")) +
    scale_shape_manual(values = c(Sample = 16, Check = 17, QC = 15, ISTD = 18)) +
    labs(x = sprintf("PC1 (%.1f%%)", 100 * ve[1]), y = sprintf("PC2 (%.1f%%)", 100 * ve[2])) + theme_qc
}

create_residual_plot <- function(condition) {
  named_file <- file.path(new_data_dir, sprintf("postSERRF_%s_named_reaggregated.csv", condition))
  fitted_file <- file.path(new_data_dir, ifelse(condition == "CTL", "7_CTL_lipid_SpATS_fitted.csv", "8_LIN_lipid_SpATS_fitted.csv"))
  named <- read_csv(named_file, show_col_types = FALSE, name_repair = "minimal")
  fitted <- read_csv(fitted_file, show_col_types = FALSE, name_repair = "minimal")
  lipid <- "MGDG(18:3/18:3)"
  obs_row <- named %>% filter(Compound_Name == lipid)
  obs <- tibble(Sample = names(obs_row)[str_detect(names(obs_row), "^S1_")], Observed = as.numeric(obs_row[1, str_detect(names(obs_row), "^S1_")]))
  d <- fitted %>% select(Sample, row, col, Fitted = all_of(lipid)) %>% inner_join(obs, by = "Sample") %>% mutate(Residual = Observed - Fitted)
  # Plot the longer field dimension horizontally in both panels.
  if (diff(range(d$col)) >= diff(range(d$row))) {
    d <- d %>% mutate(x = col, y = row); xlab <- "Field column"; ylab <- "Field row"
  } else {
    d <- d %>% mutate(x = row, y = col); xlab <- "Field row"; ylab <- "Field column"
  }
  lim <- max(abs(d$Residual), na.rm = TRUE)
  ggplot(d, aes(x, y, fill = Residual)) +
    geom_tile(colour = "white", linewidth = .15) +
    scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0, limits = c(-lim, lim), oob = squish) +
    scale_y_reverse() +
    labs(x = xlab, y = ylab, fill = "Residual") +
    theme_qc + theme(legend.position = "right")
}

# Panel A remains sourced from the same original CTL raw-acquisition file.
p_a <- create_run_order_plot(file.path(project_root, "data/raw_lipid_intensities/SetA_lipid_FLO2019Control.csv"))
p_b <- create_run_order_plot(file.path(project_root, "data/raw_lipid_intensities/SetB_lipid_FLO2022_lowP.csv"))
p_c <- create_rsd_plot(file.path(new_data_dir, "SERRF Result_CTL/QC-RSDs.csv"))
p_d <- create_rsd_plot(file.path(new_data_dir, "SERRF Result_LIN/QC-RSDs.csv"))
p_e <- create_pca_plot(file.path(new_data_dir, "3_CTL_SERRF_input.csv"), file.path(new_data_dir, "SERRF Result_CTL/normalized by - SERRF.csv"))
p_f <- create_pca_plot(file.path(new_data_dir, "4_LIN_SERRF_input.csv"), file.path(new_data_dir, "SERRF Result_LIN/normalized by - SERRF.csv"))
p_g <- create_residual_plot("CTL")
p_h <- create_residual_plot("LIN")

figure <- (p_a + p_b) / (p_c + p_d) / (p_e + p_f) / (p_g + p_h) +
  plot_layout(heights = c(1.05, .9, 1.05, 1.05)) +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold", size = 16)))
ggsave(out_file, figure, width = 16, height = 20, units = "in", dpi = 300, bg = "white")
message("Saved: ", out_file)
