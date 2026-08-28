#!/usr/bin/env Rscript

# Updated S2 QC figure from the rebuilt pre-SERRF, SERRF, and SpATS pipeline.
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
# Panels G and H read the SpATS fitted matrices written by the preprocessing pipeline (14_run_postSERRF_SpATS.R, now in _preprocessing_moved/).
# Run that script first; this one only reads its output.
new_data_dir <- Sys.getenv("NEW_DATA_DIR", "data/new_data")
repo_root <- Sys.getenv("SOLD_REPO", ".")
out_file  <- file.path(repo_root, "fig/supp/SuppFig_S2_QC_RunOrder_SERRF_PCA_SpATS.png")
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

# House styling. plot_theme comes from _common.R and is used unmodified except
# for the legend: eight panels each carrying an inside-panel legend box would
# cover the data, so legends sit under their panel here.
source("scripts/new_new_script/_common.R")

theme_qc <- plot_theme + theme(
  legend.position   = "bottom",
  legend.background = element_blank(),
  legend.direction  = "horizontal",
  legend.text       = element_text(size = 13),
  plot.margin       = margin(10, 12, 10, 12),
  plot.tag          = element_text(face = "bold", size = 20),
  plot.tag.position = c(0.01, 0.99)
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

create_run_order_plot <- function(raw_file, panel_title) {
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
  qc <- dat %>% filter(Type == "QC") %>% pull(TIC)
  qc_cv <- if (length(qc) > 1 && mean(qc) > 0) 100 * sd(qc) / mean(qc) else NA_real_
  # Reported to the console for the caption; the paper does not put it in the panel.
  message("Panel ", panel_title, " QC TIC CV: ",
          ifelse(is.finite(qc_cv), sprintf("%.1f%%", qc_cv), "not available"))
  ggplot(dat, aes(Run, TIC, colour = Type)) +
    geom_line(data = filter(dat, Type == "Sample"), aes(group = 1), linewidth = 0.3, alpha = 0.35) +
    geom_point(size = 1.6, alpha = 0.8) +
    scale_colour_manual(values = c(Sample = "#440154", Check = "#FDE725", QC = "#E64B35", Blank = "#777777", ISTD = "#00A087")) +
    scale_y_continuous(labels = scientific) +
    labs(tag = panel_title, x = "Run order", y = "TIC") + theme_qc
}

canonical_injection_name <- function(x) {
  type <- case_when(
    str_detect(x, "S1_Run") ~ "sample",
    str_detect(x, "QC_Run") ~ "qc",
    str_detect(x, "ISTD") ~ "validate",
    str_detect(x, "InjBL") ~ "blank",
    TRUE ~ NA_character_
  )
  run <- case_when(
    type == "sample" ~ str_match(x, "S1_Run(\\d+)")[, 2],
    type == "qc" ~ str_match(x, "QC_Run(\\d+)")[, 2],
    type == "validate" ~ str_match(x, "ISTD_Run(\\d+)")[, 2],
    type == "blank" ~ str_match(x, "InjBL[^_]*_Run(\\d+)")[, 2],
    TRUE ~ NA_character_
  )
  sample_id <- str_extract(x, "PI\\d+|CHECK\\d+_\\d+|Check_\\d+")
  case_when(
    type == "sample" ~ paste0("S1_", sample_id, "_Run", run),
    type == "qc" ~ paste0("QC_Run", run),
    type == "validate" ~ paste0("ISTD_Run", run),
    type == "blank" ~ paste0("InjBL_Run", run),
    TRUE ~ NA_character_
  )
}

# The pre-SERRF matrix is rebuilt from the vendor export, using the same
# canonical injection names and the same retention-time cut as the SERRF upload
# step, so panels E/F need no intermediate upload file on disk.
read_preserrf_matrix <- function(raw_file) {
  raw <- read.csv(raw_file, check.names = FALSE, stringsAsFactors = FALSE)
  raw <- raw[is.finite(raw[["row retention time"]]) & raw[["row retention time"]] >= 1, , drop = FALSE]
  peak_columns <- grep("Peak area$", names(raw), value = TRUE)
  injections <- canonical_injection_name(peak_columns)
  keep <- !is.na(injections) & !str_detect(injections, "^InjBL_") & !duplicated(injections)
  peak_columns <- peak_columns[keep]; injections <- injections[keep]
  mat <- as.matrix(raw[, peak_columns, drop = FALSE])
  storage.mode(mat) <- "numeric"
  colnames(mat) <- injections
  rownames(mat) <- as.character(raw[["row ID"]])
  list(mat = mat, meta = tibble(Injection = injections, Type = classify_type(injections)))
}

read_serrf_normalized <- function(path) {
  x <- read_csv(path, show_col_types = FALSE, name_repair = "minimal")
  ids <- as.character(x[[1]])
  mat <- as.matrix(x[, -1, drop = FALSE]); storage.mode(mat) <- "numeric"; rownames(mat) <- ids
  list(mat = mat, meta = tibble(Injection = colnames(mat), Type = classify_type(colnames(mat))))
}

create_rsd_plot <- function(path, panel_title) {
  x <- read_csv(path, show_col_types = FALSE)
  d <- tibble(Stage = factor(c("Before SERRF", "After SERRF"), levels = c("Before SERRF", "After SERRF")),
              RSD = c(median(x$QC_none, na.rm = TRUE) * 100, median(x$QC_SERRF, na.rm = TRUE) * 100),
              Pct = c(mean(x$QC_none < .30, na.rm = TRUE) * 100, mean(x$QC_SERRF < .30, na.rm = TRUE) * 100))
  message("Panel ", panel_title, " QC-RSD <30%: ", sprintf("%.1f%% to %.1f%%", d$Pct[1], d$Pct[2]))
  ggplot(d, aes(Stage, RSD, fill = Stage)) +
    geom_col(width = .62, colour = "black", linewidth = .35) +
    geom_text(aes(label = sprintf("%.1f%%", RSD)), vjust = -.35, size = 4.6, fontface = "bold") +
    scale_fill_manual(values = c("Before SERRF" = "#9E9E9E", "After SERRF" = "#009E73")) +
    labs(tag = panel_title, x = NULL, y = "Median QC-RSD (%)") +
    expand_limits(y = max(d$RSD) * 1.15) + theme_qc + theme(legend.position = "none")
}

create_pca_plot <- function(raw_file, normalized_path, panel_title) {
  raw <- read_preserrf_matrix(raw_file); serrf <- read_serrf_normalized(normalized_path)
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
  message("Panel ", panel_title, " matched injections: ", n, "; matched features: ", length(ids))
  ggplot(d, aes(PC1, PC2, colour = Stage, shape = Type)) +
    geom_point(size = 1.8, alpha = .75) +
    scale_colour_manual(values = c(Raw = "#7A7A7A", SERRF = "#009E73")) +
    scale_shape_manual(values = c(Sample = 16, Check = 17, QC = 15, ISTD = 18)) +
    labs(tag = panel_title, x = sprintf("PC1 (%.1f%%)", 100 * ve[1]),
         y = sprintf("PC2 (%.1f%%)", 100 * ve[2])) + theme_qc
}

create_residual_plot <- function(condition, panel_title) {
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
    scale_y_reverse() + coord_fixed(ratio = .55) +
    labs(tag = panel_title, x = xlab, y = ylab, fill = "Residual") +
    theme_qc + theme(legend.position = "right",
                     legend.direction = "vertical",
                     legend.title = element_text(size = 12, face = "bold"))
}

# Panel A remains sourced from the same original CTL raw-acquisition file.
p_a <- create_run_order_plot(file.path(project_root, "data/raw_lipid_intensities/SetA_lipid_FLO2019Control.csv"), "A")
p_b <- create_run_order_plot(file.path(project_root, "data/raw_lipid_intensities/SetB_lipid_FLO2022_lowP.csv"), "B")
p_c <- create_rsd_plot(file.path(new_data_dir, "SERRF Result_CTL/QC-RSDs.csv"), "C")
p_d <- create_rsd_plot(file.path(new_data_dir, "SERRF Result_LIN/QC-RSDs.csv"), "D")
p_e <- create_pca_plot(file.path(project_root, "data/raw_lipid_intensities/SetA_lipid_FLO2019Control.csv"), file.path(new_data_dir, "SERRF Result_CTL/normalized by - SERRF.csv"), "E")
p_f <- create_pca_plot(file.path(project_root, "data/raw_lipid_intensities/SetB_lipid_FLO2022_lowP.csv"), file.path(new_data_dir, "SERRF Result_LIN/normalized by - SERRF.csv"), "F")
p_g <- create_residual_plot("CTL", "G")
p_h <- create_residual_plot("LIN", "H")

figure <- (p_a + p_b) / (p_c + p_d) / (p_e + p_f) / (p_g + p_h) + plot_layout(heights = c(1.05, .9, 1.05, 1.05))
ggsave(out_file, figure, width = 20, height = 26, units = "in", dpi = 300, bg = "white")
message("Saved: ", out_file)
