# ═══════════════════════════════════════════════════════════════════════════════
# MINIMAL SCRIPT: Section 1 (QC) + Section 2 (Lipid variation) ONLY
# ═══════════════════════════════════════════════════════════════════════════════
# Outputs:
#   - Supp Table S1: Class ratio statistics with jackknife
#   - Supp Table S2: Species-level CLR jackknife
#   - Supp Table S3: Species stability by class
#   - Supp Fig S1: QC panels (A-H)
#   - Supp Fig S2: Top variance lipids
#   - Supp Table S4: Top variance lipids table
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(vroom)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(forcats)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
  library(patchwork)
  library(viridis)
  library(grid)
  library(scales)
})

# Create output directories
dir.create("fig/main", recursive = TRUE, showWarnings = FALSE)
dir.create("fig/supp", recursive = TRUE, showWarnings = FALSE)
dir.create("table/supp", recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# Plot Theme
# ─────────────────────────────────────────────────────────────────────────────
plot_theme_compact <- theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(color = "black", size = 10),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

# ─────────────────────────────────────────────────────────────────────────────
# Color Palettes
# ─────────────────────────────────────────────────────────────────────────────
class_colors <- c(
  PC = "#00441B", PA = "#1B7837", PE = "#41AB5D", LPC = "#66C2A4", LPE = "#2CA25F",
  PG = "#78C679", PS = "#C2E699", DG = "#54278F", DGDG = "#F768A1", MG = "#8941ED",
  MGDG = "#FBB4D9", SQDG = "#9D4D6C", TG = "#ED804A"
)

condition_colors <- c(Control = "#440154FF", LowInput = "#FDE725FF")

valid_classes <- names(class_colors)

# ─────────────────────────────────────────────────────────────────────────────
# Helper Functions
# ─────────────────────────────────────────────────────────────────────────────
get_lipid_class <- function(x) {
  cls <- str_extract(x, "^[A-Za-z0-9]+(?=\\()")
  cls <- ifelse(is.na(cls), NA_character_, cls)
  ifelse(cls %in% valid_classes, cls, "Other")
}

save_fig <- function(p, filename, w = 10, h = 8, dpi = 300, subdir = "main") {
  p <- p + theme(
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 15)
  )
  path <- file.path("fig", subdir, filename)
  ggsave(path, plot = p, width = w, height = h, units = "in", dpi = dpi, bg = "white")
  message("Saved: ", path)
}

# ═══════════════════════════════════════════════════════════════════════════════
# PART 1: SUPPLEMENTARY TABLES S1, S2, S3 (Jackknife Analysis)
# ═══════════════════════════════════════════════════════════════════════════════

message("\n══════════════════════════════════════════════════════════════")
message("SUPPLEMENTARY TABLES S1, S2, S3: JACKKNIFE STABILITY ANALYSIS")
message("══════════════════════════════════════════════════════════════\n")

control_file <- "data/SPATS_fitted/non_normalized_intensities/Final_subset_control_all_lipids_fitted_phenotype_non_normalized.csv"
lowinput_file <- "data/SPATS_fitted/non_normalized_intensities/Final_subset_lowinput_all_lipids_fitted_phenotype_non_normalized.csv"

control <- vroom::vroom(control_file, show_col_types = FALSE) %>% dplyr::select(-c(2, 3, 4))
lowinput <- vroom::vroom(lowinput_file, show_col_types = FALSE) %>% dplyr::select(-c(2, 3, 4))
colnames(control)[1] <- "Compound_Name"
colnames(lowinput)[1] <- "Compound_Name"

message("Control samples: ", nrow(control))
message("LowInput samples: ", nrow(lowinput))

# Valid classes for class extraction
ratio_classes <- c("TG", "DG", "MG", "PC", "PE", "DGDG", "MGDG", "SQDG", "LPC", "LPE", "PG", "PA", "PS")
class_pat <- paste0("\\b(", paste(ratio_classes, collapse = "|"), ")\\b")

long_prep_global <- function(df) {
  df_long <- df %>%
    tidyr::pivot_longer(-c(Compound_Name, Condition), names_to = "Lipid", values_to = "Intensity") %>%
    dplyr::rename(Sample = Compound_Name)

  df_long <- df_long %>%
    dplyr::group_by(Sample) %>%
    dplyr::mutate(
      TIC = sum(Intensity, na.rm = TRUE),
      rel_abund = Intensity / TIC
    ) %>%
    dplyr::ungroup()

  df_long <- df_long %>%
    dplyr::group_by(Sample) %>%
    dplyr::mutate(
      minpos = dplyr::if_else(all(is.na(rel_abund) | rel_abund <= 0), NA_real_,
                              min(rel_abund[rel_abund > 0], na.rm = TRUE)),
      eps = dplyr::if_else(is.na(minpos), 0, minpos * 0.5)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(log_rel = log10(rel_abund + eps))

  df_long %>%
    dplyr::mutate(Class = stringr::str_extract(Lipid, class_pat)) %>%
    dplyr::filter(!is.na(Class)) %>%
    dplyr::group_by(Sample, Condition, Class) %>%
    dplyr::summarise(class_log = mean(log_rel, na.rm = TRUE), .groups = "drop")
}

combined <- dplyr::bind_rows(
  control %>% dplyr::mutate(Condition = "Control"),
  lowinput %>% dplyr::mutate(Condition = "LowInput")
)
long_all <- long_prep_global(combined)

wide_log <- long_all %>%
  tidyr::pivot_wider(names_from = Class, values_from = class_log, values_fill = NA_real_)

# Ratio lists
opls_ratios <- c("MG/SQDG", "MG/MGDG", "PS/SQDG", "DGDG/MG", "PG/SQDG", "DG/MG", "PE/SQDG",
                 "LPC/PS", "DGDG/PS", "DG/PS", "PC/PS", "LPE/SQDG", "LPC/MG", "PE/PS",
                 "MGDG/PG", "PG/PS", "LPC/LPE", "MGDG/PE", "MG/PG", "SQDG/TG", "PA/PS",
                 "LPE/MGDG", "PC/SQDG")
pcorr_ratios <- c("TG/PE", "TG/MGDG", "TG/MG", "TG/DG", "PS/PC", "PS/MGDG", "PS/DG",
                  "PG/MGDG", "PG/DGDG", "PE/PA")

ratio_specs <- tibble::tibble(
  ratio_label = c(opls_ratios, pcorr_ratios),
  ratio_set = c(rep("OPLS-DA ratios", length(opls_ratios)), rep("Partial correlation pairs", length(pcorr_ratios)))
) %>%
  tidyr::separate(ratio_label, into = c("Num", "Den"), sep = "/", remove = FALSE) %>%
  dplyr::mutate(Num = stringr::str_trim(Num), Den = stringr::str_trim(Den),
                ratio_id = paste0(Num, "_", Den))

ratio_tbl <- wide_log
for (i in seq_len(nrow(ratio_specs))) {
  num <- ratio_specs$Num[i]; den <- ratio_specs$Den[i]; rid <- ratio_specs$ratio_id[i]
  if (num %in% names(ratio_tbl) && den %in% names(ratio_tbl)) {
    ratio_tbl[[rid]] <- ratio_tbl[[num]] - ratio_tbl[[den]]
  } else {
    ratio_tbl[[rid]] <- NA_real_
  }
}

ratio_long <- ratio_tbl %>%
  dplyr::select(Sample, Condition, dplyr::all_of(ratio_specs$ratio_id)) %>%
  tidyr::pivot_longer(-c(Sample, Condition), names_to = "ratio_id", values_to = "Value") %>%
  dplyr::left_join(ratio_specs %>% dplyr::select(ratio_id, ratio_label, ratio_set), by = "ratio_id") %>%
  dplyr::filter(!is.na(Value) & is.finite(Value))

# ─────────────────────────────────────────────────────────────────────────────
# SUPPLEMENTARY TABLE S1: Class Ratio Statistics
# ─────────────────────────────────────────────────────────────────────────────
message("\n── Computing ratio statistics (Wilcoxon + Jackknife) ──")

one_ratio_tests <- function(df, ratio_id_val) {
  xi <- df %>% filter(Condition == "LowInput", ratio_id == ratio_id_val) %>% pull(Value)
  x0 <- df %>% filter(Condition == "Control", ratio_id == ratio_id_val) %>% pull(Value)
  ratio_lab <- df %>% filter(ratio_id == ratio_id_val) %>% pull(ratio_label) %>% unique() %>% first()

  if (length(xi) < 3 || length(x0) < 3) {
    return(tibble(Ratio = ratio_lab, n_C = length(x0), n_LI = length(xi),
                  median_C = NA_real_, median_LI = NA_real_, effect_log10 = NA_real_,
                  p_wilcox = NA_real_, jackknife_stability = NA_real_))
  }

  wt <- wilcox.test(xi, x0, alternative = "two.sided", exact = FALSE)
  d_full <- sign(median(xi) - median(x0))
  N <- length(xi) + length(x0)
  keep <- logical(N)
  for (i in seq_len(N)) {
    xi2 <- xi; x0_2 <- x0
    if (i <= length(xi)) xi2 <- xi2[-i] else x0_2 <- x0_2[-(i - length(xi))]
    keep[i] <- sign(median(xi2) - median(x0_2)) == d_full
  }
  stab <- mean(keep)

  tibble(Ratio = ratio_lab, n_C = length(x0), n_LI = length(xi),
         median_C = round(median(x0), 4), median_LI = round(median(xi), 4),
         effect_log10 = round(median(xi) - median(x0), 4),
         p_wilcox = wt$p.value, jackknife_stability = round(stab, 3))
}

ratios_to_test <- unique(ratio_long$ratio_id)
ratio_stats <- purrr::map_dfr(ratios_to_test, function(r) one_ratio_tests(ratio_long, r)) %>%
  dplyr::mutate(
    p_adj_BH = p.adjust(p_wilcox, method = "BH"),
    effect_fc = round(10^effect_log10, 2),
    direction = dplyr::case_when(effect_log10 > 0 ~ "LI > C", effect_log10 < 0 ~ "LI < C", TRUE ~ "No change"),
    significance = dplyr::case_when(p_adj_BH < 0.001 ~ "***", p_adj_BH < 0.01 ~ "**", p_adj_BH < 0.05 ~ "*", TRUE ~ "ns")
  ) %>%
  dplyr::arrange(p_adj_BH)

ratio_stats_out <- ratio_stats %>%
  dplyr::mutate(p_wilcox = format(p_wilcox, digits = 3, scientific = TRUE),
                p_adj_BH = format(p_adj_BH, digits = 3, scientific = TRUE)) %>%
  dplyr::select(Ratio, n_C, n_LI, median_C, median_LI, effect_log10, effect_fc, direction,
                p_wilcox, p_adj_BH, significance, jackknife_stability)

write.csv(ratio_stats_out, "table/supp/SuppTable_S1_Ratio_Statistics.csv", row.names = FALSE)
message("Saved: table/supp/SuppTable_S1_Ratio_Statistics.csv")

# ─────────────────────────────────────────────────────────────────────────────
# SUPPLEMENTARY TABLE S2: Species-level Jackknife (CLR)
# ─────────────────────────────────────────────────────────────────────────────
message("\n── Computing SPECIES-level jackknife stability (CLR-transformed) ──")

control_species <- colnames(control)[-1]
lowinput_species <- colnames(lowinput)[-1]
common_species <- intersect(control_species, lowinput_species)
# Keep only true lipid species columns (e.g., "PC(16:0/18:1)"),
# matching the species universe used in SuppTable_S6 summaries.
common_species <- common_species[stringr::str_detect(common_species, "\\(")]
message(sprintf("  Found %d common lipid species", length(common_species)))

species_long <- dplyr::bind_rows(
  control %>% dplyr::mutate(Condition = "Control"),
  lowinput %>% dplyr::mutate(Condition = "LowInput")
) %>%
  tidyr::pivot_longer(-c(Compound_Name, Condition), names_to = "Lipid", values_to = "Intensity") %>%
  dplyr::rename(Sample = Compound_Name) %>%
  dplyr::filter(Lipid %in% common_species) %>%
  dplyr::group_by(Sample) %>%
  dplyr::mutate(
    minpos = min(Intensity[Intensity > 0], na.rm = TRUE),
    Intensity_nozero = dplyr::if_else(is.na(Intensity) | Intensity <= 0, minpos * 0.5, Intensity)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Sample) %>%
  dplyr::mutate(
    log_intensity = log(Intensity_nozero),
    geom_mean_log = mean(log_intensity, na.rm = TRUE),
    CLR = log_intensity - geom_mean_log
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Class = dplyr::case_when(
      stringr::str_detect(Lipid, "^TG\\(") ~ "TG", stringr::str_detect(Lipid, "^DG\\(") ~ "DG",
      stringr::str_detect(Lipid, "^MG\\(") ~ "MG", stringr::str_detect(Lipid, "^PC\\(") ~ "PC",
      stringr::str_detect(Lipid, "^PE\\(") ~ "PE", stringr::str_detect(Lipid, "^PG\\(") ~ "PG",
      stringr::str_detect(Lipid, "^PA\\(") ~ "PA", stringr::str_detect(Lipid, "^PS\\(") ~ "PS",
      stringr::str_detect(Lipid, "^DGDG\\(") ~ "DGDG", stringr::str_detect(Lipid, "^MGDG\\(") ~ "MGDG",
      stringr::str_detect(Lipid, "^SQDG\\(") ~ "SQDG", stringr::str_detect(Lipid, "^LPC\\(") ~ "LPC",
      stringr::str_detect(Lipid, "^LPE\\(") ~ "LPE", stringr::str_detect(Lipid, "^FA\\(") ~ "FA",
      stringr::str_detect(Lipid, "^Cer\\(") ~ "Cer", stringr::str_detect(Lipid, "^GalCer\\(") ~ "GalCer",
      stringr::str_detect(Lipid, "^SM\\(") ~ "SM", stringr::str_detect(Lipid, "^AEG\\(") ~ "AEG",
      TRUE ~ stringr::str_extract(Lipid, "^[A-Za-z]+")
    )
  ) %>%
  dplyr::filter(!is.na(CLR), is.finite(CLR))

one_species_tests <- function(df, lipid_name) {
  xi <- df %>% dplyr::filter(Condition == "LowInput", Lipid == lipid_name) %>% dplyr::pull(CLR)
  x0 <- df %>% dplyr::filter(Condition == "Control", Lipid == lipid_name) %>% dplyr::pull(CLR)
  lipid_class <- df %>% dplyr::filter(Lipid == lipid_name) %>% dplyr::pull(Class) %>% unique() %>% first()

  if (length(xi) < 3 || length(x0) < 3) {
    return(tibble::tibble(Lipid = lipid_name, Class = lipid_class, n_C = length(x0), n_LI = length(xi),
                          median_CLR_C = NA_real_, median_CLR_LI = NA_real_, effect_CLR = NA_real_,
                          p_wilcox = NA_real_, jackknife_stability = NA_real_))
  }

  wt <- tryCatch(wilcox.test(xi, x0, alternative = "two.sided", exact = FALSE), error = function(e) list(p.value = NA_real_))
  d_full <- sign(median(xi) - median(x0))
  N <- length(xi) + length(x0)
  keep <- logical(N)
  for (i in seq_len(N)) {
    xi2 <- xi; x0_2 <- x0
    if (i <= length(xi)) xi2 <- xi2[-i] else x0_2 <- x0_2[-(i - length(xi))]
    keep[i] <- sign(median(xi2) - median(x0_2)) == d_full
  }
  stab <- mean(keep)

  tibble::tibble(Lipid = lipid_name, Class = lipid_class, n_C = length(x0), n_LI = length(xi),
                 median_CLR_C = round(median(x0), 4), median_CLR_LI = round(median(xi), 4),
                 effect_CLR = round(median(xi) - median(x0), 4),
                 p_wilcox = wt$p.value, jackknife_stability = round(stab, 3))
}

all_species <- unique(species_long$Lipid)
message(sprintf("  Testing %d individual species...", length(all_species)))

species_stats <- purrr::map_dfr(all_species, function(sp) one_species_tests(species_long, sp)) %>%
  dplyr::mutate(
    p_adj_BH = p.adjust(p_wilcox, method = "BH"),
    effect_fc = round(exp(effect_CLR), 2),
    direction = dplyr::case_when(effect_CLR > 0 ~ "LI > C", effect_CLR < 0 ~ "LI < C", TRUE ~ "No change"),
    significance = dplyr::case_when(p_adj_BH < 0.001 ~ "***", p_adj_BH < 0.01 ~ "**", p_adj_BH < 0.05 ~ "*", TRUE ~ "ns")
  ) %>%
  dplyr::arrange(p_adj_BH)

species_stats_out <- species_stats %>%
  dplyr::mutate(p_wilcox = format(p_wilcox, digits = 3, scientific = TRUE),
                p_adj_BH = format(p_adj_BH, digits = 3, scientific = TRUE)) %>%
  dplyr::select(Lipid, Class, n_C, n_LI, median_CLR_C, median_CLR_LI, effect_CLR, effect_fc, direction,
                p_wilcox, p_adj_BH, significance, jackknife_stability)

write.csv(species_stats_out, "table/supp/SuppTable_S2_Species_Jackknife.csv", row.names = FALSE)
message("Saved: table/supp/SuppTable_S2_Species_Jackknife.csv")

# ─────────────────────────────────────────────────────────────────────────────
# SUPPLEMENTARY TABLE S3: Summary by Class
# ─────────────────────────────────────────────────────────────────────────────
class_stability_summary <- species_stats %>%
  dplyr::group_by(Class) %>%
  dplyr::summarise(
    n_species = dplyr::n(),
    n_significant = sum(p_adj_BH < 0.05, na.rm = TRUE),
    n_stable_1.0 = sum(jackknife_stability == 1, na.rm = TRUE),
    n_stable_0.95 = sum(jackknife_stability >= 0.95, na.rm = TRUE),
    mean_stability = round(mean(jackknife_stability, na.rm = TRUE), 3),
    min_stability = round(min(jackknife_stability, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_species))

write.csv(class_stability_summary, "table/supp/SuppTable_S3_Species_Stability_by_Class.csv", row.names = FALSE)
message("Saved: table/supp/SuppTable_S3_Species_Stability_by_Class.csv")

message("\n Tables S1, S2, S3 complete!")

# ═══════════════════════════════════════════════════════════════════════════════
# PART 2: SUPPLEMENTARY FIGURE S1 (QC)
# ═══════════════════════════════════════════════════════════════════════════════

message("\n══════════════════════════════════════════════════════════════")
message("SUPPLEMENTARY FIGURE S1: QC - RUN ORDER, SERRF RSD, PCA, SPATS")
message("══════════════════════════════════════════════════════════════\n")

classify_qc_injection_type <- function(x) {
  dplyr::case_when(
    str_detect(x, regex("QC", ignore_case = TRUE)) ~ "QC",
    str_detect(x, regex("CHECK", ignore_case = TRUE)) ~ "Check",
    str_detect(x, regex("ISTD", ignore_case = TRUE)) ~ "ISTD",
    str_detect(x, regex("InjBL|Blank", ignore_case = TRUE)) ~ "Blank",
    str_detect(x, regex("_PI\\d+|^S\\d+_", ignore_case = TRUE)) ~ "Sample",
    TRUE ~ "Other"
  )
}

create_run_order_plot <- function(file_candidates, condition_name) {
  file_path <- file_candidates[file.exists(file_candidates)][1]
  if (is.na(file_path)) {
    return(ggplot() + theme_void() + ggtitle(paste0(condition_name, ": Data not available")))
  }

  dat <- suppressMessages(vroom(file_path, show_col_types = FALSE, delim = ",", progress = FALSE))
  run_cols <- names(dat)[str_detect(names(dat), "Run\\d+")]
  if (length(run_cols) == 0) {
    return(ggplot() + theme_void() + ggtitle(paste0(condition_name, ": No run columns found")))
  }

  run_meta <- tibble(Injection = run_cols) %>%
    mutate(Run = as.integer(str_extract(Injection, "(?<=Run)\\d+")),
           Type = classify_qc_injection_type(Injection),
           QC_Label = dplyr::case_when(
             Type == "QC" & str_detect(Injection, regex("Pre", ignore_case = TRUE)) ~ "QC Pre",
             Type == "QC" & str_detect(Injection, regex("Post", ignore_case = TRUE)) ~ "QC Post",
             TRUE ~ NA_character_
           )) %>%
    mutate(Run = ifelse(is.na(Run), row_number(), Run)) %>%
    arrange(Run)

  run_matrix <- dat %>%
    dplyr::select(all_of(run_meta$Injection)) %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(.x)))) %>%
    as.matrix()
  mode(run_matrix) <- "numeric"
  keep_rows <- rowSums(is.finite(run_matrix), na.rm = TRUE) > 0
  run_matrix <- run_matrix[keep_rows, , drop = FALSE]

  plot_df <- run_meta %>%
    mutate(TIC = colSums(run_matrix, na.rm = TRUE)) %>%
    filter(Type %in% c("Sample", "Check", "QC", "Blank", "ISTD")) %>%
    arrange(Run)

  qc_tic <- plot_df %>% filter(Type == "QC") %>% pull(TIC)
  qc_cv <- if (length(qc_tic) >= 2 && mean(qc_tic, na.rm = TRUE) > 0) {
    100 * sd(qc_tic, na.rm = TRUE) / mean(qc_tic, na.rm = TRUE)
  } else NA_real_

  type_counts <- plot_df %>% count(Type) %>% mutate(label = paste0(Type, " n=", n)) %>% pull(label) %>% paste(collapse = " | ")
  subtitle_txt <- if (is.finite(qc_cv)) paste0(type_counts, sprintf(" | QC TIC CV=%.1f%%", qc_cv)) else type_counts

  ggplot(plot_df, aes(x = Run, y = TIC, color = Type)) +
    geom_line(data = plot_df %>% filter(Type == "Sample"), aes(group = 1), alpha = 0.35, linewidth = 0.35) +
    geom_point(size = 1.9, alpha = 0.85) +
    ggrepel::geom_text_repel(data = plot_df %>% filter(!is.na(QC_Label)), aes(label = QC_Label),
                             color = "black", size = 3, seed = 42, box.padding = 0.25) +
    scale_color_manual(values = c("Sample" = "#440154FF", "Check" = "#FDE725FF", "QC" = "#E64B35FF",
                                  "Blank" = "#7A7A7AFF", "ISTD" = "#00A087FF")) +
    scale_y_continuous(labels = scales::scientific) +
    labs(title = paste0(condition_name, " - Run Order TIC"), subtitle = subtitle_txt, x = "Run Order", y = "TIC") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9.5, color = "grey35"),
          axis.title = element_text(face = "bold"), axis.text = element_text(color = "black"),
          legend.position = "bottom", panel.grid.minor = element_blank(),
          axis.line.x = element_line(color = "black", linewidth = 0.5),
          axis.line.y = element_line(color = "black", linewidth = 0.5))
}

suppfig_s1a <- create_run_order_plot(c("data/raw_lipid_intensities/SetA_lipid_FLO2019Control.csv"), "A) Control")
suppfig_s1b <- create_run_order_plot(c("data/raw_lipid_intensities/SetB_lipid_FLO2022_lowP.csv"), "B) LowInput")
message("Created: S1A & S1B (Run Order)")

# S1C & S1D: SERRF RSD
create_qc_rsd_bar_plot <- function(rsd_path, panel_title) {
  if (!file.exists(rsd_path)) {
    return(ggplot() + theme_void() + ggtitle(paste0(panel_title, ": RSD file missing")))
  }
  rsd <- suppressMessages(vroom(rsd_path, show_col_types = FALSE))
  if (!all(c("QC_none", "QC_SERRF") %in% names(rsd))) {
    return(ggplot() + theme_void() + ggtitle(paste0(panel_title, ": Invalid RSD format")))
  }
  raw_rsd <- as.numeric(rsd$QC_none) * 100
  serrf_rsd <- as.numeric(rsd$QC_SERRF) * 100

  plot_df <- tibble(
    Stage = factor(c("Before SERRF", "After SERRF"), levels = c("Before SERRF", "After SERRF")),
    Median_RSD = c(median(raw_rsd, na.rm = TRUE), median(serrf_rsd, na.rm = TRUE)),
    Pct_lt_30 = c(mean(raw_rsd < 30, na.rm = TRUE) * 100, mean(serrf_rsd < 30, na.rm = TRUE) * 100)
  )

  ggplot(plot_df, aes(x = Stage, y = Median_RSD, fill = Stage)) +
    geom_col(width = 0.6, color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%", Median_RSD)), vjust = -0.35, size = 3.8, fontface = "bold") +
    scale_fill_manual(values = c("Before SERRF" = "#9E9E9E", "After SERRF" = "#009E73")) +
    labs(title = panel_title, subtitle = paste0("QC-RSD < 30%: ", sprintf("%.1f%%", plot_df$Pct_lt_30[1]), " -> ",
                                                 sprintf("%.1f%%", plot_df$Pct_lt_30[2])),
         x = NULL, y = "Median QC-RSD (%)") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9.5, color = "grey35"),
          axis.title = element_text(face = "bold"), axis.text = element_text(color = "black"),
          legend.position = "none", panel.grid.minor = element_blank(),
          axis.line.x = element_line(color = "black", linewidth = 0.5),
          axis.line.y = element_line(color = "black", linewidth = 0.5))
}

suppfig_s1c <- create_qc_rsd_bar_plot("data/SEERF_normalized_intensities/A_QC-RSDs.csv", "C) Control: QC-RSD After SERRF")
suppfig_s1d <- create_qc_rsd_bar_plot("data/SEERF_normalized_intensities/B_QC-RSDs.csv", "D) LowInput: QC-RSD After SERRF")
message("Created: S1C & S1D (SERRF RSD)")

# S1E & S1F: PCA Comparison (Raw vs SERRF)
load_matrix_with_meta <- function(file_path, feature_id_candidates = c("row ID", "label")) {
  dat <- suppressMessages(vroom(file_path, show_col_types = FALSE, delim = ",", progress = FALSE))
  run_cols <- names(dat)[str_detect(names(dat), "Run\\d+")]
  if (length(run_cols) == 0) {
    stop("No run columns found in: ", file_path)
  }

  feature_col <- feature_id_candidates[feature_id_candidates %in% names(dat)][1]
  if (is.na(feature_col)) {
    feature_id <- as.character(seq_len(nrow(dat)))
  } else {
    feature_id <- as.character(dat[[feature_col]])
  }

  mat <- dat %>%
    dplyr::select(all_of(run_cols)) %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(.x)))) %>%
    as.matrix()
  mode(mat) <- "numeric"

  keep_feat <- !is.na(feature_id) & nzchar(feature_id) & !duplicated(feature_id)
  mat <- mat[keep_feat, , drop = FALSE]
  feature_id <- feature_id[keep_feat]

  keep_rows <- rowSums(is.finite(mat), na.rm = TRUE) > 0
  mat <- mat[keep_rows, , drop = FALSE]
  feature_id <- feature_id[keep_rows]

  meta <- tibble(Injection = run_cols) %>%
    mutate(
      Run = as.integer(str_extract(Injection, "(?<=Run)\\d+")),
      Type = classify_qc_injection_type(Injection)
    ) %>%
    mutate(Run = ifelse(is.na(Run), row_number(), Run)) %>%
    arrange(Run)

  mat <- mat[, meta$Injection, drop = FALSE]
  list(mat = mat, feature_id = feature_id, meta = meta, file = file_path)
}

create_pca_raw_vs_serrf <- function(raw_candidates, serrf_path, panel_title) {
  raw_path <- raw_candidates[file.exists(raw_candidates)][1]
  if (is.na(raw_path) || !file.exists(serrf_path)) {
    return(ggplot() + theme_void() +
             ggtitle(paste0(panel_title, ": Raw/SERRF files not available")))
  }

  raw <- load_matrix_with_meta(raw_path, c("row ID", "label"))
  serrf <- load_matrix_with_meta(serrf_path, c("label", "row ID"))

  allowed_types <- c("Sample", "Check", "QC", "ISTD")
  raw_meta <- raw$meta %>% filter(Type %in% allowed_types)
  serrf_meta <- serrf$meta %>% filter(Type %in% allowed_types)

  matched_inj <- raw_meta %>%
    dplyr::select(Injection_raw = Injection, Run, Type) %>%
    inner_join(
      serrf_meta %>% dplyr::select(Injection_serrf = Injection, Run, Type),
      by = c("Run", "Type")
    ) %>%
    arrange(Run)

  matched_feat <- intersect(raw$feature_id, serrf$feature_id)
  if (nrow(matched_inj) < 10 || length(matched_feat) < 100) {
    return(ggplot() + theme_void() +
             ggtitle(paste0(panel_title, ": Insufficient matched data")))
  }

  raw_idx <- match(matched_feat, raw$feature_id)
  serrf_idx <- match(matched_feat, serrf$feature_id)

  raw_mat <- raw$mat[raw_idx, matched_inj$Injection_raw, drop = FALSE]
  serrf_mat <- serrf$mat[serrf_idx, matched_inj$Injection_serrf, drop = FALSE]

  raw_x <- t(log10(pmax(raw_mat, 0) + 1))
  serrf_x <- t(log10(pmax(serrf_mat, 0) + 1))
  combo_x <- rbind(raw_x, serrf_x)

  for (j in seq_len(ncol(combo_x))) {
    v <- combo_x[, j]
    v[!is.finite(v)] <- NA_real_
    med <- median(v, na.rm = TRUE)
    if (!is.finite(med)) med <- 0
    v[is.na(v)] <- med
    combo_x[, j] <- v
  }

  sds <- apply(combo_x, 2, sd)
  keep_var <- is.finite(sds) & sds > 0
  combo_x <- combo_x[, keep_var, drop = FALSE]

  if (ncol(combo_x) < 2) {
    return(ggplot() + theme_void() +
             ggtitle(paste0(panel_title, ": PCA not available")))
  }

  pca <- prcomp(combo_x, center = TRUE, scale. = TRUE)
  var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
  n_match <- nrow(matched_inj)

  pca_df <- tibble(
    Stage = rep(c("Raw", "SERRF"), each = n_match),
    Run = rep(matched_inj$Run, times = 2),
    Type = rep(matched_inj$Type, times = 2),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2]
  )

  ggplot(pca_df, aes(x = PC1, y = PC2, color = Stage, shape = Type)) +
    geom_point(size = 2.1, alpha = 0.85) +
    stat_ellipse(
      data = pca_df %>% filter(Type == "QC"),
      aes(group = Stage, color = Stage),
      type = "norm",
      level = 0.8,
      linewidth = 0.7,
      linetype = "dashed",
      show.legend = FALSE,
      na.rm = TRUE
    ) +
    scale_color_manual(values = c("Raw" = "#7A7A7A", "SERRF" = "#009E73")) +
    scale_shape_manual(values = c("Sample" = 16, "Check" = 17, "QC" = 15, "ISTD" = 18)) +
    labs(
      title = panel_title,
      subtitle = paste0(
        "Matched injections: ", n_match,
        " | matched features: ", ncol(combo_x),
        " | dashed ellipse = QC cluster"
      ),
      x = paste0("PC1 (", round(var_exp[1] * 100, 1), "%)"),
      y = paste0("PC2 (", round(var_exp[2] * 100, 1), "%)")
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9.5, color = "grey35"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5)
    )
}

suppfig_s1e <- create_pca_raw_vs_serrf(
  c(
    "data/SetA_lipid_FLO2019Control.csv",
    "data/raw_lipid_intensities/SetA_lipid_FLO2019Control.csv",
    "data/raw_lipid_intensities/A_final_lipids.csv",
    "data/raw_lipid_intensities/A_cleaned_lipids_unsummed.csv",
    "data/raw_lipid_intensities/A_cleaned_lipids.csv"
  ),
  "data/SEERF_normalized_intensities/A_normalized_SERRF.csv",
  "E) Control PCA: Raw vs SERRF"
)

suppfig_s1f <- create_pca_raw_vs_serrf(
  c(
    "data/SetB_lipid_FLO2022_lowP.csv",
    "data/raw_lipid_intensities/SetB_lipid_FLO2022_lowP.csv",
    "data/raw_lipid_intensities/B_final_lipids.csv",
    "data/raw_lipid_intensities/B_cleaned_lipids_unsummed.csv",
    "data/raw_lipid_intensities/B_cleaned_lipids.csv"
  ),
  "data/SEERF_normalized_intensities/B_normalized_SERRF.csv",
  "F) LowInput PCA: Raw vs SERRF"
)
message("Created: S1E & S1F (PCA Raw vs SERRF)")

# S1G & S1H: SpATS Residual Spatial Maps
prepare_spats_input <- function(raw_summed_path, fieldmap_path) {
  raw <- read.csv(raw_summed_path, stringsAsFactors = FALSE, check.names = FALSE) %>%
    dplyr::select(-c(2:13))

  pheno_mat <- as.data.frame(t(raw), stringsAsFactors = FALSE)
  colnames(pheno_mat) <- pheno_mat[1, ]
  pheno_mat <- pheno_mat[-1, ]

  annot_rows <- c("Class", "SubClass", "Sub_subclass")
  data_rows <- setdiff(rownames(pheno_mat), annot_rows)
  pheno_num <- pheno_mat[data_rows, , drop = FALSE]
  pheno_num[] <- lapply(pheno_num, as.numeric)

  pheno <- pheno_num %>% as_tibble(rownames = "LineRaw")

  fieldmap <- read.csv(fieldmap_path, header = FALSE, stringsAsFactors = FALSE)
  map_long <- fieldmap %>%
    mutate(row = row_number()) %>%
    pivot_longer(cols = starts_with("V"), names_to = "col", values_to = "MapRaw") %>%
    mutate(
      col = as.integer(str_remove(col, "^V")),
      MapRaw = as.character(MapRaw)
    ) %>%
    filter(!is.na(MapRaw), MapRaw != "X")

  canonical_id <- function(x) {
    x %>%
      str_to_upper() %>%
      str_replace("^PI_", "PI") %>%
      str_replace_all("[^A-Z0-9]", "")
  }

  pheno <- pheno %>% mutate(Key = canonical_id(LineRaw))
  map_long <- map_long %>% mutate(Key = canonical_id(MapRaw))

  map_long %>%
    inner_join(pheno, by = "Key") %>%
    dplyr::select(LineRaw = LineRaw, row, col, everything())
}

create_spats_residual_plot <- function(raw_summed_path, fieldmap_path, lipid_name, panel_title) {
  if (!file.exists(raw_summed_path) || !file.exists(fieldmap_path)) {
    return(ggplot() + theme_void() +
             ggtitle(paste0(panel_title, ": input file missing")))
  }

  if (!requireNamespace("SpATS", quietly = TRUE)) {
    return(ggplot() + theme_void() +
             ggtitle(paste0(panel_title, ": SpATS package not installed")))
  }

  dat <- prepare_spats_input(raw_summed_path, fieldmap_path)
  if (!lipid_name %in% names(dat)) {
    return(ggplot() + theme_void() +
             ggtitle(paste0(panel_title, ": lipid not found - ", lipid_name)))
  }

  fit_dat <- dat %>%
    dplyr::select(LineRaw, row, col, all_of(lipid_name)) %>%
    dplyr::rename(Trait = all_of(lipid_name)) %>%
    filter(is.finite(Trait))

  if (nrow(fit_dat) < 20) {
    return(ggplot() + theme_void() +
             ggtitle(paste0(panel_title, ": insufficient points for SpATS")))
  }

  fit_dat <- fit_dat %>% mutate(LineRaw = as.factor(LineRaw))
  nseg_col <- max(2, round(dplyr::n_distinct(fit_dat$col) * 0.4))
  nseg_row <- max(2, round(dplyr::n_distinct(fit_dat$row) * 0.4))

  spats_fit <- tryCatch(
    SpATS::SpATS(
      response = "Trait",
      spatial = ~ SpATS::SAP(col, row, nseg = c(nseg_col, nseg_row), degree = 3, pord = 2),
      genotype = "LineRaw",
      data = fit_dat,
      genotype.as.random = TRUE,
      control = list(tolerance = 1e-03, maxit = 500)
    ),
    error = function(e) NULL
  )

  if (is.null(spats_fit)) {
    return(ggplot() + theme_void() +
             ggtitle(paste0(panel_title, ": SpATS model failed")))
  }

  resid_df <- fit_dat %>% mutate(Residual = as.numeric(residuals(spats_fit)))

  ggplot(resid_df, aes(x = col, y = row, fill = Residual)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient2(
      low = "#3B4CC0",
      mid = "white",
      high = "#B40426",
      midpoint = 0,
      name = "Residual"
    ) +
    scale_y_reverse() +
    coord_fixed() +
    labs(
      title = panel_title,
      subtitle = paste0("Lipid: ", lipid_name, " | residuals after SpATS spatial fit"),
      x = "Field column",
      y = "Field row"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9.5, color = "grey35"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.grid = element_blank(),
      legend.position = "right",
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5)
    )
}

target_spats_lipid <- "MGDG(18:3/18:3)"

suppfig_s1g <- create_spats_residual_plot(
  "data/summed_lipid_intensities/A_final_summed_lipids.csv",
  "data/fieldmap/fieldmap_control.csv",
  target_spats_lipid,
  "G) Control SpATS Residual Map"
)

suppfig_s1h <- create_spats_residual_plot(
  "data/summed_lipid_intensities/B_final_summed_lipids.csv",
  "data/fieldmap/fieldmap_lowinput.csv",
  target_spats_lipid,
  "H) LowInput SpATS Residual Map"
)
message("Created: S1G & S1H (SpATS residual maps)")

suppfig_s1_qc <- (suppfig_s1a + suppfig_s1b) / (suppfig_s1c + suppfig_s1d) /
  (suppfig_s1e + suppfig_s1f) / (suppfig_s1g + suppfig_s1h) +
  plot_layout(heights = c(1.05, 0.95, 1, 1))

save_fig(suppfig_s1_qc, "SuppFig_S1_QC_RunOrder_SERRF_PCA_SpATS.png", w = 16, h = 20, subdir = "supp")
save_fig(suppfig_s1_qc, "SuppFig_S1_QC_RunOrder_SERRF_PCA_SpATS.pdf", w = 16, h = 20, subdir = "supp")
message(" Supplementary Figure S1 (QC) complete!")

# ═══════════════════════════════════════════════════════════════════════════════
# PART 3: SUPPLEMENTARY FIGURE S2 + TABLE S4 (Top Variance)
# ═══════════════════════════════════════════════════════════════════════════════

message("\n══════════════════════════════════════════════════════════════")
message("SUPPLEMENTARY FIGURE S2 + TABLE S4: TOP VARIANCE LIPIDS")
message("══════════════════════════════════════════════════════════════\n")

# Load raw data for TIC normalization
control_raw <- vroom(control_file, show_col_types = FALSE) %>% dplyr::select(-c(2, 3, 4))
low_raw <- vroom(lowinput_file, show_col_types = FALSE) %>% dplyr::select(-c(2, 3, 4))
colnames(control_raw)[1] <- "Compound_Name"
colnames(low_raw)[1] <- "Compound_Name"
control_raw <- control_raw %>% mutate(Condition = "Control")
low_raw <- low_raw %>% mutate(Condition = "LowInput")

# TIC normalize
tic_normalize <- function(df) {
  X <- df %>% dplyr::select(-Compound_Name, -Condition) %>% dplyr::select(where(is.numeric)) %>% as.matrix()
  rs <- rowSums(X, na.rm = TRUE)
  rs[rs == 0 | is.na(rs)] <- NA_real_
  Xpct <- sweep(X, 1, rs, FUN = "/") * 100
  Xpct[is.na(Xpct)] <- 0
  bind_cols(df %>% dplyr::select(Compound_Name, Condition), as.data.frame(Xpct))
}

control_tic <- tic_normalize(control_raw)
low_tic <- tic_normalize(low_raw)

find_high_variance_lipids <- function(df_tic, condition_name, top_n = 10) {
  lipid_cols <- setdiff(names(df_tic), c("Compound_Name", "Condition"))
  lipid_cols <- lipid_cols[grepl("\\(", lipid_cols)]

  X <- df_tic %>% dplyr::select(all_of(lipid_cols)) %>%
    mutate(across(everything(), ~suppressWarnings(as.numeric(.x))))

  X %>%
    summarise(across(everything(), ~var(.x, na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "Lipid", values_to = "Variance") %>%
    arrange(desc(Variance)) %>%
    mutate(Rank = row_number(), Condition = condition_name, Class = get_lipid_class(Lipid)) %>%
    filter(Rank <= top_n)
}

top_var_control <- find_high_variance_lipids(control_tic, "Control", top_n = 10)
top_var_low <- find_high_variance_lipids(low_tic, "LowInput", top_n = 10)

message("Top variance lipids: ", paste(head(top_var_control$Lipid, 3), collapse = ", "), ", ...")

# S2A: Control
suppfig_s2a <- ggplot(top_var_control, aes(x = reorder(Lipid, Variance), y = Variance, fill = Class)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) + coord_flip() +
  scale_fill_manual(values = class_colors, na.value = "grey70") +
  labs(title = "A) Control - Top 10 High Variance Lipids", x = NULL, y = "Variance") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(color = "black", size = 9), plot.title = element_text(face = "bold"),
        legend.position = "right", panel.grid.major.y = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

# S2B: LowInput
suppfig_s2b <- ggplot(top_var_low, aes(x = reorder(Lipid, Variance), y = Variance, fill = Class)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) + coord_flip() +
  scale_fill_manual(values = class_colors, na.value = "grey70") +
  labs(title = "B) LowInput - Top 10 High Variance Lipids", x = NULL, y = "Variance") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(color = "black", size = 9), plot.title = element_text(face = "bold"),
        legend.position = "right", panel.grid.major.y = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

# S2C: Overlap scatter
overlap_lipids <- intersect(top_var_control$Lipid, top_var_low$Lipid)
var_wide <- full_join(
  top_var_control %>% dplyr::select(Lipid, Variance, Class) %>% dplyr::rename(Var_Ctrl = Variance),
  top_var_low %>% dplyr::select(Lipid, Variance) %>% dplyr::rename(Var_Low = Variance),
  by = "Lipid"
) %>% mutate(Class = ifelse(is.na(Class), get_lipid_class(Lipid), Class), InBoth = Lipid %in% overlap_lipids)

suppfig_s2c <- ggplot(var_wide, aes(x = Var_Ctrl, y = Var_Low)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = Class, shape = InBoth), size = 3, alpha = 0.8) +
  geom_text_repel(aes(label = Lipid), size = 2.8, max.overlaps = 15) +
  scale_color_manual(values = class_colors, na.value = "grey70") +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 17)) +
  labs(title = "C) Variance Comparison", subtitle = sprintf("Overlap: %d lipids", length(overlap_lipids)),
       x = "Variance (Control)", y = "Variance (LowInput)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "right",
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

# S2D: Class distribution
combined_var <- bind_rows(top_var_control %>% mutate(Condition = "Control"),
                          top_var_low %>% mutate(Condition = "LowInput"))
class_counts <- combined_var %>% count(Condition, Class) %>%
  mutate(Condition = factor(Condition, levels = c("Control", "LowInput")))

suppfig_s2d <- ggplot(class_counts, aes(x = Class, y = n, fill = Condition)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = condition_colors) +
  labs(title = "D) Class Distribution", x = "Lipid Class", y = "Count") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top",
        panel.grid.major.x = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

# Combine
suppfig_s2_top <- suppfig_s2a + suppfig_s2b + plot_layout(guides = "collect") & theme(legend.position = "right")
suppfig_s2_bottom <- suppfig_s2c + suppfig_s2d + plot_layout(widths = c(1.2, 1))
suppfig_s2 <- suppfig_s2_top / suppfig_s2_bottom + plot_layout(heights = c(1, 1))

save_fig(suppfig_s2, "SuppFig_S2_TopVariance_Lipids.png", w = 14, h = 12, subdir = "supp")
save_fig(suppfig_s2, "SuppFig_S2_TopVariance_Lipids.pdf", w = 14, h = 12, subdir = "supp")

# Table S4
top_var_combined <- bind_rows(top_var_control %>% mutate(Condition = "Control"),
                               top_var_low %>% mutate(Condition = "LowInput")) %>%
  dplyr::select(Condition, Rank, Lipid, Class, Variance) %>% arrange(Condition, Rank)

write.csv(top_var_combined, "table/supp/SuppTable_S4_TopVariance_Lipids.csv", row.names = FALSE)
message("Saved: table/supp/SuppTable_S4_TopVariance_Lipids.csv")

message("\n Supplementary Figure S2 + Table S4 complete!")

message("\n══════════════════════════════════════════════════════════════")
message("ALL DONE! Generated:")
message("  - table/supp/SuppTable_S1_Ratio_Statistics.csv")
message("  - table/supp/SuppTable_S2_Species_Jackknife.csv")
message("  - table/supp/SuppTable_S3_Species_Stability_by_Class.csv")
message("  - fig/supp/SuppFig_S1_QC_RunOrder_SERRF_PCA_SpATS.png/pdf")
message("  - fig/supp/SuppFig_S2_TopVariance_Lipids.png/pdf")
message("  - table/supp/SuppTable_S4_TopVariance_Lipids.csv")
message("══════════════════════════════════════════════════════════════\n")
