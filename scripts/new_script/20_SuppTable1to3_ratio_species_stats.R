# ==============================================================================
# Supplementary Tables S1-S3 -- class-ratio statistics, species-level jackknife,
# and per-class stability summary.
#
# These three tables are produced by make_suppfig_s9_lipid_ratios() inside
# 21_high_variance_lipids.R. That script is a 6,000-line monolith that also
# rebuilds Figure 1, the S1 QC panel, OPLS and several supplementary figures
# under an OLDER numbering scheme (its S3 is the species counts, its S4 the
# PCA), so running it whole overwrites current figures with mislabelled ones.
#
# This script is that function lifted VERBATIM, with only its two SuppFig_S7
# writes commented out (that figure is written by 28_class_logratio_stats.R).
# The statistics are unchanged from the original.
#
# Run from the repository root.
# Output: table/supp/SuppTable_S1_Ratio_Statistics.csv
#         table/supp/SuppTable_S2_Species_Jackknife.csv
#         table/supp/SuppTable_S3_Species_Stability_by_Class.csv
# ==============================================================================

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                                                                              ║
# ║                    SoLD PAPER - FINAL FIGURES SCRIPT                         ║
# ║                                                                              ║
# ║  This script generates all main and supplementary figures for the paper.    ║
# ║                                                                              ║
# ║  FIGURES:                                                                    ║
# ║    Fig 1: Lipidomics Landscape (Overview)                                    ║
# ║           1A: Experimental Design (placeholder)                              ║
# ║           1B: Class-Level CLR Contrast                                       ║
# ║           1C: log2FC                                                         ║
# ║           1D: LION Enrichment Analysis                                       ║
# ║           1E: Chemical Shifts                                                ║
# ║    Fig 2: GWAS (TO BE ADDED)                                                 ║
# ║    Fig 3: OPLS-DA Analysis                                                   ║
# ║    Fig 4: Predictive Modeling (RF + SHAP, Plant Height)                      ║
# ║    Fig 6: Summary Model (TO BE ADDED)                                        ║
# ║                                                                              ║
# ║  SUPPLEMENTARY FIGURES:                                                      ║
# ║    Supp Fig S1: QC - Run Order, SERRF RSD, PCA, and SpATS residuals         ║
# ║    Supp Fig S2: %TIC and ALR Compositional Context Contrasts                ║
# ║    Supp Fig S3: Lipid Species Counts & Overview                              ║
# ║    Supp Fig S4: PCA of Lipids (Control & LowInput)                           ║
# ║    Supp Fig S5: Top Variance Lipids (Genotype Variation)                     ║
# ║    Supp Fig S8: Partial Correlation heatmaps (controlling for total TIC)    ║
# ║    Supp Fig S9: RF/SHAP details + Plant Height comparison                    ║
# ║    Supp Fig S10: Flowering Time combined (S11 A/B + S10 C/D)                ║
# ║    Supp Fig S11: Flowering Time SHAP Heatmap (source panel for S10)         ║
# ║    Supp Fig S12+: Additional GWAS (TO BE ADDED)                              ║
# ║                                                                              ║
# ║  SUPPLEMENTARY TABLES:                                                       ║
# ║    - Species Composition by Class (Top 10 per class, Control & LowInput)    ║
# ║    - Species Comparison (Control vs LowInput with log2FC)                    ║
# ║    - Species Change Summary (N increased/decreased per class)                ║
# ║    - Plant Height Comparison (Our Data vs Boyles top 20 SHAP)               ║
# ║    - Flowering Time Environment Overlap                                      ║
# ║                                                                              ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

# ══════════════════════════════════════════════════════════════════════════════
# 0) SETUP: PACKAGES, THEME, PALETTES
# ══════════════════════════════════════════════════════════════════════════════

if (identical(Sys.getenv("ONLY_FIG3_GWAS", ""), "1") ||
    identical(Sys.getenv("ONLY_SUPPFIG_S8_GWAS", ""), "1")) {
  # Minimal library set to avoid OpenMP-heavy packages when only GWAS is needed
  # and pin threads to avoid OpenMP SHM issues in this environment.
  Sys.setenv(OMP_NUM_THREADS = "1", KMP_USE_SHM = "0")
  suppressPackageStartupMessages({
    library(vroom)
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(ggplot2)
    library(patchwork)
    library(viridis)
    library(grid)
    library(scales)
  })
} else {
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
    library(ppcor)        # For partial correlation analysis
    library(Hmisc)
    library(corrplot)
    library(uwot)
    library(ropls)
    library(igraph)
    library(ranger)
    library(treeshap)
    library(shapviz)
    library(caret)
    library(mlr)
    library(tuneRanger)
  })
}

# Create output directories
dir.create("fig/main", recursive = TRUE, showWarnings = FALSE)
dir.create("fig/supp", recursive = TRUE, showWarnings = FALSE)
dir.create("table/supp", recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# Plot Theme (Publication Quality)
# ─────────────────────────────────────────────────────────────────────────────

# Base axis line theme - add this to all plots for consistent axis lines
axis_lines <- theme(
 axis.line.x = element_line(color = "black", linewidth = 0.5),
 axis.line.y = element_line(color = "black", linewidth = 0.5)
)

plot_theme <- theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 10)
    ),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    plot.margin = margin(15, 15, 15, 15)
  )

# For large multi-panel figures
plot_theme_large <- theme_minimal(base_size = 24) +
  theme(
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    axis.title.x = element_text(size = 24, face = "bold"),
    axis.title.y = element_text(size = 24, face = "bold"),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    plot.margin = margin(15, 15, 15, 15)
  )

# Compact theme for smaller panels
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
  PC   = "#00441B",
  PA   = "#1B7837",
  PE   = "#41AB5D",
  LPC  = "#66C2A4",
  LPE  = "#2CA25F",

  PG   = "#78C679",
  PS   = "#C2E699",
  DG   = "#54278F",
  DGDG = "#F768A1",
  MG   = "#8941ED",
  MGDG = "#FBB4D9",
  SQDG = "#9D4D6C",
  TG   = "#ED804A"
)

condition_colors <- c(
  Control = "#440154FF",
  LowInput = "#FDE725FF"
)

group_colors <- c(
  "Glycolipids"    = "#E7298A",
  "Phospholipids"  = "#1B9E77",
  "Neutral"        = "#6A51A3",
  "Storage"        = "#ED804A"
)

valid_classes <- names(class_colors)
classes_with_lyso <- unique(c(valid_classes, "LPC", "LPE"))

# ─────────────────────────────────────────────────────────────────────────────
# Helper Functions
# ─────────────────────────────────────────────────────────────────────────────
get_lipid_class <- function(x) {
  cls <- str_extract(x, "^[A-Za-z0-9]+(?=\\()")
  cls <- ifelse(is.na(cls), NA_character_, cls)
  ifelse(cls %in% valid_classes, cls, "Other")
}

parse_lipid_features <- function(lipid_name) {
  inside <- str_match(lipid_name, "\\(([^\\)]+)\\)")[, 2]
  if (is.na(inside)) {
    return(tibble(total_c = NA_real_, total_db = NA_real_, n_chains = NA_integer_))
  }
  parts <- str_split(inside, "/", simplify = TRUE)
  parts <- parts[parts != ""]
  n_chains <- length(parts)

  if (n_chains == 1 && str_detect(parts[1], "^\\d+\\s*:\\s*\\d+$")) {
    td <- str_split(parts[1], ":", simplify = TRUE)
    return(tibble(total_c = as.numeric(str_trim(td[1])),
                  total_db = as.numeric(str_trim(td[2])),
                  n_chains = 1L))
  }

  carb <- db <- rep(NA_real_, n_chains)
  for (i in seq_len(n_chains)) {
    if (str_detect(parts[i], "^\\d+\\s*:\\s*\\d+$")) {
      td <- str_split(parts[i], ":", simplify = TRUE)
      carb[i] <- as.numeric(str_trim(td[1]))
      db[i] <- as.numeric(str_trim(td[2]))
    }
  }
  tibble(total_c = sum(carb, na.rm = TRUE),
         total_db = sum(db, na.rm = TRUE),
         n_chains = as.integer(n_chains))
}

class_to_group <- function(cls) {
  case_when(
    cls %in% c("PC", "PA", "PE", "LPC", "LPE", "PG", "PS") ~ "Phospholipids",
    cls %in% c("MGDG", "DGDG", "SQDG") ~ "Glycolipids",
    cls %in% c("DG", "MG") ~ "Neutral",
    cls %in% c("TG") ~ "Storage",
    TRUE ~ "Other"
  )
}

# Save function
save_fig <- function(p, filename, w = 10, h = 8, dpi = 300, subdir = "main") {
  # Global axis text/title scaling so all saved plots are more readable.
  axis_text_size_global <- 13
  axis_title_size_global <- 15

  p <- p +
    theme(
      axis.text = element_text(size = axis_text_size_global),
      axis.text.x = element_text(size = axis_text_size_global),
      axis.text.y = element_text(size = axis_text_size_global),
      axis.title = element_text(size = axis_title_size_global),
      axis.title.x = element_text(size = axis_title_size_global),
      axis.title.y = element_text(size = axis_title_size_global)
    )

  path <- if (filename == "SuppFig_TopVariance_Lipids.png") {
    "/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database/fig/new_figures/SuppFig_TopVariance_Lipids.png"
  } else {
    file.path("fig", subdir, filename)
  }
  ggsave(path, plot = p, width = w, height = h, units = "in", dpi = dpi, bg = "white")
  message("Saved: ", path)
}

format_small_num <- function(x, digits = 2, threshold = 0.01) {
  sapply(x, function(v) {
    if (!is.finite(v)) return(NA_character_)
    if (v == 0) return(formatC(0, format = "f", digits = digits))
    if (abs(v) < threshold) {
      sign_prefix <- ifelse(v < 0, "-", "")
      return(paste0(sign_prefix, "<", formatC(threshold, format = "f", digits = 2)))
    }
    formatC(v, format = "f", digits = digits)
  })
}

# Resolve exact paths and glob patterns, then return newest existing file.
resolve_latest_existing_file <- function(path_candidates) {
  expanded <- unique(unlist(lapply(path_candidates, function(p) {
    if (grepl("[\\*\\?\\[]", p)) {
      Sys.glob(p)
    } else {
      p
    }
  })))
  expanded <- expanded[file.exists(expanded)]
  if (length(expanded) == 0) {
    return(NA_character_)
  }
  info <- file.info(expanded)
  expanded[which.max(info$mtime)]
}


make_suppfig_s9_lipid_ratios <- function() {
  message("\n══════════════════════════════════════════════════════════════")
  message("SUPPLEMENTARY FIGURE S9: LIPID CLASS RATIO PLOTS")
  message("══════════════════════════════════════════════════════════════\n")

  control_file <- "data/SPATS_fitted/non_normalized_intensities/Final_subset_control_all_lipids_fitted_phenotype_non_normalized.csv"
  lowinput_file <- "data/SPATS_fitted/non_normalized_intensities/Final_subset_lowinput_all_lipids_fitted_phenotype_non_normalized.csv"

  if (!file.exists(control_file) || !file.exists(lowinput_file)) {
    stop("Missing ratio input files: ", control_file, " or ", lowinput_file)
  }

  control <- vroom::vroom(control_file, show_col_types = FALSE) %>% dplyr::select(-c(2, 3, 4))
  lowinput <- vroom::vroom(lowinput_file, show_col_types = FALSE) %>% dplyr::select(-c(2, 3, 4))
  colnames(control)[1] <- "Compound_Name"
  colnames(lowinput)[1] <- "Compound_Name"

  # Valid classes for class extraction (avoid DG matching DGDG, etc via word-boundary).
  ratio_classes <- c("TG", "DG", "MG", "PC", "PE", "DGDG", "MGDG", "SQDG", "LPC", "LPE", "PG", "PA", "PS")
  class_pat <- paste0("\\b(", paste(ratio_classes, collapse = "|"), ")\\b")

  long_prep_global <- function(df) {
    df_long <- df %>%
      tidyr::pivot_longer(-c(Compound_Name, Condition),
                          names_to = "Lipid", values_to = "Intensity") %>%
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
        minpos = dplyr::if_else(
          all(is.na(rel_abund) | rel_abund <= 0),
          NA_real_,
          min(rel_abund[rel_abund > 0], na.rm = TRUE)
        ),
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

  # Ratio lists (OPLS-DA VIP ratios and partial-correlation pairs)
  opls_ratios <- c(
    "MG/SQDG", "MG/MGDG", "PS/SQDG", "DGDG/MG", "PG/SQDG", "DG/MG", "PE/SQDG",
    "LPC/PS", "DGDG/PS", "DG/PS", "PC/PS", "LPE/SQDG", "LPC/MG", "PE/PS",
    "MGDG/PG", "PG/PS", "LPC/LPE", "MGDG/PE", "MG/PG", "SQDG/TG", "PA/PS",
    "LPE/MGDG", "PC/SQDG"
  )
  pcorr_ratios <- c(
    "TG/PE", "TG/MGDG", "TG/MG", "TG/DG", "PS/PC", "PS/MGDG", "PS/DG",
    "PG/MGDG", "PG/DGDG", "PE/PA"
  )

  ratio_specs <- tibble::tibble(
    ratio_label = c(opls_ratios, pcorr_ratios),
    ratio_set = c(rep("OPLS-DA ratios", length(opls_ratios)), rep("Partial correlation pairs", length(pcorr_ratios)))
  ) %>%
    tidyr::separate(ratio_label, into = c("Num", "Den"), sep = "/", remove = FALSE) %>%
    dplyr::mutate(
      Num = stringr::str_trim(Num),
      Den = stringr::str_trim(Den),
      ratio_id = paste0(Num, "_", Den)
    )

  ratio_tbl <- wide_log
  for (i in seq_len(nrow(ratio_specs))) {
    num <- ratio_specs$Num[i]
    den <- ratio_specs$Den[i]
    rid <- ratio_specs$ratio_id[i]
    if (num %in% names(ratio_tbl) && den %in% names(ratio_tbl)) {
      ratio_tbl[[rid]] <- ratio_tbl[[num]] - ratio_tbl[[den]]  # log10(A/B)
    } else {
      ratio_tbl[[rid]] <- NA_real_
    }
  }

  ratio_long <- ratio_tbl %>%
    dplyr::select(Sample, Condition, dplyr::all_of(ratio_specs$ratio_id)) %>%
    tidyr::pivot_longer(-c(Sample, Condition), names_to = "ratio_id", values_to = "Value") %>%
    dplyr::left_join(ratio_specs %>% dplyr::select(ratio_id, ratio_label, ratio_set), by = "ratio_id") %>%
    dplyr::filter(!is.na(Value) & is.finite(Value))

  ratio_long <- ratio_long %>%
    dplyr::mutate(
      Condition = factor(Condition, levels = c("Control", "LowInput")),
      ratio_set = factor(ratio_set, levels = c("OPLS-DA ratios", "Partial correlation pairs"))
    )

  plot_one_set <- function(set_name, facet_levels, ncol = 5) {
    df <- ratio_long %>%
      dplyr::filter(ratio_set == set_name) %>%
      dplyr::mutate(ratio_label = factor(ratio_label, levels = facet_levels))

    stat_df <- df %>%
      dplyr::group_by(ratio_label) %>%
      dplyr::summarise(
        p_value = tryCatch(t.test(Value ~ Condition)$p.value, error = function(e) NA_real_),
        y_max = max(Value, na.rm = TRUE),
        y_min = min(Value, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        y_span = dplyr::if_else(is.finite(y_max - y_min) & (y_max - y_min) > 0, y_max - y_min, 0.30),
        y_line = y_max + 0.12 * y_span,
        y_text = y_max + 0.20 * y_span,
        star = dplyr::case_when(
          is.na(p_value) ~ "",
          p_value < 0.001 ~ "***",
          p_value < 0.01 ~ "**",
          p_value < 0.05 ~ "*",
          TRUE ~ ""
        )
      ) %>%
      dplyr::filter(star != "")

    ggplot(df, aes(x = Condition, y = Value, fill = Condition)) +
      geom_violin(width = 0.9, alpha = 0.25, trim = TRUE, color = NA) +
      geom_boxplot(width = 0.55, alpha = 0.85, color = "black", outlier.shape = NA, linewidth = 0.25) +
      geom_segment(
        data = stat_df,
        inherit.aes = FALSE,
        aes(x = 1, xend = 2, y = y_line, yend = y_line),
        color = "black",
        linewidth = 0.25
      ) +
      geom_text(
        data = stat_df,
        inherit.aes = FALSE,
        aes(x = 1.5, y = y_text, label = star),
        size = 3.2,
        color = "black"
      ) +
      scale_fill_manual(values = condition_colors) +
      labs(x = NULL, y = "log10(class ratio)") +
      plot_theme_compact +
      theme(
        legend.position = "none",
        strip.text = element_text(face = "bold", size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      facet_wrap(~ ratio_label, scales = "free_y", ncol = ncol)
  }

  # Use a consistent column count so facet panels match in size across A/B.
  p_opls <- plot_one_set("OPLS-DA ratios", opls_ratios, ncol = 5) +
    ggtitle("A) OPLS-DA Discriminatory Ratios (log10 scale)")

  p_pcorr <- plot_one_set("Partial correlation pairs", pcorr_ratios, ncol = 5) +
    ggtitle("B) Partial-Correlation-Motivated Ratios (log10 scale)")

  suppfig_s9 <- (p_opls / p_pcorr) +
    # Panel A has 5 rows (23 ratios at 5 cols); Panel B has 2 rows (10 ratios at 5 cols).
    patchwork::plot_layout(heights = c(5, 2)) +
    patchwork::plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 14, face = "bold"))

  # save_fig(suppfig_s9, "SuppFig_S7_Lipid_Ratios.png", w = 14, h = 13, subdir = "supp")   # disabled: SuppFig_S7 is written by 28_class_logratio_stats.R
  # save_fig(suppfig_s9, "SuppFig_S7_Lipid_Ratios.pdf", w = 14, h = 13, subdir = "supp")   # disabled: SuppFig_S7 is written by 28_class_logratio_stats.R
  message("✓ Supplementary Figure S9 (ratio plots) complete!")

  # ─────────────────────────────────────────────────────────────────────────────
  # SUPPLEMENTARY TABLE S1: Lipid Ratio Statistics (Wilcoxon + Jackknife)
  # ─────────────────────────────────────────────────────────────────────────────
  # Robust statistics for each lipid class ratio comparing LowInput vs Control.
  # Uses Wilcoxon rank-sum test with Hodges-Lehmann estimator for median shift
  # and jackknife leave-one-out stability to confirm population-wide effects.
  # ─────────────────────────────────────────────────────────────────────────────

  message("\n── Computing ratio statistics (Wilcoxon + Jackknife) ──")

  # Helper function: Wilcoxon test + jackknife stability
  one_ratio_tests <- function(df, ratio_id_val) {
    xi <- df %>% filter(Condition == "LowInput", ratio_id == ratio_id_val) %>% pull(Value)
    x0 <- df %>% filter(Condition == "Control",  ratio_id == ratio_id_val) %>% pull(Value)

    # Get ratio label for output
    ratio_lab <- df %>%
      filter(ratio_id == ratio_id_val) %>%
      pull(ratio_label) %>%
      unique() %>%
      first()

    if (length(xi) < 3 || length(x0) < 3) {
      return(tibble(
        Ratio = ratio_lab,
        n_C = length(x0), n_LI = length(xi),
        median_C = NA_real_, median_LI = NA_real_,
        effect_log10 = NA_real_,
        p_wilcox = NA_real_,
        jackknife_stability = NA_real_
      ))
    }

    # Wilcoxon rank-sum test
    wt <- wilcox.test(xi, x0, alternative = "two.sided", exact = FALSE)

    # Jackknife sign stability: does direction hold when removing any single sample?
    d_full <- sign(median(xi) - median(x0))
    N <- length(xi) + length(x0)
    keep <- logical(N)
    for (i in seq_len(N)) {
      xi2 <- xi; x0_2 <- x0
      if (i <= length(xi)) xi2 <- xi2[-i] else x0_2 <- x0_2[-(i - length(xi))]
      keep[i] <- sign(median(xi2) - median(x0_2)) == d_full
    }
    stab <- mean(keep)

    tibble(
      Ratio = ratio_lab,
      n_C = length(x0),
      n_LI = length(xi),
      median_C = round(median(x0), 4),
      median_LI = round(median(xi), 4),
      effect_log10 = round(median(xi) - median(x0), 4),
      p_wilcox = wt$p.value,
      jackknife_stability = round(stab, 3)
    )
  }

  # Run tests for all ratios
  ratios_to_test <- unique(ratio_long$ratio_id)

  ratio_stats <- purrr::map_dfr(ratios_to_test, function(r) {
    one_ratio_tests(ratio_long, r)
  }) %>%
    dplyr::mutate(
      p_adj_BH = p.adjust(p_wilcox, method = "BH"),
      effect_fc = round(10^effect_log10, 2),
      direction = dplyr::case_when(
        effect_log10 > 0 ~ "LI > C",
        effect_log10 < 0 ~ "LI < C",
        TRUE ~ "No change"
      ),
      significance = dplyr::case_when(
        p_adj_BH < 0.001 ~ "***",
        p_adj_BH < 0.01 ~ "**",
        p_adj_BH < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    ) %>%
    dplyr::arrange(p_adj_BH)

  # Format p-values for table
  ratio_stats_out <- ratio_stats %>%
    dplyr::mutate(
      p_wilcox = format(p_wilcox, digits = 3, scientific = TRUE),
      p_adj_BH = format(p_adj_BH, digits = 3, scientific = TRUE)
    ) %>%
    dplyr::select(
      Ratio, n_C, n_LI, median_C, median_LI,
      effect_log10, effect_fc, direction,
      p_wilcox, p_adj_BH, significance,
      jackknife_stability
    )

  # Save supplementary table
  write.csv(ratio_stats_out, "table/supp/SuppTable_S1_Ratio_Statistics.csv", row.names = FALSE)
  message("✓ Saved: table/supp/SuppTable_S1_Ratio_Statistics.csv")

  # Print summary
  n_sig <- sum(ratio_stats$p_adj_BH < 0.05, na.rm = TRUE)
  n_stable <- sum(ratio_stats$jackknife_stability == 1, na.rm = TRUE)
  message(sprintf("  → %d/%d ratios significant (BH-adj p < 0.05)", n_sig, nrow(ratio_stats)))
  message(sprintf("  → %d/%d ratios with perfect jackknife stability (1.0)", n_stable, nrow(ratio_stats)))

  # ─────────────────────────────────────────────────────────────────────────────
  # SUPPLEMENTARY TABLE S2: Individual Species Jackknife Stability
  # ─────────────────────────────────────────────────────────────────────────────
  # Jackknife stability on individual lipid species (not class aggregates)
  # to confirm species-level changes are also population-wide.
  # ─────────────────────────────────────────────────────────────────────────────

  message("\n── Computing SPECIES-level jackknife stability (CLR-transformed) ──")

  # Find COMMON species between Control and LowInput
  control_species <- colnames(control)[-1]  # exclude Compound_Name
  lowinput_species <- colnames(lowinput)[-1]
  common_species <- intersect(control_species, lowinput_species)
  # Keep only true lipid species columns (e.g., "PC(16:0/18:1)"),
  # matching the species universe used in SuppTable_S6 summaries.
  common_species <- common_species[stringr::str_detect(common_species, "\\(")]
  message(sprintf("  Found %d common lipid species between Control and LowInput", length(common_species)))

  # Prepare species-level data with CLR transformation
  # CLR removes compositional constraint: CLR(x_i) = log(x_i / geometric_mean(all_species))
  species_long <- dplyr::bind_rows(
    control %>% dplyr::mutate(Condition = "Control"),
    lowinput %>% dplyr::mutate(Condition = "LowInput")
  ) %>%
    tidyr::pivot_longer(-c(Compound_Name, Condition),
                        names_to = "Lipid", values_to = "Intensity") %>%
    dplyr::rename(Sample = Compound_Name) %>%
    dplyr::filter(Lipid %in% common_species) %>%
    # Handle zeros: replace with half the minimum positive value per sample
    dplyr::group_by(Sample) %>%
    dplyr::mutate(
      minpos = min(Intensity[Intensity > 0], na.rm = TRUE),
      Intensity_nozero = dplyr::if_else(
        is.na(Intensity) | Intensity <= 0,
        minpos * 0.5,
        Intensity
      )
    ) %>%
    dplyr::ungroup() %>%
    # CLR transformation per sample
    dplyr::group_by(Sample) %>%
    dplyr::mutate(
      log_intensity = log(Intensity_nozero),
      geom_mean_log = mean(log_intensity, na.rm = TRUE),
      CLR = log_intensity - geom_mean_log  # CLR = log(x_i / geom_mean)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      # Extract class - try multiple patterns for different lipid naming conventions
      Class = dplyr::case_when(
        stringr::str_detect(Lipid, "^TG\\(") ~ "TG",
        stringr::str_detect(Lipid, "^DG\\(") ~ "DG",
        stringr::str_detect(Lipid, "^MG\\(") ~ "MG",
        stringr::str_detect(Lipid, "^PC\\(") ~ "PC",
        stringr::str_detect(Lipid, "^PE\\(") ~ "PE",
        stringr::str_detect(Lipid, "^PG\\(") ~ "PG",
        stringr::str_detect(Lipid, "^PA\\(") ~ "PA",
        stringr::str_detect(Lipid, "^PS\\(") ~ "PS",
        stringr::str_detect(Lipid, "^DGDG\\(") ~ "DGDG",
        stringr::str_detect(Lipid, "^MGDG\\(") ~ "MGDG",
        stringr::str_detect(Lipid, "^SQDG\\(") ~ "SQDG",
        stringr::str_detect(Lipid, "^LPC\\(") ~ "LPC",
        stringr::str_detect(Lipid, "^LPE\\(") ~ "LPE",
        stringr::str_detect(Lipid, "^FA\\(") ~ "FA",
        stringr::str_detect(Lipid, "^Cer\\(") ~ "Cer",
        stringr::str_detect(Lipid, "^GalCer\\(") ~ "GalCer",
        stringr::str_detect(Lipid, "^SM\\(") ~ "SM",
        stringr::str_detect(Lipid, "^CL\\(") ~ "CL",
        stringr::str_detect(Lipid, "^AEG\\(") ~ "AEG",
        TRUE ~ stringr::str_extract(Lipid, "^[A-Za-z]+")
      )
    ) %>%
    dplyr::filter(!is.na(CLR), is.finite(CLR))

  # Helper function for species-level tests (using CLR values)
  one_species_tests <- function(df, lipid_name) {
    xi <- df %>% dplyr::filter(Condition == "LowInput", Lipid == lipid_name) %>% dplyr::pull(CLR)
    x0 <- df %>% dplyr::filter(Condition == "Control",  Lipid == lipid_name) %>% dplyr::pull(CLR)

    # Get class for output
    lipid_class <- df %>%
      dplyr::filter(Lipid == lipid_name) %>%
      dplyr::pull(Class) %>%
      unique() %>%
      first()

    if (length(xi) < 3 || length(x0) < 3) {
      return(tibble::tibble(
        Lipid = lipid_name,
        Class = lipid_class,
        n_C = length(x0), n_LI = length(xi),
        median_CLR_C = NA_real_, median_CLR_LI = NA_real_,
        effect_CLR = NA_real_,
        p_wilcox = NA_real_,
        jackknife_stability = NA_real_
      ))
    }

    # Wilcoxon rank-sum test on CLR values
    wt <- tryCatch(
      wilcox.test(xi, x0, alternative = "two.sided", exact = FALSE),
      error = function(e) list(p.value = NA_real_)
    )

    # Jackknife sign stability
    d_full <- sign(median(xi) - median(x0))
    N <- length(xi) + length(x0)
    keep <- logical(N)
    for (i in seq_len(N)) {
      xi2 <- xi; x0_2 <- x0
      if (i <= length(xi)) xi2 <- xi2[-i] else x0_2 <- x0_2[-(i - length(xi))]
      keep[i] <- sign(median(xi2) - median(x0_2)) == d_full
    }
    stab <- mean(keep)

    tibble::tibble(
      Lipid = lipid_name,
      Class = lipid_class,
      n_C = length(x0),
      n_LI = length(xi),
      median_CLR_C = round(median(x0), 4),
      median_CLR_LI = round(median(xi), 4),
      effect_CLR = round(median(xi) - median(x0), 4),  # CLR difference (ln scale)
      p_wilcox = wt$p.value,
      jackknife_stability = round(stab, 3)
    )
  }

  # Run tests for all species
  all_species <- unique(species_long$Lipid)
  message(sprintf("  Testing %d individual species...", length(all_species)))

  species_stats <- purrr::map_dfr(all_species, function(sp) {
    one_species_tests(species_long, sp)
  }, .progress = FALSE) %>%
    dplyr::mutate(
      p_adj_BH = p.adjust(p_wilcox, method = "BH"),
      # CLR is in natural log, so fold-change = exp(effect_CLR)
      effect_fc = round(exp(effect_CLR), 2),
      direction = dplyr::case_when(
        effect_CLR > 0 ~ "LI > C",
        effect_CLR < 0 ~ "LI < C",
        TRUE ~ "No change"
      ),
      significance = dplyr::case_when(
        p_adj_BH < 0.001 ~ "***",
        p_adj_BH < 0.01 ~ "**",
        p_adj_BH < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    ) %>%
    dplyr::arrange(p_adj_BH)

  # Format and save
  species_stats_out <- species_stats %>%
    dplyr::mutate(
      p_wilcox = format(p_wilcox, digits = 3, scientific = TRUE),
      p_adj_BH = format(p_adj_BH, digits = 3, scientific = TRUE)
    ) %>%
    dplyr::select(
      Lipid, Class, n_C, n_LI, median_CLR_C, median_CLR_LI,
      effect_CLR, effect_fc, direction,
      p_wilcox, p_adj_BH, significance,
      jackknife_stability
    )

  write.csv(species_stats_out, "table/supp/SuppTable_S2_Species_Jackknife.csv", row.names = FALSE)
  message("✓ Saved: table/supp/SuppTable_S2_Species_Jackknife.csv")

  # Summary statistics
  n_species_sig <- sum(species_stats$p_adj_BH < 0.05, na.rm = TRUE)
  n_species_stable <- sum(species_stats$jackknife_stability == 1, na.rm = TRUE)
  n_species_high_stable <- sum(species_stats$jackknife_stability >= 0.95, na.rm = TRUE)
  message(sprintf("  → %d/%d species significant (BH-adj p < 0.05)", n_species_sig, nrow(species_stats)))
  message(sprintf("  → %d/%d species with perfect jackknife stability (1.0)", n_species_stable, nrow(species_stats)))
  message(sprintf("  → %d/%d species with high jackknife stability (≥0.95)", n_species_high_stable, nrow(species_stats)))

  # Summary by class
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
  message("✓ Saved: table/supp/SuppTable_S3_Species_Stability_by_Class.csv")
  print(class_stability_summary)
}

# ---------------------------------------------------------------------------
make_suppfig_s9_lipid_ratios()

message("\n== Supplementary Tables S1-S3 written ==")
for (f in c("table/supp/SuppTable_S1_Ratio_Statistics.csv",
            "table/supp/SuppTable_S2_Species_Jackknife.csv",
            "table/supp/SuppTable_S3_Species_Stability_by_Class.csv")) {
  n <- nrow(read.csv(f))
  message(sprintf("  %-56s %4d rows", f, n))
}
