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

  save_fig(suppfig_s9, "SuppFig_S7_Lipid_Ratios.png", w = 14, h = 13, subdir = "supp")
  save_fig(suppfig_s9, "SuppFig_S7_Lipid_Ratios.pdf", w = 14, h = 13, subdir = "supp")
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

# Allow running only this panel without regenerating the full figure suite.
if (identical(Sys.getenv("ONLY_SUPPFIG_S9_RATIOS", ""), "1")) {
  make_suppfig_s9_lipid_ratios()
  message("ONLY_SUPPFIG_S9_RATIOS=1: done.")
  quit(save = "no", status = 0)
}

# ---- Figure 3: GWAS panels (3 columns x 4 rows) ----
make_fig3_gwas <- function() {
  message("\n══════════════════════════════════════════════════════════════")
  message("FIGURE 3: GWAS PANELS (Highest-count split, PHR1, DGAT1)")
  message("══════════════════════════════════════════════════════════════\n")

  base_dir <- "results/gwas_raw"
  highest_individual_dir <- file.path(base_dir, "highest_count_individual")
  highest_sr_dir <- file.path(base_dir, "highest_count_sum_ratio_Expansin_SORBI_3004G120000")
  phr1_dir <- file.path(base_dir, "PHR1")
  dgat1_dir <- file.path(base_dir, "DGAT1")

  stopifnot(
    dir.exists(highest_individual_dir),
    dir.exists(highest_sr_dir),
    dir.exists(phr1_dir),
    dir.exists(dgat1_dir)
  )

  list_sorted <- function(path) {
    files <- list.files(path, pattern = "\\.assoc\\.txt$", full.names = TRUE)
    if (length(files) == 0) {
      files <- list.files(path, pattern = "\\.txt$", full.names = TRUE)
    }
    sort(files)
  }

  # Highest-count split: 2 individual + 2 sum/ratio
  highest_individual_files <- list_sorted(highest_individual_dir)
  highest_sr_files <- list_sorted(highest_sr_dir)
  if (length(highest_individual_files) < 2 || length(highest_sr_files) < 2) {
    stop("Need at least 2 files each in highest_count_individual and highest_count_sum_ratio_Expansin_SORBI_3004G120000.")
  }
  highest_individual_files <- highest_individual_files[1:2]
  highest_sr_files <- highest_sr_files[1:2]

  # PHR1: first 3 (sorted) + PE(16:0_18:1)
  phr1_files_all <- list_sorted(phr1_dir)
  phr1_first3 <- phr1_files_all[1:3]
  phr1_pe <- phr1_files_all[grepl("PE\\(16:0_18:1\\)", phr1_files_all)]
  if (length(phr1_pe) != 1) {
    stop("Could not uniquely find PE(16:0_18:1) file in PHR1 folder.")
  }
  phr1_files <- c(phr1_first3, phr1_pe)

  # DGAT1: first 4 files (sorted)
  dgat1_files <- list_sorted(dgat1_dir)[1:4]

  # Marker SNP positions (chr, pos)
  marker_highest_individual <- list(chr = 10, pos = 5083051)
  marker_highest_sr <- list(chr = 4, pos = 12910117)
  marker_phr1 <- list(chr = 1, pos = 67813097)
  marker_dgat1 <- list(chr = 10, pos = 50123510)

  read_assoc <- function(f) {
    df <- vroom::vroom(
      f,
      show_col_types = FALSE,
      col_select = c(chr, ps, p_wald)
    )
    # Normalize column names
    names(df) <- tolower(names(df))
    # Expected columns: chr, ps (position), p_wald
    if (!all(c("chr", "ps", "p_wald") %in% names(df))) {
      stop("Missing required columns in ", f, ". Expected chr, ps, p_wald.")
    }
    df %>%
      dplyr::transmute(
        chr = as.integer(chr),
        pos = as.numeric(ps),
        p = as.numeric(p_wald)
      ) %>%
      dplyr::filter(is.finite(chr) & is.finite(pos) & is.finite(p) & p > 0)
  }

  prep_manhattan <- function(df) {
    # cumulative positions by chromosome
    chr_info <- df %>%
      dplyr::group_by(chr) %>%
      dplyr::summarise(chr_len = max(pos, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(chr) %>%
      dplyr::mutate(chr_start = dplyr::lag(cumsum(chr_len), default = 0))

    df2 <- df %>%
      dplyr::left_join(chr_info, by = "chr") %>%
      dplyr::mutate(pos_cum = pos + chr_start,
                    logp = -log10(p))

    axis_df <- chr_info %>%
      dplyr::mutate(center = chr_start + chr_len / 2)

    list(df = df2, axis = axis_df)
  }

  label_from_file <- function(f) {
    lab <- sub("_mod_sub_.*$", "", basename(f))
    lab <- sub("\\.txt$", "", lab)
    # display slashes in lipid names
    lab <- gsub("_", "/", lab, fixed = TRUE)
    lab
  }

  prep_panel <- function(f, marker, col_title) {
    df <- read_assoc(f)
    prep <- prep_manhattan(df)
    chr_info <- df %>%
      dplyr::group_by(chr) %>%
      dplyr::summarise(chr_len = max(pos, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(chr) %>%
      dplyr::mutate(chr_start = dplyr::lag(cumsum(chr_len), default = 0))
    vline_x <- marker$pos + chr_info$chr_start[match(marker$chr, chr_info$chr)]

    thin_for_plot <- function(df_plot, max_points = 450000L, keep_logp_ge = 5, cap_logp = 20) {
      # Drop extreme outliers above cap from plotting (do not clamp to cap)
      df_plot <- df_plot %>% dplyr::filter(logp <= cap_logp)
      n <- nrow(df_plot)
      if (n <= max_points) return(df_plot)
      hi <- which(df_plot$logp >= keep_logp_ge)
      if (length(hi) >= max_points) {
        hi <- hi[order(df_plot$logp[hi], decreasing = TRUE)]
        return(df_plot[sort(hi[seq_len(max_points)]), , drop = FALSE])
      }
      rest <- setdiff(seq_len(n), hi)
      set.seed(42)
      keep_rest <- sample(rest, max_points - length(hi))
      keep <- sort(c(hi, keep_rest))
      df_plot[keep, , drop = FALSE]
    }

    list(
      df = thin_for_plot(prep$df),
      axis = prep$axis,
      vline = vline_x,
      title = paste0(col_title, "\n", label_from_file(f))
    )
  }

  highest_panels <- list(
    prep_panel(highest_individual_files[1], marker_highest_individual, "Highest count individual"),
    prep_panel(highest_individual_files[2], marker_highest_individual, "Highest count individual"),
    prep_panel(highest_sr_files[1], marker_highest_sr, "Highest count sums and ratios"),
    prep_panel(highest_sr_files[2], marker_highest_sr, "Highest count sums and ratios")
  )
  phr1_panels <- list(
    prep_panel(phr1_files[1], marker_phr1, "PHR1"),
    prep_panel(phr1_files[2], marker_phr1, "PHR1"),
    prep_panel(phr1_files[3], marker_phr1, "PHR1"),
    prep_panel(phr1_files[4], marker_phr1, "PHR1")
  )
  dgat1_panels <- list(
    prep_panel(dgat1_files[1], marker_dgat1, "DGAT1"),
    prep_panel(dgat1_files[2], marker_dgat1, "DGAT1"),
    prep_panel(dgat1_files[3], marker_dgat1, "DGAT1"),
    prep_panel(dgat1_files[4], marker_dgat1, "DGAT1")
  )

  draw_grid <- function(outfile, device = c("png", "pdf")) {
    device <- match.arg(device)
    if (device == "png") {
      png(outfile, width = 3000, height = 3600, res = 300)
    } else {
      pdf(outfile, width = 12, height = 14)
    }
    par(mfrow = c(4, 3), mar = c(2.8, 3.2, 2.8, 1.0), oma = c(2, 2, 2, 1))
    group_labels <- list(
      list(idx = 1,  lab = "A"),  # Highest count individual (row1,col1)
      list(idx = 7,  lab = "B"),  # Highest count sums/ratios (row3,col1)
      list(idx = 2,  lab = "C"),  # PHR1 (row1,col2)
      list(idx = 3,  lab = "D")   # DGAT1 (row1,col3)
    )

    # Row-major order: each row = Highest_count, PHR1, DGAT1
    panels_ordered <- list()
    for (r in 1:4) {
      panels_ordered <- c(panels_ordered, list(highest_panels[[r]], phr1_panels[[r]], dgat1_panels[[r]]))
    }

    for (i in seq_along(panels_ordered)) {
      p <- panels_ordered[[i]]
      row <- ((i - 1) %/% 3) + 1
      col <- ((i - 1) %% 3) + 1

      show_x <- row == 4
      show_y_label <- col == 1

      # alternate chromosome colors
      pt_cols <- ifelse(p$df$chr %% 2 == 0, "grey70", "grey40")

      plot(
        p$df$pos_cum, p$df$logp,
        pch = 20, cex = 0.45, col = pt_cols,
        xlab = "", ylab = "",
        xaxt = "n", yaxt = "n",
        bty = if (show_x || show_y_label) "l" else "n"
      )
      abline(v = p$vline, lty = 3, lwd = 1.6, col = "#D95F02")
      abline(h = 7, lty = 2, lwd = 0.9, col = "grey55")
      abline(h = 5, lty = 2, lwd = 0.9, col = "#2C7FB8")
      title(main = p$title, cex.main = 0.95, font.main = 2)

      if (show_x) {
        axis(1, at = p$axis$center, labels = p$axis$chr, cex.axis = 0.9, las = 1)
      }
      # Always show y-axis ticks/labels (each panel has its own scale)
      axis(2, cex.axis = 0.9)
      if (show_y_label) {
        mtext(expression(-log[10](p)), side = 2, line = 2.3, cex = 0.95)
      }

      # Group labels: A=highest individual, B=highest sums/ratios, C=PHR1, D=DGAT1
      label_here <- group_labels[vapply(group_labels, function(x) x$idx == i, logical(1))]
      if (length(label_here) == 1) {
        usr <- par("usr")
        text(x = usr[1], y = usr[4], labels = label_here[[1]]$lab,
             adj = c(-0.2, 1.2), cex = 1.2, font = 2)
      }
    }

    dev.off()
  }

  draw_grid("fig/main/Figure3_GWAS.png", device = "png")
  message("✓ Figure 3 (GWAS panels) complete!")
}

# Allow running only this panel without regenerating the full figure suite.
if (identical(Sys.getenv("ONLY_FIG3_GWAS", ""), "1")) {
  make_fig3_gwas()
  message("ONLY_FIG3_GWAS=1: done.")
  quit(save = "no", status = 0)
}

# ---- Supplementary Figure S8: GWAS panels (2 columns x 3 rows) ----
make_suppfig_s8_gwas <- function() {
  message("\n══════════════════════════════════════════════════════════════")
  message("SUPP FIG S8: GWAS PANELS (alpha, zeax, SQDG + CDPK7/Tryptophan/PHT)")
  message("══════════════════════════════════════════════════════════════\n")

  base_dir <- "results/gwas_raw"
  cdpk7_dir <- file.path(base_dir, "CDPK7_SORBI_3001G520900")
  trp_dir <- file.path(base_dir, "Tryptophan_SORBI_3001G056400")
  pht_dir <- file.path(base_dir, "PHT_SORBI_3006G026800")
  sum_dir <- file.path(base_dir, "Sum_SQDG_APK3")
  alpha_dir <- file.path(base_dir, "alpha_carotene")
  zeax_dir <- file.path(base_dir, "zeaxanthin")

  stopifnot(
    dir.exists(cdpk7_dir), dir.exists(trp_dir), dir.exists(pht_dir),
    dir.exists(sum_dir), dir.exists(alpha_dir), dir.exists(zeax_dir)
  )

  list_sorted_txt <- function(path) {
    sort(list.files(path, pattern = "\\.txt$", full.names = TRUE))
  }

  cdpk7_files <- list_sorted_txt(cdpk7_dir)
  trp_files <- list_sorted_txt(trp_dir)
  pht_files <- list_sorted_txt(pht_dir)
  if (length(cdpk7_files) < 3 || length(trp_files) < 3 || length(pht_files) < 3) {
    stop("CDPK7, Tryptophan, and PHT folders must each have at least 3 GWAS files.")
  }
  cdpk7_files <- cdpk7_files[1:3]
  trp_files <- trp_files[1:3]
  pht_files <- pht_files[1:3]

  sum_file <- list_sorted_txt(sum_dir)
  alpha_file <- list_sorted_txt(alpha_dir)
  zeax_file <- list_sorted_txt(zeax_dir)
  if (length(sum_file) != 1 || length(alpha_file) != 1 || length(zeax_file) != 1) {
    stop("Expected exactly one GWAS file in Sum_SQDG_APK3, alpha_carotene, zeaxanthin folders.")
  }

  marker_sum <- list(chr = 5, pos = 67974522)
  marker_alpha <- list(chr = 6, pos = 55355499)
  marker_zeax <- list(chr = 2, pos = 3107262)
  marker_cdpk7 <- list(chr = 1, pos = 78622695)
  marker_trp <- list(chr = 1, pos = 4237343)
  marker_pht <- list(chr = 6, pos = 4906326)

  read_assoc <- function(f) {
    df <- vroom::vroom(
      f,
      show_col_types = FALSE,
      col_select = c(chr, ps, p_wald)
    )
    names(df) <- tolower(names(df))
    df %>%
      dplyr::transmute(
        chr = as.integer(chr),
        pos = as.numeric(ps),
        p = as.numeric(p_wald)
      ) %>%
      dplyr::filter(is.finite(chr) & is.finite(pos) & is.finite(p) & p > 0)
  }

  prep_manhattan <- function(df) {
    chr_info <- df %>%
      dplyr::group_by(chr) %>%
      dplyr::summarise(chr_len = max(pos, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(chr) %>%
      dplyr::mutate(chr_start = dplyr::lag(cumsum(chr_len), default = 0))

    df2 <- df %>%
      dplyr::left_join(chr_info, by = "chr") %>%
      dplyr::mutate(pos_cum = pos + chr_start,
                    logp = -log10(p))

    axis_df <- chr_info %>%
      dplyr::mutate(center = chr_start + chr_len / 2)

    list(df = df2, axis = axis_df, chr_info = chr_info)
  }

  label_from_file <- function(f) {
    lab <- sub("_mod_sub_.*$", "", basename(f))
    lab <- sub("\\.txt$", "", lab)
    lab <- gsub("_", "/", lab, fixed = TRUE)
    lab
  }

  prep_panel <- function(f, marker, col_title) {
    df <- read_assoc(f)
    prep <- prep_manhattan(df)
    vline_x <- marker$pos + prep$chr_info$chr_start[match(marker$chr, prep$chr_info$chr)]
    list(
      df = prep$df,
      axis = prep$axis,
      vline = vline_x,
      title = paste0(col_title, "\n", label_from_file(f))
    )
  }

  # Top row single-GWAS panels (A/B/C)
  panel_A <- prep_panel(alpha_file[1], marker_alpha, "alpha_carotene")
  panel_B <- prep_panel(zeax_file[1], marker_zeax, "zeaxanthin")
  panel_C <- prep_panel(sum_file[1], marker_sum, "Sum_SQDG_APK3")

  # Three separate GWAS per folder (rows 2-4; D/E/F topics)
  cdpk7_panels <- lapply(cdpk7_files, function(f) prep_panel(f, marker_cdpk7, "CDPK7_SORBI_3001G520900"))
  trp_panels   <- lapply(trp_files,   function(f) prep_panel(f, marker_trp, "Tryptophan_SORBI_3001G056400"))
  pht_panels   <- lapply(pht_files,   function(f) prep_panel(f, marker_pht, "PHT_SORBI_3006G026800"))

  out_png <- file.path("fig/supp", "SuppFig_S8_GWAS.png")
  png(out_png, width = 3300, height = 3000, res = 300)
  par(mfrow = c(4, 3), mar = c(3.0, 3.6, 2.8, 1.0), oma = c(2, 2, 2, 1))

  panels_ordered <- list(
    panel_A, panel_B, panel_C,
    cdpk7_panels[[1]], trp_panels[[1]], pht_panels[[1]],
    cdpk7_panels[[2]], trp_panels[[2]], pht_panels[[2]],
    cdpk7_panels[[3]], trp_panels[[3]], pht_panels[[3]]
  )
  y_cap <- 20

  for (i in seq_along(panels_ordered)) {
    p <- panels_ordered[[i]]
    row <- ((i - 1) %/% 3) + 1
    col <- ((i - 1) %% 3) + 1

    show_x <- row == 4
    show_y_label <- col == 1

    # alternate chromosome colors
    pt_cols <- ifelse(p$df$chr %% 2 == 0, "grey70", "grey40")

    keep <- is.finite(p$df$logp) & (p$df$logp <= y_cap)
    x_plot <- p$df$pos_cum[keep]
    y_plot <- p$df$logp[keep]
    col_plot <- pt_cols[keep]

    y_top <- suppressWarnings(max(y_plot, na.rm = TRUE))
    if (!is.finite(y_top)) y_top <- 7.2
    y_top <- max(7.2, y_top * 1.03)
    plot(
      x_plot, y_plot,
      pch = 20, cex = 0.45, col = col_plot,
      xlab = "", ylab = "",
      xaxt = "n", yaxt = "n",
      ylim = c(0, y_top),
      bty = if (show_x || show_y_label) "l" else "n"
    )
    # marker SNP (only where provided)
    if (is.finite(p$vline)) {
      abline(v = p$vline, lty = 3, lwd = 1.6, col = "#D95F02")
    }
    # thresholds: genome-wide (7) and suggestive (5)
    abline(h = 7, lty = 2, lwd = 0.9, col = "grey55")
    abline(h = 5, lty = 2, lwd = 0.9, col = "#2C7FB8")
    title(main = p$title, cex.main = 0.95, font.main = 2)

    if (show_x) {
      axis(1, at = p$axis$center, labels = p$axis$chr, cex.axis = 1.0, las = 1)
    }
    axis(2, cex.axis = 1.0)
    if (show_y_label) {
      mtext(expression(-log[10](p)), side = 2, line = 2.4, cex = 1.05)
    }

    # Topic labels only: A/B/C on row1, D/E/F on row2
    panel_label <- ""
    if (row == 1 && col == 1) panel_label <- "A"
    if (row == 1 && col == 2) panel_label <- "B"
    if (row == 1 && col == 3) panel_label <- "C"
    if (row == 2 && col == 1) panel_label <- "D"
    if (row == 2 && col == 2) panel_label <- "E"
    if (row == 2 && col == 3) panel_label <- "F"
    if (panel_label != "") {
      usr <- par("usr")
      text(x = usr[1], y = usr[4], labels = panel_label,
           adj = c(-0.2, 1.2), cex = 1.2, font = 2)
    }
  }

  dev.off()
  message("✓ Supp Fig S8 (GWAS panels) complete!")
}

# Allow running only this panel without regenerating the full figure suite.
if (identical(Sys.getenv("ONLY_SUPPFIG_S8_GWAS", ""), "1")) {
  make_suppfig_s8_gwas()
  message("ONLY_SUPPFIG_S8_GWAS=1: done.")
  quit(save = "no", status = 0)
}

# Read image from candidate paths (supports glob patterns, picks newest by mtime;
# PNG/JPEG with macOS sips fallback).
read_image_raster_any <- function(path_candidates) {
  path <- resolve_latest_existing_file(path_candidates)
  if (is.na(path)) {
    return(list(raster = NULL, path = NA_character_, err = "Image file not found"))
  }

  raster <- NULL
  err <- NULL

  if (requireNamespace("png", quietly = TRUE)) {
    raster <- tryCatch(png::readPNG(path), error = function(e) NULL)
  }

  if (is.null(raster) && requireNamespace("jpeg", quietly = TRUE)) {
    raster <- tryCatch(jpeg::readJPEG(path), error = function(e) NULL)
  }

  if (is.null(raster) && nzchar(Sys.which("sips")) && requireNamespace("png", quietly = TRUE)) {
    tmp_png <- tempfile(fileext = ".png")
    conv_status <- suppressWarnings(
      tryCatch(
        system2("sips", c("-s", "format", "png", path, "--out", tmp_png),
                stdout = TRUE, stderr = TRUE),
        error = function(e) NULL
      )
    )
    if (!is.null(conv_status) && file.exists(tmp_png)) {
      raster <- tryCatch(png::readPNG(tmp_png), error = function(e) NULL)
    }
  }

  if (is.null(raster)) {
    err <- "Unsupported image format or missing png/jpeg reader packages"
  }

  list(raster = raster, path = path, err = err)
}

crop_raster_grid <- function(raster, row_idx, col_idx, n_rows = 1, n_cols = 1) {
  d <- dim(raster)
  if (length(d) < 2) return(NULL)

  h <- d[1]
  w <- d[2]

  r_breaks <- round(seq(0, h, length.out = n_rows + 1))
  c_breaks <- round(seq(0, w, length.out = n_cols + 1))

  r1 <- max(1, r_breaks[row_idx] + 1)
  r2 <- min(h, r_breaks[row_idx + 1])
  c1 <- max(1, c_breaks[col_idx] + 1)
  c2 <- min(w, c_breaks[col_idx + 1])

  if (length(d) == 2) {
    raster[r1:r2, c1:c2, drop = FALSE]
  } else {
    raster[r1:r2, c1:c2, , drop = FALSE]
  }
}

crop_raster_col_ratio <- function(raster, left_ratio = 0.7) {
  d <- dim(raster)
  if (length(d) < 2) return(list(left = NULL, right = NULL))

  h <- d[1]
  w <- d[2]
  c_cut <- max(1, min(w - 1, round(w * left_ratio)))

  if (length(d) == 2) {
    left <- raster[1:h, 1:c_cut, drop = FALSE]
    right <- raster[1:h, (c_cut + 1):w, drop = FALSE]
  } else {
    left <- raster[1:h, 1:c_cut, , drop = FALSE]
    right <- raster[1:h, (c_cut + 1):w, , drop = FALSE]
  }

  list(left = left, right = right)
}

panel_from_raster <- function(raster, fallback_label = "Panel unavailable") {
  if (is.null(raster)) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = fallback_label, size = 4.2, hjust = 0.5) +
        theme_void() +
        theme(panel.border = element_rect(color = "grey70", fill = NA, linewidth = 1))
    )
  }

  ggplot() +
    annotation_custom(
      grob = grid::rasterGrob(raster, interpolate = TRUE),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    ) +
    coord_cartesian(expand = FALSE) +
    theme_void()
}

# Bootstrap SHAP stability for defensible feature ranking
compute_shap_stability_bootstrap <- function(X, y,
                                             mtry, min_node_size, sample_fraction,
                                             n_boot = 100, top_k = 20,
                                             num_trees = 500, seed = 123,
                                             num_threads = max(1, parallel::detectCores() - 1)) {

  X <- as.data.frame(X)
  feature_names <- colnames(X)
  n <- nrow(X)
  p <- length(feature_names)

  if (n < 20 || p < 2) {
    message("  SHAP stability skipped: insufficient rows/features")
    return(NULL)
  }

  rank_mat <- matrix(NA_real_, nrow = n_boot, ncol = p, dimnames = list(NULL, feature_names))
  top_mat <- matrix(NA_real_, nrow = n_boot, ncol = p, dimnames = list(NULL, feature_names))
  shap_mat <- matrix(NA_real_, nrow = n_boot, ncol = p, dimnames = list(NULL, feature_names))

  successful <- 0L

  for (b in seq_len(n_boot)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    Xb <- X[idx, , drop = FALSE]
    yb <- y[idx]

    rf_b <- tryCatch(
      ranger(
        x = Xb,
        y = yb,
        num.trees = num_trees,
        mtry = mtry,
        min.node.size = min_node_size,
        sample.fraction = sample_fraction,
        importance = "none",
        num.threads = num_threads,
        seed = seed + b
      ),
      error = function(e) NULL
    )
    if (is.null(rf_b)) next

    ts_b <- tryCatch({
      um_b <- unify(rf_b, Xb)
      treeshap(um_b, Xb)
    }, error = function(e) NULL)
    if (is.null(ts_b)) next

    imp <- colMeans(abs(ts_b$shaps), na.rm = TRUE)
    if (length(imp) == 0 || all(!is.finite(imp))) next
    imp <- imp[feature_names]

    ord <- order(imp, decreasing = TRUE, na.last = NA)
    ranked_names <- names(imp)[ord]

    ranks <- setNames(rep(NA_real_, p), feature_names)
    ranks[ranked_names] <- seq_along(ranked_names)
    rank_mat[b, ] <- ranks

    top_flags <- setNames(rep(0, p), feature_names)
    top_names <- ranked_names[seq_len(min(top_k, length(ranked_names)))]
    top_flags[top_names] <- 1
    top_mat[b, ] <- top_flags

    shap_mat[b, names(imp)] <- imp
    successful <- successful + 1L

    if (b %% 10 == 0 || b == n_boot) {
      message("  SHAP stability bootstrap: ", b, "/", n_boot, " (successful: ", successful, ")")
    }
  }

  if (successful == 0) {
    message("  SHAP stability failed: no successful bootstrap runs")
    return(NULL)
  }

  qfun <- function(x, prob) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NA_real_)
    unname(quantile(x, probs = prob, type = 7))
  }

  sel_freq <- colMeans(top_mat, na.rm = TRUE)
  sel_freq[is.nan(sel_freq)] <- NA_real_

  mean_shap <- colMeans(shap_mat, na.rm = TRUE)
  mean_shap[is.nan(mean_shap)] <- NA_real_

  stability <- tibble(
    Lipid = feature_names,
    SelectionFreqTopK = sel_freq,
    SelectionCI_Low = apply(top_mat, 2, qfun, prob = 0.025),
    SelectionCI_High = apply(top_mat, 2, qfun, prob = 0.975),
    MedianRank = apply(rank_mat, 2, function(x) {
      x <- x[is.finite(x)]
      if (length(x) == 0) NA_real_ else median(x)
    }),
    RankCI_Low = apply(rank_mat, 2, qfun, prob = 0.025),
    RankCI_High = apply(rank_mat, 2, qfun, prob = 0.975),
    MeanAbsSHAP_Boot = mean_shap,
    SHAP_CI_Low = apply(shap_mat, 2, qfun, prob = 0.025),
    SHAP_CI_High = apply(shap_mat, 2, qfun, prob = 0.975),
    SuccessfulBootstraps = successful,
    RequestedBootstraps = n_boot
  ) %>%
    arrange(desc(SelectionFreqTopK), MedianRank, desc(MeanAbsSHAP_Boot))

  stability
}

compute_regression_metrics <- function(obs, pred, train_mean = NULL) {
  obs <- as.numeric(obs)
  pred <- as.numeric(pred)
  keep <- is.finite(obs) & is.finite(pred)
  obs <- obs[keep]
  pred <- pred[keep]

  if (length(obs) == 0) {
    return(list(RMSE = NA_real_, MAE = NA_real_, R2 = NA_real_, R2_Cor = NA_real_))
  }

  if (is.null(train_mean) || !is.finite(train_mean)) {
    train_mean <- mean(obs, na.rm = TRUE)
  }

  sse <- sum((obs - pred)^2)
  sst <- sum((obs - train_mean)^2)

  list(
    RMSE = sqrt(mean((obs - pred)^2)),
    MAE = mean(abs(obs - pred)),
    R2 = if (is.finite(sst) && sst > 0) 1 - (sse / sst) else NA_real_,
    R2_Cor = if (length(obs) > 1) suppressWarnings(cor(obs, pred)^2) else NA_real_
  )
}

residualize_train_test_by_pc <- function(train_X, test_X, train_y, test_y, train_pc, test_pc) {
  train_pc <- as.data.frame(train_pc)
  test_pc <- as.data.frame(test_pc)

  if (ncol(train_pc) == 0) {
    return(list(
      train_X = as.matrix(train_X),
      test_X = as.matrix(test_X),
      train_y = as.numeric(train_y),
      test_y = as.numeric(test_y)
    ))
  }

  df_train_y <- data.frame(Phenotype = as.numeric(train_y), train_pc)
  lm_fit <- lm(Phenotype ~ ., data = df_train_y)
  train_y_adj <- resid(lm_fit)
  test_y_adj <- as.numeric(test_y) - predict(lm_fit, newdata = test_pc)

  X_train <- as.matrix(train_X)
  X_test <- as.matrix(test_X)
  Xd_train <- model.matrix(~ ., data = train_pc)
  Xd_test <- model.matrix(~ ., data = test_pc)
  coef_mat <- qr.solve(Xd_train, X_train)

  train_X_adj <- X_train - (Xd_train %*% coef_mat)
  test_X_adj <- X_test - (Xd_test %*% coef_mat)

  list(train_X = train_X_adj, test_X = test_X_adj, train_y = train_y_adj, test_y = test_y_adj)
}

residualize_full_by_pc <- function(X, y, pc_df = NULL) {
  X <- as.matrix(X)
  y <- as.numeric(y)

  if (is.null(pc_df)) {
    return(list(X = X, y = y))
  }

  pc_df <- as.data.frame(pc_df)
  if (ncol(pc_df) == 0) {
    return(list(X = X, y = y))
  }

  df_y <- data.frame(Phenotype = y, pc_df)
  lm_fit <- lm(Phenotype ~ ., data = df_y)
  y_adj <- resid(lm_fit)

  Xd <- model.matrix(~ ., data = pc_df)
  coef_mat <- qr.solve(Xd, X)
  X_adj <- X - (Xd %*% coef_mat)

  list(X = X_adj, y = y_adj)
}

run_repeated_rf_eval <- function(X, y, pc_df = NULL,
                                 train_prop = 0.8, n_repeats = 10, seed_base = 42,
                                 tune_num_trees = 500, fit_num_trees = 1000,
                                 num_threads = max(1, parallel::detectCores() - 1)) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- nrow(X)

  if (!is.null(pc_df)) {
    pc_df <- as.data.frame(pc_df)
    if (nrow(pc_df) != n) {
      stop("PC matrix rows must match X rows for repeated RF evaluation.")
    }
  }

  metrics_list <- vector("list", n_repeats)
  split_info <- vector("list", n_repeats)
  pars_list <- vector("list", n_repeats)
  successful <- 0L

  for (i in seq_len(n_repeats)) {
    set.seed(seed_base + i)

    if (!is.null(pc_df) && ncol(pc_df) > 0 && nrow(pc_df) > 10) {
      n_centers <- min(5, nrow(pc_df) - 1)
      if (n_centers >= 2) {
        km <- tryCatch(kmeans(pc_df, centers = n_centers, nstart = 20), error = function(e) NULL)
        clusters <- if (!is.null(km)) km$cluster else rep(1, n)
      } else {
        clusters <- rep(1, n)
      }
    } else {
      clusters <- rep(1, n)
    }

    train_idx <- tryCatch(
      createDataPartition(factor(clusters), p = train_prop, list = FALSE),
      error = function(e) integer(0)
    )

    if (length(train_idx) < 5 || length(train_idx) > (n - 5)) {
      fallback_n <- max(5, floor(train_prop * n))
      fallback_n <- min(fallback_n, n - 5)
      train_idx <- sort(sample.int(n, size = fallback_n))
    }
    test_idx <- setdiff(seq_len(n), train_idx)

    train_X_raw <- X[train_idx, , drop = FALSE]
    test_X_raw <- X[test_idx, , drop = FALSE]
    train_y_raw <- y[train_idx]
    test_y_raw <- y[test_idx]

    if (!is.null(pc_df) && ncol(pc_df) > 0) {
      train_pc <- pc_df[train_idx, , drop = FALSE]
      test_pc <- pc_df[test_idx, , drop = FALSE]
      adj <- residualize_train_test_by_pc(
        train_X = train_X_raw, test_X = test_X_raw,
        train_y = train_y_raw, test_y = test_y_raw,
        train_pc = train_pc, test_pc = test_pc
      )
    } else {
      adj <- list(train_X = train_X_raw, test_X = test_X_raw, train_y = train_y_raw, test_y = test_y_raw)
    }

    train_df <- data.frame(as.data.frame(adj$train_X), Phenotype = adj$train_y)
    regr_task <- makeRegrTask(data = train_df, target = "Phenotype")

    tune_res <- tryCatch(
      tuneRanger(
        task = regr_task,
        measure = list(rmse),
        num.trees = tune_num_trees,
        tune.parameters = c("mtry", "min.node.size", "sample.fraction"),
        num.threads = num_threads,
        show.info = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(tune_res)) next

    pars <- tune_res$recommended.pars
    rf_fit <- tryCatch(
      ranger(
        x = as.data.frame(adj$train_X),
        y = adj$train_y,
        num.trees = fit_num_trees,
        mtry = pars$mtry,
        min.node.size = pars$min.node.size,
        sample.fraction = pars$sample.fraction,
        importance = "none",
        num.threads = num_threads,
        seed = seed_base + i
      ),
      error = function(e) NULL
    )
    if (is.null(rf_fit)) next

    pred <- predict(rf_fit, data = as.data.frame(adj$test_X))$predictions
    mets <- compute_regression_metrics(adj$test_y, pred, train_mean = mean(adj$train_y))

    successful <- successful + 1L
    metrics_list[[i]] <- tibble(
      Repeat = i,
      RMSE = mets$RMSE,
      MAE = mets$MAE,
      R2 = mets$R2,
      R2_Cor = mets$R2_Cor,
      mtry = as.numeric(pars$mtry),
      min.node.size = as.numeric(pars$min.node.size),
      sample.fraction = as.numeric(pars$sample.fraction),
      N_Train = length(adj$train_y),
      N_Test = length(adj$test_y)
    )
    split_info[[i]] <- list(
      obs = adj$test_y,
      pred = pred,
      train_y = adj$train_y,
      n_train = length(adj$train_y),
      n_test = length(adj$test_y),
      split_id = i
    )

    if (i %% 2 == 0 || i == n_repeats) {
      message("  RF outer split ", i, "/", n_repeats, " (successful: ", successful, ")")
    }
  }

  metrics_df <- bind_rows(metrics_list)
  if (nrow(metrics_df) == 0) {
    stop("Repeated RF evaluation failed: no successful train/test runs.")
  }

  summary_tbl <- metrics_df %>%
    summarise(
      RMSE_Mean = mean(RMSE, na.rm = TRUE),
      RMSE_SD = sd(RMSE, na.rm = TRUE),
      MAE_Mean = mean(MAE, na.rm = TRUE),
      MAE_SD = sd(MAE, na.rm = TRUE),
      R2_Mean = mean(R2, na.rm = TRUE),
      R2_SD = sd(R2, na.rm = TRUE),
      R2_Cor_Mean = mean(R2_Cor, na.rm = TRUE),
      R2_Cor_SD = sd(R2_Cor, na.rm = TRUE)
    )

  rep_row <- which.min(abs(metrics_df$RMSE - median(metrics_df$RMSE, na.rm = TRUE)))
  rep_repeat <- metrics_df$Repeat[rep_row]
  rep_info <- split_info[[rep_repeat]]

  median_mtry <- max(1, min(ncol(X), round(median(metrics_df$mtry, na.rm = TRUE))))
  median_min_node <- max(1, round(median(metrics_df$min.node.size, na.rm = TRUE)))
  median_sample_fraction <- median(metrics_df$sample.fraction, na.rm = TRUE)

  list(
    metrics = metrics_df,
    summary = as.list(summary_tbl[1, ]),
    representative = rep_info,
    median_params = list(
      mtry = median_mtry,
      min.node.size = median_min_node,
      sample.fraction = median_sample_fraction
    ),
    successful = nrow(metrics_df)
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# 1) LOAD DATA
# ══════════════════════════════════════════════════════════════════════════════

message("\n══════════════════════════════════════════════════════════════")
message("Loading data...")
message("══════════════════════════════════════════════════════════════\n")

# Read raw intensity data (SPATS-fitted)
control_raw <- vroom("data/SPATS_fitted/non_normalized_intensities/Final_subset_control_all_lipids_fitted_phenotype_non_normalized.csv",
                     show_col_types = FALSE) %>%
  dplyr::select(-c(2, 3, 4)) %>%
  dplyr::rename(Compound_Name = 1) %>%
  dplyr::mutate(Condition = "Control")

low_raw <- vroom("data/SPATS_fitted/non_normalized_intensities/Final_subset_lowinput_all_lipids_fitted_phenotype_non_normalized.csv",
                 show_col_types = FALSE) %>%
  dplyr::select(-c(2, 3, 4)) %>%
  dplyr::rename(Compound_Name = 1) %>%
  dplyr::mutate(Condition = "LowInput")

message("Control samples: ", nrow(control_raw))
message("LowInput samples: ", nrow(low_raw))

# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
#
#                         FIGURE 1: LIPIDOMICS LANDSCAPE
#
# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════

message("\n══════════════════════════════════════════════════════════════")
message("FIGURE 1: LIPIDOMICS LANDSCAPE")
message("══════════════════════════════════════════════════════════════\n")

# ─────────────────────────────────────────────────────────────────────────────
# Fig 1B: TIC Composition Stacked Bar
# ─────────────────────────────────────────────────────────────────────────────

# Function to compute class %TIC
classsum_tic_long <- function(df, label,
                              drop_classes = c("DGTS"),
                              valid_classes = classes_with_lyso) {

  class_levels <- unique(valid_classes)
  class_levels <- class_levels[order(nchar(class_levels), decreasing = TRUE)]
  class_pat <- paste0("^(", paste(class_levels, collapse = "|"), ")(?=\\()")

  df %>%
    dplyr::select(-any_of("Condition")) %>%
    pivot_longer(-Compound_Name, names_to = "Lipid", values_to = "Intensity") %>%
    mutate(Class = str_extract(Lipid, class_pat)) %>%
    filter(!is.na(Class), !Class %in% drop_classes) %>%
    group_by(Sample = Compound_Name, Class) %>%
    summarise(sum_int = sum(as.numeric(Intensity), na.rm = TRUE), .groups = "drop") %>%
    mutate(Condition = label) %>%
    group_by(Sample, Condition) %>%
    mutate(pct = 100 * sum_int / sum(sum_int, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(nchar(Sample) >= 5)
}

# Compute class %TIC for both conditions
ctrl_pct <- classsum_tic_long(control_raw, "Control")
low_pct <- classsum_tic_long(low_raw, "LowInput")
pct_tbl <- bind_rows(ctrl_pct, low_pct)

# Mean %TIC per class per condition
pct_mean <- pct_tbl %>%
  group_by(Condition, Class) %>%
  summarise(pct_mean = mean(pct, na.rm = TRUE), .groups = "drop")

# Compute fold changes
fold_tbl <- pct_mean %>%
  pivot_wider(names_from = Condition, values_from = pct_mean, values_fill = 0) %>%
  mutate(
    eps = 1e-8,
    FoldChange = (LowInput + eps) / (Control + eps),
    log2FC = log2(FoldChange),
    Direction = ifelse(log2FC > 0, "Increased", "Decreased")
  )

# Create legend labels
legend_tbl <- pct_mean %>%
  pivot_wider(names_from = Condition, values_from = pct_mean, values_fill = 0) %>%
  mutate(
    legend_label = sprintf("%s (C %.1f%%, L %.1f%%)", Class, Control, LowInput)
  )
label_map <- setNames(legend_tbl$legend_label, legend_tbl$Class)

# Order classes
classes_all <- sort(unique(pct_mean$Class))
pct_mean <- pct_mean %>% mutate(Class = factor(Class, levels = classes_all))

# Fig 1B (%TIC) moved to Supplementary Figure S2A
suppfig_s2a_tic <- ggplot(pct_mean, aes(x = pct_mean, y = Condition, fill = Class)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(pct_mean >= 3, sprintf("%.1f%%", pct_mean), "")),
            position = position_stack(vjust = 0.5),
            colour = "white", size = 3.5, fontface = "bold") +
  scale_fill_manual(
    values = class_colors[classes_all],
    breaks = classes_all,
    labels = label_map[classes_all],
    name = NULL
  ) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 100, 25),
                     labels = function(x) paste0(x, "%")) +
  labs(x = "% of TIC (mean across samples)", y = NULL,
       title = "A) Class %TIC Composition") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.4, "cm"),
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.major.y = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5)
  )

# Helper: closure + zero replacement for compositional transforms
close_replace_rows <- function(mat, delta_frac = 0.5) {
  mat <- as.matrix(mat)
  out <- matrix(NA_real_, nrow = nrow(mat), ncol = ncol(mat))
  colnames(out) <- colnames(mat)
  rownames(out) <- rownames(mat)

  for (i in seq_len(nrow(mat))) {
    x <- mat[i, ]
    x[!is.finite(x)] <- 0
    x[x < 0] <- 0
    if (sum(x) <= 0) next

    if (any(x == 0)) {
      pos <- x[x > 0]
      if (length(pos) > 0) {
        delta <- min(pos) * delta_frac
        x[x == 0] <- delta
      }
    }

    out[i, ] <- x / sum(x)
  }
  out
}

# Helper: bootstrap descriptive contrast (LowInput - Control)
bootstrap_context_contrast <- function(mat, cond_vec, n_boot = 500, seed = 2026) {
  mat <- as.matrix(mat)
  keep <- complete.cases(mat) & !is.na(cond_vec)
  mat <- mat[keep, , drop = FALSE]
  cond_vec <- cond_vec[keep]

  idx_ctrl <- which(cond_vec == "Control")
  idx_low <- which(cond_vec == "LowInput")
  if (length(idx_ctrl) < 5 || length(idx_low) < 5) return(NULL)

  est <- colMeans(mat[idx_low, , drop = FALSE], na.rm = TRUE) -
    colMeans(mat[idx_ctrl, , drop = FALSE], na.rm = TRUE)

  boot <- matrix(NA_real_, nrow = n_boot, ncol = ncol(mat))
  colnames(boot) <- colnames(mat)

  set.seed(seed)
  for (b in seq_len(n_boot)) {
    ctrl_b <- sample(idx_ctrl, size = length(idx_ctrl), replace = TRUE)
    low_b <- sample(idx_low, size = length(idx_low), replace = TRUE)
    boot[b, ] <- colMeans(mat[low_b, , drop = FALSE], na.rm = TRUE) -
      colMeans(mat[ctrl_b, , drop = FALSE], na.rm = TRUE)
  }

  tibble(
    Feature = colnames(mat),
    Effect = est,
    CI_Low = apply(boot, 2, function(x) quantile(x, 0.025, na.rm = TRUE)),
    CI_High = apply(boot, 2, function(x) quantile(x, 0.975, na.rm = TRUE))
  )
}

# Build class-by-sample composition matrix (%TIC -> composition)
class_comp_wide <- pct_tbl %>%
  dplyr::select(Sample, Condition, Class, pct) %>%
  pivot_wider(names_from = Class, values_from = pct, values_fill = 0)

class_cols_comp <- setdiff(names(class_comp_wide), c("Sample", "Condition"))
class_comp_raw <- as.matrix(class_comp_wide[, class_cols_comp, drop = FALSE]) / 100
rownames(class_comp_raw) <- class_comp_wide$Sample
class_comp <- close_replace_rows(class_comp_raw, delta_frac = 0.5)

# CLR contrast (main Figure 1B)
clr_class <- log(class_comp) - rowMeans(log(class_comp), na.rm = TRUE)
clr_tbl <- bootstrap_context_contrast(clr_class, class_comp_wide$Condition, n_boot = 500, seed = 1101) %>%
  dplyr::rename(Class = Feature) %>%
  mutate(
    AbsEffect = abs(Effect),
    Class = factor(Class, levels = rev(Class[order(AbsEffect, decreasing = TRUE)]))
  )

# ALR contrast (Supplementary Figure S2B)
alr_ref_class <- if ("TG" %in% colnames(class_comp)) "TG" else {
  names(sort(colMeans(class_comp, na.rm = TRUE), decreasing = TRUE))[1]
}
alr_num <- setdiff(colnames(class_comp), alr_ref_class)
alr_class <- log(sweep(class_comp[, alr_num, drop = FALSE], 1, class_comp[, alr_ref_class], "/"))
alr_tbl <- bootstrap_context_contrast(alr_class, class_comp_wide$Condition, n_boot = 500, seed = 1201) %>%
  mutate(
    Class = str_replace(Feature, paste0("_vs_", alr_ref_class, "$"), ""),
    AbsEffect = abs(Effect),
    Class = factor(Class, levels = rev(Class[order(AbsEffect, decreasing = TRUE)]))
  )

write.csv(clr_tbl, "table/supp/SuppTable_S2A_Class_CLR_Contrast.csv", row.names = FALSE)
write.csv(alr_tbl, "table/supp/SuppTable_S2B_Class_ALR_Contrast.csv", row.names = FALSE)

# Main Figure 1B: CLR effect-size ranked contrast
fig1b <- ggplot(clr_tbl, aes(x = Class, y = Effect, fill = Class)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.5) +
  geom_col(width = 0.72, color = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = CI_Low, ymax = CI_High), width = 0.2, linewidth = 0.6) +
  coord_flip() +
  scale_fill_manual(values = class_colors, na.value = "grey70") +
  scale_y_continuous(labels = function(x) format_small_num(x, digits = 2, threshold = 0.01)) +
  labs(
    x = NULL,
    y = "CLR contrast (LowInput - Control)",
    title = "Class-Level Compositional Contrast (CLR)",
    subtitle = "Ordered by |effect size|"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    panel.grid.major.y = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5)
  )

# Supplementary Figure S2B: ALR sensitivity contrast
suppfig_s2b_alr <- ggplot(alr_tbl, aes(x = Class, y = Effect, color = Class)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.5) +
  geom_errorbar(aes(ymin = CI_Low, ymax = CI_High), width = 0.2, linewidth = 0.7) +
  geom_point(size = 2.2) +
  coord_flip() +
  scale_color_manual(values = class_colors, na.value = "grey70") +
  scale_y_continuous(labels = function(x) format_small_num(x, digits = 2, threshold = 0.01)) +
  labs(
    x = NULL,
    y = "ALR contrast (LowInput - Control)",
    title = paste0("B) Class ALR Contrast (Ref = ", alr_ref_class, ")")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.major.y = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5)
  )

message("Created: Fig 1B (CLR effect-size ranked contrast)")

# ─────────────────────────────────────────────────────────────────────────────
# Fig 1C: log2FC Bar Chart
# ─────────────────────────────────────────────────────────────────────────────

fig1c <- ggplot(fold_tbl, aes(x = reorder(Class, log2FC), y = log2FC, fill = Direction)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  scale_fill_manual(values = c("Increased" = "#D55E00", "Decreased" = "#0072B2")) +
  labs(x = "Lipid Class", y = expression(log[2]*"(LowInput / Control)"),
       title = "Class-Level Changes Under Low-Input") +
  coord_flip() +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.major.y = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5)
  )

message("Created: Fig 1C (log2FC)")

# ─────────────────────────────────────────────────────────────────────────────
# Fig 1D: LION Enrichment Analysis
# ─────────────────────────────────────────────────────────────────────────────

# Read LION enrichment data
lion_file <- "table/Linex2/LION-enrichment.csv"

if (file.exists(lion_file)) {
  lion_df <- vroom(lion_file, show_col_types = FALSE)

  # Filter for significant (p < 0.05)
  lion_sig <- lion_df %>%
    dplyr::filter(`p-value` < 0.05) %>%
    dplyr::rename(
      p_value = `p-value`,
      FDR_q = `FDR q-value`,
      Description = Discription  # Fix typo in source data
    ) %>%
    dplyr::mutate(
      logP = -log10(p_value),
      logQ = -log10(FDR_q),
      Description = forcats::fct_reorder(Description, logQ)
    )

  # Create LION enrichment plot
  fig1d <- ggplot(lion_sig, aes(x = logQ, y = Description)) +
    # FDR = 0.05 threshold line
    geom_vline(xintercept = -log10(0.05), color = "red", linetype = "dashed", linewidth = 0.8) +
    # Facet by up/down regulation
    facet_grid(Regulated ~ ., scales = "free_y", space = "free_y") +
    # Points sized by annotated lipids, colored by -log10(p)
    geom_point(aes(size = Annotated, color = logP), alpha = 0.8) +
    scale_color_viridis_c(
      name = expression(-log[10](p)),
      option = "C"
    ) +
    scale_size_continuous(name = "Annotated\nlipids", range = c(2, 8)) +
    labs(
      x = expression(-log[10](FDR~q~value)),
      y = NULL,
      title = "LION Enrichment: LowInput vs Control"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      axis.title.x = element_text(face = "bold", size = 11),
      axis.text.y = element_text(size = 9, color = "black"),
      axis.text.x = element_text(size = 10, color = "black"),
      strip.text = element_text(face = "bold", size = 10),
      panel.grid.major.y = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5)
    )

  message("Created: Fig 1D (LION Enrichment)")
} else {
  # Placeholder if LION file not found
  fig1d <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "LION Enrichment\n(table/Linex2/LION-enrichment.csv not found)",
             size = 4, hjust = 0.5) +
    theme_void() +
    theme(panel.border = element_rect(color = "grey70", fill = NA, linewidth = 1))
  message("Warning: LION file not found, using placeholder")
}

# ─────────────────────────────────────────────────────────────────────────────
# TIC Normalization & Lipid Processing (used for multiple figures)
# ─────────────────────────────────────────────────────────────────────────────

# TIC normalize function
tic_normalize <- function(df) {
  X <- df %>%
    dplyr::select(-Compound_Name, -Condition) %>%
    dplyr::select(where(is.numeric)) %>%
    as.matrix()

  rs <- rowSums(X, na.rm = TRUE)
  rs[rs == 0 | is.na(rs)] <- NA_real_
  Xpct <- sweep(X, 1, rs, FUN = "/") * 100
  Xpct[is.na(Xpct)] <- 0

  bind_cols(df %>% dplyr::select(Compound_Name, Condition), as.data.frame(Xpct))
}

# Subset to valid lipid classes
subset_lipids <- function(df, keep_classes = valid_classes, drop_classes = c("DGTS")) {
  lipid_cols <- names(df)[grepl("\\(", names(df))]
  keep_pat <- paste0("^(", paste(keep_classes, collapse = "|"), ")\\(")
  keep_cols <- lipid_cols[grepl(keep_pat, lipid_cols)]

  if (length(drop_classes) > 0) {
    drop_pat <- paste0("^(", paste(drop_classes, collapse = "|"), ")\\(")
    keep_cols <- keep_cols[!grepl(drop_pat, keep_cols)]
  }

  df %>% dplyr::select(Compound_Name, Condition, all_of(keep_cols))
}

# Process data for PCA (will be used in supplementary)
control_sub <- subset_lipids(control_raw)
low_sub <- subset_lipids(low_raw)

control_tic <- tic_normalize(control_sub)
low_tic <- tic_normalize(low_sub)

# Align columns
shared_cols <- intersect(
  setdiff(names(control_tic), c("Compound_Name", "Condition")),
  setdiff(names(low_tic), c("Compound_Name", "Condition"))
)

control_tic <- control_tic %>% dplyr::select(Compound_Name, Condition, all_of(shared_cols))
low_tic <- low_tic %>% dplyr::select(Compound_Name, Condition, all_of(shared_cols))

# ─────────────────────────────────────────────────────────────────────────────
# Fig 1E: Chemical Shifts (Double Bond vs Total Carbon)
# ─────────────────────────────────────────────────────────────────────────────

calculate_class_chem <- function(df_tic, condition_name) {
  lipid_cols <- setdiff(names(df_tic), c("Compound_Name", "Condition"))
  lipid_cols <- lipid_cols[grepl("\\(", lipid_cols)]

  X <- df_tic %>%
    dplyr::select(all_of(lipid_cols)) %>%
    mutate(across(everything(), ~suppressWarnings(as.numeric(.x))))

  mean_pct <- colMeans(as.matrix(X), na.rm = TRUE)

  tibble(Lipid = lipid_cols) %>%
    mutate(
      Class = get_lipid_class(Lipid),
      MeanPctTIC = mean_pct[Lipid],
      Condition = condition_name
    ) %>%
    bind_cols(purrr::map_dfr(.$Lipid, parse_lipid_features)) %>%
    filter(!is.na(total_c), Class %in% valid_classes) %>%
    group_by(Class, Condition) %>%
    summarise(
      WeightedMeanC = weighted.mean(total_c, MeanPctTIC, na.rm = TRUE),
      WeightedMeanDB = weighted.mean(total_db, MeanPctTIC, na.rm = TRUE),
      n_species = n(),
      .groups = "drop"
    )
}

chem_control <- calculate_class_chem(control_tic, "Control")
chem_low <- calculate_class_chem(low_tic, "LowInput")
chem_summary <- bind_rows(chem_control, chem_low)

fig1e <- ggplot(chem_summary, aes(x = WeightedMeanC, y = WeightedMeanDB, color = Class, shape = Condition)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_path(aes(group = Class), arrow = arrow(type = "closed", length = unit(0.12, "inches")),
            linewidth = 0.8, alpha = 0.7) +
  scale_color_manual(values = class_colors) +
  scale_shape_manual(values = c("Control" = 16, "LowInput" = 17)) +
  labs(
    x = "Weighted Mean Total Carbons",
    y = "Weighted Mean Double Bonds",
    title = "Chemical Shifts Between Conditions"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5)
  )

message("Created: Fig 1E (Chemical shifts)")

# ─────────────────────────────────────────────────────────────────────────────
# Assemble Figure 1
# ─────────────────────────────────────────────────────────────────────────────

# Fig 1A image panel (working pipeline)
fig1a_candidates <- c(
  "fig/main/individual_figs/Fig1A_working_pipeline.png",
  "fig/main/individual_figs/Fig1A_working_pipeline.jpeg",
  "fig/main/individual_figs/Fig1A_working_pipeline.jpg"
)
fig1a_path <- fig1a_candidates[file.exists(fig1a_candidates)][1]
fig1a_raster <- NULL
fig1a_err <- NULL

if (!is.na(fig1a_path)) {
  # 1) Try PNG reader first (works for real PNG files)
  if (requireNamespace("png", quietly = TRUE)) {
    fig1a_raster <- tryCatch(png::readPNG(fig1a_path), error = function(e) NULL)
  }

  # 2) Try JPEG reader if PNG reader failed
  if (is.null(fig1a_raster) && requireNamespace("jpeg", quietly = TRUE)) {
    fig1a_raster <- tryCatch(jpeg::readJPEG(fig1a_path), error = function(e) NULL)
  }

  # 3) Last fallback on macOS: convert to temporary PNG via `sips`
  if (is.null(fig1a_raster) && nzchar(Sys.which("sips")) && requireNamespace("png", quietly = TRUE)) {
    tmp_png <- tempfile(fileext = ".png")
    conv_status <- suppressWarnings(
      tryCatch(
        system2("sips", c("-s", "format", "png", fig1a_path, "--out", tmp_png),
                stdout = TRUE, stderr = TRUE),
        error = function(e) NULL
      )
    )
    if (!is.null(conv_status) && file.exists(tmp_png)) {
      fig1a_raster <- tryCatch(png::readPNG(tmp_png), error = function(e) NULL)
    }
  }

  if (is.null(fig1a_raster)) {
    fig1a_err <- "Could not decode image (check png/jpeg package availability)"
  }
} else {
  fig1a_err <- paste0("Image file not found. Checked: ", paste(fig1a_candidates, collapse = ", "))
}

if (!is.null(fig1a_raster)) {
  fig1a_panel <- ggplot() +
    annotation_custom(
      grob = grid::rasterGrob(fig1a_raster, interpolate = TRUE),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    ) +
    coord_cartesian(expand = FALSE) +
    theme_void()
} else {
  fig1a_panel <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = paste0("Fig 1A image not loaded:\n", fig1a_path, "\n(", fig1a_err, ")"),
             size = 4.2, hjust = 0.5) +
    theme_void() +
    theme(panel.border = element_rect(color = "grey70", fill = NA, linewidth = 1))
}

# Combine panels
fig1_row1 <- fig1a_panel
fig1_row2 <- fig1b + fig1e + plot_layout(widths = c(1.05, 1))
fig1_row3 <- fig1d

figure1 <- (fig1_row1 / fig1_row2 / fig1_row3) +
  plot_layout(heights = c(0.85, 1, 1.15)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold"))

save_fig(figure1, "Figure1_Lipidomics_Landscape.png", w = 14, h = 16)
save_fig(figure1, "Figure1_Lipidomics_Landscape.pdf", w = 14, h = 16)

message("\n✓ Figure 1 complete!\n")

# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
#
#     SUPPLEMENTARY FIGURE S2: COMPOSITIONAL CONTEXT CONTRASTS
#
# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════

message("\n══════════════════════════════════════════════════════════════")
message("SUPPLEMENTARY FIGURE S2: COMPOSITIONAL CONTEXT CONTRASTS")
message("══════════════════════════════════════════════════════════════\n")

suppfig_s2_comp <- suppfig_s2a_tic + suppfig_s2b_alr +
  plot_layout(widths = c(1.2, 1))

save_fig(suppfig_s2_comp, "SuppFig_S2_Compositional_Contrasts.png", w = 15, h = 7, subdir = "supp")
save_fig(suppfig_s2_comp, "SuppFig_S2_Compositional_Contrasts.pdf", w = 15, h = 7, subdir = "supp")

message("✓ Supplementary Figure S2 (Compositional) complete!\n")


# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
#
#     SUPPLEMENTARY TABLE: SPECIES COMPOSITION BY CLASS
#
# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════

message("\n══════════════════════════════════════════════════════════════")
message("SUPPLEMENTARY TABLE: SPECIES COMPOSITION BY CLASS")
message("══════════════════════════════════════════════════════════════\n")

# Function to get top species within each class
get_species_composition <- function(df_tic, condition_name, top_n_per_class = 10) {

  lipid_cols <- setdiff(names(df_tic), c("Compound_Name", "Condition"))
  lipid_cols <- lipid_cols[grepl("\\(", lipid_cols)]

  # Calculate mean % TIC for each lipid
  X <- df_tic %>%
    dplyr::select(all_of(lipid_cols)) %>%
    mutate(across(everything(), ~suppressWarnings(as.numeric(.x))))

  mean_pct <- colMeans(as.matrix(X), na.rm = TRUE)
  sd_pct <- apply(as.matrix(X), 2, sd, na.rm = TRUE)

  # Build species table
  species_tbl <- tibble(
    Lipid = lipid_cols,
    Class = get_lipid_class(lipid_cols),
    Mean_PctTIC = mean_pct[lipid_cols],
    SD_PctTIC = sd_pct[lipid_cols],
    Condition = condition_name
  ) %>%
    filter(!is.na(Class), Class != "Other") %>%
    # Parse carbon and double bond info
    bind_cols(purrr::map_dfr(.$Lipid, parse_lipid_features)) %>%
    arrange(Class, desc(Mean_PctTIC))

  # Get top N per class
  top_species <- species_tbl %>%
    group_by(Class) %>%
    slice_head(n = top_n_per_class) %>%
    mutate(Rank_in_Class = row_number()) %>%
    ungroup()

  top_species
}

# Get species composition for both conditions
species_ctrl <- get_species_composition(control_tic, "Control", top_n_per_class = 10)
species_low <- get_species_composition(low_tic, "LowInput", top_n_per_class = 10)

message("Top species per class:")
message("  Control: ", nrow(species_ctrl), " entries across ", n_distinct(species_ctrl$Class), " classes")
message("  LowInput: ", nrow(species_low), " entries across ", n_distinct(species_low$Class), " classes")

# Combine into a single table
species_combined <- bind_rows(species_ctrl, species_low) %>%
  dplyr::select(Condition, Class, Rank_in_Class, Lipid, Mean_PctTIC, SD_PctTIC, total_c, total_db) %>%
  dplyr::rename(
    `Mean %TIC` = Mean_PctTIC,
    `SD %TIC` = SD_PctTIC,
    `Total Carbons` = total_c,
    `Double Bonds` = total_db
  ) %>%
  mutate(
    `Mean %TIC` = round(`Mean %TIC`, 4),
    `SD %TIC` = round(`SD %TIC`, 4)
  ) %>%
  arrange(Class, Condition, Rank_in_Class)

# Save combined table
write.csv(species_combined, "table/supp/SuppTable_Species_Composition_by_Class.csv", row.names = FALSE)
message("✓ Saved: table/supp/SuppTable_Species_Composition_by_Class.csv")

# Create a wide format comparing Control vs LowInput for the same species
species_comparison <- species_ctrl %>%
  dplyr::select(Class, Lipid, Mean_PctTIC, total_c, total_db) %>%
  dplyr::rename(Control_MeanPct = Mean_PctTIC) %>%
  full_join(
    species_low %>%
      dplyr::select(Class, Lipid, Mean_PctTIC) %>%
      dplyr::rename(LowInput_MeanPct = Mean_PctTIC),
    by = c("Class", "Lipid")
  ) %>%
  mutate(
    Control_MeanPct = replace_na(Control_MeanPct, 0),
    LowInput_MeanPct = replace_na(LowInput_MeanPct, 0),
    log2FC = log2((LowInput_MeanPct + 0.001) / (Control_MeanPct + 0.001)),
    Direction = case_when(
      log2FC > 0.5 ~ "Increased in LowInput",
      log2FC < -0.5 ~ "Decreased in LowInput",
      TRUE ~ "No change"
    )
  ) %>%
  arrange(Class, desc(abs(log2FC)))

write.csv(species_comparison, "table/supp/SuppTable_Species_Comparison_Ctrl_vs_Low.csv", row.names = FALSE)
message("✓ Saved: table/supp/SuppTable_Species_Comparison_Ctrl_vs_Low.csv")

# Summary by class: how many species increased/decreased
class_change_summary <- species_comparison %>%
  group_by(Class, Direction) %>%
  summarise(N_Species = n(), .groups = "drop") %>%
  pivot_wider(names_from = Direction, values_from = N_Species, values_fill = 0)

write.csv(class_change_summary, "table/supp/SuppTable_Species_Change_Summary.csv", row.names = FALSE)
message("✓ Saved: table/supp/SuppTable_Species_Change_Summary.csv")

message("\n✓ Species Composition Tables complete!\n")


# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
#
#     SUPPLEMENTARY FIGURE S8: PARTIAL CORRELATION ANALYSIS
#
# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════

message("\n══════════════════════════════════════════════════════════════")
message("SUPPLEMENTARY FIGURE S8: PARTIAL CORRELATION ANALYSIS")
message("══════════════════════════════════════════════════════════════\n")

# Function to compute partial correlations controlling for total lipid content
compute_partial_correlations <- function(df_tic, condition_name, top_classes = NULL) {

  lipid_cols <- setdiff(names(df_tic), c("Compound_Name", "Condition"))
  lipid_cols <- lipid_cols[grepl("\\(", lipid_cols)]

  # Get class percentages (as in the correlation heatmap)
  X <- df_tic %>%
    dplyr::select(all_of(lipid_cols)) %>%
    mutate(across(everything(), ~suppressWarnings(as.numeric(.x)))) %>%
    as.matrix()

  # Assign classes directly from lipid names (includes LPC/LPE when present)
  classes <- str_extract(lipid_cols, "^[A-Za-z0-9]+(?=\\()")

  # Calculate class-level sums
  unique_classes <- sort(unique(classes[!is.na(classes)]))
  if (!is.null(top_classes)) {
    unique_classes <- intersect(unique_classes, top_classes)
  }

  class_mat <- sapply(unique_classes, function(cl) {
    idx <- which(classes == cl)
    if (length(idx) > 0) rowSums(X[, idx, drop = FALSE], na.rm = TRUE) else rep(0, nrow(X))
  })
  colnames(class_mat) <- unique_classes

  # Convert to proportions (% TIC)
  total_tic <- rowSums(class_mat, na.rm = TRUE)
  total_tic[total_tic == 0] <- NA
  class_pct <- sweep(class_mat, 1, total_tic, "/") * 100
  class_pct[is.na(class_pct)] <- 0

  # Remove samples with NA
  complete_idx <- complete.cases(class_pct)
  class_pct <- class_pct[complete_idx, ]
  total_tic_clean <- total_tic[complete_idx]

  if (nrow(class_pct) < 10) {
    message("  Not enough complete samples for partial correlation")
    return(NULL)
  }

  n_classes <- ncol(class_pct)

  # Regular Pearson correlation
  regular_cor <- cor(class_pct, use = "pairwise.complete.obs", method = "pearson")

  # Partial correlation controlling for total TIC
  # Using ppcor package
  partial_cor_mat <- matrix(NA, n_classes, n_classes)
  partial_pval_mat <- matrix(NA, n_classes, n_classes)
  rownames(partial_cor_mat) <- colnames(partial_cor_mat) <- unique_classes
  rownames(partial_pval_mat) <- colnames(partial_pval_mat) <- unique_classes

  for (i in 1:(n_classes - 1)) {
    for (j in (i + 1):n_classes) {
      tryCatch({
        # Partial correlation: x = class i, y = class j, z = total TIC
        pc_result <- pcor.test(class_pct[, i], class_pct[, j], total_tic_clean, method = "pearson")
        partial_cor_mat[i, j] <- pc_result$estimate
        partial_cor_mat[j, i] <- pc_result$estimate
        partial_pval_mat[i, j] <- pc_result$p.value
        partial_pval_mat[j, i] <- pc_result$p.value
      }, error = function(e) {
        partial_cor_mat[i, j] <<- NA
        partial_cor_mat[j, i] <<- NA
      })
    }
  }
  diag(partial_cor_mat) <- 1
  diag(partial_pval_mat) <- 0

  list(
    regular = regular_cor,
    partial = partial_cor_mat,
    partial_pval = partial_pval_mat,
    condition = condition_name,
    n_samples = nrow(class_pct)
  )
}

# Compute for both conditions
message("Computing partial correlations (controlling for total TIC)...")
control_sub_corr <- subset_lipids(
  control_raw,
  keep_classes = classes_with_lyso,
  drop_classes = c("DGTS")
)
low_sub_corr <- subset_lipids(
  low_raw,
  keep_classes = classes_with_lyso,
  drop_classes = c("DGTS")
)

control_tic_corr <- tic_normalize(control_sub_corr)
low_tic_corr <- tic_normalize(low_sub_corr)
shared_cols_corr <- intersect(
  setdiff(names(control_tic_corr), c("Compound_Name", "Condition")),
  setdiff(names(low_tic_corr), c("Compound_Name", "Condition"))
)
control_tic_corr <- control_tic_corr %>% dplyr::select(Compound_Name, Condition, all_of(shared_cols_corr))
low_tic_corr <- low_tic_corr %>% dplyr::select(Compound_Name, Condition, all_of(shared_cols_corr))

pcor_control <- compute_partial_correlations(control_tic_corr, "Control")
pcor_lowinput <- compute_partial_correlations(low_tic_corr, "LowInput")

# Align classes between conditions so panel ordering is consistent
if (!is.null(pcor_control) && !is.null(pcor_lowinput)) {
  common_classes <- intersect(colnames(pcor_control$partial), colnames(pcor_lowinput$partial))
} else {
  common_classes <- NULL
}

# Identify class pairs with significant sign flip between Control and LowInput
get_signflip_pair_keys <- function(pcor_control, pcor_lowinput, classes,
                                   alpha = 0.05) {
  if (is.null(pcor_control) || is.null(pcor_lowinput) || length(classes) == 0) {
    return(character(0))
  }

  pairs <- expand.grid(Class1 = classes, Class2 = classes, stringsAsFactors = FALSE) %>%
    filter(Class1 != Class2) %>%
    mutate(
      pair_key = ifelse(Class1 < Class2,
                        paste(Class1, Class2, sep = "||"),
                        paste(Class2, Class1, sep = "||")),
      r_ctrl = pmap_dbl(list(Class1, Class2), ~ pcor_control$partial[..1, ..2]),
      p_ctrl = pmap_dbl(list(Class1, Class2), ~ pcor_control$partial_pval[..1, ..2]),
      r_low  = pmap_dbl(list(Class1, Class2), ~ pcor_lowinput$partial[..1, ..2]),
      p_low  = pmap_dbl(list(Class1, Class2), ~ pcor_lowinput$partial_pval[..1, ..2]),
      sig_ctrl = !is.na(p_ctrl) & p_ctrl < alpha,
      sig_low  = !is.na(p_low) & p_low < alpha,
      flip_sig = sig_ctrl & sig_low & !is.na(r_ctrl) & !is.na(r_low) & (r_ctrl * r_low < 0)
    )

  unique(pairs$pair_key[pairs$flip_sig])
}

# Function to create partial-correlation heatmap with sign-flip highlighting
create_partial_correlation_plot <- function(pcor_result, title_prefix, class_order = NULL,
                                            flip_pair_keys = character(0),
                                            alpha = 0.05, sig_only = TRUE) {

  if (is.null(pcor_result)) {
    return(ggplot() + theme_void() + ggtitle(paste(title_prefix, "- Not available")))
  }

  classes <- colnames(pcor_result$partial)
  if (!is.null(class_order)) {
    classes <- intersect(class_order, classes)
  }
  n_samples <- pcor_result$n_samples

  # Build data frame for plotting
  plot_data <- expand.grid(Class1 = classes, Class2 = classes, stringsAsFactors = FALSE) %>%
    mutate(
      idx1 = match(Class1, classes),
      idx2 = match(Class2, classes)
    )

  # Partial correlation only
  plot_data <- plot_data %>%
    rowwise() %>%
    mutate(
      r = ifelse(idx1 == idx2, 1, pcor_result$partial[Class1, Class2]),
      p_value = ifelse(idx1 == idx2, 0, pcor_result$partial_pval[Class1, Class2])
    ) %>%
    ungroup() %>%
    mutate(
      IsDiag = idx1 == idx2,
      IsUpper = idx1 < idx2,
      pair_key = ifelse(Class1 < Class2,
                        paste(Class1, Class2, sep = "||"),
                        paste(Class2, Class1, sep = "||")),
      MarkSig = idx1 != idx2 & !is.na(p_value) & p_value < alpha,
      MarkFlip = MarkSig & pair_key %in% flip_pair_keys,
      r_plot = r,
      sig_stars = case_when(
        is.na(p_value) ~ "",
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ ""
      ),
      label = case_when(
        IsDiag ~ sprintf("%.2f", r),
        IsUpper ~ sig_stars,
        TRUE ~ sprintf("%.2f", r)
      ),
      label_color = ifelse(IsDiag, "white", "black"),
      # Factor for ordering
      Class1 = factor(Class1, levels = classes),
      Class2 = factor(Class2, levels = rev(classes))
    )

  # Create the plot
  p <- ggplot(plot_data, aes(x = Class1, y = Class2, fill = r_plot)) +
    geom_tile(color = "white", linewidth = 0.5) +
    # Diagonal black tiles
    geom_tile(
      data = plot_data %>% filter(IsDiag),
      aes(x = Class1, y = Class2),
      inherit.aes = FALSE,
      fill = "black",
      color = "white",
      linewidth = 0.5
    ) +
    # Significant sign-flip pairs between conditions: pink square outline
    geom_point(data = plot_data %>% filter(MarkFlip),
               aes(x = Class1, y = Class2),
               shape = 0, size = 8.6, stroke = 1.55, color = "#C51B7D", fill = NA, inherit.aes = FALSE) +
    geom_text(aes(label = label, color = label_color), size = 3) +
    scale_color_identity(guide = "none") +
    scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                         midpoint = 0, limits = c(-1, 1),
                         name = "Partial r") +
    labs(
      title = paste0(title_prefix, " (n = ", n_samples, ")"),
      subtitle = "Partial correlation (controlling for total TIC); all values shown",
      x = NULL, y = NULL,
      caption = paste0(
        "Upper triangle stars: * p<0.05, ** p<0.01, *** p<0.001. ",
        "Pink square: sign flip vs other condition, significant in both (p<", alpha, ")"
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
      axis.text.y = element_text(color = "black", size = 10),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey30"),
      plot.caption = element_text(size = 9, color = "grey40", hjust = 0.5),
      legend.position = "right",
      panel.grid = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5)
    ) +
    coord_fixed()

  p
}

# Create comparison plots (partial correlations only)
flip_pair_keys <- get_signflip_pair_keys(
  pcor_control, pcor_lowinput,
  classes = common_classes,
  alpha = 0.05
)

suppfig_pcor_ctrl <- create_partial_correlation_plot(
  pcor_control, "A) Control",
  class_order = common_classes,
  flip_pair_keys = flip_pair_keys,
  alpha = 0.05,
  sig_only = FALSE
)
suppfig_pcor_low <- create_partial_correlation_plot(
  pcor_lowinput, "B) LowInput",
  class_order = common_classes,
  flip_pair_keys = flip_pair_keys,
  alpha = 0.05,
  sig_only = FALSE
)

# Combine into single figure
suppfig_partial_cor <- suppfig_pcor_ctrl + suppfig_pcor_low +
  plot_layout(ncol = 2, guides = "collect")

save_fig(suppfig_partial_cor, "SuppFig_S6_Partial_Correlations.png", w = 16, h = 8, subdir = "supp")
save_fig(suppfig_partial_cor, "SuppFig_S6_Partial_Correlations.pdf", w = 16, h = 8, subdir = "supp")

message("\n✓ Partial Correlation Analysis complete!\n")


# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
#
#     SUPPLEMENTARY FIGURE S1: QC - RUN ORDER, SERRF RSD, AND PCA
#
# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════

message("\n══════════════════════════════════════════════════════════════")
message("SUPPLEMENTARY FIGURE S1: QC - RUN ORDER, SERRF RSD, PCA, AND SPATS RESIDUALS")
message("══════════════════════════════════════════════════════════════\n")

# ─────────────────────────────────────────────────────────────────────────────
# S1A & S1B: Run Order Stability Plots (Control & LowInput)
# Shows intensity stability across LC-MS runs for QC assessment
# ─────────────────────────────────────────────────────────────────────────────

classify_qc_injection_type <- function(x) {
  case_when(
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
    message("  File not found for ", condition_name, ": ",
            paste(file_candidates, collapse = " | "))
    return(ggplot() + theme_void() +
             ggtitle(paste0(condition_name, ": Data not available")))
  }

  dat <- suppressMessages(vroom(file_path, show_col_types = FALSE, delim = ",", progress = FALSE))
  run_cols <- names(dat)[str_detect(names(dat), "Run\\d+")]
  if (length(run_cols) == 0) {
    message("  No run columns detected in: ", file_path)
    return(ggplot() + theme_void() +
             ggtitle(paste0(condition_name, ": No run columns found")))
  }

  run_meta <- tibble(Injection = run_cols) %>%
    mutate(
      Run = as.integer(str_extract(Injection, "(?<=Run)\\d+")),
      Type = classify_qc_injection_type(Injection),
      QC_Label = case_when(
        Type == "QC" & str_detect(Injection, regex("Pre", ignore_case = TRUE)) ~ "QC Pre",
        Type == "QC" & str_detect(Injection, regex("Post", ignore_case = TRUE)) ~ "QC Post",
        TRUE ~ NA_character_
      )
    ) %>%
    mutate(Run = ifelse(is.na(Run), row_number(), Run)) %>%
    arrange(Run)

  run_matrix <- dat %>%
    dplyr::select(all_of(run_meta$Injection)) %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(.x)))) %>%
    as.matrix()
  mode(run_matrix) <- "numeric"
  keep_rows <- rowSums(is.finite(run_matrix), na.rm = TRUE) > 0
  run_matrix <- run_matrix[keep_rows, , drop = FALSE]
  run_matrix <- run_matrix[, run_meta$Injection, drop = FALSE]

  plot_df <- run_meta %>%
    mutate(TIC = colSums(run_matrix, na.rm = TRUE)) %>%
    filter(Type %in% c("Sample", "Check", "QC", "Blank", "ISTD")) %>%
    arrange(Run)

  qc_tic <- plot_df %>% filter(Type == "QC") %>% pull(TIC)
  qc_cv <- if (length(qc_tic) >= 2 && mean(qc_tic, na.rm = TRUE) > 0) {
    100 * sd(qc_tic, na.rm = TRUE) / mean(qc_tic, na.rm = TRUE)
  } else {
    NA_real_
  }

  type_counts <- plot_df %>%
    count(Type) %>%
    mutate(label = paste0(Type, " n=", n)) %>%
    pull(label) %>%
    paste(collapse = " | ")

  subtitle_txt <- if (is.finite(qc_cv)) {
    paste0(type_counts, sprintf(" | QC TIC CV=%.1f%%", qc_cv))
  } else {
    type_counts
  }

  p <- ggplot(plot_df, aes(x = Run, y = TIC, color = Type)) +
    geom_line(data = plot_df %>% filter(Type == "Sample"),
              aes(group = 1), alpha = 0.35, linewidth = 0.35) +
    geom_point(size = 1.9, alpha = 0.85) +
    ggrepel::geom_text_repel(
      data = plot_df %>% filter(!is.na(QC_Label)),
      aes(label = QC_Label),
      color = "black",
      size = 3,
      seed = 42,
      box.padding = 0.25,
      point.padding = 0.2,
      min.segment.length = 0
    ) +
    scale_color_manual(
      values = c(
        "Sample" = "#440154FF",
        "Check" = "#FDE725FF",
        "QC" = "#E64B35FF",
        "Blank" = "#7A7A7AFF",
        "ISTD" = "#00A087FF"
      ),
      breaks = c("Sample", "Check", "QC", "Blank", "ISTD"),
      name = "Injection Type"
    ) +
    scale_y_continuous(labels = scales::scientific) +
    labs(
      title = paste0(condition_name, " - Run Order TIC from Raw File"),
      subtitle = subtitle_txt,
      x = "Run Order",
      y = "Total Ion Current (sum of all features)"
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

  p
}

# Create run order plots for both conditions
suppfig_s1a <- create_run_order_plot(
  c(
    "data/SetA_lipid_FLO2019Control.csv",
    "data/raw_lipid_intensities/SetA_lipid_FLO2019Control.csv",
    "data/raw_lipid_intensities/A_final_lipids.csv",
    "data/raw_lipid_intensities/A_cleaned_lipids_unsummed.csv",
    "data/raw_lipid_intensities/A_cleaned_lipids.csv"
  ),
  "A) Control"
)

suppfig_s1b <- create_run_order_plot(
  c(
    "data/SetB_lipid_FLO2022_lowP.csv",
    "data/raw_lipid_intensities/SetB_lipid_FLO2022_lowP.csv",
    "data/raw_lipid_intensities/B_final_lipids.csv",
    "data/raw_lipid_intensities/B_cleaned_lipids_unsummed.csv",
    "data/raw_lipid_intensities/B_cleaned_lipids.csv"
  ),
  "B) LowInput"
)

message("Created: S1A & S1B (Run Order Stability)")

# ─────────────────────────────────────────────────────────────────────────────
# S1C & S1D: SERRF QC-RSD Improvement
# ─────────────────────────────────────────────────────────────────────────────

create_qc_rsd_bar_plot <- function(rsd_path, panel_title) {
  if (!file.exists(rsd_path)) {
    message("  RSD file not found: ", rsd_path)
    return(ggplot() + theme_void() +
             ggtitle(paste0(panel_title, ": RSD file missing")))
  }

  rsd <- suppressMessages(vroom(rsd_path, show_col_types = FALSE))
  req_cols <- c("QC_none", "QC_SERRF")
  if (!all(req_cols %in% names(rsd))) {
    message("  Required RSD columns missing in: ", rsd_path)
    return(ggplot() + theme_void() +
             ggtitle(paste0(panel_title, ": Invalid RSD format")))
  }

  raw_rsd <- as.numeric(rsd$QC_none) * 100
  serrf_rsd <- as.numeric(rsd$QC_SERRF) * 100

  plot_df <- tibble(
    Stage = factor(c("Before SERRF", "After SERRF"), levels = c("Before SERRF", "After SERRF")),
    Median_RSD = c(median(raw_rsd, na.rm = TRUE), median(serrf_rsd, na.rm = TRUE)),
    Pct_lt_30 = c(mean(raw_rsd < 30, na.rm = TRUE) * 100, mean(serrf_rsd < 30, na.rm = TRUE) * 100)
  )

  improve_pct <- if (is.finite(plot_df$Median_RSD[1]) && plot_df$Median_RSD[1] > 0) {
    100 * (plot_df$Median_RSD[1] - plot_df$Median_RSD[2]) / plot_df$Median_RSD[1]
  } else {
    NA_real_
  }

  p <- ggplot(plot_df, aes(x = Stage, y = Median_RSD, fill = Stage)) +
    geom_col(width = 0.6, color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%", Median_RSD)),
              vjust = -0.35, size = 3.8, fontface = "bold") +
    scale_fill_manual(values = c("Before SERRF" = "#9E9E9E", "After SERRF" = "#009E73")) +
    labs(
      title = panel_title,
      subtitle = paste0(
        "QC-RSD < 30%: ",
        sprintf("%.1f%%", plot_df$Pct_lt_30[1]), " → ",
        sprintf("%.1f%%", plot_df$Pct_lt_30[2]),
        ifelse(is.finite(improve_pct), paste0(" | Median drop: ", sprintf("%.1f%%", improve_pct)), "")
      ),
      x = NULL,
      y = "Median QC-RSD (%)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9.5, color = "grey35"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5)
    )

  p
}

suppfig_s1c <- create_qc_rsd_bar_plot(
  "data/SEERF_normalized_intensities/A_QC-RSDs.csv",
  "C) Control: QC-RSD Improvement After SERRF"
)

suppfig_s1d <- create_qc_rsd_bar_plot(
  "data/SEERF_normalized_intensities/B_QC-RSDs.csv",
  "D) LowInput: QC-RSD Improvement After SERRF"
)

message("Created: S1C & S1D (SERRF QC-RSD Improvement)")

# ─────────────────────────────────────────────────────────────────────────────
# S1E & S1F: PCA Comparison (Raw vs SERRF)
# ─────────────────────────────────────────────────────────────────────────────

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

  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Stage, shape = Type)) +
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

  p
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

# ─────────────────────────────────────────────────────────────────────────────
# S1G & S1H: SpATS Residual Spatial Maps for One Lipid
# ─────────────────────────────────────────────────────────────────────────────

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

  pheno <- pheno_num %>%
    as_tibble(rownames = "LineRaw")

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

  pheno <- pheno %>%
    mutate(Key = canonical_id(LineRaw))

  map_long <- map_long %>%
    mutate(Key = canonical_id(MapRaw))

  filtered <- map_long %>%
    inner_join(pheno, by = "Key") %>%
    dplyr::select(LineRaw = LineRaw, row, col, everything())

  filtered
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

  fit_dat <- fit_dat %>%
    mutate(LineRaw = as.factor(LineRaw))

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

  resid_df <- fit_dat %>%
    mutate(Residual = as.numeric(residuals(spats_fit)))

  p <- ggplot(resid_df, aes(x = col, y = row, fill = Residual)) +
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

  p
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

# ─────────────────────────────────────────────────────────────────────────────
# Combine S1 panels
# ─────────────────────────────────────────────────────────────────────────────

suppfig_s1_qc <- (suppfig_s1a + suppfig_s1b) /
  (suppfig_s1c + suppfig_s1d) /
  (suppfig_s1e + suppfig_s1f) /
  (suppfig_s1g + suppfig_s1h) +
  plot_layout(heights = c(1.05, 0.95, 1, 1))

save_fig(suppfig_s1_qc, "SuppFig_S1_QC_RunOrder_SERRF_PCA_SpATS.png", w = 16, h = 20, subdir = "supp")
save_fig(suppfig_s1_qc, "SuppFig_S1_QC_RunOrder_SERRF_PCA_SpATS.pdf", w = 16, h = 20, subdir = "supp")
# Backward-compatible filename used in previous drafts
save_fig(suppfig_s1_qc, "SuppFig_S1_QC_RunOrder_SERRF_PCA.png", w = 16, h = 20, subdir = "supp")
save_fig(suppfig_s1_qc, "SuppFig_S1_QC_RunOrder_SERRF_PCA.pdf", w = 16, h = 20, subdir = "supp")
# Backward-compatible filename used in previous drafts
save_fig(suppfig_s1_qc, "SuppFig_S1_QC_RunOrder_SpATS.png", w = 16, h = 20, subdir = "supp")
save_fig(suppfig_s1_qc, "SuppFig_S1_QC_RunOrder_SpATS.pdf", w = 16, h = 20, subdir = "supp")

message("✓ Supplementary Figure S1 (QC) complete!\n")


# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
#
#       SUPPLEMENTARY FIGURE S3: LIPID SPECIES COUNTS & OVERVIEW
#
# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════

message("\n══════════════════════════════════════════════════════════════")
message("SUPPLEMENTARY FIGURE S3: LIPID SPECIES COUNTS & OVERVIEW")
message("══════════════════════════════════════════════════════════════\n")

# Load lipid class info
lipid_class_info <- vroom("data/lipid_class/final_lipid_classes.csv",
                          show_col_types = FALSE) %>%
  transmute(Lipid = Lipids, SuperClass = Class)

message("Loaded lipid class info: ", nrow(lipid_class_info), " entries")
message("SuperClasses: ", paste(unique(lipid_class_info$SuperClass), collapse = ", "))

# Get lipid columns from each condition
get_lipid_names <- function(df) {

  cols <- names(df)
  # Keep columns with parentheses (lipid species notation)
  lipid_cols <- cols[grepl("\\(", cols)]
  lipid_cols
}

ctrl_lipids <- get_lipid_names(control_raw)
low_lipids <- get_lipid_names(low_raw)

# Count totals and overlaps
n_ctrl <- length(ctrl_lipids)
n_low <- length(low_lipids)
common_lipids <- intersect(ctrl_lipids, low_lipids)
n_common <- length(common_lipids)
ctrl_only <- setdiff(ctrl_lipids, low_lipids)
low_only <- setdiff(low_lipids, ctrl_lipids)

message("Control lipids: ", n_ctrl)
message("LowInput lipids: ", n_low)
message("Common: ", n_common)
message("Control-only: ", length(ctrl_only))
message("LowInput-only: ", length(low_only))

# ─────────────────────────────────────────────────────────────────────────────
# S2A: Total Species Count (Venn-style bar)
# ─────────────────────────────────────────────────────────────────────────────
total_counts <- tibble(
  Category = c("Control Total", "LowInput Total", "Common", "Control Only", "LowInput Only"),
  Count = c(n_ctrl, n_low, n_common, length(ctrl_only), length(low_only)),
  Type = c("Total", "Total", "Overlap", "Unique", "Unique")
)

total_counts$Category <- factor(total_counts$Category,
                                 levels = c("Control Total", "LowInput Total",
                                           "Common", "Control Only", "LowInput Only"))

suppfig_s2a_counts <- ggplot(total_counts, aes(x = Category, y = Count, fill = Type)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = Count), vjust = -0.3, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("Total" = "#440154FF", "Overlap" = "#21908CFF", "Unique" = "#FDE725FF")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "A) Lipid Species Counts",
    subtitle = sprintf("Total unique species: %d", length(unique(c(ctrl_lipids, low_lipids)))),
    x = NULL, y = "Number of Species"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0),
    plot.subtitle = element_text(color = "grey40"),
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5)
  )

# ─────────────────────────────────────────────────────────────────────────────
# S2B: Species Count by Lipid Class (using get_lipid_class)
# ─────────────────────────────────────────────────────────────────────────────

# Count by class for each condition
count_by_class <- function(lipids, condition_name) {
  tibble(Lipid = lipids) %>%
    mutate(Class = get_lipid_class(Lipid)) %>%
    count(Class, name = "Count") %>%
    mutate(Condition = condition_name)
}

class_ctrl <- count_by_class(ctrl_lipids, "Control")
class_low <- count_by_class(low_lipids, "LowInput")
class_counts <- bind_rows(class_ctrl, class_low) %>%
  filter(!is.na(Class), Class != "Other")

# Order classes by total count
class_order <- class_counts %>%
  group_by(Class) %>%
  summarise(Total = sum(Count)) %>%
  arrange(desc(Total)) %>%
  pull(Class)

class_counts$Class <- factor(class_counts$Class, levels = rev(class_order))

suppfig_s2b_counts <- ggplot(class_counts, aes(x = Class, y = Count, fill = Condition)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7,
           color = "black", linewidth = 0.3) +
  geom_text(aes(label = Count),
            position = position_dodge(width = 0.8),
            hjust = -0.2, size = 3) +
  coord_flip() +
  scale_fill_manual(values = condition_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "B) Species Count by Lipid Class",
    x = NULL, y = "Number of Species"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0),
    legend.position = "top",
    panel.grid.major.y = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5)
  )

# ─────────────────────────────────────────────────────────────────────────────
# S2C: Common vs Unique by Class
# ─────────────────────────────────────────────────────────────────────────────

# For each class, count common and unique species
overlap_by_class <- function(ctrl_lipids, low_lipids) {
  all_lipids <- unique(c(ctrl_lipids, low_lipids))

  tibble(Lipid = all_lipids) %>%
    mutate(
      Class = get_lipid_class(Lipid),
      InControl = Lipid %in% ctrl_lipids,
      InLowInput = Lipid %in% low_lipids,
      Status = case_when(
        InControl & InLowInput ~ "Common",
        InControl & !InLowInput ~ "Control Only",
        !InControl & InLowInput ~ "LowInput Only"
      )
    ) %>%
    filter(!is.na(Class), Class != "Other") %>%
    count(Class, Status) %>%
    mutate(Status = factor(Status, levels = c("Common", "Control Only", "LowInput Only")))
}

overlap_class <- overlap_by_class(ctrl_lipids, low_lipids)
overlap_class$Class <- factor(overlap_class$Class, levels = rev(class_order))

suppfig_s2c_counts <- ggplot(overlap_class, aes(x = Class, y = n, fill = Status)) +
  geom_col(position = "stack", width = 0.7, color = "black", linewidth = 0.2) +
  coord_flip() +
  scale_fill_manual(values = c("Common" = "#21908CFF",
                                "Control Only" = "#440154FF",
                                "LowInput Only" = "#FDE725FF")) +
  labs(
    title = "C) Common vs Unique Species by Class",
    x = NULL, y = "Number of Species"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0),
    legend.position = "top",
    panel.grid.major.y = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5)
  )

# ─────────────────────────────────────────────────────────────────────────────
# S2D: Species Count by SuperClass (from lipid_class_info)
# ─────────────────────────────────────────────────────────────────────────────

# Match lipids to SuperClass
count_by_superclass <- function(lipids, condition_name, class_info) {
  tibble(Lipid = lipids) %>%
    left_join(class_info, by = "Lipid") %>%
    mutate(SuperClass = ifelse(is.na(SuperClass), "Unclassified", SuperClass)) %>%
    count(SuperClass, name = "Count") %>%
    mutate(Condition = condition_name)
}

super_ctrl <- count_by_superclass(ctrl_lipids, "Control", lipid_class_info)
super_low <- count_by_superclass(low_lipids, "LowInput", lipid_class_info)
super_counts <- bind_rows(super_ctrl, super_low)

# Order by total count
super_order <- super_counts %>%
  group_by(SuperClass) %>%
  summarise(Total = sum(Count)) %>%
  arrange(desc(Total)) %>%
  pull(SuperClass)

super_counts$SuperClass <- factor(super_counts$SuperClass, levels = rev(super_order))

# Define SuperClass colors
superclass_colors <- c(
  "Glycerolipid" = "#ED804A",
  "Glycerophospholipid" = "#00441B",
  "Sterol" = "#7570B3",
  "Terpenoid" = "#E7298A",
  "Fatty acid" = "#66A61E",
  "Sphingolipid" = "#D95F02",
  "Steroid" = "#A6761D",
  "Prenol" = "#1B9E77",
  "Ether lipid" = "#E6AB02",
  "Tetrapyrrole" = "#666666",
  "Flavonoid" = "#F768A1",
  "Betaine lipid" = "#8941ED",
  "Unclassified" = "grey70"
)

suppfig_s2d_counts <- ggplot(super_counts, aes(x = SuperClass, y = Count, fill = Condition)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7,
           color = "black", linewidth = 0.3) +
  geom_text(aes(label = Count),
            position = position_dodge(width = 0.8),
            hjust = -0.2, size = 3) +
  coord_flip() +
  scale_fill_manual(values = condition_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "D) Species Count by SuperClass",
    subtitle = "Classification from lipid database",
    x = NULL, y = "Number of Species"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0),
    plot.subtitle = element_text(color = "grey40"),
    legend.position = "top",
    panel.grid.major.y = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5)
  )

# ─────────────────────────────────────────────────────────────────────────────
# Combine S2 panels
# ─────────────────────────────────────────────────────────────────────────────
suppfig_s3_counts_top <- suppfig_s2a_counts + suppfig_s2b_counts + plot_layout(widths = c(1, 1.3))
suppfig_s3_counts_bottom <- suppfig_s2c_counts + suppfig_s2d_counts + plot_layout(widths = c(1, 1.3))

suppfig_s3_counts <- suppfig_s3_counts_top / suppfig_s3_counts_bottom +
  plot_layout(heights = c(1, 1))

save_fig(suppfig_s3_counts, "SuppFig_S3_Lipid_Species_Counts.png", w = 14, h = 12, subdir = "supp")
save_fig(suppfig_s3_counts, "SuppFig_S3_Lipid_Species_Counts.pdf", w = 14, h = 12, subdir = "supp")

# Save summary tables
species_summary <- tibble(
  Metric = c("Control Total", "LowInput Total", "Common", "Control Only", "LowInput Only", "All Unique"),
  Count = c(n_ctrl, n_low, n_common, length(ctrl_only), length(low_only),
            length(unique(c(ctrl_lipids, low_lipids))))
)
write.csv(species_summary, "table/supp/SuppTable_S3a_Species_Summary.csv", row.names = FALSE)

class_summary <- class_counts %>%
  pivot_wider(names_from = Condition, values_from = Count, values_fill = 0) %>%
  arrange(desc(Control + LowInput))
write.csv(class_summary, "table/supp/SuppTable_S3b_Species_by_Class.csv", row.names = FALSE)

super_summary <- super_counts %>%
  pivot_wider(names_from = Condition, values_from = Count, values_fill = 0) %>%
  arrange(desc(Control + LowInput))
write.csv(super_summary, "table/supp/SuppTable_S3c_Species_by_SuperClass.csv", row.names = FALSE)

message("✓ Supplementary Figure S3 (Lipid Species Counts) complete!\n")


# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
#
#       SUPPLEMENTARY FIGURE S4: PCA OF LIPIDS (Control & LowInput)
#
# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════

message("\n══════════════════════════════════════════════════════════════")
message("SUPPLEMENTARY FIGURE S4: PCA OF LIPIDS")
message("══════════════════════════════════════════════════════════════\n")

# Run PCA (lipids as points, samples as features) using CLR
run_lipid_pca <- function(df_tic, pseudocount = 1e-6) {
  X <- df_tic %>%
    dplyr::select(-Compound_Name, -Condition) %>%
    as.matrix()

  # CLR transformation
  P <- X / pmax(rowSums(X, na.rm = TRUE), 1e-12)
  P <- P + pseudocount
  P <- P / rowSums(P, na.rm = TRUE)

  L <- log(P)
  X_clr <- L - rowMeans(L, na.rm = TRUE)

  # Transpose: lipids as rows, samples as columns
  V <- t(X_clr)
  V <- t(scale(t(V), center = TRUE, scale = TRUE))
  V[is.na(V)] <- 0

  pc <- prcomp(V, center = FALSE, scale. = FALSE)

  out <- as.data.frame(pc$x[, 1:2])
  colnames(out) <- c("PC1", "PC2")
  out$Lipid <- rownames(pc$x)
  out$Class <- get_lipid_class(out$Lipid)
  out$Group <- class_to_group(out$Class)

  ve <- (pc$sdev^2) / sum(pc$sdev^2)
  attr(out, "pve") <- ve[1:2]
  out
}

pca_control <- run_lipid_pca(control_tic)
pca_low <- run_lipid_pca(low_tic)

# Plot function for PCA
plot_lipid_pca <- function(pca_df, title = "", show_legend = TRUE) {
  pve <- attr(pca_df, "pve")
  xlab <- sprintf("PC1 (%.1f%%)", 100 * pve[1])
  ylab <- sprintf("PC2 (%.1f%%)", 100 * pve[2])

  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Class)) +
    geom_point(size = 2.5, alpha = 0.8) +
    stat_ellipse(aes(group = Group, color = NULL),
                 linetype = "dashed", linewidth = 0.8, color = "grey50",
                 level = 0.95, type = "norm") +
    scale_color_manual(values = class_colors) +
    labs(title = title, x = xlab, y = ylab) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = if (show_legend) "right" else "none",
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5)
    )
  p
}

suppfig_s3_control <- plot_lipid_pca(pca_control, "A) Control: Lipid PCA", show_legend = FALSE)
suppfig_s3_low <- plot_lipid_pca(pca_low, "B) LowInput: Lipid PCA", show_legend = TRUE)

suppfig_s4_pca <- suppfig_s3_control + suppfig_s3_low +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

save_fig(suppfig_s4_pca, "SuppFig_S4_PCA_Lipids.png", w = 14, h = 7, subdir = "supp")
save_fig(suppfig_s4_pca, "SuppFig_S4_PCA_Lipids.pdf", w = 14, h = 7, subdir = "supp")

message("✓ Supplementary Figure S4 (PCA) complete!\n")


# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
#
#       SUPPLEMENTARY FIGURE S5: TOP VARIANCE LIPID SPECIES
#
# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════

message("\n══════════════════════════════════════════════════════════════")
message("SUPPLEMENTARY FIGURE S5: TOP VARIANCE LIPID SPECIES (GENOTYPE VARIATION)")
message("══════════════════════════════════════════════════════════════\n")

# Function to find high variance lipids across genotypes
find_high_variance_lipids <- function(df_tic, condition_name, top_n = 10) {
  lipid_cols <- setdiff(names(df_tic), c("Compound_Name", "Condition"))
  lipid_cols <- lipid_cols[grepl("\\(", lipid_cols)]

  X <- df_tic %>%
    dplyr::select(all_of(lipid_cols)) %>%
    mutate(across(everything(), ~suppressWarnings(as.numeric(.x))))

  X %>%
    summarise(across(everything(), ~var(.x, na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "Lipid", values_to = "Variance") %>%
    arrange(desc(Variance)) %>%
    mutate(
      Rank = row_number(),
      Condition = condition_name,
      Class = get_lipid_class(Lipid)
    ) %>%
    filter(Rank <= top_n)
}

# Find top 10 high variance lipids for each condition
top_var_control <- find_high_variance_lipids(control_tic, "Control", top_n = 10)
top_var_low <- find_high_variance_lipids(low_tic, "LowInput", top_n = 10)

message("Top variance lipids identified:")
message("  Control: ", paste(head(top_var_control$Lipid, 3), collapse = ", "), ", ...")
message("  LowInput: ", paste(head(top_var_low$Lipid, 3), collapse = ", "), ", ...")

# ─────────────────────────────────────────────────────────────────────────────
# S3A: Control - Top 10 High Variance Lipids
# ─────────────────────────────────────────────────────────────────────────────
suppfig_s3a_var <- ggplot(top_var_control,
                          aes(x = reorder(Lipid, Variance), y = Variance, fill = Class)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  coord_flip() +
  scale_fill_manual(values = class_colors, na.value = "grey70") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = NULL,
    y = "Variance"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(color = "black", size = 9),
    axis.text.x = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0),
    plot.subtitle = element_text(hjust = 0, color = "grey40"),
    legend.position = "right",
    panel.grid.major.y = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5)
  )

# ─────────────────────────────────────────────────────────────────────────────
# S3B: LowInput - Top 10 High Variance Lipids
# ─────────────────────────────────────────────────────────────────────────────
suppfig_s3b_var <- ggplot(top_var_low,
                          aes(x = reorder(Lipid, Variance), y = Variance, fill = Class)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  coord_flip() +
  scale_fill_manual(values = class_colors, na.value = "grey70") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = NULL,
    y = "Variance"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(color = "black", size = 9),
    axis.text.x = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0),
    plot.subtitle = element_text(hjust = 0, color = "grey40"),
    legend.position = "right",
    panel.grid.major.y = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5)
  )

# ─────────────────────────────────────────────────────────────────────────────
# S3C: Overlap comparison (Venn-style summary)
# ─────────────────────────────────────────────────────────────────────────────
# Find which lipids are in both top 10 lists
overlap_lipids <- intersect(top_var_control$Lipid, top_var_low$Lipid)
ctrl_only <- setdiff(top_var_control$Lipid, top_var_low$Lipid)
low_only <- setdiff(top_var_low$Lipid, top_var_control$Lipid)

# Combine for comparison plot
combined_var <- bind_rows(
  top_var_control %>% mutate(Condition = "Control"),
  top_var_low %>% mutate(Condition = "LowInput")
)

# Create comparison scatter
var_wide <- full_join(
  top_var_control %>% dplyr::select(Lipid, Variance, Class) %>% dplyr::rename(Var_Ctrl = Variance),
  top_var_low %>% dplyr::select(Lipid, Variance) %>% dplyr::rename(Var_Low = Variance),
  by = "Lipid"
) %>%
  mutate(
    Class = ifelse(is.na(Class), get_lipid_class(Lipid), Class),
    InBoth = Lipid %in% overlap_lipids
  )

suppfig_s3c_var <- ggplot(var_wide, aes(x = Var_Ctrl, y = Var_Low)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = Class, shape = InBoth), size = 3, alpha = 0.8) +
  geom_text_repel(aes(label = Lipid), size = 2.8, max.overlaps = 15) +
  scale_color_manual(values = class_colors, na.value = "grey70") +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 17),
                     labels = c("TRUE" = "In both top 10", "FALSE" = "Unique"),
                     name = "Overlap") +
  labs(
    x = "Variance (Control)",
    y = "Variance (LowInput)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0),
    plot.subtitle = element_text(hjust = 0, color = "grey40"),
    legend.position = "right",
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5)
  )

# ─────────────────────────────────────────────────────────────────────────────
# S3D: Class distribution of high variance lipids
# ─────────────────────────────────────────────────────────────────────────────
class_counts <- combined_var %>%
  count(Condition, Class) %>%
  mutate(Condition = factor(Condition, levels = c("Control", "LowInput")))

suppfig_s3d_var <- ggplot(class_counts, aes(x = Class, y = n, fill = Condition)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = condition_colors) +
  labs(
    x = "Lipid Class",
    y = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0),
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5)
  )

# ─────────────────────────────────────────────────────────────────────────────
# Combine S3 panels
# ─────────────────────────────────────────────────────────────────────────────
suppfig_s5_var_top <- suppfig_s3a_var + suppfig_s3b_var + plot_layout(guides = "collect") &
  theme(legend.position = "right")
suppfig_s5_var_bottom <- suppfig_s3c_var + suppfig_s3d_var + plot_layout(widths = c(1.2, 1))

suppfig_s5_var <- suppfig_s5_var_top / suppfig_s5_var_bottom +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold", size = 16)))

save_fig(suppfig_s5_var, "SuppFig_TopVariance_Lipids.png", w = 14, h = 12, subdir = "supp")

# Save table of top variance lipids
top_var_combined <- bind_rows(
  top_var_control %>% mutate(Condition = "Control"),
  top_var_low %>% mutate(Condition = "LowInput")
) %>%
  dplyr::select(Condition, Rank, Lipid, Class, Variance) %>%
  arrange(Condition, Rank)

write.csv(top_var_combined, "table/supp/SuppTable_S4_TopVariance_Lipids.csv", row.names = FALSE)

message("✓ Supplementary Figure S2 (Top Variance Lipids) complete!\n")

# TEMPORARY STOP - Remove this when ready for next sections
if (identical(Sys.getenv("STOP_AFTER_SECTION2", ""), "1")) {
  message("STOP_AFTER_SECTION2=1: Stopping after Section 2.")
  quit(save = "no", status = 0)
}

# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
#
#                       FIGURE 3: OPLS-DA ANALYSIS
#
# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════

message("\n══════════════════════════════════════════════════════════════")
message("FIGURE 3: OPLS-DA ANALYSIS")
message("══════════════════════════════════════════════════════════════\n")

message("Running OPLS-DA...")

# Prepare data for OPLS-DA (using class ratios like in script 14)
long_prep_global <- function(df) {
  class_pat <- paste0("\\b(", paste(classes_with_lyso, collapse = "|"), ")\\b")

  df %>%
    pivot_longer(-c(Compound_Name, Condition), names_to = "Lipid", values_to = "Intensity") %>%
    dplyr::rename(Sample = Compound_Name) %>%
    group_by(Sample) %>%
    mutate(TIC = sum(Intensity, na.rm = TRUE),
           rel_abund = Intensity / TIC) %>%
    ungroup() %>%
    group_by(Sample) %>%
    mutate(minpos = ifelse(all(is.na(rel_abund) | rel_abund <= 0), NA_real_,
                           min(rel_abund[rel_abund > 0], na.rm = TRUE)),
           eps = ifelse(is.na(minpos), 0, minpos * 0.5)) %>%
    ungroup() %>%
    mutate(log_rel = log10(rel_abund + eps)) %>%
    mutate(Class = str_extract(Lipid, class_pat)) %>%
    filter(!is.na(Class)) %>%
    group_by(Sample, Condition, Class) %>%
    summarise(class_log = mean(log_rel, na.rm = TRUE), .groups = "drop")
}

# Combine and process
combined <- bind_rows(
  control_raw %>% dplyr::mutate(Condition = "Control"),
  low_raw %>% dplyr::mutate(Condition = "LowInput")
)

long_all <- long_prep_global(combined) %>%
  filter(!is.na(class_log), nchar(Sample) >= 5)

# Pivot to wide and create ratios
wide_log <- long_all %>%
  pivot_wider(names_from = Class, values_from = class_log, values_fill = NA_real_)

# Create all pairwise ratios (log-ratio = subtraction in log space)
long_classes <- wide_log %>%
  pivot_longer(cols = -c(Sample, Condition), names_to = "Class", values_to = "Abundance")

ratios_long <- long_classes %>%
  inner_join(long_classes, by = c("Sample", "Condition"), suffix = c(".num", ".den")) %>%
  filter(Class.num < Class.den) %>%
  transmute(Sample, Condition,
            RatioName = paste(Class.num, Class.den, sep = "/"),
            Ratio = Abundance.num - Abundance.den)  # Subtraction = log-ratio

# Filter complete cases
full_n <- ratios_long %>% group_by(Sample) %>% tally() %>% pull(n) %>% max()
ratios_clean <- ratios_long %>% group_by(Sample) %>% filter(n() == full_n) %>% ungroup()

# Create matrix for OPLS-DA
wide_ratios <- ratios_clean %>%
  pivot_wider(id_cols = c(Sample, Condition), names_from = RatioName, values_from = Ratio)

meta <- wide_ratios %>% dplyr::select(Sample, Condition)
lipid_mat <- wide_ratios %>% dplyr::select(-Sample, -Condition) %>% as.matrix()
class_vec <- factor(meta$Condition)

# Fit OPLS-DA
set.seed(123)
opls_mod <- opls(lipid_mat, class_vec, predI = 1, orthoI = NA, scaleC = "standard",
                 fig.pdfC = "none", info.txtC = "none")

# Extract scores
scores_df <- data.frame(
  Sample = rownames(opls_mod@scoreMN),
  t1 = opls_mod@scoreMN[, 1],
  to1 = opls_mod@orthoScoreMN[, 1],
  Condition = class_vec
)

# ─────────────────────────────────────────────────────────────────────────────
# S4A: OPLS-DA scores plot
# ─────────────────────────────────────────────────────────────────────────────
suppfig_s5a <- ggplot(scores_df, aes(x = t1, y = to1, color = Condition)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 2.5, alpha = 0.7) +
  stat_ellipse(aes(fill = Condition), geom = "polygon", type = "norm",
               level = 0.95, alpha = 0.15, colour = NA) +
  scale_color_manual(values = condition_colors) +
  scale_fill_manual(values = condition_colors) +
  labs(
    x = paste0("t1 (", round(opls_mod@summaryDF$`R2X(cum)` * 100, 1), "% of X variance)"),
    y = "to1 (orthogonal)",
    title = "A) OPLS-DA Scores Plot"
  ) +
  coord_fixed() +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "grey70", linewidth = 0.3),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5)
  )

# ─────────────────────────────────────────────────────────────────────────────
# S4B: VIP scores
# ─────────────────────────────────────────────────────────────────────────────
vip_vec <- slot(opls_mod, "vipVn")
vip_df <- tibble(Ratio = names(vip_vec), VIP = as.numeric(vip_vec)) %>%
  filter(VIP > 1) %>%
  arrange(desc(VIP)) %>%
  slice_head(n = 25)

suppfig_s5b <- ggplot(vip_df, aes(x = reorder(Ratio, VIP), y = VIP)) +
  geom_col(fill = "#440154FF", width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(x = NULL, y = "VIP Score", title = "B) Top Discriminatory Class Ratios (VIP > 1)") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0),
    panel.grid.major.y = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5)
  )

# ─────────────────────────────────────────────────────────────────────────────
# S4C: Permutation test
# ─────────────────────────────────────────────────────────────────────────────
message("Running OPLS-DA permutation test (n=200)...")
perm_res <- opls(lipid_mat, class_vec, predI = 1, permI = 200, scaleC = "standard",
                 fig.pdfC = "none", info.txtC = "none")

perm_mat <- as.data.frame(perm_res@suppLs$permMN)
colnames(perm_mat) <- c("R2Xcum", "R2Ycum", "Q2cum", "RMSEE", "pre", "ort", "unused")
perm_only <- perm_mat[-1, ]

obsR2Y <- opls_mod@summaryDF$`R2Y(cum)`
obsQ2 <- opls_mod@summaryDF$`Q2(cum)`

pR2Y <- (sum(perm_only$R2Ycum >= obsR2Y) + 1) / (nrow(perm_only) + 1)
pQ2 <- (sum(perm_only$Q2cum >= obsQ2) + 1) / (nrow(perm_only) + 1)

perm_plot_df <- perm_only %>%
  dplyr::select(R2Ycum, Q2cum) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric, R2Ycum = "R²Y", Q2cum = "Q²"))

# Use ggbreak for zigzag axis break
if (!requireNamespace("ggbreak", quietly = TRUE)) {
  install.packages("ggbreak")
}
library(ggbreak)

# Calculate break range: from max permutation value to min observed value
perm_max_val <- max(perm_plot_df$value, na.rm = TRUE)
obs_min_val <- min(obsR2Y, obsQ2)
break_start <- perm_max_val + 0.05
break_end <- obs_min_val - 0.05

suppfig_s5c <- ggplot(perm_plot_df, aes(x = value, fill = metric)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity", color = "white") +
  geom_vline(xintercept = obsR2Y, linetype = "dashed", color = "#440154FF", linewidth = 1) +
  geom_vline(xintercept = obsQ2, linetype = "dashed", color = "#FDE725FF", linewidth = 1) +
  # Place annotations to the LEFT of the lines so they're visible
  annotate("text", x = obsR2Y - 0.01, y = Inf,
           label = sprintf("R²Y = %.3f\np = %.3f", obsR2Y, pR2Y),
           vjust = 1.5, hjust = 1, size = 3.5, color = "#440154FF", fontface = "bold") +
  annotate("text", x = obsQ2 - 0.01, y = Inf,
           label = sprintf("Q² = %.3f\np = %.3f", obsQ2, pQ2),
           vjust = 3.5, hjust = 1, size = 3.5, color = "#FDE725FF", fontface = "bold") +
  scale_fill_manual(values = c("R²Y" = "#440154FF", "Q²" = "#FDE725FF")) +
  # Add zigzag break to skip empty space between histogram and observed values
  scale_x_break(c(break_start, break_end), scales = 0.5) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  labs(x = "Metric Value", y = "Frequency", title = "C) Permutation Test (n=200)") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5)
  )

# ─────────────────────────────────────────────────────────────────────────────
# S4D: Model summary metrics
# ─────────────────────────────────────────────────────────────────────────────
model_metrics <- tibble(
  Metric = c("R²X (cumulative)", "R²Y (cumulative)", "Q² (cumulative)",
             "Predictive components", "Orthogonal components"),
  Value = c(
    round(opls_mod@summaryDF$`R2X(cum)`, 3),
    round(opls_mod@summaryDF$`R2Y(cum)`, 3),
    round(opls_mod@summaryDF$`Q2(cum)`, 3),
    opls_mod@summaryDF$pre,
    opls_mod@summaryDF$ort
  )
)

suppfig_s5d <- ggplot(model_metrics, aes(x = Metric, y = as.numeric(Value))) +
  geom_col(fill = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E"), width = 0.6) +
  geom_text(aes(label = Value), vjust = -0.5, size = 4, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Value", title = "D) Model Quality Metrics") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0),
    panel.grid.major.x = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5)
  )

# # ─────────────────────────────────────────────────────────────────────────────
# # Combine S4 panels
# # ─────────────────────────────────────────────────────────────────────────────
suppfig_s5_top <- suppfig_s5a + suppfig_s5b + plot_layout(widths = c(1, 1.2))

# ggbreak plots need special handling with aplot for combining
if (!requireNamespace("aplot", quietly = TRUE)) {
  install.packages("aplot")
}
library(aplot)

# Use aplot::plot_list for ggbreak compatibility
suppfig_s5_bottom <- aplot::plot_list(suppfig_s5c, suppfig_s5d, widths = c(1.2, 1))

suppfig_s5 <- suppfig_s5_top / suppfig_s5_bottom +
  plot_layout(heights = c(1, 0.8))

save_fig(suppfig_s5, "Figure2_OPLS_DA.png", w = 14, h = 12)
save_fig(suppfig_s5, "Figure2_OPLS_DA.pdf", w = 14, h = 12)

# Save VIP table
write.csv(vip_df, "table/supp/SuppTable_S6_OPLS_VIP_ratios.csv", row.names = FALSE)

message("✓ Figure 3 (OPLS-DA) complete!\n")


# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
#
#          ADDITIONAL ANALYSIS: LIPID-ENZYME NETWORK ANALYSIS
#
# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════

message("\n══════════════════════════════════════════════════════════════")
message("ADDITIONAL ANALYSIS: LIPID-ENZYME NETWORK")
message("══════════════════════════════════════════════════════════════\n")

# Network plotting function (from script 22)
plot_gene_lipid_network <- function(df, lipid_names, gene_text, p_value,
                                    use_tic = FALSE, r_min = 0.0,
                                    out_png = NULL, width = 7.5, height = 7.5) {

  # Standardize lipid names
  standardize_name <- function(x) {
    x %>% str_replace_all("\\s+", "") %>% str_replace_all("_", "/")
  }

  # Find matching columns
  cn <- colnames(df)
  cn_std <- standardize_name(cn)
  target_std <- standardize_name(lipid_names)

  hit_idx <- match(target_std, cn_std)
  if (any(is.na(hit_idx))) {
    missing <- lipid_names[is.na(hit_idx)]
    warning("Missing lipids: ", paste(missing, collapse = ", "))
    hit_idx <- hit_idx[!is.na(hit_idx)]
  }

  if (length(hit_idx) < 2) {
    message("Not enough lipids found for network")
    return(NULL)
  }

  lipid_cols <- cn[hit_idx]

  # Extract and optionally TIC-normalize
  X <- df %>%
    dplyr::select(all_of(lipid_cols)) %>%
    mutate(across(everything(), ~suppressWarnings(as.numeric(.x)))) %>%
    as.matrix()

  if (use_tic) {
    rs <- rowSums(X, na.rm = TRUE)
    rs[rs == 0 | is.na(rs)] <- NA_real_
    X <- sweep(X, 1, rs, FUN = "/") * 100
    X[is.na(X)] <- 0
  }

  R <- cor(X, use = "pairwise.complete.obs", method = "pearson")

  # Build nodes
  nodes <- data.frame(
    name = colnames(R),
    class = str_extract(colnames(R), "^[A-Za-z]+"),
    stringsAsFactors = FALSE
  )

  n <- nrow(nodes)
  theta <- seq(pi/2, pi/2 - 2*pi, length.out = n + 1)[1:n]
  nodes$x <- cos(theta) * 0.55 + 2.25
  nodes$y <- sin(theta) * 0.55
  nodes$fill <- unname(class_colors[nodes$class])
  nodes$fill[is.na(nodes$fill)] <- "grey70"

  # Build edges
  comb <- t(combn(colnames(R), 2))
  edges <- data.frame(
    from = comb[, 1], to = comb[, 2],
    r = mapply(function(a, b) R[a, b], comb[, 1], comb[, 2]),
    stringsAsFactors = FALSE
  ) %>%
    filter(abs(r) >= r_min)

  eplot <- edges %>%
    left_join(nodes %>% dplyr::select(name, x, y), by = c("from" = "name")) %>%
    left_join(nodes %>% dplyr::select(name, x, y), by = c("to" = "name"), suffix = c("_from", "_to")) %>%
    mutate(x1 = x_from, y1 = y_from, x2 = x_to, y2 = y_to,
           xm = (x1 + x2)/2, ym = (y1 + y2)/2,
           lab = sprintf("%.2f", r),
           w = scales::rescale(abs(r), to = c(0.4, 1.4)))

  p_txt <- paste0("P = ", format(p_value, scientific = TRUE))

  p <- ggplot() +
    annotate("segment", x = 0.70, y = 0, xend = 1.65, yend = 0,
             linewidth = 0.9, arrow = arrow(length = unit(0.2, "inches"), type = "closed")) +
    annotate("text", x = 1.18, y = 0.23, label = p_txt, color = "red3",
             fontface = "italic", size = 4) +
    annotate("text", x = -0.15, y = 0, label = gene_text, hjust = 0, size = 4.5) +
    geom_segment(data = eplot, aes(x = x1, y = y1, xend = x2, yend = y2, linewidth = w),
                 color = "grey30", lineend = "round") +
    geom_text(data = eplot, aes(x = xm, y = ym, label = lab), color = "red3", size = 3.5) +
    geom_point(data = nodes, aes(x = x, y = y), shape = 21, size = 8,
               stroke = 1, fill = nodes$fill, color = "grey20") +
    geom_text(data = nodes, aes(x = x, y = y, label = name), vjust = -1.3, size = 3.5) +
    scale_linewidth_identity() +
    coord_equal(xlim = c(-0.2, 3.2), ylim = c(-1.1, 1.1), expand = FALSE) +
    theme_void(base_size = 12)

  if (!is.null(out_png)) {
    ggsave(out_png, p, width = width, height = height, dpi = 300, bg = "white")
    message("Saved: ", out_png)
  }

  p
}

# Create network panels for key genes
# PHR1 network
lipids_phr1 <- c("PC(16:0/20:3)", "PC(16:0/22:5)", "PC(16:0/22:6)",
                 "PC(18:1/20:1)", "PC(18:1/24:1)", "PC(18:2/20:0)", "PE(16:0/18:1)")

fig4_phr1_ctrl <- plot_gene_lipid_network(
  df = control_raw, lipid_names = lipids_phr1,
  gene_text = "SORBI_3001G384300\n(PHR1)", p_value = 1.5e-9
)

fig4_phr1_low <- plot_gene_lipid_network(
  df = low_raw, lipid_names = lipids_phr1,
  gene_text = "SORBI_3001G384300\n(PHR1)", p_value = 6.59e-11
)

# DGAT1 network
lipids_dgat1 <- c("TG(18:1/18:3/22:0)", "TG(18:1/20:3/22:0)", "TG(18:2/18:2/18:4)",
                  "TG(18:2/20:3/22:0)", "TG(18:3/18:3/18:3)")

fig4_dgat1_low <- plot_gene_lipid_network(
  df = low_raw, lipid_names = lipids_dgat1,
  gene_text = "SORBI_3010G170000\n(DGAT1)", p_value = 1.5e-9
)

# LCAT4 network
lipids_lcat4 <- c("DG(18:2/18:3)", "DGDG(16:0/18:2)", "MG(18:2)")

fig4_lcat4_low <- plot_gene_lipid_network(
  df = low_raw, lipid_names = lipids_lcat4,
  gene_text = "SORBI_3006G214500\n(LCAT4)", p_value = 4.11e-10
)

# Combine network panels
if (!is.null(fig4_phr1_ctrl) && !is.null(fig4_phr1_low)) {
  figure4_top <- fig4_phr1_ctrl + fig4_phr1_low +
    plot_annotation(title = "PHR1: Phospholipid Network",
                    subtitle = "Control (left) vs LowInput (right)")
}

message("\n✓ Figure 4 networks created!\n")


# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
#
#                    FIGURE 4: PREDICTIVE MODELING (RF + SHAP)
#
# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════

message("\n══════════════════════════════════════════════════════════════")
message("FIGURE 4: PREDICTIVE MODELING (RF + SHAP)")
message("══════════════════════════════════════════════════════════════\n")

# ─────────────────────────────────────────────────────────────────────────────
# Configuration: Choose phenotype to predict
# ─────────────────────────────────────────────────────────────────────────────
# Options: "FloweringTime", "PlantHeight", "GrainYield", etc.
PHENOTYPE_NAME <- "PlantHeight"
RUN_RF_FRESH <- TRUE  # Set to FALSE to use pre-computed results if available
RUN_SHAP_STABILITY <- FALSE
SAVE_S9B_STABILITY_FIG <- FALSE  # Keep stability table, do not save standalone S9b figure
SAVE_LEGACY_RF_FIG <- FALSE  # Do not overwrite final Figure 4 with the older 4-panel diagnostic composite
SHAP_STABILITY_N_BOOT <- 100
SHAP_STABILITY_TOP_K <- 20
SHAP_STABILITY_NUM_TREES <- 500

# ─────────────────────────────────────────────────────────────────────────────
# Load and prepare data
# ─────────────────────────────────────────────────────────────────────────────
message("Preparing data for Random Forest...")

# Load lipids (use Control condition - can change to low_raw for LowInput)
lipids_rf <- control_raw %>%
  dplyr::select(-Condition) %>%
  dplyr::rename(Line = Compound_Name) %>%
  filter(nchar(Line) >= 5)

# Log10 transform + median center
pseudocount_rf <- 1
lipid_cols_rf <- setdiff(names(lipids_rf), "Line")

lipids_rf <- lipids_rf %>%
  mutate(across(all_of(lipid_cols_rf), ~ {
    x <- suppressWarnings(as.numeric(.x))
    x <- log10(x + pseudocount_rf)
    x - median(x, na.rm = TRUE)
  }))

# Load phenotype data
pheno_file <- "data/phenotypes/control_field_phenotypes.csv"
if (file.exists(pheno_file)) {
  pheno <- vroom(pheno_file, show_col_types = FALSE) %>%
    dplyr::select(1, 2) %>%  # Column 2 = Plant Height (Our Data)
    dplyr::rename(Line = 1, Phenotype = 2) %>%
    filter(!is.na(Phenotype))

  message("Loaded phenotype: ", PHENOTYPE_NAME)
  message("  Samples with phenotype: ", nrow(pheno))

} else {
  message("Phenotype file not found: ", pheno_file)
  message("Skipping RF/SHAP analysis.")
  RUN_RF_FRESH <- FALSE
}

# ─────────────────────────────────────────────────────────────────────────────
# Run RF + SHAP Pipeline
# ─────────────────────────────────────────────────────────────────────────────
if (RUN_RF_FRESH && exists("pheno")) {

  # Merge lipids and phenotype
  dat_rf <- lipids_rf %>%
    inner_join(pheno, by = "Line")

  message("Samples after merge: ", nrow(dat_rf))

  # Load population structure PCs (if available)
  pc_file <- "table/PCA_SAP.csv"
  PC <- NULL
  if (file.exists(pc_file)) {
    PC <- vroom(pc_file, show_col_types = FALSE) %>%
      dplyr::select(1:4) %>%
      dplyr::rename(Line = 1)

    # Align samples
    common_lines <- intersect(dat_rf$Line, PC$Line)
    dat_rf <- dat_rf %>% filter(Line %in% common_lines) %>% arrange(Line)
    PC <- PC %>% filter(Line %in% common_lines) %>% arrange(Line)
    message("Loaded PCs and aligned to phenotype/lipid samples")
  } else {
    message("No PC file found - proceeding without PC residualization")
  }

  pheno_raw <- dat_rf$Phenotype
  lipid_mat_rf <- dat_rf %>% dplyr::select(-Line, -Phenotype) %>% as.matrix()
  pc_mat <- if (!is.null(PC)) PC %>% dplyr::select(-Line) %>% as.data.frame() else NULL

  # ─────────────────────────────────────────────────────────────────────────────
  # Repeated outer split evaluation (leakage-safe residualization inside split)
  # ─────────────────────────────────────────────────────────────────────────────
  RF_OUTER_REPEATS <- 10
  message("Running repeated outer train/test evaluation (n = ", RF_OUTER_REPEATS, ")...")

  rf_eval <- run_repeated_rf_eval(
    X = lipid_mat_rf,
    y = pheno_raw,
    pc_df = pc_mat,
    train_prop = 0.8,
    n_repeats = RF_OUTER_REPEATS,
    seed_base = 42,
    tune_num_trees = 500,
    fit_num_trees = 1000,
    num_threads = max(1, parallel::detectCores() - 1)
  )

  rmse_test <- rf_eval$summary$RMSE_Mean
  rmse_test_sd <- rf_eval$summary$RMSE_SD
  r2_test <- rf_eval$summary$R2_Mean
  r2_test_sd <- rf_eval$summary$R2_SD

  preds_test <- rf_eval$representative$pred
  test_y <- rf_eval$representative$obs
  train_y <- rf_eval$representative$train_y

  message("  Repeated-split RMSE (mean ± SD): ", round(rmse_test, 4), " ± ", round(rmse_test_sd, 4))
  message("  Repeated-split R² (mean ± SD): ", round(r2_test, 4), " ± ", round(r2_test_sd, 4))

  tune_res <- list(recommended.pars = rf_eval$median_params)
  message("Using median tuned parameters across repeats for final SHAP model:")
  message("  mtry: ", tune_res$recommended.pars$mtry)
  message("  min.node.size: ", tune_res$recommended.pars$min.node.size)
  message("  sample.fraction: ", round(tune_res$recommended.pars$sample.fraction, 4))

  # Residualize all data only for SHAP interpretation model (no test evaluation here)
  final_adj <- residualize_full_by_pc(lipid_mat_rf, pheno_raw, pc_df = pc_mat)
  lipid_mat_adj <- final_adj$X
  pheno_resid <- final_adj$y

  if (!is.null(pc_mat)) {
    message("Applied PC residualization on full data for SHAP model")
  }

  # Regr task on full data for supplementary CV panel
  regr_task <- makeRegrTask(
    data = data.frame(as.data.frame(lipid_mat_adj), Phenotype = pheno_resid),
    target = "Phenotype"
  )

  # ─────────────────────────────────────────────────────────────────────────────
  # Fit Final Model on ALL data (for SHAP interpretation only)
  # ─────────────────────────────────────────────────────────────────────────────
  message("Fitting final model on ALL data for SHAP analysis...")

  rf_final <- ranger(
    x = as.data.frame(lipid_mat_adj),
    y = pheno_resid,
    num.trees = 1000,
    mtry = tune_res$recommended.pars$mtry,
    min.node.size = tune_res$recommended.pars$min.node.size,
    sample.fraction = tune_res$recommended.pars$sample.fraction,
    importance = "none",
    num.threads = parallel::detectCores() - 1,
    seed = 42
  )

  # ─────────────────────────────────────────────────────────────────────────────
  # TreeSHAP Analysis
  # ─────────────────────────────────────────────────────────────────────────────
  message("Computing TreeSHAP values (this may take a moment)...")

  final_X_df <- as.data.frame(lipid_mat_adj)
  um <- unify(rf_final, final_X_df)
  ts <- treeshap(um, final_X_df)

  # Global importance: mean |SHAP|
  shap_rank <- tibble(
    Lipid = colnames(ts$shaps),
    MeanAbsSHAP = colMeans(abs(ts$shaps))
  ) %>%
    arrange(desc(MeanAbsSHAP))

  # Save SHAP results
  write.csv(shap_rank, paste0("table/supp/SuppTable_lipid_SHAP_", PHENOTYPE_NAME, ".csv"),
            row.names = FALSE)

  message("✓ SHAP values computed and saved")

  # ─────────────────────────────────────────────────────────────────────────────
  # SHAP Stability Analysis (bootstrap)
  # ─────────────────────────────────────────────────────────────────────────────
  shap_stability <- NULL
  if (RUN_SHAP_STABILITY) {
    message("Running bootstrap SHAP stability analysis...")

    shap_stability <- compute_shap_stability_bootstrap(
      X = final_X_df,
      y = pheno_resid,
      mtry = tune_res$recommended.pars$mtry,
      min_node_size = tune_res$recommended.pars$min.node.size,
      sample_fraction = tune_res$recommended.pars$sample.fraction,
      n_boot = SHAP_STABILITY_N_BOOT,
      top_k = SHAP_STABILITY_TOP_K,
      num_trees = SHAP_STABILITY_NUM_TREES,
      seed = 2026
    )

    if (!is.null(shap_stability)) {
      shap_stability <- shap_stability %>%
        mutate(Class = get_lipid_class(Lipid))

      write.csv(
        shap_stability,
        paste0("table/supp/SuppTable_SHAP_Stability_", PHENOTYPE_NAME, ".csv"),
        row.names = FALSE
      )

      stable_top <- shap_stability %>%
        slice_head(n = 20) %>%
        mutate(Lipid = factor(Lipid, levels = rev(Lipid)))

      p_stab_a <- ggplot(stable_top, aes(x = Lipid, y = SelectionFreqTopK, fill = Class)) +
        geom_col(width = 0.7, color = "black", linewidth = 0.2) +
        geom_errorbar(aes(ymin = SelectionCI_Low, ymax = SelectionCI_High), width = 0.2) +
        coord_flip() +
        scale_fill_manual(values = class_colors, na.value = "grey70") +
        scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
        labs(
          title = paste0("A) SHAP Stability: Top-", SHAP_STABILITY_TOP_K, " Selection Frequency"),
          x = NULL, y = "Bootstrap Selection Frequency"
        ) +
        theme_minimal(base_size = 11) +
        theme(
          axis.text = element_text(color = "black", size = 9),
          axis.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0),
          legend.position = "right",
          panel.grid.major.y = element_blank(),
          axis.line.x = element_line(color = "black", linewidth = 0.5),
          axis.line.y = element_line(color = "black", linewidth = 0.5)
        )

      p_stab_b <- ggplot(stable_top, aes(y = Lipid, x = MedianRank, color = Class)) +
        geom_segment(aes(x = RankCI_Low, xend = RankCI_High, yend = Lipid), linewidth = 0.8, alpha = 0.8) +
        geom_point(size = 2) +
        scale_color_manual(values = class_colors, na.value = "grey70") +
        labs(
          title = "B) Rank Stability (Median and 95% CI)",
          x = "Rank (lower is more stable)", y = NULL
        ) +
        theme_minimal(base_size = 11) +
        theme(
          axis.text = element_text(color = "black", size = 9),
          axis.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0),
          legend.position = "none",
          panel.grid.major.y = element_blank(),
          axis.line.x = element_line(color = "black", linewidth = 0.5),
          axis.line.y = element_line(color = "black", linewidth = 0.5)
        )

      suppfig_shap_stability <- p_stab_a + p_stab_b + plot_layout(widths = c(1.2, 1))
      if (isTRUE(SAVE_S9B_STABILITY_FIG)) {
        save_fig(suppfig_shap_stability, paste0("SuppFig_S9b_SHAP_Stability_", PHENOTYPE_NAME, ".png"),
                 w = 14, h = 8, subdir = "supp")
        save_fig(suppfig_shap_stability, paste0("SuppFig_S9b_SHAP_Stability_", PHENOTYPE_NAME, ".pdf"),
                 w = 14, h = 8, subdir = "supp")
        message("✓ SHAP stability analysis complete (table + supplementary figure)")
      } else {
        message("✓ SHAP stability analysis complete (table only; S9b figure disabled)")
      }
    }
  }

  # ─────────────────────────────────────────────────────────────────────────────
  # Fig 5A: Observed vs Predicted
  # ─────────────────────────────────────────────────────────────────────────────
  df_pred <- tibble(Observed = test_y, Predicted = preds_test)

  fig5a <- ggplot(df_pred, aes(x = Observed, y = Predicted)) +
    geom_point(alpha = 0.6, size = 2, color = "#440154FF") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
    geom_smooth(method = "lm", se = FALSE, color = "#FDE725FF", linewidth = 1) +
    annotate("text", x = -Inf, y = Inf,
             label = sprintf("R² (mean±SD) = %.3f±%.3f\nRMSE (mean±SD) = %.3f±%.3f",
                             r2_test, r2_test_sd, rmse_test, rmse_test_sd),
             hjust = -0.1, vjust = 1.3, size = 4, fontface = "bold") +
    labs(title = paste0("A) Observed vs Predicted (", PHENOTYPE_NAME, ")"),
         x = "Observed (PC-residualized)",
         y = "RF Predicted") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5)
    )

  # ─────────────────────────────────────────────────────────────────────────────
  # Fig 5B: Top SHAP Importance Bar
  # ─────────────────────────────────────────────────────────────────────────────
  top_n_shap <- 20
  top_shap <- shap_rank %>% slice_head(n = top_n_shap)

  # Extract lipid class for coloring
  top_shap <- top_shap %>%
    mutate(Class = get_lipid_class(Lipid))

  fig5b <- ggplot(top_shap, aes(x = reorder(Lipid, MeanAbsSHAP), y = MeanAbsSHAP, fill = Class)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.2) +
    coord_flip() +
    scale_fill_manual(values = class_colors, na.value = "grey70") +
    labs(title = paste0("B) Top ", top_n_shap, " Lipids by SHAP Importance"),
         x = NULL, y = "Mean |SHAP|") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text = element_text(color = "black", size = 9),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0),
      legend.position = "right",
      panel.grid.major.y = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5)
    )

  # ─────────────────────────────────────────────────────────────────────────────
  # Fig 5C: SHAP Beeswarm (for top features)
  # ─────────────────────────────────────────────────────────────────────────────
  # Select top features for beeswarm
  top_features <- shap_rank %>% slice_head(n = 15) %>% pull(Lipid)

  shap_long <- as_tibble(ts$shaps) %>%
    mutate(Sample = row_number()) %>%
    pivot_longer(-Sample, names_to = "Lipid", values_to = "SHAP") %>%
    filter(Lipid %in% top_features) %>%
    left_join(
      as_tibble(final_X_df) %>%
        mutate(Sample = row_number()) %>%
        pivot_longer(-Sample, names_to = "Lipid", values_to = "FeatureValue"),
      by = c("Sample", "Lipid")
    ) %>%
    mutate(Lipid = factor(Lipid, levels = rev(top_features)))

  fig5c <- ggplot(shap_long, aes(x = SHAP, y = Lipid, color = FeatureValue)) +
    geom_jitter(height = 0.25, size = 1, alpha = 0.6) +
    scale_color_gradient2(low = "blue", mid = "grey90", high = "red",
                          midpoint = 0, name = "Feature\nValue") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    labs(title = "C) SHAP Value Distribution", x = "SHAP Value", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text = element_text(color = "black", size = 9),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0),
      legend.position = "right",
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5)
    )

  # ─────────────────────────────────────────────────────────────────────────────
  # Fig 5D: Cumulative SHAP importance (Pareto)
  # ─────────────────────────────────────────────────────────────────────────────
  shap_cum <- shap_rank %>%
    mutate(
      cum_imp = cumsum(MeanAbsSHAP) / sum(MeanAbsSHAP),
      rank = row_number()
    )

  cutoff_80 <- which(shap_cum$cum_imp >= 0.80)[1]

  fig5d <- ggplot(shap_cum, aes(x = rank, y = cum_imp)) +
    geom_line(linewidth = 1, color = "#440154FF") +
    geom_point(size = 0.5, color = "#440154FF") +
    geom_hline(yintercept = 0.80, linetype = "dashed", color = "red") +
    geom_vline(xintercept = cutoff_80, linetype = "dashed", color = "red") +
    annotate("text", x = cutoff_80 + 5, y = 0.82,
             label = sprintf("80%% at rank %d", cutoff_80),
             hjust = 0, color = "red", size = 3.5) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(title = "D) Cumulative SHAP Importance",
         x = "Feature Rank", y = "Cumulative Importance") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5)
    )

  # ─────────────────────────────────────────────────────────────────────────────
  # Combine Figure 5
  # ─────────────────────────────────────────────────────────────────────────────
  fig5_top <- fig5a + fig5d + plot_layout(widths = c(1, 1))
  fig5_bottom <- fig5b + fig5c + plot_layout(widths = c(1, 1.2))

  figure5 <- fig5_top / fig5_bottom + plot_layout(heights = c(1, 1.3))

  if (isTRUE(SAVE_LEGACY_RF_FIG)) {
    save_fig(figure5, paste0("Legacy_Figure4_RF_SHAP_", PHENOTYPE_NAME, ".png"),
             w = 14, h = 12, subdir = "supp")
    save_fig(figure5, paste0("Legacy_Figure4_RF_SHAP_", PHENOTYPE_NAME, ".pdf"),
             w = 14, h = 12, subdir = "supp")
    message("\n✓ Legacy RF composite saved (supp).")
  } else {
    message("\n✓ Legacy RF composite generated in memory (not saved; final Figure 4 is assembled later).")
  }

} else {
  # Load pre-computed results if available
  shap_file <- paste0("table/supp/SuppTable_lipid_SHAP_", PHENOTYPE_NAME, ".csv")

  if (file.exists(shap_file)) {
    message("Loading pre-computed SHAP results...")
    shap_rank <- vroom(shap_file, show_col_types = FALSE)

    top_shap <- shap_rank %>%
      slice_head(n = 20) %>%
      mutate(Class = get_lipid_class(Lipid))

    fig5_simple <- ggplot(top_shap, aes(x = reorder(Lipid, MeanAbsSHAP), y = MeanAbsSHAP, fill = Class)) +
      geom_col(width = 0.7, color = "black", linewidth = 0.2) +
      coord_flip() +
      scale_fill_manual(values = class_colors, na.value = "grey70") +
      labs(title = paste0("Top 20 Lipids Predicting ", PHENOTYPE_NAME),
           x = NULL, y = "Mean |SHAP|") +
      theme_minimal(base_size = 12) +
      theme(
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "right",
        panel.grid.major.y = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5)
      )

    if (isTRUE(SAVE_LEGACY_RF_FIG)) {
      save_fig(fig5_simple, paste0("Legacy_Figure4_SHAP_", PHENOTYPE_NAME, ".png"),
               w = 10, h = 8, subdir = "supp")
      message("✓ Legacy SHAP-only plot saved (supp).")
    } else {
      message("✓ Legacy SHAP-only plot generated in memory (not saved).")
    }

  } else {
    message("No pre-computed SHAP results found.")
    message("Set RUN_RF_FRESH <- TRUE and ensure phenotype file exists to run full pipeline.")
  }
}


# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
#
#     SUPPLEMENTARY TABLE S1: MULTI-PHENOTYPE SHAP COMPARISON
#
# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════

message("\n══════════════════════════════════════════════════════════════")
message("SUPPLEMENTARY TABLE S1: MULTI-PHENOTYPE SHAP COMPARISON")
message("══════════════════════════════════════════════════════════════\n")

# Define all phenotypes to analyze
# Load Wei Jialu flowering time data with multiple environments
FT_WeiJialu_file <- "data/phenotypes/Flowering_Time_Plant_Height_WeiJialu.csv"
if (file.exists(FT_WeiJialu_file)) {
  FT_WeiJialu <- read.csv(FT_WeiJialu_file) %>% dplyr::select(-4)
  message("Wei Jialu environments: ", paste(unique(FT_WeiJialu$Env), collapse = ", "))

  # Select the different environments for analysis (save as temp files for the pipeline)
  FT_WeiJialu_MI20 <- FT_WeiJialu %>% dplyr::filter(Env == "MI20") %>% dplyr::select(-2)
  FT_WeiJialu_NE20 <- FT_WeiJialu %>% dplyr::filter(Env == "NE20") %>% dplyr::select(-2)
  FT_WeiJialu_IA21_2 <- FT_WeiJialu %>% dplyr::filter(Env == "IA21_2") %>% dplyr::select(-2)
  FT_WeiJialu_IA21_1 <- FT_WeiJialu %>% dplyr::filter(Env == "IA21_1") %>% dplyr::select(-2)

  # Save as temp files so the pipeline can read them
  dir.create("data/phenotypes/temp", showWarnings = FALSE, recursive = TRUE)
  write.csv(FT_WeiJialu_MI20, "data/phenotypes/temp/FT_WeiJialu_MI20.csv", row.names = FALSE)
  write.csv(FT_WeiJialu_NE20, "data/phenotypes/temp/FT_WeiJialu_NE20.csv", row.names = FALSE)
  write.csv(FT_WeiJialu_IA21_2, "data/phenotypes/temp/FT_WeiJialu_IA21_2.csv", row.names = FALSE)
  write.csv(FT_WeiJialu_IA21_1, "data/phenotypes/temp/FT_WeiJialu_IA21_1.csv", row.names = FALSE)
}

PHENOTYPE_CONFIG <- list(
  # Our field data - Plant Height
  PH_ours = list(file = "data/phenotypes/control_field_phenotypes.csv", col = 2, name = "Plant Height (Our Data)"),

  # Others' research - Plant Height for comparison
  PH_Boyles = list(file = "data/phenotypes/plantheight_diameter_SAP.csv", col = 2, name = "Plant Height (Boyles et al.)"),

  # Flowering Time - Our data (DTA = Days to Anthesis)
  FT_ours = list(file = "data/phenotypes/control_field_phenotypes.csv", col = 3, name = "Flowering Time DTA (Our Data)"),

  # Flowering Time - Wei Jialu across environments (GDD = Growing Degree Days)
  FT_MI20 = list(file = "data/phenotypes/temp/FT_WeiJialu_MI20.csv", col = 2, name = "Flowering Time MI20 GDD (Wei et al.)"),
  FT_NE20 = list(file = "data/phenotypes/temp/FT_WeiJialu_NE20.csv", col = 2, name = "Flowering Time NE20 GDD (Wei et al.)"),
  FT_IA21_1 = list(file = "data/phenotypes/temp/FT_WeiJialu_IA21_1.csv", col = 2, name = "Flowering Time IA21_1 GDD (Wei et al.)"),
  FT_IA21_2 = list(file = "data/phenotypes/temp/FT_WeiJialu_IA21_2.csv", col = 2, name = "Flowering Time IA21_2 GDD (Wei et al.)")
)


# Function to run RF + SHAP for a phenotype with leakage-safe evaluation
run_rf_shap_proper <- function(lipids_df, pheno_config, pc_df = NULL, n_outer_repeats = 10,
                               run_stability = TRUE, stability_n_boot = 100,
                               stability_top_k = 20, stability_num_trees = 500) {

  # Load phenotype
  if (!file.exists(pheno_config$file)) {
    message("    File not found: ", pheno_config$file)
    return(NULL)
  }

  pheno <- vroom(pheno_config$file, delim = ",", show_col_types = FALSE) %>%
    dplyr::select(1, pheno_config$col) %>%
    dplyr::rename(Line = 1, Phenotype = 2) %>%
    dplyr::filter(!is.na(Phenotype)) %>%
    dplyr::mutate(Phenotype = as.numeric(Phenotype)) %>%
    dplyr::filter(!is.na(Phenotype))

  # Merge with lipids
  dat <- lipids_df %>%
    inner_join(pheno, by = "Line")

  if (nrow(dat) < 50) {
    message("    Too few samples after merge: ", nrow(dat))
    return(NULL)
  }

  message("    Samples: ", nrow(dat))

  # Align PCs if available
  pc_mat <- NULL
  if (!is.null(pc_df)) {
    common_lines <- intersect(dat$Line, pc_df$Line)
    if (length(common_lines) >= 50) {
      dat <- dat %>% filter(Line %in% common_lines) %>% arrange(Line)
      pc_sub <- pc_df %>% filter(Line %in% common_lines) %>% arrange(Line)
      pc_mat <- pc_sub %>% dplyr::select(-Line) %>% as.data.frame()
      message("    PC-aligned samples: ", nrow(dat))
    }
  }

  y_raw <- dat$Phenotype
  X_raw <- dat %>% dplyr::select(-Line, -Phenotype) %>% as.matrix()

  rf_eval <- run_repeated_rf_eval(
    X = X_raw,
    y = y_raw,
    pc_df = pc_mat,
    train_prop = 0.8,
    n_repeats = n_outer_repeats,
    seed_base = 42,
    tune_num_trees = 500,
    fit_num_trees = 1000,
    num_threads = max(1, parallel::detectCores() - 1)
  )

  test_r2 <- rf_eval$summary$R2_Mean
  test_r2_sd <- rf_eval$summary$R2_SD
  test_rmse <- rf_eval$summary$RMSE_Mean
  test_rmse_sd <- rf_eval$summary$RMSE_SD

  message("    Repeated-split R² (mean ± SD): ", round(test_r2, 4), " ± ", round(test_r2_sd, 4))
  message("    Repeated-split RMSE (mean ± SD): ", round(test_rmse, 4), " ± ", round(test_rmse_sd, 4))

  tune_res <- list(recommended.pars = rf_eval$median_params)
  message("    Median tuned mtry: ", tune_res$recommended.pars$mtry,
          " | min.node.size: ", tune_res$recommended.pars$min.node.size)

  # Fit on ALL data for SHAP interpretation
  full_adj <- residualize_full_by_pc(X_raw, y_raw, pc_df = pc_mat)
  X <- full_adj$X
  y <- full_adj$y

  rf_final <- ranger(
    x = as.data.frame(X),
    y = y,
    num.trees = 1000,
    mtry = tune_res$recommended.pars$mtry,
    min.node.size = tune_res$recommended.pars$min.node.size,
    sample.fraction = tune_res$recommended.pars$sample.fraction,
    importance = "none",
    num.threads = parallel::detectCores() - 1,
    seed = 42
  )

  # TreeSHAP on full data
  X_df <- as.data.frame(X)
  um <- unify(rf_final, X_df)
  ts <- treeshap(um, X_df)

  # Global importance
  shap_df <- tibble(
    Lipid = colnames(ts$shaps),
    MeanAbsSHAP = colMeans(abs(ts$shaps))
  ) %>%
    arrange(desc(MeanAbsSHAP)) %>%
    mutate(
      Rank = row_number(),
      Class = get_lipid_class(Lipid)
    )

  shap_stability <- NULL
  if (run_stability) {
    message("    Running SHAP stability bootstrap...")
    shap_stability <- compute_shap_stability_bootstrap(
      X = X_df,
      y = y,
      mtry = tune_res$recommended.pars$mtry,
      min_node_size = tune_res$recommended.pars$min.node.size,
      sample_fraction = tune_res$recommended.pars$sample.fraction,
      n_boot = stability_n_boot,
      top_k = stability_top_k,
      num_trees = stability_num_trees,
      seed = 2026
    )

    if (!is.null(shap_stability)) {
      shap_stability <- shap_stability %>%
        mutate(
          Class = get_lipid_class(Lipid),
          Rank = row_number()
        )
    }
  }

  list(
    shap = shap_df,
    stability = shap_stability,
    r2 = test_r2,
    r2_sd = test_r2_sd,
    rmse = test_rmse,
    rmse_sd = test_rmse_sd,
    n_samples = length(y),
    n_train = rf_eval$representative$n_train,
    n_test = rf_eval$representative$n_test,
    mtry = tune_res$recommended.pars$mtry,
    min_node_size = tune_res$recommended.pars$min.node.size
  )
}

# Run analysis for all phenotypes
RUN_MULTI_PHENO <- TRUE  # Set to FALSE to skip
SAVE_FT_STABILITY_FIG <- FALSE  # Keep FT stability table, but do not output FT stability figure panel
SAVE_PH_COMPARISON_FIG <- FALSE  # Keep table, but do not output standalone S9a
SAVE_FT_SOURCE_HEATMAP_FIG <- FALSE  # Do not output standalone S11 unless explicitly needed

if (RUN_MULTI_PHENO) {
  MULTI_SHAP_STABILITY_N_BOOT <- 100
  MULTI_SHAP_STABILITY_TOP_K <- 20
  MULTI_SHAP_STABILITY_NUM_TREES <- 500

  # Load PC file if available
  pc_file <- "table/PCA_SAP.csv"
  PC <- NULL
  if (file.exists(pc_file)) {
    PC <- vroom(pc_file, show_col_types = FALSE) %>%
      dplyr::select(1:4) %>%
      dplyr::rename(Line = 1)
    message("Loaded population structure PCs")
  }

  # Store results
  all_shap_results <- list()
  all_stability_results <- list()
  model_summary <- list()

  for (pheno_id in names(PHENOTYPE_CONFIG)) {
    config <- PHENOTYPE_CONFIG[[pheno_id]]
    message("\nProcessing: ", config$name)

    tryCatch({
      result <- run_rf_shap_proper(
        lipids_rf, config, PC,
        run_stability = FALSE,
        stability_n_boot = MULTI_SHAP_STABILITY_N_BOOT,
        stability_top_k = MULTI_SHAP_STABILITY_TOP_K,
        stability_num_trees = MULTI_SHAP_STABILITY_NUM_TREES
      )

      if (!is.null(result)) {
        all_shap_results[[pheno_id]] <- result$shap %>%
          mutate(Phenotype = config$name, PhenotypeID = pheno_id)

        if (!is.null(result$stability)) {
          all_stability_results[[pheno_id]] <- result$stability %>%
            mutate(Phenotype = config$name, PhenotypeID = pheno_id)

          write.csv(
            all_stability_results[[pheno_id]],
            paste0("table/supp/SuppTable_SHAP_Stability_", pheno_id, ".csv"),
            row.names = FALSE
          )
        }

        model_summary[[pheno_id]] <- tibble(
          PhenotypeID = pheno_id,
          Phenotype = config$name,
          N_Samples = result$n_samples,
          N_Train = result$n_train,
          N_Test = result$n_test,
          Test_R2 = round(result$r2, 4),
          Test_R2_SD = round(result$r2_sd, 4),
          Test_RMSE = round(result$rmse, 4),
          Test_RMSE_SD = round(result$rmse_sd, 4),
          Tuned_mtry = result$mtry,
          Tuned_min_node_size = result$min_node_size
        )

        message("    ✓ Test R²: ", round(result$r2, 4))
      }
    }, error = function(e) {
      message("    ERROR: ", e$message)
    })
  }

  # ─────────────────────────────────────────────────────────────────────────────
  # Create Supplementary Table S1: Top 20 SHAP lipids per phenotype
  # ─────────────────────────────────────────────────────────────────────────────
  if (length(all_shap_results) > 0) {

    # Combine all results
    combined_shap <- bind_rows(all_shap_results)

    # Top 20 per phenotype
    top20_per_pheno <- combined_shap %>%
      filter(Rank <= 20) %>%
      dplyr::select(Phenotype, Rank, Lipid, Class, MeanAbsSHAP) %>%
      arrange(Phenotype, Rank)

    # Save as wide format table (one column per phenotype)
    top20_wide <- top20_per_pheno %>%
      dplyr::select(Phenotype, Rank, Lipid) %>%
      pivot_wider(names_from = Phenotype, values_from = Lipid)

    # Save tables
    write.csv(top20_per_pheno, "table/supp/SuppTable_S1_Top20_SHAP_all_phenotypes_long.csv", row.names = FALSE)
    write.csv(top20_wide, "table/supp/SuppTable_S1_Top20_SHAP_all_phenotypes_wide.csv", row.names = FALSE)

    # Model summary table
    model_df <- bind_rows(model_summary)
    write.csv(model_df, "table/supp/SuppTable_S1_RF_model_summary.csv", row.names = FALSE)

    message("\n✓ Supplementary Table S1 saved:")
    message("  - table/supp/SuppTable_S1_Top20_SHAP_all_phenotypes_long.csv")
    message("  - table/supp/SuppTable_S1_Top20_SHAP_all_phenotypes_wide.csv")
    message("  - table/supp/SuppTable_S1_RF_model_summary.csv")

    # ─────────────────────────────────────────────────────────────────────────────
    # Plant Height Comparison: Our data vs Boyles et al.
    # ─────────────────────────────────────────────────────────────────────────────
    if ("PH_ours" %in% names(all_shap_results) && "PH_Boyles" %in% names(all_shap_results)) {

      message("\n─── PLANT HEIGHT COMPARISON: Our Data vs Boyles et al. ───")

      ph_ours_top <- all_shap_results$PH_ours %>% filter(Rank <= 20)
      ph_boyles_top <- all_shap_results$PH_Boyles %>% filter(Rank <= 20)

      ph_ours_lipids <- ph_ours_top %>% pull(Lipid)
      ph_boyles_lipids <- ph_boyles_top %>% pull(Lipid)

      overlap_ph <- intersect(ph_ours_lipids, ph_boyles_lipids)
      only_ours <- setdiff(ph_ours_lipids, ph_boyles_lipids)
      only_boyles <- setdiff(ph_boyles_lipids, ph_ours_lipids)

      message("  Top 20 overlap: ", length(overlap_ph), " lipids")
      if (length(overlap_ph) > 0) {
        message("  Overlapping lipids: ", paste(overlap_ph, collapse = ", "))
      }

      # Create comparison table
      ph_comparison <- tibble(
        Rank = 1:20,
        PH_OurData = ph_ours_lipids,
        PH_OurData_SHAP = ph_ours_top$MeanAbsSHAP,
        PH_Boyles = ph_boyles_lipids,
        PH_Boyles_SHAP = ph_boyles_top$MeanAbsSHAP
      )

      write.csv(ph_comparison, "table/supp/SuppTable_PlantHeight_comparison.csv", row.names = FALSE)

      # ─────────────────────────────────────────────────────────────────────────────
      # SUPPLEMENTARY FIGURE: Plant Height SHAP Comparison
      # ─────────────────────────────────────────────────────────────────────────────
      # Create Venn-style bar showing overlap
      ph_venn_data <- tibble(
        Category = c("Our Data Only", "Both", "Boyles Only"),
        Count = c(length(only_ours), length(overlap_ph), length(only_boyles)),
        Type = c("Unique", "Overlap", "Unique")
      ) %>%
        mutate(Category = factor(Category, levels = c("Our Data Only", "Both", "Boyles Only")))

      suppfig_ph_venn <- ggplot(ph_venn_data, aes(x = Category, y = Count, fill = Type)) +
        geom_col(width = 0.6, color = "black", linewidth = 0.3) +
        geom_text(aes(label = Count), vjust = -0.3, size = 5, fontface = "bold") +
        scale_fill_manual(values = c("Overlap" = "#21908CFF", "Unique" = "#FDE725FF")) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
        labs(title = "A) Top 20 SHAP Lipids: Plant Height Comparison",
             subtitle = "Our Data vs Boyles et al.",
             x = NULL, y = "Number of Lipids") +
        theme_minimal(base_size = 12) +
        theme(
          legend.position = "none",
          axis.text = element_text(color = "black", size = 11),
          axis.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
          panel.grid.major.x = element_blank(),
          axis.line.x = element_line(color = "black", linewidth = 0.5),
          axis.line.y = element_line(color = "black", linewidth = 0.5)
        )

      # Scatter plot comparing SHAP values for common lipids
      common_lipids_all <- union(ph_ours_lipids, ph_boyles_lipids)
      ph_scatter_data <- tibble(Lipid = common_lipids_all) %>%
        left_join(all_shap_results$PH_ours %>% dplyr::select(Lipid, MeanAbsSHAP) %>%
                    dplyr::rename(SHAP_Ours = MeanAbsSHAP), by = "Lipid") %>%
        left_join(all_shap_results$PH_Boyles %>% dplyr::select(Lipid, MeanAbsSHAP) %>%
                    dplyr::rename(SHAP_Boyles = MeanAbsSHAP), by = "Lipid") %>%
        mutate(
          SHAP_Ours = replace_na(SHAP_Ours, 0),
          SHAP_Boyles = replace_na(SHAP_Boyles, 0),
          InBoth = Lipid %in% overlap_ph,
          Class = get_lipid_class(Lipid)
        )

      suppfig_ph_scatter <- ggplot(ph_scatter_data, aes(x = SHAP_Ours, y = SHAP_Boyles)) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
        geom_point(aes(color = Class, shape = InBoth), size = 3, alpha = 0.7) +
        geom_text_repel(data = ph_scatter_data %>% filter(InBoth),
                        aes(label = Lipid), size = 2.5, max.overlaps = 10) +
        scale_color_manual(values = class_colors, na.value = "grey70") +
        scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 17),
                           labels = c("TRUE" = "In both top 20", "FALSE" = "Unique")) +
        labs(title = "B) SHAP Importance Correlation",
             subtitle = sprintf("Top 20 lipids (r = %.2f for overlap)",
                                cor(ph_scatter_data$SHAP_Ours[ph_scatter_data$InBoth],
                                    ph_scatter_data$SHAP_Boyles[ph_scatter_data$InBoth],
                                    use = "complete.obs")),
             x = "Mean |SHAP| - Our Data", y = "Mean |SHAP| - Boyles et al.") +
        theme_minimal(base_size = 12) +
        theme(
          legend.position = "right",
          axis.text = element_text(color = "black"),
          axis.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
          axis.line.x = element_line(color = "black", linewidth = 0.5),
          axis.line.y = element_line(color = "black", linewidth = 0.5)
        )

      suppfig_ph_comparison <- suppfig_ph_venn + suppfig_ph_scatter +
        plot_layout(widths = c(1, 1.3))

      if (isTRUE(SAVE_PH_COMPARISON_FIG)) {
        save_fig(suppfig_ph_comparison, "SuppFig_S9a_PlantHeight_SHAP_Comparison.png",
                 w = 14, h = 6, subdir = "supp")
        save_fig(suppfig_ph_comparison, "SuppFig_S9a_PlantHeight_SHAP_Comparison.pdf",
                 w = 14, h = 6, subdir = "supp")
        message("  ✓ Plant Height comparison figure saved")
      } else {
        message("  ✓ Plant Height comparison table saved (S9a figure output disabled)")
      }
    }

    # ─────────────────────────────────────────────────────────────────────────────
    # Flowering Time SHAP Stability: Ours vs Other Environments
    # ─────────────────────────────────────────────────────────────────────────────
    ft_stability_ids <- c("FT_ours", "FT_MI20", "FT_NE20", "FT_IA21_1", "FT_IA21_2")
    ft_stability_available <- ft_stability_ids[ft_stability_ids %in% names(all_stability_results)]

    if ("FT_ours" %in% ft_stability_available && length(setdiff(ft_stability_available, "FT_ours")) >= 1) {
      message("\n─── FLOWERING TIME SHAP STABILITY (OURS VS OTHERS) ───")

      quant_safe <- function(x, p) {
        x <- x[is.finite(x)]
        if (length(x) == 0) return(NA_real_)
        as.numeric(quantile(x, probs = p, na.rm = TRUE, type = 7))
      }

      make_stability_panels <- function(stab_df, title_a, title_b, top_n = 20) {
        plot_df <- stab_df %>%
          arrange(desc(SelectionFreqTopK), MedianRank, desc(MeanAbsSHAP_Boot)) %>%
          slice_head(n = top_n) %>%
          mutate(Lipid = factor(Lipid, levels = rev(Lipid)))

        p_a <- ggplot(plot_df, aes(x = Lipid, y = SelectionFreqTopK, fill = Class)) +
          geom_col(width = 0.7, color = "black", linewidth = 0.2) +
          geom_errorbar(aes(ymin = SelectionCI_Low, ymax = SelectionCI_High), width = 0.2) +
          coord_flip() +
          scale_fill_manual(values = class_colors, na.value = "grey70") +
          scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
          labs(title = title_a, x = NULL, y = "Bootstrap Selection Frequency") +
          theme_minimal(base_size = 11) +
          theme(
            axis.text = element_text(color = "black", size = 9),
            axis.title = element_text(face = "bold"),
            plot.title = element_text(face = "bold", hjust = 0),
            legend.position = "right",
            panel.grid.major.y = element_blank(),
            axis.line.x = element_line(color = "black", linewidth = 0.5),
            axis.line.y = element_line(color = "black", linewidth = 0.5)
          )

        p_b <- ggplot(plot_df, aes(y = Lipid, x = MedianRank, color = Class)) +
          geom_segment(aes(x = RankCI_Low, xend = RankCI_High, yend = Lipid), linewidth = 0.8, alpha = 0.8) +
          geom_point(size = 2) +
          scale_color_manual(values = class_colors, na.value = "grey70") +
          labs(title = title_b, x = "Rank (lower is more stable)", y = NULL) +
          theme_minimal(base_size = 11) +
          theme(
            axis.text = element_text(color = "black", size = 9),
            axis.title = element_text(face = "bold"),
            plot.title = element_text(face = "bold", hjust = 0),
            legend.position = "none",
            panel.grid.major.y = element_blank(),
            axis.line.x = element_line(color = "black", linewidth = 0.5),
            axis.line.y = element_line(color = "black", linewidth = 0.5)
          )

        list(a = p_a, b = p_b)
      }

      ft_ours_stability <- all_stability_results$FT_ours %>%
        arrange(desc(SelectionFreqTopK), MedianRank, desc(MeanAbsSHAP_Boot))

      ft_other_ids <- setdiff(ft_stability_available, "FT_ours")
      ft_others_stability <- bind_rows(lapply(ft_other_ids, function(id) {
        all_stability_results[[id]] %>%
          mutate(EnvironmentID = id)
      })) %>%
        group_by(Lipid) %>%
        summarise(
          SelFreqMean = mean(SelectionFreqTopK, na.rm = TRUE),
          SelFreqLow = quant_safe(SelectionFreqTopK, 0.025),
          SelFreqHigh = quant_safe(SelectionFreqTopK, 0.975),
          RankMean = mean(MedianRank, na.rm = TRUE),
          RankLow = quant_safe(MedianRank, 0.025),
          RankHigh = quant_safe(MedianRank, 0.975),
          SHAPMean = mean(MeanAbsSHAP_Boot, na.rm = TRUE),
          SHAPLow = quant_safe(MeanAbsSHAP_Boot, 0.025),
          SHAPHigh = quant_safe(MeanAbsSHAP_Boot, 0.975),
          N_Environments = dplyr::n_distinct(EnvironmentID),
          .groups = "drop"
        ) %>%
        mutate(
          SelectionFreqTopK = SelFreqMean,
          SelectionCI_Low = SelFreqLow,
          SelectionCI_High = SelFreqHigh,
          MedianRank = RankMean,
          RankCI_Low = RankLow,
          RankCI_High = RankHigh,
          MeanAbsSHAP_Boot = SHAPMean,
          SHAP_CI_Low = SHAPLow,
          SHAP_CI_High = SHAPHigh,
          Class = get_lipid_class(Lipid),
          Rank = row_number()
        ) %>%
        dplyr::select(
          Lipid, SelectionFreqTopK, SelectionCI_Low, SelectionCI_High,
          MedianRank, RankCI_Low, RankCI_High, MeanAbsSHAP_Boot,
          SHAP_CI_Low, SHAP_CI_High, N_Environments, Class, Rank
        ) %>%
        arrange(desc(SelectionFreqTopK), MedianRank, desc(MeanAbsSHAP_Boot))

      write.csv(ft_others_stability,
                "table/supp/SuppTable_FT_SHAP_Stability_Others_Aggregated.csv",
                row.names = FALSE)

      if (isTRUE(SAVE_FT_STABILITY_FIG)) {
        ft_panels_ours <- make_stability_panels(
          ft_ours_stability,
          "A) FT Ours: Selection Frequency",
          "B) FT Ours: Rank Stability",
          top_n = MULTI_SHAP_STABILITY_TOP_K
        )
        ft_panels_others <- make_stability_panels(
          ft_others_stability,
          "C) FT Others: Selection Frequency",
          "D) FT Others: Rank Stability",
          top_n = MULTI_SHAP_STABILITY_TOP_K
        )

        suppfig_ft_stability <- (ft_panels_ours$a + ft_panels_ours$b) /
          (ft_panels_others$a + ft_panels_others$b) +
          plot_layout(widths = c(1.2, 1), heights = c(1, 1))

        save_fig(suppfig_ft_stability, "SuppFig_S10b_FloweringTime_SHAP_Stability.png",
                 w = 14, h = 12, subdir = "supp")
        save_fig(suppfig_ft_stability, "SuppFig_S10b_FloweringTime_SHAP_Stability.pdf",
                 w = 14, h = 12, subdir = "supp")

        message("  ✓ Saved S10b (Flowering Time SHAP stability: ours vs others)")
      } else {
        message("  ✓ FT stability tables saved (figure output disabled by SAVE_FT_STABILITY_FIG = FALSE)")
      }
    }

    # ─────────────────────────────────────────────────────────────────────────────
    # Flowering Time Across Environments Comparison
    # ─────────────────────────────────────────────────────────────────────────────
    ft_envs <- c("FT_ours", "FT_MI20", "FT_NE20", "FT_IA21_1", "FT_IA21_2")
    ft_available <- ft_envs[ft_envs %in% names(all_shap_results)]

    if (length(ft_available) >= 2) {
      message("\n─── FLOWERING TIME ACROSS ENVIRONMENTS ───")
      message("  Environments available: ", paste(ft_available, collapse = ", "))

      # Get top 20 lipids for each environment
      ft_top_lists <- lapply(ft_available, function(env) {
        all_shap_results[[env]] %>% filter(Rank <= 20) %>% pull(Lipid)
      })
      names(ft_top_lists) <- ft_available

      # Count how many environments each lipid appears in
      all_ft_lipids <- unique(unlist(ft_top_lists))
      ft_overlap_counts <- sapply(all_ft_lipids, function(lip) {
        sum(sapply(ft_top_lists, function(x) lip %in% x))
      })

      ft_overlap_df <- tibble(
        Lipid = all_ft_lipids,
        N_Environments = ft_overlap_counts,
        Class = get_lipid_class(all_ft_lipids)
      ) %>%
        arrange(desc(N_Environments), Lipid)

      # Save table
      write.csv(ft_overlap_df, "table/supp/SuppTable_FT_Environment_Overlap.csv", row.names = FALSE)

      # Create bar chart showing lipids by number of environments
      ft_bar_data <- ft_overlap_df %>%
        count(N_Environments) %>%
        mutate(N_Environments = factor(N_Environments))

      suppfig_ft_bar <- ggplot(ft_bar_data, aes(x = N_Environments, y = n)) +
        geom_col(fill = "#440154FF", width = 0.6, color = "black", linewidth = 0.3) +
        geom_text(aes(label = n), vjust = -0.3, size = 4, fontface = "bold") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
        labs(title = "A) Lipid Consistency Across Environments",
             subtitle = sprintf("Top 20 SHAP lipids across %d environments", length(ft_available)),
             x = "Number of Environments", y = "Number of Lipids") +
        theme_minimal(base_size = 12) +
        theme(
          axis.text = element_text(color = "black", size = 11),
          axis.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
          panel.grid.major.x = element_blank(),
          axis.line.x = element_line(color = "black", linewidth = 0.5),
          axis.line.y = element_line(color = "black", linewidth = 0.5)
        )

      # Build lipid-wise prevalence panel (requested replacement for old panel D)
      n_ft_env <- length(ft_available)
      ft_mean_shap <- bind_rows(lapply(ft_available, function(env) {
        all_shap_results[[env]] %>%
          dplyr::select(Lipid, MeanAbsSHAP)
      })) %>%
        group_by(Lipid) %>%
        summarise(MeanAbsSHAP_Across = mean(MeanAbsSHAP, na.rm = TRUE), .groups = "drop")

      ft_prevalence <- ft_overlap_df %>%
        left_join(ft_mean_shap, by = "Lipid") %>%
        mutate(MeanAbsSHAP_Across = replace_na(MeanAbsSHAP_Across, 0)) %>%
        filter(N_Environments >= 2) %>%
        arrange(desc(N_Environments), desc(MeanAbsSHAP_Across), Lipid) %>%
        mutate(Lipid = factor(Lipid, levels = rev(unique(Lipid))))

      # Safety fallback: if no lipid appears in >=2 environments, keep all rows.
      if (nrow(ft_prevalence) == 0) {
        ft_prevalence <- ft_overlap_df %>%
          left_join(ft_mean_shap, by = "Lipid") %>%
          mutate(MeanAbsSHAP_Across = replace_na(MeanAbsSHAP_Across, 0)) %>%
          arrange(desc(N_Environments), desc(MeanAbsSHAP_Across), Lipid) %>%
          mutate(Lipid = factor(Lipid, levels = rev(unique(Lipid))))
      }

      suppfig_ft_prevalence <- ggplot(ft_prevalence,
                                      aes(x = N_Environments, y = Lipid, fill = Class)) +
        geom_col(width = 0.75, color = "black", linewidth = 0.25) +
        scale_fill_manual(values = class_colors, na.value = "grey70") +
        scale_x_continuous(
          breaks = seq_len(n_ft_env),
          limits = c(0, n_ft_env + 0.1),
          minor_breaks = NULL
        ) +
        labs(title = "D) Lipid Environment Prevalence",
             subtitle = "Lipids shown if present in >=2 environments (top-20 SHAP lists)",
             x = "Number of Environments", y = "Lipid") +
        theme_minimal(base_size = 11) +
        theme(
          axis.text.y = element_text(color = "black", size = 7),
          axis.text.x = element_text(color = "black", size = 10),
          axis.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
          legend.position = "right",
          panel.grid.major.y = element_blank(),
          axis.line.x = element_line(color = "black", linewidth = 0.5),
          axis.line.y = element_line(color = "black", linewidth = 0.5)
        )

      message("  Lipids in all ", length(ft_available), " environments: ",
              sum(ft_overlap_counts == n_ft_env))

      # ─────────────────────────────────────────────────────────────────────────────
      # S11: Detailed SHAP Heatmap Across FT Environments
      # ─────────────────────────────────────────────────────────────────────────────
      message("\n─── CREATING DETAILED FT ENVIRONMENT HEATMAP ───")

      # Collect SHAP values for all lipids across environments
      ft_shap_combined <- bind_rows(lapply(ft_available, function(env) {
        all_shap_results[[env]] %>%
          mutate(Environment = gsub("FT_", "", env))
      }))

      # Get union of top 25 lipids from each environment
      top_lipids_per_env <- lapply(ft_available, function(env) {
        all_shap_results[[env]] %>% filter(Rank <= 25) %>% pull(Lipid)
      })
      union_top_lipids <- unique(unlist(top_lipids_per_env))

      # Create matrix for heatmap
      ft_heatmap_data <- ft_shap_combined %>%
        filter(Lipid %in% union_top_lipids) %>%
        dplyr::select(Lipid, Environment, MeanAbsSHAP) %>%
        pivot_wider(names_from = Environment, values_from = MeanAbsSHAP, values_fill = 0)

      # Calculate row ordering by mean SHAP across environments
      ft_heatmap_data <- ft_heatmap_data %>%
        mutate(Mean_Across_Envs = rowMeans(dplyr::select(., -Lipid), na.rm = TRUE)) %>%
        arrange(desc(Mean_Across_Envs))

      # Limit to top 40 for readability
      ft_heatmap_data <- ft_heatmap_data %>% slice_head(n = 40)

      # Convert to long format for ggplot
      ft_heatmap_long <- ft_heatmap_data %>%
        dplyr::select(-Mean_Across_Envs) %>%
        pivot_longer(-Lipid, names_to = "Environment", values_to = "SHAP") %>%
        mutate(
          Lipid = factor(Lipid, levels = rev(ft_heatmap_data$Lipid)),
          # Add label suffix: ours=DTA, others=GDD
          Environment = case_when(
            Environment == "ours" ~ "Ours (DTA)",
            TRUE ~ paste0(Environment, " (GDD)")
          ),
          Environment = factor(Environment, levels = c("Ours (DTA)", "MI20 (GDD)", "NE20 (GDD)", "IA21_1 (GDD)", "IA21_2 (GDD)")),
          Class = get_lipid_class(Lipid)
        )

      # Create the heatmap
      suppfig_ft_heatmap <- ggplot(ft_heatmap_long, aes(x = Environment, y = Lipid, fill = SHAP)) +
        geom_tile(color = "white", linewidth = 0.3) +
        scale_fill_viridis_c(option = "plasma", name = "Mean |SHAP|",
                             limits = c(0, max(ft_heatmap_long$SHAP, na.rm = TRUE))) +
        labs(title = "A) Lipid SHAP Values Across Flowering Time Environments",
             subtitle = "Top lipids ranked by mean importance across environments",
             x = "Environment", y = NULL) +
        theme_minimal(base_size = 11) +
        theme(
          axis.text.y = element_text(color = "black", size = 8),
          axis.text.x = element_text(color = "black", size = 11, face = "bold",
                                     angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10),
          legend.position = "right",
          panel.grid = element_blank(),
          axis.line.x = element_line(color = "black", linewidth = 0.5),
          axis.line.y = element_line(color = "black", linewidth = 0.5)
        )

      # Create class annotation bar
      class_annotation <- ft_heatmap_data %>%
        mutate(
          Lipid = factor(Lipid, levels = rev(Lipid)),
          Class = get_lipid_class(Lipid)
        )

      suppfig_ft_class_bar <- ggplot(class_annotation, aes(x = 1, y = Lipid, fill = Class)) +
        geom_tile(color = "white", linewidth = 0.3) +
        scale_fill_manual(values = class_colors, na.value = "grey70") +
        theme_void() +
        theme(
          legend.position = "none",
          axis.text.y = element_blank(),
          plot.margin = margin(0, 0, 0, 0)
        ) +
        coord_cartesian(xlim = c(0.5, 1.5))

      # Create environment R² annotation
      ft_r2_data <- model_df %>%
        filter(PhenotypeID %in% ft_available) %>%
        mutate(
          Environment = gsub("FT_", "", PhenotypeID),
          Environment = case_when(
            Environment == "ours" ~ "Ours (DTA)",
            TRUE ~ paste0(Environment, " (GDD)")
          ),
          Environment = factor(Environment, levels = c("Ours (DTA)", "MI20 (GDD)", "NE20 (GDD)", "IA21_1 (GDD)", "IA21_2 (GDD)"))
        )

      suppfig_ft_r2_bar <- ggplot(ft_r2_data, aes(x = Environment, y = Test_R2, fill = Test_R2)) +
        geom_col(width = 0.7, color = "black", linewidth = 0.3) +
        geom_text(aes(label = format_small_num(Test_R2, digits = 3, threshold = 0.01)),
                  vjust = -0.3, size = 3.5, fontface = "bold") +
        scale_fill_viridis_c(option = "viridis", limits = c(0, max(ft_r2_data$Test_R2, na.rm = TRUE) * 1.2)) +
        scale_y_continuous(
          expand = expansion(mult = c(0, 0.2)),
          labels = function(x) format_small_num(x, digits = 2, threshold = 0.01)
        ) +
        labs(title = "B) Model R² by Environment", x = NULL, y = "Test R²") +
        theme_minimal(base_size = 10) +
        theme(
          legend.position = "none",
          axis.text.x = element_text(color = "black", size = 10, face = "bold",
                                     angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(color = "black"),
          axis.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
          panel.grid.major.x = element_blank(),
          axis.line.x = element_line(color = "black", linewidth = 0.5),
          axis.line.y = element_line(color = "black", linewidth = 0.5)
        )

      # Optional standalone source heatmap (S11) - disabled by default
      suppfig_ft_detailed <- (suppfig_ft_class_bar | suppfig_ft_heatmap | suppfig_ft_r2_bar) +
        plot_layout(widths = c(0.05, 1, 0.45))

      # Final S11 is a single 2x2 figure: A heatmap, B R², C count-by-env, D lipid-wise prevalence
      suppfig_ft_envs <- (suppfig_ft_heatmap + suppfig_ft_r2_bar) /
        (suppfig_ft_bar + suppfig_ft_prevalence) +
        plot_annotation(tag_levels = "A") &
        theme(plot.tag = element_text(size = 14, face = "bold"))

      save_fig(suppfig_ft_envs, "SuppFig_S11_FloweringTime_Environment_Comparison.png",
               w = 14, h = 10, subdir = "supp")
      save_fig(suppfig_ft_envs, "SuppFig_S11_FloweringTime_Environment_Comparison.pdf",
               w = 14, h = 10, subdir = "supp")

      if (isTRUE(SAVE_FT_SOURCE_HEATMAP_FIG)) {
        save_fig(suppfig_ft_detailed, "SuppFig_S11_FloweringTime_SHAP_Heatmap.png",
                 w = 12.5, h = 10, subdir = "supp")
        save_fig(suppfig_ft_detailed, "SuppFig_S11_FloweringTime_SHAP_Heatmap.pdf",
                 w = 12.5, h = 10, subdir = "supp")
        message("  ✓ Saved S11 and standalone FT SHAP heatmap")
      } else {
        message("  ✓ Saved S11 only (standalone FT SHAP heatmap output disabled)")
      }

      # Also create a summary table with SHAP values across environments
      # Get environment columns (exclude Lipid and Mean_Across_Envs)
      env_cols <- setdiff(names(ft_heatmap_data), c("Lipid", "Mean_Across_Envs"))

      ft_shap_table <- ft_heatmap_data %>%
        mutate(
          Class = get_lipid_class(Lipid),
          N_Environments = rowSums(dplyr::select(., all_of(env_cols)) > 0)
        ) %>%
        dplyr::select(Lipid, Class, Mean_Across_Envs, N_Environments, everything()) %>%
        arrange(desc(Mean_Across_Envs))

      write.csv(ft_shap_table, "table/supp/SuppTable_FT_SHAP_Across_Environments.csv", row.names = FALSE)
      message("  ✓ FT SHAP table saved")

      # Print top consistent lipids
      message("\n─── TOP 10 LIPIDS BY MEAN SHAP ACROSS FT ENVIRONMENTS ───")
      print(ft_shap_table %>%
              slice_head(n = 10) %>%
              dplyr::select(Lipid, Class, Mean_Across_Envs, N_Environments))
    }

    # Print model summary
    message("\n─── MODEL PERFORMANCE SUMMARY ───")
    print(as.data.frame(model_df))
  }

} else {
  message("Multi-phenotype analysis skipped (RUN_MULTI_PHENO = FALSE)")
}


# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
#
#     SUPPLEMENTARY FIGURE S10: RF/SHAP ADDITIONAL DETAILS
#
# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════

message("\n══════════════════════════════════════════════════════════════")
message("SUPPLEMENTARY FIGURE S10: RF/SHAP ADDITIONAL DETAILS")
message("══════════════════════════════════════════════════════════════\n")

# Only create if RF was run
if (exists("shap_rank") && exists("rf_final")) {

  # ─────────────────────────────────────────────────────────────────────────────
  # S7A: Phenotype distribution (raw vs residualized)
  # ─────────────────────────────────────────────────────────────────────────────
  if (exists("dat_rf") && exists("pheno_resid")) {
    df_pheno_dist <- tibble(
      Raw = dat_rf$Phenotype,
      `PC-Residualized` = pheno_resid
    ) %>%
      pivot_longer(everything(), names_to = "Type", values_to = "Value")

    suppfig_s5a <- ggplot(df_pheno_dist, aes(x = Value, fill = Type)) +
      geom_density(alpha = 0.6, color = NA) +
      scale_fill_manual(values = c("Raw" = "#440154FF", "PC-Residualized" = "#FDE725FF")) +
      labs(x = "Phenotypic value", y = "Density") +
      theme_minimal(base_size = 12) +
      theme(
        legend.position = "top",
        legend.title = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5)
      )
  } else {
    suppfig_s5a <- ggplot() + theme_void() + ggtitle("A) Phenotype Distribution (not available)")
  }

  # ─────────────────────────────────────────────────────────────────────────────
  # S7B: Class composition of top SHAP lipids
  # ─────────────────────────────────────────────────────────────────────────────
  top50_classes <- shap_rank %>%
    slice_head(n = 50) %>%
    mutate(Class = get_lipid_class(Lipid)) %>%
    count(Class) %>%
    filter(!is.na(Class))

  suppfig_s5b <- ggplot(top50_classes, aes(x = reorder(Class, n), y = n, fill = Class)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.3) +
    coord_flip() +
    scale_fill_manual(values = class_colors, na.value = "grey70") +
    labs(title = "B) Class Composition of Top 50 SHAP Lipids",
         x = NULL, y = "Count") +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0),
      panel.grid.major.y = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5)
    )

  # ─────────────────────────────────────────────────────────────────────────────
  # S7C: Feature correlation with SHAP importance
  # ─────────────────────────────────────────────────────────────────────────────
  # Compute correlation of each lipid with phenotype
  if (exists("lipid_mat_adj") && exists("pheno_resid")) {
    cor_with_pheno <- sapply(1:ncol(lipid_mat_adj), function(i) {
      cor(lipid_mat_adj[, i], pheno_resid, use = "pairwise.complete.obs")
    })

    cor_shap_df <- tibble(
      Lipid = colnames(lipid_mat_adj),
      Correlation = cor_with_pheno
    ) %>%
      inner_join(shap_rank, by = "Lipid") %>%
      mutate(
        Class = get_lipid_class(Lipid),
        IsTop20 = Lipid %in% (shap_rank %>% slice_head(n = 20) %>% pull(Lipid))
      )

    suppfig_s5c <- ggplot(cor_shap_df, aes(x = abs(Correlation), y = MeanAbsSHAP)) +
      geom_point(aes(color = Class, size = IsTop20), alpha = 0.6) +
      geom_smooth(method = "lm", se = FALSE, color = "grey40", linetype = "dashed") +
      scale_color_manual(values = class_colors, na.value = "grey70") +
      scale_size_manual(values = c("TRUE" = 3, "FALSE" = 1.5), guide = "none") +
      labs(title = "C) Correlation vs SHAP Importance",
           x = "|Correlation with Phenotype|", y = "Mean |SHAP|") +
      theme_minimal(base_size = 12) +
      theme(
        legend.position = "right",
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5)
      )
  } else {
    suppfig_s5c <- ggplot() + theme_void() + ggtitle("C) Correlation vs SHAP (not available)")
  }

  # ─────────────────────────────────────────────────────────────────────────────
  # S7D: Cross-validation metrics
  # ─────────────────────────────────────────────────────────────────────────────
  # Run 5-fold CV
  message("Running 5-fold cross-validation...")

  rf_learner <- makeLearner(
    "regr.ranger",
    par.vals = list(
      num.trees = 500,
      mtry = tune_res$recommended.pars$mtry,
      min.node.size = tune_res$recommended.pars$min.node.size,
      sample.fraction = tune_res$recommended.pars$sample.fraction
    ),
    predict.type = "response"
  )

  cv_res <- resample(
    learner = rf_learner,
    task = regr_task,
    resampling = makeResampleDesc("CV", iters = 5L),
    measures = list(rmse, mae, rsq),
    show.info = FALSE
  )

  cv_df <- as.data.frame(cv_res$measures.test) %>%
    mutate(Fold = row_number()) %>%
    pivot_longer(-Fold, names_to = "Metric", values_to = "Value") %>%
    mutate(Metric = recode(Metric,
                           "rmse.test.rmse" = "RMSE",
                           "mae.test.mean" = "MAE",
                           "rsq.test.mean" = "R²"))

  cv_means <- cv_df %>% group_by(Metric) %>% summarise(Mean = mean(Value), .groups = "drop")
  mid_fold <- mean(cv_df$Fold, na.rm = TRUE)

  suppfig_s5d <- ggplot(cv_df, aes(x = Fold, y = Value, color = Metric)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    geom_hline(data = cv_means, aes(yintercept = Mean, color = Metric),
               linetype = "dashed", linewidth = 0.8) +
    geom_text(
      data = cv_means,
      aes(
        x = mid_fold,
        y = Mean,
        label = paste0(tolower(Metric), " mean = ", sprintf("%.3f", Mean)),
        color = Metric
      ),
      angle = 30, hjust = 0.5, vjust = -0.5, size = 3.3
    ) +
    scale_color_manual(values = c("MAE" = "#440154FF", "RMSE" = "#21908CFF", "R²" = "#7AD151FF")) +
    labs(title = "5 Fold CV Metrics by Iteration (RF)",
         x = "CV iteration", y = "Metric value") +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = c(0.95, 0.95),
      legend.justification = c("right", "top"),
      legend.background = element_rect(fill = "white", color = "grey70", linewidth = 0.4),
      legend.title = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5)
    )

  # ─────────────────────────────────────────────────────────────────────────────
  # Combine Supplementary Figure S10
  # ─────────────────────────────────────────────────────────────────────────────
  suppfig_s5_top <- suppfig_s5a + suppfig_s5b + plot_layout(widths = c(1, 1))
  suppfig_s5_bottom <- suppfig_s5c + suppfig_s5d + plot_layout(widths = c(1.2, 1))

  suppfig_s5 <- (suppfig_s5_top / suppfig_s5_bottom) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 15, face = "bold"))

  save_fig(suppfig_s5, paste0("SuppFig_S10_RF_SHAP_details_", PHENOTYPE_NAME, ".png"),
           w = 14, h = 10, subdir = "supp")
  save_fig(suppfig_s5, paste0("SuppFig_S10_RF_SHAP_details_", PHENOTYPE_NAME, ".pdf"),
           w = 14, h = 10, subdir = "supp")

  message("✓ Supplementary Figure S10 (RF/SHAP Details) complete!")

} else {
  message("RF/SHAP not run - skipping Supplementary Figure S10")
}


# ══════════════════════════════════════════════════════════════════════════════
# REASSEMBLE RF/SHAP FIGURES (FIG4, S10, S11)
# ─────────────────────────────────────────────────────────────────────────────

message("\n══════════════════════════════════════════════════════════════")
message("REASSEMBLING RF/SHAP FIGURES (FIG4, S10, S11)")
message("══════════════════════════════════════════════════════════════\n")

load_panel_from_candidates <- function(path_candidates, missing_label) {
  img <- read_image_raster_any(path_candidates)
  if (!is.na(img$path)) {
    message("  Using panel source: ", img$path)
  }
  panel_from_raster(img$raster, fallback_label = missing_label)
}

remove_matching_files <- function(glob_patterns, label) {
  files <- unique(unlist(lapply(glob_patterns, Sys.glob)))
  files <- files[file.exists(files)]
  if (length(files) == 0) return(invisible(0L))
  removed <- file.remove(files)
  message("  Removed ", sum(removed), " ", label, " file(s)")
  invisible(sum(removed))
}

build_stability_panels_from_table <- function(stability_csv, top_n = 20) {
  if (is.na(stability_csv) || !file.exists(stability_csv)) return(NULL)

  stab <- tryCatch(vroom::vroom(stability_csv, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(stab) || nrow(stab) == 0) return(NULL)

  req_cols <- c("Lipid", "SelectionFreqTopK", "SelectionCI_Low", "SelectionCI_High",
                "MedianRank", "RankCI_Low", "RankCI_High")
  if (!all(req_cols %in% names(stab))) return(NULL)

  if (!("Class" %in% names(stab))) {
    stab <- stab %>% mutate(Class = get_lipid_class(Lipid))
  }

  stable_top <- stab %>%
    mutate(
      Class = ifelse(is.na(Class) | Class == "", get_lipid_class(Lipid), as.character(Class))
    ) %>%
    arrange(desc(SelectionFreqTopK), MedianRank, desc(MeanAbsSHAP_Boot)) %>%
    slice_head(n = top_n) %>%
    mutate(Lipid = factor(Lipid, levels = rev(unique(Lipid))))

  if (nrow(stable_top) == 0) return(NULL)

  p_stab_b <- ggplot(stable_top, aes(x = Lipid, y = SelectionFreqTopK, fill = Class)) +
    geom_col(width = 0.72, color = "black", linewidth = 0.2) +
    geom_errorbar(aes(ymin = SelectionCI_Low, ymax = SelectionCI_High), width = 0.2, linewidth = 0.35) +
    coord_flip() +
    scale_fill_manual(values = class_colors, na.value = "grey70") +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
    labs(
      title = paste0("B) SHAP Stability: Top-", top_n, " Selection Frequency"),
      x = NULL, y = "Bootstrap Selection Frequency"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text = element_text(color = "black", size = 8.5),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0),
      legend.position = "none",
      panel.grid.major.y = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5)
    )

  p_stab_c <- ggplot(stable_top, aes(y = Lipid, x = MedianRank, color = Class)) +
    geom_segment(aes(x = RankCI_Low, xend = RankCI_High, yend = Lipid), linewidth = 0.8, alpha = 0.8) +
    geom_point(size = 2) +
    scale_color_manual(values = class_colors, na.value = "grey70") +
    labs(
      title = "C) Rank Stability (Median and 95% CI)",
      x = "Rank (lower is more stable)", y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text = element_text(color = "black", size = 8.5),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0),
      legend.position = "none",
      panel.grid.major.y = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5)
    )

  list(panel_b = p_stab_b, panel_c = p_stab_c)
}

build_top_shap_panel_from_table <- function(shap_csv, top_n = 20) {
  if (is.na(shap_csv) || !file.exists(shap_csv)) return(NULL)

  shap_df <- tryCatch(vroom::vroom(shap_csv, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(shap_df) || nrow(shap_df) == 0) return(NULL)
  if (!all(c("Lipid", "MeanAbsSHAP") %in% names(shap_df))) return(NULL)

  top_df <- shap_df %>%
    mutate(Class = get_lipid_class(Lipid)) %>%
    arrange(desc(MeanAbsSHAP)) %>%
    slice_head(n = top_n)

  if (nrow(top_df) == 0) return(NULL)

  ggplot(top_df, aes(x = reorder(Lipid, MeanAbsSHAP), y = MeanAbsSHAP, fill = Class)) +
    geom_col(width = 0.72, color = "black", linewidth = 0.2) +
    coord_flip() +
    scale_fill_manual(values = class_colors, na.value = "grey70") +
    labs(
      title = paste0("B) Top ", top_n, " SHAP Lipids (Plant Height)"),
      x = NULL,
      y = "Mean |SHAP|"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text = element_text(color = "black", size = 8.5),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0),
      legend.position = "right",
      panel.grid.major.y = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5)
    )
}

# New function: Build SHAP panel with error bars from stability data
build_top_shap_stability_panel <- function(stability_csv, top_n = 20) {
  if (is.na(stability_csv) || !file.exists(stability_csv)) return(NULL)

  stab_df <- tryCatch(vroom::vroom(stability_csv, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(stab_df) || nrow(stab_df) == 0) return(NULL)

  # Check for required columns
  required_cols <- c("Lipid", "MeanAbsSHAP_Boot", "SHAP_CI_Low", "SHAP_CI_High")
  if (!all(required_cols %in% names(stab_df))) return(NULL)

  # Get Class column if exists, otherwise derive it
  if (!"Class" %in% names(stab_df)) {
    stab_df <- stab_df %>% mutate(Class = get_lipid_class(Lipid))
  }

  top_df <- stab_df %>%
    arrange(desc(MeanAbsSHAP_Boot)) %>%
    slice_head(n = top_n) %>%
    mutate(Lipid = factor(Lipid, levels = rev(Lipid)))  # For coord_flip ordering

  if (nrow(top_df) == 0) return(NULL)

  ggplot(top_df, aes(x = Lipid, y = MeanAbsSHAP_Boot, fill = Class)) +
    geom_col(width = 0.72, color = "black", linewidth = 0.2) +
    geom_errorbar(aes(ymin = SHAP_CI_Low, ymax = SHAP_CI_High),
                  width = 0.25, linewidth = 0.4, color = "black") +
    coord_flip() +
    scale_fill_manual(values = class_colors, na.value = "grey70") +
    labs(
      title = paste0("B) Top ", top_n, " SHAP Lipids (Plant Height)"),
      subtitle = "Bootstrap stability (95% CI)",
      x = NULL,
      y = "Mean |SHAP|"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text = element_text(color = "black", size = 8.5),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 9, color = "grey40"),
      legend.position = "right",
      panel.grid.major.y = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5)
    )
}

build_obs_pred_panel_from_objects <- function() {
  if (!exists("df_pred") || !is.data.frame(df_pred)) return(NULL)
  if (!all(c("Observed", "Predicted") %in% names(df_pred))) return(NULL)
  if (nrow(df_pred) == 0) return(NULL)

  txt_bits <- c()
  if (exists("rmse_test") && is.finite(rmse_test)) {
    txt_bits <- c(txt_bits, paste0("RMSE = ", sprintf("%.3f", rmse_test)))
  }
  if (exists("r2_test") && is.finite(r2_test)) {
    txt_bits <- c(txt_bits, paste0("R² = ", sprintf("%.3f", r2_test)))
  }
  metric_txt <- if (length(txt_bits) > 0) paste(txt_bits, collapse = "\n") else NULL

  p <- ggplot(df_pred, aes(x = Observed, y = Predicted)) +
    geom_point(alpha = 0.7, size = 1.8, color = "#555555") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60") +
    labs(
      title = "Observed vs Predicted Plant Height",
      x = "Observed residual phenotype",
      y = "RF predicted residual phenotype"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5)
    )

  if (!is.null(metric_txt)) {
    p <- p + annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.2,
                      label = metric_txt, size = 3.3, fontface = "bold")
  }
  p
}

# ─────────────────────────────────────────────────────────────────────────────
# Figure 4 (main): Plant Height prediction + SHAP + external reproducibility
# ─────────────────────────────────────────────────────────────────────────────

fig4_panel_a <- build_obs_pred_panel_from_objects()
if (is.null(fig4_panel_a)) {
  fig4_panel_a <- load_panel_from_candidates(
    c(
      "fig/main/individual_figs/Fig4d_obs_vs_pred_PlantHeight_rra*.png",
      "fig/main/All_figures/Fig4d_obs_vs_pred_PlantHeight_rra*.png"
    ),
    "Figure 4A (Observed vs Predicted) not found"
  )
}

# Try stability CSV first (has error bars), fallback to regular SHAP CSV
fig4_stability_csv <- resolve_latest_existing_file(
  c("table/supp/SuppTable_SHAP_Stability_PlantHeight.csv")
)
fig4_shap_csv <- resolve_latest_existing_file(
  c("table/supp/SuppTable_lipid_SHAP_PlantHeight.csv")
)

# Prefer stability version with error bars
fig4_panel_b <- build_top_shap_stability_panel(fig4_stability_csv, top_n = 20)
if (is.null(fig4_panel_b)) {
  # Fallback to regular SHAP without error bars
  fig4_panel_b <- build_top_shap_panel_from_table(fig4_shap_csv, top_n = 20)
}
if (is.null(fig4_panel_b)) {
  fig4_panel_b <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "Figure 4B (Top SHAP lipids) not available", size = 4.2, hjust = 0.5) +
    theme_void() +
    theme(panel.border = element_rect(color = "grey70", fill = NA, linewidth = 1))
}

ph_comp_csv <- resolve_latest_existing_file(c("table/supp/SuppTable_PlantHeight_comparison.csv"))
if (!is.na(ph_comp_csv) && file.exists(ph_comp_csv)) {
  message("  Figure 4 C/D source: ", ph_comp_csv)

  ph_comp <- tryCatch(
    vroom::vroom(ph_comp_csv, delim = ",", show_col_types = FALSE),
    error = function(e) NULL
  )

  if (!is.null(ph_comp) &&
      all(c("PH_OurData", "PH_OurData_SHAP", "PH_Boyles", "PH_Boyles_SHAP") %in% names(ph_comp))) {

    ours_tbl <- ph_comp %>%
      transmute(
        Lipid = as.character(PH_OurData),
        SHAP_Ours = as.numeric(PH_OurData_SHAP)
      ) %>%
      filter(!is.na(Lipid), Lipid != "") %>%
      distinct(Lipid, .keep_all = TRUE)

    boyles_tbl <- ph_comp %>%
      transmute(
        Lipid = as.character(PH_Boyles),
        SHAP_Boyles = as.numeric(PH_Boyles_SHAP)
      ) %>%
      filter(!is.na(Lipid), Lipid != "") %>%
      distinct(Lipid, .keep_all = TRUE)

    ph_ours_lipids <- ours_tbl$Lipid
    ph_boyles_lipids <- boyles_tbl$Lipid
    overlap_ph <- intersect(ph_ours_lipids, ph_boyles_lipids)
    only_ours <- setdiff(ph_ours_lipids, ph_boyles_lipids)
    only_boyles <- setdiff(ph_boyles_lipids, ph_ours_lipids)

    ph_venn_data <- tibble(
      Category = c("Our Data Only", "Both", "Boyles Only"),
      Count = c(length(only_ours), length(overlap_ph), length(only_boyles)),
      Type = c("Unique", "Overlap", "Unique")
    ) %>%
      mutate(Category = factor(Category, levels = c("Our Data Only", "Both", "Boyles Only")))

    fig4_panel_c <- ggplot(ph_venn_data, aes(x = Category, y = Count, fill = Type)) +
      geom_col(width = 0.6, color = "black", linewidth = 0.3) +
      geom_text(aes(label = Count), vjust = -0.3, size = 4.5, fontface = "bold") +
      scale_fill_manual(values = c("Overlap" = "#21908CFF", "Unique" = "#FDE725FF")) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(title = "C) Top 20 SHAP Lipids: Overlap",
           subtitle = "Plant Height (Our Data vs Boyles et al.)",
           x = NULL, y = "Number of Lipids") +
      theme_minimal(base_size = 11) +
      theme(
        legend.position = "none",
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0),
        plot.subtitle = element_text(hjust = 0, color = "grey40"),
        panel.grid.major.x = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5)
      )

    ph_scatter_data <- tibble(Lipid = union(ph_ours_lipids, ph_boyles_lipids)) %>%
      left_join(ours_tbl, by = "Lipid") %>%
      left_join(boyles_tbl, by = "Lipid") %>%
      mutate(
        SHAP_Ours = replace_na(SHAP_Ours, 0),
        SHAP_Boyles = replace_na(SHAP_Boyles, 0),
        InBoth = Lipid %in% overlap_ph,
        Class = get_lipid_class(Lipid)
      )

    overlap_r <- suppressWarnings(
      cor(
        ph_scatter_data$SHAP_Ours[ph_scatter_data$InBoth],
        ph_scatter_data$SHAP_Boyles[ph_scatter_data$InBoth],
        use = "complete.obs"
      )
    )
    if (!is.finite(overlap_r)) overlap_r <- NA_real_

    fig4_panel_d <- ggplot(ph_scatter_data, aes(x = SHAP_Ours, y = SHAP_Boyles)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
      geom_point(aes(color = Class, shape = InBoth), size = 2.8, alpha = 0.75) +
      scale_color_manual(values = class_colors, na.value = "grey70") +
      scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 17),
                         labels = c("TRUE" = "In both top 20", "FALSE" = "Unique")) +
      labs(
        title = "D) SHAP Importance Correlation",
        subtitle = ifelse(is.na(overlap_r),
                          "Top 20 lipids (overlap correlation not available)",
                          paste0("Top 20 lipids (r = ", sprintf("%.2f", overlap_r), " for overlap)")),
        x = "Mean |SHAP| - Our Data",
        y = "Mean |SHAP| - Boyles et al."
      ) +
      theme_minimal(base_size = 11) +
      theme(
        legend.position = "right",
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0),
        plot.subtitle = element_text(hjust = 0, color = "grey40"),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5)
      )
  } else {
    fig4_panel_c <- panel_from_raster(NULL, "Figure 4C (PH overlap) not available")
    fig4_panel_d <- panel_from_raster(NULL, "Figure 4D (PH SHAP correlation) not available")
  }
} else {
  fig4_panel_c <- panel_from_raster(NULL, "Figure 4C (PH overlap) not available")
  fig4_panel_d <- panel_from_raster(NULL, "Figure 4D (PH SHAP correlation) not available")
}

figure4_rf_shap <- (fig4_panel_a + fig4_panel_b) /
  (fig4_panel_c + fig4_panel_d) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 15, face = "bold"))

save_fig(figure4_rf_shap, "Figure4_RF_SHAP_PlantHeight.png", w = 14, h = 12, subdir = "main")
save_fig(figure4_rf_shap, "Figure4_RF_SHAP_PlantHeight.pdf", w = 14, h = 12, subdir = "main")
message("  ✓ Figure 4 remapped to main-story Plant Height panels (A-D)")

# Remove standalone S9b per request
remove_matching_files(
  c(
    "fig/supp/SuppFig_S9b_SHAP_Stability_PlantHeight*.png",
    "fig/supp/SuppFig_S9b_SHAP_Stability_PlantHeight*.pdf",
    "fig/supp/SuppFig_S9b_SHAP_Stability*.png",
    "fig/supp/SuppFig_S9b_SHAP_Stability*.pdf",
    "fig/supp/all_figures/SuppFig_S9b_SHAP_Stability_PlantHeight*.png",
    "fig/supp/all_figures/SuppFig_S9b_SHAP_Stability_PlantHeight*.pdf"
  ),
  "standalone S9b"
)

# Keep S10 as diagnostics-only (A-D) generated above; do not merge with S9a in script 26
message("  ✓ Keeping SuppFig S10 as diagnostics-only (no S10 + S9a merge)")

# ─────────────────────────────────────────────────────────────────────────────
# Supp Fig S11: built directly in multi-phenotype section (A-D)
# ─────────────────────────────────────────────────────────────────────────────

message("  ✓ SuppFig S11 is generated directly (A-D); skipping legacy S11-based remap")

# Remove no-longer-needed standalone figures
remove_matching_files(
  c(
    "fig/supp/SuppFig_S9a_PlantHeight_SHAP_Comparison*.png",
    "fig/supp/SuppFig_S9a_PlantHeight_SHAP_Comparison*.pdf",
    "fig/supp/all_figures/SuppFig_S9a_PlantHeight_SHAP_Comparison*.png",
    "fig/supp/all_figures/SuppFig_S9a_PlantHeight_SHAP_Comparison*.pdf"
  ),
  "standalone S9a"
)
remove_matching_files(
  c(
    "fig/supp/SuppFig_S11_FloweringTime_SHAP_Heatmap*.png",
    "fig/supp/SuppFig_S11_FloweringTime_SHAP_Heatmap*.pdf",
    "fig/supp/all_figures/SuppFig_S11_FloweringTime_SHAP_Heatmap*.png",
    "fig/supp/all_figures/SuppFig_S11_FloweringTime_SHAP_Heatmap*.pdf"
  ),
  "standalone S11"
)

# Remove standalone S10b per request
remove_matching_files(
  c(
    "fig/supp/SuppFig_S10b_FloweringTime_SHAP_Stability*.png",
    "fig/supp/SuppFig_S10b_FloweringTime_SHAP_Stability*.pdf",
    "fig/supp/all_figures/SuppFig_S10b_FloweringTime_SHAP_Stability*.png",
    "fig/supp/all_figures/SuppFig_S10b_FloweringTime_SHAP_Stability*.pdf"
  ),
  "standalone S10b"
)


# ══════════════════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ══════════════════════════════════════════════════════════════════════════════

message("\n══════════════════════════════════════════════════════════════")
message("FIGURE GENERATION COMPLETE")
message("══════════════════════════════════════════════════════════════\n")

message("MAIN FIGURES:")
message("  Fig 1: Lipidomics Landscape (Pipeline image, CLR, Chemical shifts, LION Enrichment)")
message("  Fig 3: OPLS-DA Analysis")
message("  Fig 4: RF + SHAP Predictive Modeling (Plant Height)")
message("")
message("SUPPLEMENTARY FIGURES:")
message("  S1: QC - Run Order, SERRF RSD, PCA, and SpATS residuals")
message("  S2: %TIC and ALR Compositional Context Contrasts")
message("  S3: Lipid Species Counts & Overview")
message("  S4: PCA of Lipids (Control + LowInput)")
message("  S5: Top Variance Lipids (Genotype Variation)")
message("  S8: Partial Correlation Heatmaps (controlling for total TIC)")
message("  S10: RF/SHAP diagnostics for Plant Height (A-D)")
message("  S11: Flowering Time consistency across 5 environments (A-D)")
message("")
message("OUTPUT LOCATIONS:")
message("  Main figures: fig/main/")
message("  Supplementary: fig/supp/")
message("  Tables: table/supp/")
message("")
message("STILL TO DO:")
message("  - GWAS figure (Fig 2): finalize panels and supplementary figures")
message("  - Fig 6: Summary model (conceptual)")
message("")
message("══════════════════════════════════════════════════════════════\n")
