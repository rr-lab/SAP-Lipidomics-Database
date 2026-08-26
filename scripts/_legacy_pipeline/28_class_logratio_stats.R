# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4: RATIO-BASED MULTIVARIATE ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════
#
# This script generates ONLY the figures and tables for Section 4:
# "Ratio-based multivariate analysis reveals coordinated lipid class transitions"
#
# OUTPUTS:
#   - Figure 2 (Main): OPLS-DA Analysis (4 panels: A-D)
#   - Supplementary Figure S6: Partial Correlations
#   - Supplementary Figure S7: Lipid Ratios
#   - Supplementary Table S7: OPLS-DA VIP ratios
#
# ══════════════════════════════════════════════════════════════════════════════

message("\n")
message("══════════════════════════════════════════════════════════════")
message("SECTION 4: RATIO-BASED MULTIVARIATE ANALYSIS")
message("══════════════════════════════════════════════════════════════")
message("\n")

# ──────────────────────────────────────────────────────────────────────────────
# 0) SETUP: PACKAGES, THEME, PALETTES
# ──────────────────────────────────────────────────────────────────────────────

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
  library(patchwork)
  library(viridis)
  library(scales)
  library(ppcor)        # For partial correlation analysis
  library(ropls)        # For OPLS-DA
})

# Create output directories
dir.create("fig/main", recursive = TRUE, showWarnings = FALSE)
dir.create("fig/supp", recursive = TRUE, showWarnings = FALSE)
dir.create("table/supp", recursive = TRUE, showWarnings = FALSE)

# ──────────────────────────────────────────────────────────────────────────────
# Plot Theme (Publication Quality)
# ──────────────────────────────────────────────────────────────────────────────

plot_theme <- theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
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

# ──────────────────────────────────────────────────────────────────────────────
# Color Palettes
# ──────────────────────────────────────────────────────────────────────────────

class_colors <- c(
  PC   = "#00441B", PA   = "#1B7837", PE   = "#41AB5D", LPC  = "#66C2A4", LPE  = "#2CA25F",
  PG   = "#78C679", PS   = "#C2E699", DG   = "#54278F", DGDG = "#F768A1", MG   = "#8941ED",
  MGDG = "#FBB4D9", SQDG = "#9D4D6C", TG   = "#ED804A"
)

condition_colors <- c(Control = "#440154FF", LowInput = "#FDE725FF")

valid_classes <- names(class_colors)
classes_with_lyso <- unique(c(valid_classes, "LPC", "LPE"))

# ──────────────────────────────────────────────────────────────────────────────
# Helper Functions
# ──────────────────────────────────────────────────────────────────────────────

get_lipid_class <- function(x) {
  cls <- str_extract(x, "^[A-Za-z0-9]+(?=\\()")
  cls <- ifelse(is.na(cls), NA_character_, cls)
  ifelse(cls %in% valid_classes, cls, "Other")
}

# Save function with global axis scaling
save_fig <- function(p, filename, w = 10, h = 8, dpi = 300, subdir = "main") {
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

  path <- file.path("fig", subdir, filename)
  ggsave(path, plot = p, width = w, height = h, units = "in", dpi = dpi, bg = "white")
  message("Saved: ", path)
}

# Restrict to genotypes present in both Control and LowInput
restrict_to_condition_pairs <- function(df) {
  shared_ids <- df %>%
    dplyr::distinct(Compound_Name, Condition) %>%
    dplyr::count(Compound_Name, name = "n_cond") %>%
    dplyr::filter(n_cond == 2) %>%
    dplyr::pull(Compound_Name)

  df %>% dplyr::filter(Compound_Name %in% shared_ids)
}

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

# ══════════════════════════════════════════════════════════════════════════════
# 1) LOAD DATA
# ══════════════════════════════════════════════════════════════════════════════

message("Loading data...")

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
#                    FIGURE 2 (MAIN): OPLS-DA ANALYSIS
#
# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════

message("\n══════════════════════════════════════════════════════════════")
message("FIGURE 2 (MAIN): OPLS-DA ANALYSIS")
message("══════════════════════════════════════════════════════════════\n")

message("Running OPLS-DA...")

# Prepare data for OPLS-DA (using class ratios)
long_prep_global <- function(df) {
  class_pat <- paste0("^(", paste(classes_with_lyso, collapse = "|"), ")(?=\\()")

  df %>%
    pivot_longer(-c(Compound_Name, Condition), names_to = "Lipid", values_to = "Intensity") %>%
    dplyr::rename(Sample = Compound_Name) %>%
    mutate(Class = str_extract(Lipid, class_pat)) %>%
    filter(!is.na(Class)) %>%
    group_by(Sample, Condition, Class) %>%
    summarise(class_total = sum(Intensity, na.rm = TRUE), .groups = "drop") %>%
    group_by(Sample) %>%
    mutate(
      class_total = pmax(class_total, 0),
      minpos = ifelse(
        all(is.na(class_total) | class_total <= 0),
        NA_real_,
        min(class_total[class_total > 0], na.rm = TRUE)
      ),
      eps = ifelse(is.na(minpos), 0, minpos * 0.5),
      class_log = log10(class_total + eps)
    ) %>%
    ungroup() %>%
    dplyr::select(Sample, Condition, Class, class_log)
}

# Combine and process
combined <- bind_rows(
  control_raw %>% dplyr::mutate(Condition = "Control"),
  low_raw %>% dplyr::mutate(Condition = "LowInput")
)

combined_paired <- restrict_to_condition_pairs(combined)
message("Matched genotype intersection for OPLS-DA: ", nrow(combined_paired) / 2, " pairs")

long_all <- long_prep_global(combined_paired) %>%
  filter(!is.na(class_log), nchar(Sample) >= 5)

# Pivot to wide and create ratios
wide_log <- long_all %>%
  pivot_wider(names_from = Class, values_from = class_log, values_fill = NA_real_)

# Create all pairwise ratios (log-ratio = subtraction in log space)
long_classes <- wide_log %>%
  pivot_longer(cols = -c(Sample, Condition), names_to = "Class", values_to = "Abundance")

ratios_long <- long_classes %>%
  inner_join(
    long_classes,
    by = c("Sample", "Condition"),
    suffix = c(".num", ".den"),
    relationship = "many-to-many"
  ) %>%
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
# Fig 2A: OPLS-DA scores plot
# ─────────────────────────────────────────────────────────────────────────────
fig2a <- ggplot(scores_df, aes(x = t1, y = to1, color = Condition)) +
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
# Fig 2B: VIP scores
# ─────────────────────────────────────────────────────────────────────────────
vip_vec <- slot(opls_mod, "vipVn")
vip_df <- tibble(Ratio = names(vip_vec), VIP = as.numeric(vip_vec)) %>%
  filter(VIP > 1) %>%
  arrange(desc(VIP)) %>%
  slice_head(n = 25)

fig2b <- ggplot(vip_df, aes(x = reorder(Ratio, VIP), y = VIP)) +
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
# Fig 2C: Permutation test
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
  pivot_longer(everything(), names_to = "metric_key", values_to = "value") %>%
  mutate(
    metric_key = recode(metric_key, R2Ycum = "R2Y", Q2cum = "Q2"),
    metric_label = recode(metric_key, R2Y = "R\u00b2Y", Q2 = "Q\u00b2")
  )

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

fig2c <- ggplot(perm_plot_df, aes(x = value, fill = metric_key)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity", color = "white") +
  geom_vline(xintercept = obsR2Y, linetype = "dashed", color = "#440154FF", linewidth = 1) +
  geom_vline(xintercept = obsQ2, linetype = "dashed", color = "#FDE725FF", linewidth = 1) +
  annotate("text", x = obsR2Y - 0.01, y = Inf,
           label = sprintf("R\u00b2Y = %.3f\np = %.3f", obsR2Y, pR2Y),
           vjust = 1.5, hjust = 1, size = 3.5, color = "#440154FF", fontface = "bold") +
  annotate("text", x = obsQ2 - 0.01, y = Inf,
           label = sprintf("Q\u00b2 = %.3f\np = %.3f", obsQ2, pQ2),
           vjust = 3.5, hjust = 1, size = 3.5, color = "#FDE725FF", fontface = "bold") +
  scale_fill_manual(
    values = c("R2Y" = "#440154FF", "Q2" = "#FDE725FF"),
    labels = c("R2Y" = "R\u00b2Y", "Q2" = "Q\u00b2")
  ) +
  # scale_x_break() disabled: ggbreak 0.1.2 is incompatible with the S7 object
  # model in ggplot2 4.x (as.grob has no method for ScaleContinuousPosition).
  # The axis break is cosmetic and Figure2_OPLS_DA is not used in main.tex.
  # scale_x_break(c(break_start, break_end), scales = 0.5) +
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
# Fig 2D: Model summary metrics
# ─────────────────────────────────────────────────────────────────────────────
model_metrics <- tibble(
  Metric = c("R\u00b2X (cumulative)", "R\u00b2Y (cumulative)", "Q\u00b2 (cumulative)",
             "Predictive components", "Orthogonal components"),
  Value = c(
    round(opls_mod@summaryDF$`R2X(cum)`, 3),
    round(opls_mod@summaryDF$`R2Y(cum)`, 3),
    round(opls_mod@summaryDF$`Q2(cum)`, 3),
    opls_mod@summaryDF$pre,
    opls_mod@summaryDF$ort
  )
)

fig2d <- ggplot(model_metrics, aes(x = Metric, y = as.numeric(Value))) +
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

# ─────────────────────────────────────────────────────────────────────────────
# Combine Figure 2 panels
# ─────────────────────────────────────────────────────────────────────────────
fig2_top <- fig2a + fig2b + plot_layout(widths = c(1, 1.2))

# ggbreak plots need special handling with aplot for combining
if (!requireNamespace("aplot", quietly = TRUE)) {
  install.packages("aplot")
}
library(aplot)

# The VIP table does not depend on the figure, so write it first. This keeps
# the table reproducible even if the ggbreak/aplot composition fails.
write.csv(vip_df, "table/supp/SuppTable_S7_OPLS_VIP_ratios.csv", row.names = FALSE)
message("Saved: table/supp/SuppTable_S7_OPLS_VIP_ratios.csv")

# ggbreak's scale_x_break is incompatible with the S7-based ggplot2 4.x object
# model (as.grob has no method for ScaleContinuousPosition). Figure2_OPLS_DA is
# not used in main.tex, so a failure here must not stop the script before the
# supplementary figures below are written.
tryCatch({
  fig2_bottom <- fig2c + fig2d + plot_layout(widths = c(1.2, 1))
  fig2_combined <- fig2_top / fig2_bottom +
    plot_layout(heights = c(1, 0.8))
  save_fig(fig2_combined, "Figure2_OPLS_DA.png", w = 14, h = 12)
  save_fig(fig2_combined, "Figure2_OPLS_DA.pdf", w = 14, h = 12)
}, error = function(e) {
  message("SKIPPED Figure2_OPLS_DA (ggbreak/ggplot2 incompatibility): ",
          conditionMessage(e))
  message("  -> not used in main.tex; continuing to the supplementary figures.")
})

message("\n\u2713 Figure 2 (OPLS-DA) complete!\n")


# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
#
#     SUPPLEMENTARY FIGURE S6: CLR-SPACE CLASS CORRELATIONS
#
# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════

message("\n══════════════════════════════════════════════════════════════")
message("SUPPLEMENTARY FIGURE S6: CLR-SPACE CLASS CORRELATIONS")
message("══════════════════════════════════════════════════════════════\n")

# NOTE: We use CLR-transformed class abundances instead of partial correlations
# on %TIC. Partial correlations controlling for TIC while using %TIC data is
# problematic due to compositional closure (sum = 100%). CLR transformation
# properly handles compositional data and avoids spurious correlations.

# ──────────────────────────────────────────────────────────────────────────────
# Function to compute class-level CLR and correlations
# ──────────────────────────────────────────────────────────────────────────────
compute_clr_correlations <- function(df_raw, condition_name, valid_classes) {

 lipid_cols <- setdiff(names(df_raw), c("Compound_Name", "Condition"))
 lipid_cols <- lipid_cols[grepl("\\(", lipid_cols)]

 # Get raw intensities
 X <- df_raw %>%
   dplyr::select(all_of(lipid_cols)) %>%
   mutate(across(everything(), ~suppressWarnings(as.numeric(.x)))) %>%
   as.matrix()

 # Assign classes from lipid names
 classes <- str_extract(lipid_cols, "^[A-Za-z0-9]+(?=\\()")

 # Calculate class-level sums (raw intensities)
 unique_classes <- sort(unique(classes[!is.na(classes)]))
 unique_classes <- intersect(unique_classes, valid_classes)

 class_mat <- sapply(unique_classes, function(cl) {
   idx <- which(classes == cl)
   if (length(idx) > 0) rowSums(X[, idx, drop = FALSE], na.rm = TRUE) else rep(0, nrow(X))
 })
 colnames(class_mat) <- unique_classes

 # Convert to compositions (proportions that sum to 1)
 total <- rowSums(class_mat, na.rm = TRUE)
 total[total == 0] <- NA
 class_comp <- sweep(class_mat, 1, total, "/")

 # Zero replacement for CLR (multiplicative replacement)
 delta_frac <- 0.5
 for (i in seq_len(nrow(class_comp))) {
   x <- class_comp[i, ]
   if (any(x == 0, na.rm = TRUE)) {
     pos_vals <- x[x > 0 & !is.na(x)]
     if (length(pos_vals) > 0) {
       delta <- min(pos_vals) * delta_frac
       x[x == 0] <- delta
       class_comp[i, ] <- x / sum(x, na.rm = TRUE)  # Re-close
     }
   }
 }

 # Remove samples with NA
 complete_idx <- complete.cases(class_comp)
 class_comp <- class_comp[complete_idx, , drop = FALSE]

 if (nrow(class_comp) < 10) {
   message("  Not enough complete samples for CLR correlation")
   return(NULL)
 }

 # CLR transformation: log(x_i) - mean(log(x))
 log_comp <- log(class_comp)
 clr_mat <- log_comp - rowMeans(log_comp, na.rm = TRUE)

 # Compute Pearson correlation on CLR-transformed data
 cor_mat <- cor(clr_mat, use = "pairwise.complete.obs", method = "pearson")

 # Compute p-values using cor.test
 n_classes <- ncol(clr_mat)
 pval_mat <- matrix(NA, n_classes, n_classes)
 rownames(pval_mat) <- colnames(pval_mat) <- unique_classes

 for (i in 1:(n_classes - 1)) {
   for (j in (i + 1):n_classes) {
     tryCatch({
       ct <- cor.test(clr_mat[, i], clr_mat[, j], method = "pearson")
       pval_mat[i, j] <- ct$p.value
       pval_mat[j, i] <- ct$p.value
     }, error = function(e) {
       pval_mat[i, j] <<- NA
       pval_mat[j, i] <<- NA
     })
   }
 }
 diag(pval_mat) <- 0

 list(
   cor_mat = cor_mat,
   pval_mat = pval_mat,
   clr_mat = clr_mat,
   condition = condition_name,
   n_samples = nrow(clr_mat),
   classes = unique_classes
 )
}

# ──────────────────────────────────────────────────────────────────────────────
# Compute CLR correlations for both conditions
# ──────────────────────────────────────────────────────────────────────────────
message("Computing CLR-space correlations (compositional-safe)...")

clr_control <- compute_clr_correlations(control_raw, "Control", classes_with_lyso)
clr_lowinput <- compute_clr_correlations(low_raw, "LowInput", classes_with_lyso)

# Align classes between conditions
if (!is.null(clr_control) && !is.null(clr_lowinput)) {
 common_classes <- intersect(clr_control$classes, clr_lowinput$classes)
 # Subset to common classes
 clr_control$cor_mat <- clr_control$cor_mat[common_classes, common_classes]
 clr_control$pval_mat <- clr_control$pval_mat[common_classes, common_classes]
 clr_lowinput$cor_mat <- clr_lowinput$cor_mat[common_classes, common_classes]
 clr_lowinput$pval_mat <- clr_lowinput$pval_mat[common_classes, common_classes]
} else {
 common_classes <- NULL
}

# Compute delta (change in correlation)
if (!is.null(clr_control) && !is.null(clr_lowinput)) {
 cor_delta <- clr_lowinput$cor_mat - clr_control$cor_mat
} else {
 cor_delta <- NULL
}

# ──────────────────────────────────────────────────────────────────────────────
# Identify sign flips (significant in both, opposite signs)
# ──────────────────────────────────────────────────────────────────────────────
get_signflip_pair_keys_clr <- function(clr_ctrl, clr_low, classes, alpha = 0.05) {
 if (is.null(clr_ctrl) || is.null(clr_low) || length(classes) == 0) {
   return(character(0))
 }

 pairs <- expand.grid(Class1 = classes, Class2 = classes, stringsAsFactors = FALSE) %>%
   filter(Class1 != Class2) %>%
   mutate(
     pair_key = ifelse(Class1 < Class2,
                       paste(Class1, Class2, sep = "||"),
                       paste(Class2, Class1, sep = "||")),
     r_ctrl = pmap_dbl(list(Class1, Class2), ~ clr_ctrl$cor_mat[..1, ..2]),
     p_ctrl = pmap_dbl(list(Class1, Class2), ~ clr_ctrl$pval_mat[..1, ..2]),
     r_low  = pmap_dbl(list(Class1, Class2), ~ clr_low$cor_mat[..1, ..2]),
     p_low  = pmap_dbl(list(Class1, Class2), ~ clr_low$pval_mat[..1, ..2]),
     sig_ctrl = !is.na(p_ctrl) & p_ctrl < alpha,
     sig_low  = !is.na(p_low) & p_low < alpha,
     flip_sig = sig_ctrl & sig_low & !is.na(r_ctrl) & !is.na(r_low) & (r_ctrl * r_low < 0)
   )

 unique(pairs$pair_key[pairs$flip_sig])
}

flip_pair_keys <- get_signflip_pair_keys_clr(clr_control, clr_lowinput, common_classes, alpha = 0.05)

# ──────────────────────────────────────────────────────────────────────────────
# Function to create CLR correlation heatmap
# ──────────────────────────────────────────────────────────────────────────────
create_clr_correlation_plot <- function(clr_result, title_prefix, class_order = NULL,
                                        flip_pair_keys = character(0), alpha = 0.05,
                                        is_delta = FALSE) {

 if (is.null(clr_result) && !is_delta) {
   return(ggplot() + theme_void() + ggtitle(paste(title_prefix, "- Not available")))
 }

 if (is_delta) {
   # For delta plot, clr_result is the delta matrix directly
   cor_mat <- clr_result
   pval_mat <- NULL
   n_samples <- NA
   classes <- rownames(cor_mat)
 } else {
   cor_mat <- clr_result$cor_mat
   pval_mat <- clr_result$pval_mat
   n_samples <- clr_result$n_samples
   classes <- rownames(cor_mat)
 }

 if (!is.null(class_order)) {
   classes <- intersect(class_order, classes)
   cor_mat <- cor_mat[classes, classes, drop = FALSE]
   if (!is.null(pval_mat)) pval_mat <- pval_mat[classes, classes, drop = FALSE]
 }

 # Build data frame for plotting
 plot_data <- expand.grid(Class1 = classes, Class2 = classes, stringsAsFactors = FALSE) %>%
   mutate(idx1 = match(Class1, classes), idx2 = match(Class2, classes))

 plot_data <- plot_data %>%
   rowwise() %>%
   mutate(
     r = cor_mat[Class1, Class2],
     p_value = if (!is.null(pval_mat)) pval_mat[Class1, Class2] else NA_real_
   ) %>%
   ungroup() %>%
   mutate(
     IsDiag = idx1 == idx2,
     IsUpper = idx1 < idx2,
     pair_key = ifelse(Class1 < Class2,
                       paste(Class1, Class2, sep = "||"),
                       paste(Class2, Class1, sep = "||")),
     MarkSig = idx1 != idx2 & !is.na(p_value) & p_value < alpha,
     MarkFlip = !is_delta & MarkSig & pair_key %in% flip_pair_keys,
     r_plot = r,
     sig_stars = case_when(
       is.na(p_value) ~ "",
       p_value < 0.001 ~ "***",
       p_value < 0.01 ~ "**",
       p_value < 0.05 ~ "*",
       TRUE ~ ""
     ),
     label = case_when(
       IsDiag & !is_delta ~ sprintf("%.2f", r),
       IsDiag & is_delta ~ "",
       IsUpper & !is_delta ~ sig_stars,
       IsUpper & is_delta ~ "",
       TRUE ~ sprintf("%.2f", r)
     ),
     label_color = ifelse(IsDiag, "white", "black"),
     Class1 = factor(Class1, levels = classes),
     Class2 = factor(Class2, levels = rev(classes))
   )

 # Color scale limits
 if (is_delta) {
   fill_limits <- c(-1, 1)
   legend_name <- "\u0394r"
   subtitle_text <- "Change in correlation (LowInput - Control)"
 } else {
   fill_limits <- c(-1, 1)
   legend_name <- "CLR r"
   subtitle_text <- "CLR-space Pearson correlation"
 }

 # Create the plot
 p <- ggplot(plot_data, aes(x = Class1, y = Class2, fill = r_plot)) +
   geom_tile(color = "white", linewidth = 0.5)

 if (!is_delta) {
   p <- p + geom_tile(
     data = plot_data %>% filter(IsDiag),
     aes(x = Class1, y = Class2),
     inherit.aes = FALSE,
     fill = "black",
     color = "white",
     linewidth = 0.5
   )
 }

 if (!is_delta && length(flip_pair_keys) > 0) {
   p <- p + geom_point(
     data = plot_data %>% filter(MarkFlip),
     aes(x = Class1, y = Class2),
     shape = 0, size = 8.6, stroke = 1.55, color = "#C51B7D", fill = NA, inherit.aes = FALSE
   )
 }

 p <- p +
   geom_text(aes(label = label, color = label_color), size = 3) +
   scale_color_identity(guide = "none") +
   scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                        midpoint = 0, limits = fill_limits,
                        name = legend_name) +
   theme_minimal(base_size = 12) +
   theme(
     axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
     axis.text.y = element_text(color = "black", size = 10),
     plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
     plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey30"),
     plot.caption = element_text(size = 8, color = "grey40", hjust = 0.5),
     legend.position = "right",
     panel.grid = element_blank(),
     axis.line.x = element_line(color = "black", linewidth = 0.5),
     axis.line.y = element_line(color = "black", linewidth = 0.5)
   ) +
   coord_fixed()

 if (is_delta) {
   p <- p + labs(
     title = title_prefix,
     subtitle = subtitle_text,
     x = NULL, y = NULL
   )
 } else {
   p <- p + labs(
     title = paste0(title_prefix, " (n = ", n_samples, ")"),
     subtitle = subtitle_text,
     x = NULL, y = NULL,
     caption = paste0(
       "Upper triangle: * p<0.05, ** p<0.01, *** p<0.001. ",
       "Pink square: sign flip, significant in both conditions"
     )
   )
 }

 p
}

# ──────────────────────────────────────────────────────────────────────────────
# Create the three-panel figure
# ──────────────────────────────────────────────────────────────────────────────
suppfig_clr_ctrl <- create_clr_correlation_plot(
 clr_control, "A) Control",
 class_order = common_classes,
 flip_pair_keys = flip_pair_keys,
 alpha = 0.05
)

suppfig_clr_low <- create_clr_correlation_plot(
 clr_lowinput, "B) LowInput",
 class_order = common_classes,
 flip_pair_keys = flip_pair_keys,
 alpha = 0.05
)

suppfig_clr_delta <- create_clr_correlation_plot(
 cor_delta, "C) \u0394 Correlation",
 class_order = common_classes,
 is_delta = TRUE
)

# Combine into single figure (3 panels)
suppfig_clr_cor <- suppfig_clr_ctrl + suppfig_clr_low + suppfig_clr_delta +
 plot_layout(ncol = 3, guides = "collect")

save_fig(suppfig_clr_cor, "SuppFig_S6_CLR_Correlations.png", w = 18, h = 7, subdir = "supp")
save_fig(suppfig_clr_cor, "SuppFig_S6_CLR_Correlations.pdf", w = 18, h = 7, subdir = "supp")

# ---------------------------------------------------------------------------
# Class-level CLR correlations and the CTL->LIN change, as a table.
# The sign reversals quoted in the Results (PC-PS, PS-TG, LPC-PC) are otherwise
# only readable off the heatmap, with no numbers a reader can check.
# ---------------------------------------------------------------------------
if (!is.null(clr_control) && !is.null(clr_lowinput)) {
  shared <- intersect(rownames(clr_control$cor_mat), rownames(clr_lowinput$cor_mat))
  pairs  <- t(combn(shared, 2))

  delta_r <- data.frame(
    Class_A = pairs[, 1],
    Class_B = pairs[, 2],
    r_CTL   = clr_control$cor_mat[pairs],
    r_LIN   = clr_lowinput$cor_mat[pairs],
    p_CTL   = clr_control$pval_mat[pairs],
    p_LIN   = clr_lowinput$pval_mat[pairs],
    stringsAsFactors = FALSE
  )
  delta_r$delta_r       <- delta_r$r_LIN - delta_r$r_CTL
  delta_r$sign_reversal <- sign(delta_r$r_CTL) != sign(delta_r$r_LIN)
  delta_r$q_CTL <- p.adjust(delta_r$p_CTL, "BH")
  delta_r$q_LIN <- p.adjust(delta_r$p_LIN, "BH")
  delta_r <- delta_r[order(-abs(delta_r$delta_r)), ]

  write.csv(delta_r, "table/supp/SuppTable_S5C_Class_CLR_Correlation_Delta.csv",
            row.names = FALSE)
  message("Saved: table/supp/SuppTable_S5C_Class_CLR_Correlation_Delta.csv  (",
          nrow(delta_r), " class pairs, ",
          sum(delta_r$sign_reversal), " sign reversals)")
  message("  n samples: CTL ", clr_control$n_samples, ", LIN ", clr_lowinput$n_samples)

  key <- subset(delta_r, (Class_A == "PC"  & Class_B == "PS") |
                         (Class_A == "PS"  & Class_B == "PC") |
                         (Class_A == "PS"  & Class_B == "TG") |
                         (Class_A == "TG"  & Class_B == "PS") |
                         (Class_A == "LPC" & Class_B == "PC") |
                         (Class_A == "PC"  & Class_B == "LPC"))
  if (nrow(key)) {
    message("\n-- the three reversals quoted in the Results --")
    print(key[, c("Class_A","Class_B","r_CTL","r_LIN","delta_r",
                  "q_CTL","q_LIN","sign_reversal")], row.names = FALSE)
  }
}

message("\n\u2713 Supplementary Figure S6 (CLR Correlations) complete!\n")


# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
#
#     SUPPLEMENTARY FIGURE S7: LIPID CLASS RATIO PLOTS
#
# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════

message("\n══════════════════════════════════════════════════════════════")
message("SUPPLEMENTARY FIGURE S7: LIPID CLASS RATIO PLOTS")
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

# Valid classes for class extraction
ratio_classes <- c("TG", "DG", "MG", "PC", "PE", "DGDG", "MGDG", "SQDG", "LPC", "LPE", "PG", "PA", "PS")
class_pat <- paste0("^(", paste(ratio_classes, collapse = "|"), ")(?=\\()")

long_prep_ratio <- function(df) {
  df %>%
    tidyr::pivot_longer(-c(Compound_Name, Condition),
                        names_to = "Lipid", values_to = "Intensity") %>%
    dplyr::rename(Sample = Compound_Name) %>%
    dplyr::mutate(Class = stringr::str_extract(Lipid, class_pat)) %>%
    dplyr::filter(!is.na(Class)) %>%
    dplyr::group_by(Sample, Condition, Class) %>%
    dplyr::summarise(class_total = sum(Intensity, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(Sample) %>%
    dplyr::mutate(
      class_total = pmax(class_total, 0),
      minpos = dplyr::if_else(
        all(is.na(class_total) | class_total <= 0),
        NA_real_,
        min(class_total[class_total > 0], na.rm = TRUE)
      ),
      eps = dplyr::if_else(is.na(minpos), 0, minpos * 0.5),
      class_log = log10(class_total + eps)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(Sample, Condition, Class, class_log)
}

combined_ratio <- dplyr::bind_rows(
  control %>% dplyr::mutate(Condition = "Control"),
  lowinput %>% dplyr::mutate(Condition = "LowInput")
)
combined_ratio_paired <- restrict_to_condition_pairs(combined_ratio)
message("Matched genotype intersection for ratio plots: ", nrow(combined_ratio_paired) / 2, " pairs")

long_all_ratio <- long_prep_ratio(combined_ratio_paired)

wide_log_ratio <- long_all_ratio %>%
  tidyr::pivot_wider(names_from = Class, values_from = class_log, values_fill = NA_real_)

# Ratio lists (OPLS-DA VIP ratios and correlation-derived pairs)
opls_ratios <- c(
  "MG/SQDG", "MG/MGDG", "PS/SQDG", "DGDG/MG", "PG/SQDG", "DG/MG", "PE/SQDG",
  "LPC/PS", "DGDG/PS", "DG/PS", "PC/PS", "LPE/SQDG", "LPC/MG", "PE/PS",
  "MGDG/PG", "PG/PS", "LPC/LPE", "MGDG/PE", "MG/PG", "SQDG/TG", "PA/PS",
  "LPE/MGDG", "PC/SQDG"
)
pcorr_ratios <- c(
  "TG/DG", "PC/LPC", "PC/PS"
)

ratio_specs <- tibble::tibble(
  ratio_label = c(opls_ratios, pcorr_ratios),
  ratio_set = c(rep("OPLS-DA ratios", length(opls_ratios)), rep("Correlation-derived pairs", length(pcorr_ratios)))
) %>%
  tidyr::separate(ratio_label, into = c("Num", "Den"), sep = "/", remove = FALSE) %>%
  dplyr::mutate(
    Num = stringr::str_trim(Num),
    Den = stringr::str_trim(Den),
    ratio_id = paste0(Num, "_", Den)
  )

ratio_tbl <- wide_log_ratio
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

ratio_long_df <- ratio_tbl %>%
  dplyr::select(Sample, Condition, dplyr::all_of(ratio_specs$ratio_id)) %>%
  tidyr::pivot_longer(-c(Sample, Condition), names_to = "ratio_id", values_to = "Value") %>%
  dplyr::left_join(
    ratio_specs %>% dplyr::select(ratio_id, ratio_label, ratio_set),
    by = "ratio_id",
    relationship = "many-to-many"
  ) %>%
  dplyr::filter(!is.na(Value) & is.finite(Value))

ratio_long_df <- ratio_long_df %>%
  dplyr::mutate(
    Condition = factor(Condition, levels = c("Control", "LowInput")),
    ratio_set = factor(ratio_set, levels = c("OPLS-DA ratios", "Correlation-derived pairs"))
  )

paired_ratio_pvalue <- function(df_ratio) {
  paired_df <- df_ratio %>%
    dplyr::select(Sample, Condition, Value) %>%
    tidyr::pivot_wider(names_from = Condition, values_from = Value) %>%
    dplyr::filter(!is.na(Control), !is.na(LowInput))

  if (nrow(paired_df) < 3) return(NA_real_)

  tryCatch(
    wilcox.test(paired_df$LowInput, paired_df$Control, paired = TRUE, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
}

plot_one_set <- function(set_name, facet_levels, ncol = 5) {
  df <- ratio_long_df %>%
    dplyr::filter(ratio_set == set_name) %>%
    dplyr::mutate(ratio_label = factor(ratio_label, levels = facet_levels))

  stat_df <- df %>%
    dplyr::group_by(ratio_label) %>%
    dplyr::summarise(
      p_value = paired_ratio_pvalue(pick(Sample, Condition, Value)),
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

p_pcorr <- plot_one_set("Correlation-derived pairs", pcorr_ratios, ncol = 5) +
  ggtitle("B) Correlation-Derived Ratios (log10 scale)")

suppfig_s7 <- (p_opls / p_pcorr) +
  patchwork::plot_layout(heights = c(5, 2)) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

save_fig(suppfig_s7, "SuppFig_S7_Lipid_Ratios.png", w = 14, h = 13, subdir = "supp")
save_fig(suppfig_s7, "SuppFig_S7_Lipid_Ratios.pdf", w = 14, h = 13, subdir = "supp")

message("\n\u2713 Supplementary Figure S7 (Lipid Ratios) complete!\n")


# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════════════════

message("\n")
message("══════════════════════════════════════════════════════════════")
message("SECTION 4 COMPLETE!")
message("══════════════════════════════════════════════════════════════")
message("\n")
message("Generated outputs:")
message("  - fig/main/Figure2_OPLS_DA.png")
message("  - fig/main/Figure2_OPLS_DA.pdf")
message("  - fig/supp/SuppFig_S6_CLR_Correlations.png")
message("  - fig/supp/SuppFig_S6_CLR_Correlations.pdf")
message("  - fig/supp/SuppFig_S7_Lipid_Ratios.png")
message("  - fig/supp/SuppFig_S7_Lipid_Ratios.pdf")
message("  - table/supp/SuppTable_S7_OPLS_VIP_ratios.csv")
message("\n")
