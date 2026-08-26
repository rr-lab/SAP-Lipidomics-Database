# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3: LIPIDOME OVERVIEW
# ═══════════════════════════════════════════════════════════════════════════════
# Outputs:
#   - Figure 1: Lipidomics Landscape (main figure)
#   - Supp Fig S3: Compositional Contrasts (%TIC + ALR)
#   - Supp Fig S4: Species Counts (A-D panels)
#   - Supp Fig S5: PCA of Lipids
#   - Supp Fig S5A: Sample PCA colored by K.Cluster
#   - Supp Fig S5B: Sample PCA colored by Original_Race
#   - Supp Fig S5C: 3D sample PCA (PC1-PC2-PC3) by K.Cluster
#   - Supp Fig S5D: 3D sample PCA (PC1-PC2-PC3) by Original_Race
#   - Supp Table S5A: Class CLR Contrast
#   - Supp Table S5B: Class ALR Contrast
#   - Supp Table S6a: Species Summary
#   - Supp Table S6b: Species by Class
#   - Supp Table S6c: Species by SuperClass
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
# Color Palettes
# ─────────────────────────────────────────────────────────────────────────────
class_colors <- c(
  PC = "#00441B", PA = "#1B7837", PE = "#41AB5D", LPC = "#66C2A4", LPE = "#2CA25F",
  PG = "#78C679", PS = "#C2E699", DG = "#54278F", DGDG = "#F768A1", MG = "#8941ED",
  MGDG = "#FBB4D9", SQDG = "#9D4D6C", TG = "#ED804A"
)
# The expanded GWAS class set adds chemically distinct classes not used in LINEX.
expanded_class_colors <- c(
  class_colors,
  AEG = "#4E79A7", Cer = "#7A7A7A", FA = "#C69214", GalCer = "#8C6D31", SM = "#17A2A4"
)
expanded_class_order <- c(
  "MG", "DG", "TG", "MGDG", "DGDG", "SQDG", "PA", "PC", "PE", "PG", "PS", "LPC", "LPE",
  "AEG", "Cer", "GalCer", "SM", "FA"
)
condition_colors <- c(Control = "#440154FF", LowInput = "#FDE725FF")
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

class_to_group <- function(cls) {
  dplyr::case_when(
    cls %in% c("PC", "PA", "PE", "LPC", "LPE", "PG", "PS") ~ "Phospholipids",
    cls %in% c("MGDG", "DGDG", "SQDG") ~ "Glycolipids",
    cls %in% c("DG", "MG") ~ "Neutral",
    cls %in% c("TG") ~ "Storage",
    TRUE ~ "Other"
  )
}

parse_lipid_features <- function(lipid_name) {
  inside <- str_match(lipid_name, "\\(([^\\)]+)\\)")[, 2]
  if (is.na(inside)) return(tibble(total_c = NA_real_, total_db = NA_real_, n_chains = NA_integer_))
  parts <- str_split(inside, "/", simplify = TRUE)
  parts <- parts[parts != ""]
  n_chains <- length(parts)
  if (n_chains == 1 && str_detect(parts[1], "^\\d+\\s*:\\s*\\d+$")) {
    td <- str_split(parts[1], ":", simplify = TRUE)
    return(tibble(total_c = as.numeric(str_trim(td[1])), total_db = as.numeric(str_trim(td[2])), n_chains = 1L))
  }
  carb <- db <- rep(NA_real_, n_chains)
  for (i in seq_len(n_chains)) {
    if (str_detect(parts[i], "^\\d+\\s*:\\s*\\d+$")) {
      td <- str_split(parts[i], ":", simplify = TRUE)
      carb[i] <- as.numeric(str_trim(td[1]))
      db[i] <- as.numeric(str_trim(td[2]))
    }
  }
  tibble(total_c = sum(carb, na.rm = TRUE), total_db = sum(db, na.rm = TRUE), n_chains = as.integer(n_chains))
}

format_small_num <- function(x, digits = 2, threshold = 0.01) {
  sapply(x, function(v) {
    if (!is.finite(v)) return(NA_character_)
    if (v == 0) return(formatC(0, format = "f", digits = digits))
    if (abs(v) < threshold) return(paste0(ifelse(v < 0, "-", ""), "<", formatC(threshold, format = "f", digits = 2)))
    formatC(v, format = "f", digits = digits)
  })
}

save_fig <- function(p, filename, w = 10, h = 8, dpi = 300, subdir = "main") {
  p <- p + theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15))
  path <- file.path("fig", subdir, filename)
  ggsave(path, plot = p, width = w, height = h, units = "in", dpi = dpi, bg = "white")
  message("Saved: ", path)
}

# ─────────────────────────────────────────────────────────────────────────────
# Load Data
# ─────────────────────────────────────────────────────────────────────────────
message("\n══════════════════════════════════════════════════════════════")
message("Loading data...")
message("══════════════════════════════════════════════════════════════\n")

control_file <- "data/SPATS_fitted/non_normalized_intensities/Final_subset_control_all_lipids_fitted_phenotype_non_normalized.csv"
lowinput_file <- "data/SPATS_fitted/non_normalized_intensities/Final_subset_lowinput_all_lipids_fitted_phenotype_non_normalized.csv"

control_raw <- vroom(control_file, show_col_types = FALSE) %>% dplyr::select(-c(2, 3, 4))
low_raw <- vroom(lowinput_file, show_col_types = FALSE) %>% dplyr::select(-c(2, 3, 4))
# ─────────────────────────────────────────────────────────────────────────────
# Artefact exclusion
# ─────────────────────────────────────────────────────────────────────────────
# Two annotated features fail a basic plausibility check and are dropped before
# any analysis. Both survive only in the non-focused / expanded class sets; the
# 13 focused classes are unaffected.
#
#   Phytosphingosine   CTL median 1.6e5 (0.035% of TIC) -> LIN median 9.5e7
#                      (16.9% of TIC), the 2nd most abundant feature in the LIN
#                      lipidome. A free sphingoid long-chain base cannot be ~17%
#                      of a leaf lipidome; this is a co-eluting contaminant,
#                      misannotation, or batch/standard difference between the
#                      2019 (CTL) and 2022 (LIN) runs, not biology.
#
#   SM(d18:1/17:0)     The only member of the "SM" class. Plants do not
#                      synthesise sphingomyelin, and an odd-chain C17 species is
#                      a standard internal-standard chemotype. Its ~19x
#                      CTL-vs-LIN difference tracks batch, not genotype.
#
# Neither feature produced a GWAS candidate, so no genetic result depends on them.
artefact_features <- c("Phytosphingosine", "SM(d18:1/17:0)")

drop_artefacts <- function(df) {
  hit <- intersect(artefact_features, names(df))
  if (length(hit)) message("  Dropping artefact feature(s): ", paste(hit, collapse = ", "))
  dplyr::select(df, -dplyr::any_of(artefact_features))
}

control_raw <- drop_artefacts(control_raw)
low_raw     <- drop_artefacts(low_raw)

colnames(control_raw)[1] <- "Compound_Name"
colnames(low_raw)[1] <- "Compound_Name"
control_raw <- control_raw %>% mutate(Condition = "Control")
low_raw <- low_raw %>% mutate(Condition = "LowInput")

message("Control samples: ", nrow(control_raw))
message("LowInput samples: ", nrow(low_raw))

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 1: LIPIDOMICS LANDSCAPE
# ═══════════════════════════════════════════════════════════════════════════════

message("\n══════════════════════════════════════════════════════════════")
message("FIGURE 1: LIPIDOMICS LANDSCAPE")
message("══════════════════════════════════════════════════════════════\n")

# ─────────────────────────────────────────────────────────────────────────────
# Compute %TIC by class
# ─────────────────────────────────────────────────────────────────────────────
classsum_tic_long <- function(df, label, drop_classes = c("DGTS"), valid_classes = classes_with_lyso) {
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

ctrl_pct <- classsum_tic_long(control_raw, "Control")
low_pct <- classsum_tic_long(low_raw, "LowInput")
pct_tbl <- bind_rows(ctrl_pct, low_pct)

pct_mean <- pct_tbl %>%
  group_by(Condition, Class) %>%
  summarise(pct_mean = mean(pct, na.rm = TRUE), .groups = "drop")

# Fold changes
fold_tbl <- pct_mean %>%
  pivot_wider(names_from = Condition, values_from = pct_mean, values_fill = 0) %>%
  mutate(eps = 1e-8, FoldChange = (LowInput + eps) / (Control + eps),
         log2FC = log2(FoldChange), Direction = ifelse(log2FC > 0, "Increased", "Decreased"))

# Legend labels
legend_tbl <- pct_mean %>%
  pivot_wider(names_from = Condition, values_from = pct_mean, values_fill = 0) %>%
  mutate(legend_label = sprintf("%s (C %.1f%%, L %.1f%%)", Class, Control, LowInput))
label_map <- setNames(legend_tbl$legend_label, legend_tbl$Class)

classes_all <- sort(unique(pct_mean$Class))
pct_mean <- pct_mean %>% mutate(Class = factor(Class, levels = classes_all))

# ─────────────────────────────────────────────────────────────────────────────
# Supp Fig S3A: %TIC Composition
# ─────────────────────────────────────────────────────────────────────────────
suppfig_s3a_tic <- ggplot(pct_mean, aes(x = pct_mean, y = Condition, fill = Class)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(pct_mean >= 3, sprintf("%.1f%%", pct_mean), "")),
            position = position_stack(vjust = 0.5), colour = "white", size = 3.5, fontface = "bold") +
  scale_fill_manual(values = class_colors[classes_all], breaks = classes_all, labels = label_map[classes_all], name = NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 100, 25), labels = function(x) paste0(x, "%")) +
  labs(x = "% of TIC (mean across samples)", y = NULL, title = "A) Class %TIC Composition") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right", legend.text = element_text(size = 9),
        axis.text = element_text(color = "black", size = 11), axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5), panel.grid.major.y = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

# ─────────────────────────────────────────────────────────────────────────────
# CLR and ALR contrasts
# ─────────────────────────────────────────────────────────────────────────────
close_replace_rows <- function(mat, delta_frac = 0.5) {
  mat <- as.matrix(mat)
  out <- matrix(NA_real_, nrow = nrow(mat), ncol = ncol(mat))
  colnames(out) <- colnames(mat); rownames(out) <- rownames(mat)
  for (i in seq_len(nrow(mat))) {
    x <- mat[i, ]; x[!is.finite(x)] <- 0; x[x < 0] <- 0
    if (sum(x) <= 0) next
    if (any(x == 0)) { pos <- x[x > 0]; if (length(pos) > 0) { delta <- min(pos) * delta_frac; x[x == 0] <- delta }}
    out[i, ] <- x / sum(x)
  }
  out
}

bootstrap_context_contrast <- function(mat, cond_vec, n_boot = 500, seed = 2026) {
  mat <- as.matrix(mat); keep <- complete.cases(mat) & !is.na(cond_vec)
  mat <- mat[keep, , drop = FALSE]; cond_vec <- cond_vec[keep]
  idx_ctrl <- which(cond_vec == "Control"); idx_low <- which(cond_vec == "LowInput")
  if (length(idx_ctrl) < 5 || length(idx_low) < 5) return(NULL)
  est <- colMeans(mat[idx_low, , drop = FALSE], na.rm = TRUE) - colMeans(mat[idx_ctrl, , drop = FALSE], na.rm = TRUE)
  boot <- matrix(NA_real_, nrow = n_boot, ncol = ncol(mat)); colnames(boot) <- colnames(mat)
  set.seed(seed)
  for (b in seq_len(n_boot)) {
    ctrl_b <- sample(idx_ctrl, size = length(idx_ctrl), replace = TRUE)
    low_b <- sample(idx_low, size = length(idx_low), replace = TRUE)
    boot[b, ] <- colMeans(mat[low_b, , drop = FALSE], na.rm = TRUE) - colMeans(mat[ctrl_b, , drop = FALSE], na.rm = TRUE)
  }
  tibble(Feature = colnames(mat), Effect = est,
         CI_Low = apply(boot, 2, function(x) quantile(x, 0.025, na.rm = TRUE)),
         CI_High = apply(boot, 2, function(x) quantile(x, 0.975, na.rm = TRUE)))
}

class_comp_wide <- pct_tbl %>%
  dplyr::select(Sample, Condition, Class, pct) %>%
  pivot_wider(names_from = Class, values_from = pct, values_fill = 0)

class_cols_comp <- setdiff(names(class_comp_wide), c("Sample", "Condition"))
class_comp_raw <- as.matrix(class_comp_wide[, class_cols_comp, drop = FALSE]) / 100
rownames(class_comp_raw) <- class_comp_wide$Sample
class_comp <- close_replace_rows(class_comp_raw, delta_frac = 0.5)

# CLR
clr_class <- log(class_comp) - rowMeans(log(class_comp), na.rm = TRUE)
clr_tbl <- bootstrap_context_contrast(clr_class, class_comp_wide$Condition, n_boot = 500, seed = 1101) %>%
  dplyr::rename(Class = Feature) %>%
  mutate(AbsEffect = abs(Effect), Class = factor(Class, levels = rev(Class[order(AbsEffect, decreasing = TRUE)])))

# ALR
alr_ref_class <- if ("TG" %in% colnames(class_comp)) "TG" else names(sort(colMeans(class_comp, na.rm = TRUE), decreasing = TRUE))[1]
alr_num <- setdiff(colnames(class_comp), alr_ref_class)
alr_class <- log(sweep(class_comp[, alr_num, drop = FALSE], 1, class_comp[, alr_ref_class], "/"))
alr_tbl <- bootstrap_context_contrast(alr_class, class_comp_wide$Condition, n_boot = 500, seed = 1201) %>%
  mutate(Class = str_replace(Feature, paste0("_vs_", alr_ref_class, "$"), ""),
         AbsEffect = abs(Effect), Class = factor(Class, levels = rev(Class[order(AbsEffect, decreasing = TRUE)])))

# Save tables with NEW numbering (S5)
write.csv(clr_tbl, "table/supp/SuppTable_S5A_Class_CLR_Contrast.csv", row.names = FALSE)
write.csv(alr_tbl, "table/supp/SuppTable_S5B_Class_ALR_Contrast.csv", row.names = FALSE)
message("Saved: table/supp/SuppTable_S5A_Class_CLR_Contrast.csv")
message("Saved: table/supp/SuppTable_S5B_Class_ALR_Contrast.csv")

# ─────────────────────────────────────────────────────────────────────────────
# Fig 1B: CLR contrast
# ─────────────────────────────────────────────────────────────────────────────
fig1b <- ggplot(clr_tbl, aes(x = Class, y = Effect, fill = Class)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.5) +
  geom_col(width = 0.72, color = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = CI_Low, ymax = CI_High), width = 0.2, linewidth = 0.6) +
  coord_flip() +
  scale_fill_manual(values = class_colors, na.value = "grey70") +
  scale_y_continuous(labels = function(x) format_small_num(x, digits = 2, threshold = 0.01)) +
  labs(x = NULL, y = "CLR contrast (LowInput - Control)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", axis.text = element_text(color = "black", size = 11),
        axis.title = element_text(face = "bold"),
        panel.grid.major.y = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

message("Created: Fig 1B (CLR contrast)")

# ─────────────────────────────────────────────────────────────────────────────
# Supp Fig S3B: ALR contrast
# ─────────────────────────────────────────────────────────────────────────────
suppfig_s3b_alr <- ggplot(alr_tbl, aes(x = Class, y = Effect, color = Class)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.5) +
  geom_errorbar(aes(ymin = CI_Low, ymax = CI_High), width = 0.2, linewidth = 0.7) +
  geom_point(size = 2.2) + coord_flip() +
  scale_color_manual(values = class_colors, na.value = "grey70") +
  labs(x = NULL, y = "ALR contrast (LowInput - Control)", title = paste0("B) Class ALR Contrast (Ref = ", alr_ref_class, ")")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(face = "bold"), plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid.major.y = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

# ─────────────────────────────────────────────────────────────────────────────
# TIC Normalization for chemical shifts and PCA
# ─────────────────────────────────────────────────────────────────────────────
tic_normalize <- function(df) {
  X <- df %>% dplyr::select(-Compound_Name, -Condition) %>% dplyr::select(where(is.numeric)) %>% as.matrix()
  rs <- rowSums(X, na.rm = TRUE); rs[rs == 0 | is.na(rs)] <- NA_real_
  Xpct <- sweep(X, 1, rs, FUN = "/") * 100; Xpct[is.na(Xpct)] <- 0
  bind_cols(df %>% dplyr::select(Compound_Name, Condition), as.data.frame(Xpct))
}

subset_lipids <- function(df, keep_classes = valid_classes, drop_classes = c("DGTS")) {
  lipid_cols <- names(df)[grepl("\\(", names(df))]
  keep_pat <- paste0("^(", paste(keep_classes, collapse = "|"), ")\\(")
  keep_cols <- lipid_cols[grepl(keep_pat, lipid_cols)]
  if (length(drop_classes) > 0) { drop_pat <- paste0("^(", paste(drop_classes, collapse = "|"), ")\\("); keep_cols <- keep_cols[!grepl(drop_pat, keep_cols)] }
  df %>% dplyr::select(Compound_Name, Condition, all_of(keep_cols))
}

control_sub <- subset_lipids(control_raw)
low_sub <- subset_lipids(low_raw)
control_tic <- tic_normalize(control_sub)
low_tic <- tic_normalize(low_sub)

shared_cols <- intersect(setdiff(names(control_tic), c("Compound_Name", "Condition")),
                         setdiff(names(low_tic), c("Compound_Name", "Condition")))
control_tic <- control_tic %>% dplyr::select(Compound_Name, Condition, all_of(shared_cols))
low_tic <- low_tic %>% dplyr::select(Compound_Name, Condition, all_of(shared_cols))

# ─────────────────────────────────────────────────────────────────────────────
# Fig 1C: Chemical Shifts
# ─────────────────────────────────────────────────────────────────────────────
calculate_class_chem <- function(df_tic, condition_name) {
  lipid_cols <- setdiff(names(df_tic), c("Compound_Name", "Condition"))
  lipid_cols <- lipid_cols[grepl("\\(", lipid_cols)]
  X <- df_tic %>% dplyr::select(all_of(lipid_cols)) %>% mutate(across(everything(), ~suppressWarnings(as.numeric(.x))))
  mean_pct <- colMeans(as.matrix(X), na.rm = TRUE)
  tibble(Lipid = lipid_cols) %>%
    mutate(Class = get_lipid_class(Lipid), MeanPctTIC = mean_pct[Lipid], Condition = condition_name) %>%
    bind_cols(purrr::map_dfr(.$Lipid, parse_lipid_features)) %>%
    filter(!is.na(total_c), Class %in% valid_classes) %>%
    group_by(Class, Condition) %>%
    summarise(WeightedMeanC = weighted.mean(total_c, MeanPctTIC, na.rm = TRUE),
              WeightedMeanDB = weighted.mean(total_db, MeanPctTIC, na.rm = TRUE), n_species = n(), .groups = "drop")
}

chem_control <- calculate_class_chem(control_tic, "Control")
chem_low <- calculate_class_chem(low_tic, "LowInput")
chem_summary <- bind_rows(chem_control, chem_low)

fig1c <- ggplot(chem_summary, aes(x = WeightedMeanC, y = WeightedMeanDB, color = Class, shape = Condition)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_path(aes(group = Class), arrow = arrow(type = "closed", length = unit(0.12, "inches")), linewidth = 0.8, alpha = 0.7) +
  scale_color_manual(values = class_colors) +
  scale_shape_manual(values = c("Control" = 16, "LowInput" = 17)) +
  labs(x = "Weighted Mean Total Carbons", y = "Weighted Mean Double Bonds") +
  theme_minimal(base_size = 12) +
  theme(axis.text = element_text(color = "black"), axis.title = element_text(face = "bold"),
        legend.position = "right",
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

message("Created: Fig 1C (Chemical shifts)")

# ─────────────────────────────────────────────────────────────────────────────
# Fig 1D: LION Enrichment
# ─────────────────────────────────────────────────────────────────────────────
lion_file <- "table/Linex2/LION-enrichment.csv"
if (file.exists(lion_file)) {
  lion_df <- vroom(lion_file, show_col_types = FALSE)
  lion_sig <- lion_df %>%
    dplyr::filter(`p-value` < 0.05) %>%
    dplyr::rename(p_value = `p-value`, FDR_q = `FDR q-value`, Description = Discription) %>%
    dplyr::mutate(logP = -log10(p_value), logQ = -log10(FDR_q), Description = forcats::fct_reorder(Description, logQ))

  fig1d <- ggplot(lion_sig, aes(x = logQ, y = Description)) +
    geom_vline(xintercept = -log10(0.05), color = "red", linetype = "dashed", linewidth = 0.8) +
    facet_grid(Regulated ~ ., scales = "free_y", space = "free_y") +
    geom_point(aes(size = Annotated, color = logP), alpha = 0.8) +
    scale_color_viridis_c(name = expression(-log[10](p)), option = "C") +
    scale_size_continuous(name = "Annotated\nlipids", range = c(2, 8)) +
    labs(x = expression(-log[10](FDR~q~value)), y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.title.x = element_text(face = "bold", size = 11),
          axis.text.y = element_text(size = 9, color = "black"),
          strip.text = element_text(face = "bold", size = 10),
          legend.position = "right",
          axis.line.x = element_line(color = "black", linewidth = 0.5),
          axis.line.y = element_line(color = "black", linewidth = 0.5))
  message("Created: Fig 1D (LION Enrichment)")
} else {
  fig1d <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "LION Enrichment\n(file not found)", size = 4) +
    theme_void() + theme(panel.border = element_rect(color = "grey70", fill = NA))
  message("Warning: LION file not found")
}

# ─────────────────────────────────────────────────────────────────────────────
# Assemble Figure 1
# ─────────────────────────────────────────────────────────────────────────────
# The workflow schematic that used to be panel A here has moved to Figure 1
# (17_figure1_analytical_comparability.R). The Results cite the workflow in
# their opening paragraph, before Figure 1, so keeping it here created a
# forward reference -- and when the file was missing this script silently
# emitted a "Fig 1A (placeholder)" box into the compiled manuscript.
fig1_row1 <- fig1b + fig1c + plot_layout(widths = c(1.05, 1))
fig1_row2 <- fig1d

figure1 <- (fig1_row1 / fig1_row2) +
  plot_layout(heights = c(1, 1.15)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 16, face = "bold"))

save_fig(figure1, "Figure1_Lipidomics_Landscape.png", w = 14, h = 13)
save_fig(figure1, "Figure1_Lipidomics_Landscape.pdf", w = 14, h = 13)
message("\n Figure 1 complete!")

# ═══════════════════════════════════════════════════════════════════════════════
# SUPP FIG S3: COMPOSITIONAL CONTRASTS
# ═══════════════════════════════════════════════════════════════════════════════
message("\n══════════════════════════════════════════════════════════════")
message("SUPPLEMENTARY FIGURE S3: COMPOSITIONAL CONTRASTS")
message("══════════════════════════════════════════════════════════════\n")

suppfig_s3_comp <- suppfig_s3a_tic + suppfig_s3b_alr + plot_layout(widths = c(1.2, 1))
save_fig(suppfig_s3_comp, "SuppFig_S3_Compositional_Contrasts.png", w = 15, h = 7, subdir = "supp")
save_fig(suppfig_s3_comp, "SuppFig_S3_Compositional_Contrasts.pdf", w = 15, h = 7, subdir = "supp")
message(" Supp Fig S3 (Compositional) complete!")

# ─────────────────────────────────────────────────────────────────────────────
# SUPP FIG S8: NON-FOCUSED LIPID-CLASS CONTEXT
# ─────────────────────────────────────────────────────────────────────────────
message("\n══════════════════════════════════════════════════════════════")
message("SUPPLEMENTARY FIGURE S8: NON-FOCUSED LIPID-CLASS CONTEXT")
message("══════════════════════════════════════════════════════════════\n")

# This panel intentionally excludes the 13 classes used in the focused
# composition, ratio-network, and LINEX analyses. The condition-specific class
# annotation files use "::" to join alternative identifications, so aliases are
# resolved before class-level aggregation.
focused_class_order <- c("MG", "DG", "TG", "MGDG", "DGDG", "SQDG", "PA", "PC", "PE", "PG", "PS", "LPC", "LPE")
nonfocused_class_order <- c(
  "Ether lipid", "Sphingolipid", "Fatty acid", "Fatty acid amide", "Glycerolipid",
  "Glycerophospholipid", "Betaine lipid", "Sterol", "Prenol", "Terpenoid"
)
nonfocused_class_colors <- c(
  "Ether lipid" = "#4E79A7", Sphingolipid = "#7A7A7A", "Fatty acid" = "#C69214",
  "Fatty acid amide" = "#B07AA1", Glycerolipid = "#F28E2B",
  Glycerophospholipid = "#59A14F", "Betaine lipid" = "#76B7B2", Sterol = "#9C755F",
  Prenol = "#EDC948", Terpenoid = "#BAB0AC"
)

normalize_annotation_class <- function(x) {
  x <- str_replace(x, "_$", "")
  recode(
    x,
    "Fatty acid derivative" = "Fatty acid",
    "N-acylethanolamine" = "Fatty acid amide",
    "Steroid" = "Sterol",
    .default = x
  )
}

annotation_files <- c(
  "data/lipid_class/lipid_classes.csv",
  "data/lipid_class/lipid_classes_control.csv",
  "data/lipid_class/lipid_classes_lowinput.csv",
  "data/lipid_class/final_lipid_classes.csv"
)
lipid_annotation <- map_dfr(annotation_files, function(path) {
  vroom(path, show_col_types = FALSE) %>%
    transmute(AnnotationClass = normalize_annotation_class(Class), Alias = str_split(Lipids, "::")) %>%
    unnest(Alias) %>%
    filter(!is.na(AnnotationClass))
}) %>%
  distinct(Alias, .keep_all = TRUE)

prefix_to_annotation_class <- function(lipid) {
  code <- str_extract(lipid, "^[A-Za-z0-9]+(?=\\()")
  case_when(
    code %in% focused_class_order ~ "Focused",
    code == "AEG" ~ "Ether lipid",
    code %in% c("Cer", "GalCer", "SM") ~ "Sphingolipid",
    code == "NAE" ~ "Fatty acid amide",
    code == "CL" ~ "Glycerophospholipid",
    code == "DGTS" ~ "Betaine lipid",
    code == "FA" ~ "Fatty acid",
    TRUE ~ NA_character_
  )
}

summarize_nonfocused_classes <- function(df, condition_label) {
  df %>%
    select(-any_of("Condition")) %>%
    pivot_longer(-Compound_Name, names_to = "Lipid", values_to = "Intensity") %>%
    left_join(lipid_annotation, by = c("Lipid" = "Alias")) %>%
    mutate(
      PrefixClass = prefix_to_annotation_class(Lipid),
      Class = coalesce(PrefixClass, AnnotationClass),
      Class = na_if(Class, "Focused")
    ) %>%
    filter(Class %in% nonfocused_class_order) %>%
    group_by(Sample = Compound_Name, Class) %>%
    summarise(sum_int = sum(as.numeric(Intensity), na.rm = TRUE), .groups = "drop") %>%
    mutate(Condition = condition_label) %>%
    group_by(Sample, Condition) %>%
    mutate(pct = 100 * sum_int / sum(sum_int, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(nchar(Sample) >= 5)
}

ctrl_pct_nonfocused <- summarize_nonfocused_classes(control_raw, "Control")
low_pct_nonfocused <- summarize_nonfocused_classes(low_raw, "LowInput")
pct_tbl_nonfocused <- bind_rows(ctrl_pct_nonfocused, low_pct_nonfocused) %>%
  mutate(Class = factor(Class, levels = nonfocused_class_order))

pct_mean_nonfocused <- pct_tbl_nonfocused %>%
  group_by(Condition, Class) %>%
  summarise(pct_mean = mean(pct, na.rm = TRUE), .groups = "drop") %>%
  complete(Condition, Class = factor(nonfocused_class_order, levels = nonfocused_class_order), fill = list(pct_mean = 0))

nonfocused_label_tbl <- pct_mean_nonfocused %>%
  pivot_wider(names_from = Condition, values_from = pct_mean, values_fill = 0) %>%
  mutate(legend_label = sprintf("%s (CTL %.2f%%, LIN %.2f%%)", Class, Control, LowInput))
nonfocused_label_map <- setNames(nonfocused_label_tbl$legend_label, nonfocused_label_tbl$Class)

suppfig_s8a_nonfocused_tic <- ggplot(pct_mean_nonfocused, aes(x = pct_mean, y = Condition, fill = Class)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.25) +
  geom_text(aes(label = ifelse(pct_mean >= 5, sprintf("%.1f%%", pct_mean), "")),
            position = position_stack(vjust = 0.5), colour = "white", size = 3.2, fontface = "bold") +
  scale_fill_manual(values = nonfocused_class_colors[nonfocused_class_order],
                    breaks = nonfocused_class_order,
                    labels = nonfocused_label_map[nonfocused_class_order], name = NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 100, 25), labels = function(x) paste0(x, "%")) +
  labs(x = "Relative abundance (% of the non-focused lipid pool)", y = NULL,
       title = "A) Non-Focused Lipid Classes") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right", legend.text = element_text(size = 8),
        axis.text = element_text(color = "black", size = 11), axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5), panel.grid.major.y = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

class_comp_wide_nonfocused <- pct_tbl_nonfocused %>%
  select(Sample, Condition, Class, pct) %>%
  pivot_wider(names_from = Class, values_from = pct, values_fill = 0) %>%
  select(Sample, Condition, any_of(nonfocused_class_order))

nonfocused_matrix <- as.matrix(class_comp_wide_nonfocused[, nonfocused_class_order, drop = FALSE]) / 100
rownames(nonfocused_matrix) <- class_comp_wide_nonfocused$Sample
nonfocused_comp <- close_replace_rows(nonfocused_matrix, delta_frac = 0.5)
nonfocused_clr <- log(nonfocused_comp) - rowMeans(log(nonfocused_comp), na.rm = TRUE)

nonfocused_clr_tbl <- bootstrap_context_contrast(
  nonfocused_clr, class_comp_wide_nonfocused$Condition, n_boot = 500, seed = 1301
) %>%
  rename(Class = Feature) %>%
  mutate(Class = factor(Class, levels = rev(nonfocused_class_order)))

suppfig_s8b_nonfocused_clr <- ggplot(nonfocused_clr_tbl, aes(x = Class, y = Effect, fill = Class)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.5) +
  geom_col(width = 0.72, color = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = CI_Low, ymax = CI_High), width = 0.2, linewidth = 0.6) +
  coord_flip() +
  scale_fill_manual(values = nonfocused_class_colors[nonfocused_class_order], guide = "none") +
  scale_y_continuous(labels = function(x) format_small_num(x, digits = 2, threshold = 0.01)) +
  labs(x = NULL, y = "CLR contrast (LIN - CTL)",
       title = "B) Non-Focused Class CLR Contrast") +
  theme_minimal(base_size = 12) +
  theme(axis.text = element_text(color = "black", size = 10), axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5), panel.grid.major.y = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

suppfig_s8_nonfocused_classes <- suppfig_s8a_nonfocused_tic + suppfig_s8b_nonfocused_clr +
  plot_layout(widths = c(1.35, 1))
save_fig(suppfig_s8_nonfocused_classes, "SuppFig_S8_NonFocused_Class_Composition.png", w = 17, h = 8, subdir = "supp")
write.csv(nonfocused_clr_tbl, "table/supp/NonFocused_Lipid_Class_CLR_Contrast.csv", row.names = FALSE)
message(" Supp Fig S8 (Non-Focused Lipid-Class Context) complete!")

# Final S8 uses all shared lipid features as the CLR denominator, avoiding a
# denominator defined only by the non-focused feature subset.
shared_lipids_s8 <- intersect(
  setdiff(names(control_raw), c("Compound_Name", "Condition")),
  setdiff(names(low_raw), c("Compound_Name", "Condition"))
)
global_s8_wide <- bind_rows(
  control_raw %>% select(Compound_Name, all_of(shared_lipids_s8)) %>% mutate(Condition = "Control"),
  low_raw %>% select(Compound_Name, all_of(shared_lipids_s8)) %>% mutate(Condition = "LowInput")
)
global_s8_matrix <- as.matrix(global_s8_wide[, shared_lipids_s8, drop = FALSE])
rownames(global_s8_matrix) <- global_s8_wide$Compound_Name
global_s8_comp <- close_replace_rows(global_s8_matrix, delta_frac = 0.5)
global_s8_clr <- log(global_s8_comp) - rowMeans(log(global_s8_comp), na.rm = TRUE)

# Assign every shared feature to its broad annotation class. Standard lipid
# prefixes provide a consistent fallback when a class file stores an alias.
prefix_to_broad_class <- function(lipid) {
  code <- str_extract(lipid, "^[A-Za-z0-9]+(?=\\()")
  case_when(
    code %in% c("MG", "DG", "TG", "MGDG", "DGDG", "SQDG") ~ "Glycerolipid",
    code %in% c("PA", "PC", "PE", "PG", "PS", "LPC", "LPE", "CL") ~ "Glycerophospholipid",
    code == "AEG" ~ "Ether lipid",
    code %in% c("Cer", "GalCer", "SM") ~ "Sphingolipid",
    code == "NAE" ~ "Fatty acid amide",
    code == "DGTS" ~ "Betaine lipid",
    code == "FA" ~ "Fatty acid",
    TRUE ~ NA_character_
  )
}

s8_feature_class <- tibble(Lipid = shared_lipids_s8) %>%
  left_join(lipid_annotation, by = c("Lipid" = "Alias")) %>%
  mutate(
    FeatureClass = coalesce(prefix_to_broad_class(Lipid), AnnotationClass),
    IsFocused = str_extract(Lipid, "^[A-Za-z0-9]+(?=\\()") %in% focused_class_order
  ) %>%
  filter(!is.na(FeatureClass))

s8_class_order <- c(
  "Ether lipid", "Sphingolipid", "Fatty acid", "Fatty acid amide", "Glycerolipid",
  "Glycerophospholipid", "Betaine lipid", "Sterol", "Prenol", "Terpenoid"
)
s8_class_colors <- nonfocused_class_colors[s8_class_order]

s8_class_totals <- global_s8_wide %>%
  select(Compound_Name, Condition, all_of(s8_feature_class$Lipid)) %>%
  pivot_longer(-c(Compound_Name, Condition), names_to = "Lipid", values_to = "Intensity") %>%
  left_join(s8_feature_class, by = "Lipid") %>%
  group_by(Compound_Name, Condition, FeatureClass) %>%
  summarise(Intensity = sum(as.numeric(Intensity), na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = FeatureClass, values_from = Intensity, values_fill = 0) %>%
  select(Compound_Name, Condition, any_of(s8_class_order))

# DGTS is detected only in LIN and therefore has no shared feature for this CLR.
s8_class_order_present <- intersect(s8_class_order, names(s8_class_totals))
s8_class_matrix <- as.matrix(s8_class_totals[, s8_class_order_present, drop = FALSE])
rownames(s8_class_matrix) <- s8_class_totals$Compound_Name
s8_class_comp <- close_replace_rows(s8_class_matrix, delta_frac = 0.5)
s8_class_clr <- log(s8_class_comp) - rowMeans(log(s8_class_comp), na.rm = TRUE)
s8_class_clr_tbl <- bootstrap_context_contrast(
  s8_class_clr, s8_class_totals$Condition, n_boot = 500, seed = 1303
) %>%
  rename(Class = Feature) %>%
  mutate(Class = factor(Class, levels = rev(s8_class_order_present)))

suppfig_s8a_class_clr <- ggplot(s8_class_clr_tbl, aes(x = Class, y = Effect, fill = Class)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.5) +
  geom_col(width = 0.72, color = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = CI_Low, ymax = CI_High), width = 0.2, linewidth = 0.55) +
  coord_flip() +
  scale_fill_manual(values = s8_class_colors[s8_class_order_present], guide = "none") +
  labs(x = NULL, y = "Class CLR contrast (LIN - CTL)",
       title = "A) Global CLR Contrast Across Annotated Lipid Classes") +
  theme_minimal(base_size = 12) +
  theme(axis.text = element_text(color = "black", size = 10), axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5), panel.grid.major.y = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

s8_global_clr_tbl <- bootstrap_context_contrast(
  global_s8_clr[, s8_feature_class$Lipid[!s8_feature_class$IsFocused], drop = FALSE], global_s8_wide$Condition, n_boot = 500, seed = 1302
) %>%
  rename(Lipid = Feature) %>%
  left_join(s8_feature_class %>% select(Lipid, FeatureClass), by = "Lipid") %>%
  mutate(AbsEffect = abs(Effect)) %>%
  arrange(desc(AbsEffect))

# Display the 20 largest global CLR shifts; retain all annotated features in the table.
s8_display_tbl <- s8_global_clr_tbl %>%
  slice_head(n = 20) %>%
  mutate(Lipid = factor(Lipid, levels = rev(Lipid)))

suppfig_s8b_global_clr <- ggplot(s8_display_tbl, aes(x = Lipid, y = Effect, fill = FeatureClass)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.5) +
  geom_col(width = 0.72, color = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = CI_Low, ymax = CI_High), width = 0.2, linewidth = 0.55) +
  coord_flip() +
  scale_fill_manual(values = nonfocused_class_colors, name = "Annotated class") +
  labs(x = NULL, y = "Global CLR contrast (LIN - CTL)",
       title = "B) Largest Global CLR Shifts in Individual Features") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8),
        axis.text = element_text(color = "black", size = 9), axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5), panel.grid.major.y = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

suppfig_s8_feature_validation <- suppfig_s8a_class_clr + suppfig_s8b_global_clr +
  plot_layout(widths = c(0.9, 1.6))
save_fig(suppfig_s8_feature_validation, "SuppFig_S8_NonFocused_Lipid_Class_Context.png", w = 20, h = 10, subdir = "supp")
write.csv(s8_global_clr_tbl, "table/supp/NonFocused_Feature_Global_CLR_Contrast.csv", row.names = FALSE)
write.csv(s8_class_clr_tbl, "table/supp/AllAnnotatedClass_Global_CLR_Contrast.csv", row.names = FALSE)
message(" Supp Fig S8 (Global CLR Class and Feature Contrasts) complete!")

# ═══════════════════════════════════════════════════════════════════════════════
# SUPP FIG S4: SPECIES COUNTS
# ═══════════════════════════════════════════════════════════════════════════════
message("\n══════════════════════════════════════════════════════════════")
message("SUPPLEMENTARY FIGURE S4: SPECIES COUNTS")
message("══════════════════════════════════════════════════════════════\n")

# Load lipid class info
lipid_class_file <- "data/lipid_class/final_lipid_classes.csv"
if (file.exists(lipid_class_file)) {
  lipid_class_info <- vroom(lipid_class_file, show_col_types = FALSE) %>%
    transmute(Lipid = Lipids, SuperClass = Class)
  message("Loaded lipid class info: ", nrow(lipid_class_info), " entries")
} else {
  lipid_class_info <- tibble(Lipid = character(), SuperClass = character())
  message("Warning: lipid class file not found")
}

ctrl_lipids <- names(control_raw)[grepl("\\(", names(control_raw))]
low_lipids <- names(low_raw)[grepl("\\(", names(low_raw))]
n_ctrl <- length(ctrl_lipids); n_low <- length(low_lipids)
common_lipids <- intersect(ctrl_lipids, low_lipids); n_common <- length(common_lipids)
ctrl_only <- setdiff(ctrl_lipids, low_lipids); low_only <- setdiff(low_lipids, ctrl_lipids)

message("Control lipids: ", n_ctrl, " | LowInput: ", n_low, " | Common: ", n_common)

# S4A: Total counts
total_counts <- tibble(
  Category = c("Control Total", "LowInput Total", "Common", "Control Only", "LowInput Only"),
  Count = c(n_ctrl, n_low, n_common, length(ctrl_only), length(low_only)),
  Type = c("Total", "Total", "Overlap", "Unique", "Unique")
)
total_counts$Category <- factor(total_counts$Category, levels = c("Control Total", "LowInput Total", "Common", "Control Only", "LowInput Only"))

suppfig_s4a <- ggplot(total_counts, aes(x = Category, y = Count, fill = Type)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = Count), vjust = -0.3, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("Total" = "#440154FF", "Overlap" = "#21908CFF", "Unique" = "#FDE725FF")) +
  labs(title = "A) Lipid Species Counts", x = NULL, y = "Number of Species") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top",
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

# S4B: By class
count_by_class <- function(lipids, condition_name) {
  tibble(Lipid = lipids) %>% mutate(Class = get_lipid_class(Lipid)) %>% count(Class, name = "Count") %>% mutate(Condition = condition_name)
}
class_ctrl <- count_by_class(ctrl_lipids, "Control")
class_low <- count_by_class(low_lipids, "LowInput")
class_counts <- bind_rows(class_ctrl, class_low) %>% filter(!is.na(Class), Class != "Other")
class_order <- class_counts %>% group_by(Class) %>% summarise(Total = sum(Count)) %>% arrange(desc(Total)) %>% pull(Class)
class_counts$Class <- factor(class_counts$Class, levels = rev(class_order))

suppfig_s4b <- ggplot(class_counts, aes(x = Class, y = Count, fill = Condition)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = Count), position = position_dodge(width = 0.8), hjust = -0.2, size = 3) +
  coord_flip() + scale_fill_manual(values = condition_colors) +
  labs(title = "B) Species by Lipid Class", x = NULL, y = "Number of Species") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

# S4D: By SuperClass
count_by_superclass <- function(lipids, condition_name, class_info) {
  tibble(Lipid = lipids) %>% left_join(class_info, by = "Lipid") %>%
    mutate(SuperClass = ifelse(is.na(SuperClass), "Unclassified", SuperClass)) %>%
    count(SuperClass, name = "Count") %>% mutate(Condition = condition_name)
}
super_ctrl <- count_by_superclass(ctrl_lipids, "Control", lipid_class_info)
super_low <- count_by_superclass(low_lipids, "LowInput", lipid_class_info)
super_counts <- bind_rows(super_ctrl, super_low)
super_order <- super_counts %>% group_by(SuperClass) %>% summarise(Total = sum(Count)) %>% arrange(desc(Total)) %>% pull(SuperClass)
super_counts$SuperClass <- factor(super_counts$SuperClass, levels = rev(super_order))

suppfig_s4d <- ggplot(super_counts, aes(x = SuperClass, y = Count, fill = Condition)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = Count), position = position_dodge(width = 0.8), hjust = -0.2, size = 3) +
  coord_flip() + scale_fill_manual(values = condition_colors) +
  labs(title = "D) Species by SuperClass", x = NULL, y = "Number of Species") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

# Combine
suppfig_s4_top <- suppfig_s4a + suppfig_s4b + plot_layout(widths = c(1, 1.3))
suppfig_s4_bottom <- suppfig_s4d + plot_spacer() + plot_layout(widths = c(1, 0.3))
suppfig_s4 <- suppfig_s4_top / suppfig_s4_bottom + plot_layout(heights = c(1, 1))

save_fig(suppfig_s4, "SuppFig_S4_Lipid_Species_Counts.png", w = 14, h = 12, subdir = "supp")
save_fig(suppfig_s4, "SuppFig_S4_Lipid_Species_Counts.pdf", w = 14, h = 12, subdir = "supp")

# Tables S6
species_summary <- tibble(Metric = c("Control Total", "LowInput Total", "Common", "Control Only", "LowInput Only", "All Unique"),
                          Count = c(n_ctrl, n_low, n_common, length(ctrl_only), length(low_only), length(unique(c(ctrl_lipids, low_lipids)))))
write.csv(species_summary, "table/supp/SuppTable_S6a_Species_Summary.csv", row.names = FALSE)

class_summary <- class_counts %>% pivot_wider(names_from = Condition, values_from = Count, values_fill = 0) %>% arrange(desc(Control + LowInput))
write.csv(class_summary, "table/supp/SuppTable_S6b_Species_by_Class.csv", row.names = FALSE)

super_summary <- super_counts %>% pivot_wider(names_from = Condition, values_from = Count, values_fill = 0) %>% arrange(desc(Control + LowInput))
write.csv(super_summary, "table/supp/SuppTable_S6c_Species_by_SuperClass.csv", row.names = FALSE)

message("Saved: table/supp/SuppTable_S6a,b,c")
message(" Supp Fig S4 (Species Counts) complete!")

# ═══════════════════════════════════════════════════════════════════════════════
# SUPP FIG S5: PCA OF LIPIDS
# ═══════════════════════════════════════════════════════════════════════════════
message("\n══════════════════════════════════════════════════════════════")
message("SUPPLEMENTARY FIGURE S5: PCA OF LIPIDS")
message("══════════════════════════════════════════════════════════════\n")

run_lipid_pca <- function(df_tic, pseudocount = 1e-6) {
  X <- df_tic %>% dplyr::select(-Compound_Name, -Condition) %>% as.matrix()
  P <- X / pmax(rowSums(X, na.rm = TRUE), 1e-12); P <- P + pseudocount; P <- P / rowSums(P, na.rm = TRUE)
  L <- log(P); X_clr <- L - rowMeans(L, na.rm = TRUE)
  V <- t(X_clr); V <- t(scale(t(V), center = TRUE, scale = TRUE)); V[is.na(V)] <- 0
  pc <- prcomp(V, center = FALSE, scale. = FALSE)
  out <- as.data.frame(pc$x[, 1:2]); colnames(out) <- c("PC1", "PC2")
  out$Lipid <- rownames(pc$x); out$Class <- get_lipid_class(out$Lipid); out$Group <- class_to_group(out$Class)
  ve <- (pc$sdev^2) / sum(pc$sdev^2); attr(out, "pve") <- ve[1:2]
  out
}

pca_control <- run_lipid_pca(control_tic)
pca_low <- run_lipid_pca(low_tic)

plot_lipid_pca <- function(pca_df, title = "", show_legend = TRUE) {
  pve <- attr(pca_df, "pve")
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Class)) +
    geom_point(size = 2.5, alpha = 0.8) +
    stat_ellipse(aes(group = Group, color = NULL), linetype = "dashed", linewidth = 0.8, color = "grey50", level = 0.95) +
    scale_color_manual(values = class_colors) +
    labs(title = title, x = sprintf("PC1 (%.1f%%)", 100 * pve[1]), y = sprintf("PC2 (%.1f%%)", 100 * pve[2])) +
    theme_minimal(base_size = 12) +
    theme(axis.text = element_text(color = "black"), axis.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = if (show_legend) "right" else "none",
          axis.line.x = element_line(color = "black", linewidth = 0.5),
          axis.line.y = element_line(color = "black", linewidth = 0.5))
}

suppfig_s5_control <- plot_lipid_pca(pca_control, "A) Control: Lipid PCA", show_legend = FALSE)
suppfig_s5_low <- plot_lipid_pca(pca_low, "B) LowInput: Lipid PCA", show_legend = TRUE)

suppfig_s5_pca <- suppfig_s5_control + suppfig_s5_low + plot_layout(guides = "collect") & theme(legend.position = "bottom")

save_fig(suppfig_s5_pca, "SuppFig_S5_PCA_Lipids.png", w = 14, h = 7, subdir = "supp")
save_fig(suppfig_s5_pca, "SuppFig_S5_PCA_Lipids.pdf", w = 14, h = 7, subdir = "supp")
message(" Supp Fig S5 (PCA) complete!")

# ═══════════════════════════════════════════════════════════════════════════════
# SUPP FIG S5A/B: SAMPLE PCA COLORED BY GEOLOCATION META
# ═══════════════════════════════════════════════════════════════════════════════
message("\n══════════════════════════════════════════════════════════════")
message("SUPPLEMENTARY FIGURE S5A/B: SAMPLE PCA WITH GEOLOCATION META")
message("══════════════════════════════════════════════════════════════\n")

run_sample_pca <- function(df_tic, pseudocount = 1e-6) {
  X <- df_tic %>% dplyr::select(-Compound_Name, -Condition) %>% as.matrix()
  P <- X / pmax(rowSums(X, na.rm = TRUE), 1e-12)
  P <- P + pseudocount
  P <- P / rowSums(P, na.rm = TRUE)
  L <- log(P)
  X_clr <- L - rowMeans(L, na.rm = TRUE)
  V <- scale(X_clr, center = TRUE, scale = TRUE)
  V[is.na(V)] <- 0
  pc <- prcomp(V, center = FALSE, scale. = FALSE)
  out <- as.data.frame(pc$x[, 1:4, drop = FALSE])
  colnames(out) <- c("PC1", "PC2", "PC3", "PC4")
  out$Compound_Name <- df_tic$Compound_Name
  out$Condition <- df_tic$Condition
  ve <- (pc$sdev^2) / sum(pc$sdev^2)
  attr(out, "pve") <- ve
  out
}

geoloc_candidates <- c("data/SAP_geoloc.csv", "/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database/data/SAP_geoloc.csv")
geoloc_file <- geoloc_candidates[file.exists(geoloc_candidates)][1]

if (!is.na(geoloc_file)) {
  geo <- vroom(geoloc_file, show_col_types = FALSE)
  names(geo)[1] <- "Taxa"
  geo <- geo %>%
    transmute(
      Compound_Name = as.character(Taxa),
      K.Cluster = as.character(`K.Cluster`),
      Original_Race = as.character(Original_Race)
    ) %>%
    mutate(
      K.Cluster = ifelse(is.na(K.Cluster) | trimws(K.Cluster) == "", "Unknown", paste0("K", K.Cluster)),
      Original_Race = ifelse(is.na(Original_Race) | trimws(Original_Race) == "", "Unknown", Original_Race)
    )

  pca_samp_control <- run_sample_pca(control_tic) %>% left_join(geo, by = "Compound_Name")
  pca_samp_low <- run_sample_pca(low_tic) %>% left_join(geo, by = "Compound_Name")

  kcluster_levels <- sort(unique(c(pca_samp_control$K.Cluster, pca_samp_low$K.Cluster)))
  kcluster_levels <- c(setdiff(kcluster_levels, "Unknown"), intersect("Unknown", kcluster_levels))
  kcluster_palette <- setNames(c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "grey60")[seq_along(kcluster_levels)], kcluster_levels)

  race_levels <- sort(unique(c(pca_samp_control$Original_Race, pca_samp_low$Original_Race)))
  race_levels <- c(setdiff(race_levels, "Unknown"), intersect("Unknown", race_levels))
  race_palette <- setNames(c(scales::hue_pal(l = 60, c = 100)(max(length(race_levels) - 1, 1)), "grey60")[seq_along(race_levels)], race_levels)

  plot_sample_pca_meta <- function(df, color_var, title, palette, x_pc = 2, y_pc = 3, legend_ncol = 2, show_legend = TRUE) {
    pve <- attr(df, "pve")
    x_name <- paste0("PC", x_pc)
    y_name <- paste0("PC", y_pc)
    ggplot(df, aes(x = PC1, y = PC2, color = .data[[color_var]])) +
      aes(x = .data[[x_name]], y = .data[[y_name]]) +
      geom_point(size = 2.2, alpha = 0.85) +
      scale_color_manual(values = palette, drop = FALSE) +
      labs(
        title = title,
        x = sprintf("%s (%.1f%%)", x_name, 100 * pve[x_pc]),
        y = sprintf("%s (%.1f%%)", y_name, 100 * pve[y_pc]),
        color = color_var
      ) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = if (show_legend) "right" else "none",
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size = 8),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5)
      ) +
      guides(color = guide_legend(ncol = legend_ncol, override.aes = list(size = 3, alpha = 1)))
  }

  # K.Cluster: show requested axis pairs (PC2 vs PC3, PC3 vs PC4)
  s5a_a <- plot_sample_pca_meta(
    pca_samp_control, "K.Cluster", "A) Control: PC2 vs PC3 (K.Cluster)",
    palette = kcluster_palette, x_pc = 2, y_pc = 3, legend_ncol = 1, show_legend = FALSE
  )
  s5a_b <- plot_sample_pca_meta(
    pca_samp_low, "K.Cluster", "B) LowInput: PC2 vs PC3 (K.Cluster)",
    palette = kcluster_palette, x_pc = 2, y_pc = 3, legend_ncol = 1, show_legend = TRUE
  )
  s5a_c <- plot_sample_pca_meta(
    pca_samp_control, "K.Cluster", "C) Control: PC3 vs PC4 (K.Cluster)",
    palette = kcluster_palette, x_pc = 3, y_pc = 4, legend_ncol = 1, show_legend = FALSE
  )
  s5a_d <- plot_sample_pca_meta(
    pca_samp_low, "K.Cluster", "D) LowInput: PC3 vs PC4 (K.Cluster)",
    palette = kcluster_palette, x_pc = 3, y_pc = 4, legend_ncol = 1, show_legend = TRUE
  )
  suppfig_s5a_cluster <- (s5a_a + s5a_b) / (s5a_c + s5a_d) + plot_layout(guides = "collect") & theme(legend.position = "right")

  # Original_Race: same requested axis pairs
  s5b_a <- plot_sample_pca_meta(
    pca_samp_control, "Original_Race", "A) Control: PC2 vs PC3 (Original_Race)",
    palette = race_palette, x_pc = 2, y_pc = 3, legend_ncol = 2, show_legend = FALSE
  )
  s5b_b <- plot_sample_pca_meta(
    pca_samp_low, "Original_Race", "B) LowInput: PC2 vs PC3 (Original_Race)",
    palette = race_palette, x_pc = 2, y_pc = 3, legend_ncol = 2, show_legend = TRUE
  )
  s5b_c <- plot_sample_pca_meta(
    pca_samp_control, "Original_Race", "C) Control: PC3 vs PC4 (Original_Race)",
    palette = race_palette, x_pc = 3, y_pc = 4, legend_ncol = 2, show_legend = FALSE
  )
  s5b_d <- plot_sample_pca_meta(
    pca_samp_low, "Original_Race", "D) LowInput: PC3 vs PC4 (Original_Race)",
    palette = race_palette, x_pc = 3, y_pc = 4, legend_ncol = 2, show_legend = TRUE
  )
  suppfig_s5b_race <- (s5b_a + s5b_b) / (s5b_c + s5b_d) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

  save_fig(suppfig_s5a_cluster, "SuppFig_S5A_PCA_by_KCluster.png", w = 14, h = 12, subdir = "supp")
  save_fig(suppfig_s5a_cluster, "SuppFig_S5A_PCA_by_KCluster.pdf", w = 14, h = 12, subdir = "supp")
  save_fig(suppfig_s5b_race, "SuppFig_S5B_PCA_by_Original_Race.png", w = 14, h = 13, subdir = "supp")
  save_fig(suppfig_s5b_race, "SuppFig_S5B_PCA_by_Original_Race.pdf", w = 14, h = 13, subdir = "supp")

  # Optional interactive 3D PCA (PC1-PC2-PC3), useful when 2D panels show weak separation
  if (requireNamespace("plotly", quietly = TRUE) && requireNamespace("htmlwidgets", quietly = TRUE)) {
    pca_samp_all <- bind_rows(pca_samp_control, pca_samp_low)
    pve3d <- attr(pca_samp_control, "pve")

    make_plot3d <- function(df, color_var, title_text, palette_named) {
      df_local <- df
      df_local[[color_var]] <- factor(df_local[[color_var]], levels = names(palette_named))
      plotly::plot_ly(
        data = df_local,
        x = ~PC1, y = ~PC2, z = ~PC3,
        color = as.formula(paste0("~`", color_var, "`")),
        colors = unname(palette_named),
        symbol = ~Condition,
        symbols = c("circle", "triangle-up"),
        type = "scatter3d",
        mode = "markers",
        marker = list(size = 3.5, opacity = 0.8),
        text = ~paste0(
          "Sample: ", Compound_Name,
          "<br>Condition: ", Condition,
          "<br>", color_var, ": ", .data[[color_var]]
        ),
        hoverinfo = "text"
      ) %>%
        plotly::layout(
          title = list(text = title_text),
          scene = list(
            xaxis = list(title = sprintf("PC1 (%.1f%%)", 100 * pve3d[1])),
            yaxis = list(title = sprintf("PC2 (%.1f%%)", 100 * pve3d[2])),
            zaxis = list(title = sprintf("PC3 (%.1f%%)", 100 * pve3d[3]))
          ),
          legend = list(orientation = "v")
        )
    }

    fig_s5c_3d_cluster <- make_plot3d(
      pca_samp_all, "K.Cluster",
      "SuppFig S5C: 3D Sample PCA (PC1-PC2-PC3) colored by K.Cluster",
      palette_named = kcluster_palette
    )
    fig_s5d_3d_race <- make_plot3d(
      pca_samp_all, "Original_Race",
      "SuppFig S5D: 3D Sample PCA (PC1-PC2-PC3) colored by Original_Race",
      palette_named = race_palette
    )

    htmlwidgets::saveWidget(
      fig_s5c_3d_cluster,
      file = "fig/supp/SuppFig_S5C_PCA3D_by_KCluster.html",
      selfcontained = TRUE
    )
    htmlwidgets::saveWidget(
      fig_s5d_3d_race,
      file = "fig/supp/SuppFig_S5D_PCA3D_by_Original_Race.html",
      selfcontained = TRUE
    )
    message("Saved: fig/supp/SuppFig_S5C_PCA3D_by_KCluster.html")
    message("Saved: fig/supp/SuppFig_S5D_PCA3D_by_Original_Race.html")
  } else {
    message("Skipping 3D PCA html: plotly/htmlwidgets not installed.")
  }
  message(" Supp Fig S5A/B (sample PCA by geolocation metadata) complete!")
} else {
  message("Skipping geolocation PCA plots: SAP_geoloc.csv not found.")
}

# ═══════════════════════════════════════════════════════════════════════════════
message("\n══════════════════════════════════════════════════════════════")
message("ALL DONE! Generated:")
message("  - fig/main/Figure1_Lipidomics_Landscape.png/pdf")
message("  - fig/supp/SuppFig_S3_Compositional_Contrasts.png/pdf")
message("  - fig/supp/SuppFig_S8_NonFocused_Lipid_Class_Context.png")
message("  - fig/supp/SuppFig_S4_Lipid_Species_Counts.png/pdf")
message("  - fig/supp/SuppFig_S5_PCA_Lipids.png/pdf")
message("  - fig/supp/SuppFig_S5A_PCA_by_KCluster.png/pdf")
message("  - fig/supp/SuppFig_S5B_PCA_by_Original_Race.png/pdf")
message("  - fig/supp/SuppFig_S5C_PCA3D_by_KCluster.html")
message("  - fig/supp/SuppFig_S5D_PCA3D_by_Original_Race.html")
message("  - table/supp/SuppTable_S5A_Class_CLR_Contrast.csv")
message("  - table/supp/SuppTable_S5B_Class_ALR_Contrast.csv")
message("  - table/supp/NonFocused_Lipid_Class_CLR_Contrast.csv")
message("  - table/supp/SuppTable_S6a_Species_Summary.csv")
message("  - table/supp/SuppTable_S6b_Species_by_Class.csv")
message("  - table/supp/SuppTable_S6c_Species_by_SuperClass.csv")
message("══════════════════════════════════════════════════════════════\n")
