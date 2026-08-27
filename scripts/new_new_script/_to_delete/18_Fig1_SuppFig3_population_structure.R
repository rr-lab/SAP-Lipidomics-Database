# ==============================================================================
# Figure 1        -- Lipid-class variation is largely independent of botanical
#                    race and of genome-wide genetic cluster, except for PS.
# Supp Figure S3  -- Class-composition PCA of each trial, coloured by botanical
#                    race and by genetic cluster. Four panels; a PS-gradient
#                    overlay was dropped because in CTL it tracks overall class
#                    composition (Spearman rho 0.56 on PC1), which is a separate
#                    observation from the population-structure question here.
# Supp Table      -- Kruskal-Wallis tests of every class against both groupings.
#
# Everything here is computed WITHIN a trial. There is no CTL-LIN contrast and
# no compositional (CLR/ALR) transform anywhere in this script.
#
# plot_theme is used unmodified. The only deviations, both to stop a legend
# landing on top of the data, are that the repeated legend is switched off in
# the second panel of each pair, and that the PCA panels keep the legend in the
# theme's own inside-panel position only where there is empty space for it.
#
# The published test table was produced by a Python script that read a temporary
# join file which no longer exists, so the tests are rebuilt here from repository
# data. Numbers printed at the end should be checked against the manuscript.
#
# Inputs
#   data/SPATS_fitted/non_normalized_intensities/Final_subset_{control,lowinput}_*.csv
#   data/SAP_geoloc.csv                       (Taxa, K.Cluster, Original_Race)
#
# Outputs
#   fig/main/Figure1_Population_Structure.png
#   fig/supp/SuppFig_S3_Class_PCA_Structure.png
#   table/supp/SuppTable_S25_Population_Structure_Lipid_Tests.csv
#   table/supp/SuppTable_S25b_Population_Structure_PC_Tests.csv
# ==============================================================================
source("scripts/new_new_script/_common.R")
suppressPackageStartupMessages({ library(tibble); library(tidyr) })

geoloc_csv <- Sys.getenv("SAP_GEOLOC", file.path(DATA_ROOT, "SAP_geoloc.csv"))
stopifnot(file.exists(CTL_CSV), file.exists(LIN_CSV), file.exists(geoloc_csv))

PURE        <- c("Bicolor", "Caudatum", "Durra", "Guinea", "Kafir")
RACE_ORDER  <- c(PURE, "Mixed")
CLASS_ORDER <- names(class_colors)
MIN_PER_GROUP <- 5

race_colors <- c(Bicolor = "#E69F00", Caudatum = "#0072B2", Durra = "#009E73",
                 Guinea = "#CC79A7", Kafir = "#D55E00", Mixed = "#999999")
k_colors <- c("1" = "#E69F00", "2" = "#56B4E9", "3" = "#009E73",
              "4" = "#F0E442", "5" = "#0072B2", "6" = "#CC79A7")

# ---- class composition as % of total ion current -----------------------------
# TotalLipid is the summed signal over every annotated feature, so the 13 focal
# classes do not add up to 100%; the remainder is sterols, carotenoids and the
# other non-focal categories.
class_percent_tic <- function(path, condition) {
  x <- read_trial(path)
  geno <- x[[1]]
  m <- as.matrix(x[, -1, drop = FALSE])
  storage.mode(m) <- "numeric"
  m[!is.finite(m)] <- 0
  cls <- lipid_class(colnames(m))
  total <- rowSums(m)
  out <- tibble(LineRaw = geno, Condition = condition, TotalLipid = total)
  for (cl in CLASS_ORDER) {
    j <- which(cls == cl)
    out[[cl]] <- if (length(j)) 100 * rowSums(m[, j, drop = FALSE]) / total else 0
  }
  out
}

race_group <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x %in% c("", "NA")] <- NA_character_
  x[grepl("verticilliflorum", x, ignore.case = TRUE)] <- NA_character_
  ifelse(is.na(x), NA_character_, ifelse(x %in% PURE, x, "Mixed"))
}

geo <- vroom(geoloc_csv, show_col_types = FALSE)
names(geo)[1] <- "Taxa"                       # the file carries a UTF-8 BOM
geo <- geo %>%
  dplyr::select(Taxa, K.Cluster, Original_Race) %>%
  distinct(Taxa, .keep_all = TRUE)

dat <- bind_rows(class_percent_tic(CTL_CSV, "CTL"),
                 class_percent_tic(LIN_CSV, "LIN")) %>%
  left_join(geo, by = c("LineRaw" = "Taxa")) %>%
  mutate(Condition        = factor(Condition, c("CTL", "LIN")),
         RaceGroup        = factor(race_group(Original_Race), RACE_ORDER),
         KCluster         = factor(K.Cluster, levels = 1:6),
         TotalLipid_log10 = log10(TotalLipid))

TRAITS <- c("TotalLipid", CLASS_ORDER)
trait_column <- function(tr) if (tr == "TotalLipid") "TotalLipid_log10" else tr

# ---- Kruskal-Wallis with epsilon-squared -------------------------------------
# epsilon^2 = (H - k + 1) / (n - k), the proportion of rank variance explained,
# truncated at zero. Groups with fewer than five accessions are dropped first.
kw_eps2 <- function(values, groups) {
  ok <- !is.na(values) & !is.na(groups)
  values <- values[ok]; groups <- droplevels(factor(groups[ok]))
  big <- names(which(table(groups) >= MIN_PER_GROUP))
  ok <- groups %in% big
  values <- values[ok]; groups <- droplevels(factor(groups[ok]))
  k <- nlevels(groups); n <- length(values)
  if (k < 3L) return(NULL)
  kt <- kruskal.test(values, groups)
  H <- unname(kt$statistic)
  tibble(n = n, k_groups = k, H = H, p = kt$p.value,
         epsilon2 = max(0, (H - k + 1) / (n - k)))
}

run_tests <- function(df, group_col, traits, label = group_col) {
  rows <- list()
  for (cond in levels(df$Condition)) {
    sub <- df %>% filter(Condition == cond, !is.na(.data[[group_col]]))
    for (tr in traits) {
      r <- kw_eps2(sub[[trait_column(tr)]], sub[[group_col]])
      if (is.null(r)) next
      rows[[length(rows) + 1L]] <- bind_cols(
        tibble(Condition = cond, Grouping = label, Trait = tr), r)
    }
  }
  bind_rows(rows) %>%
    group_by(Condition, Grouping) %>%
    mutate(q_BH = p.adjust(p, method = "BH")) %>%
    ungroup()
}

# "K.Cluster" is the label the published supplementary table already uses.
tests <- bind_rows(run_tests(dat, "RaceGroup", TRAITS, "RaceGroup"),
                   run_tests(dat, "KCluster", TRAITS, "K.Cluster")) %>%
  mutate(Condition = factor(Condition, c("CTL", "LIN")),
         sig = !is.na(q_BH) & q_BH < 0.05)

save_table(tests, "SuppTable_S25_Population_Structure_Lipid_Tests.csv")

# ---- Figure 1 ----------------------------------------------------------------
ord <- tests %>% group_by(Trait) %>%
  summarise(m = max(epsilon2, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(m)) %>% pull(Trait)
tests <- tests %>% mutate(Trait = factor(Trait, levels = rev(ord)))

eps2_panel <- function(grouping, show_legend = TRUE) {
  d <- tests %>% filter(Grouping == grouping)
  nsig <- sum(d$sig, na.rm = TRUE)
  ggplot(d, aes(Trait, epsilon2, fill = Condition)) +
    geom_col(position = position_dodge(.78), width = .7,
             colour = "black", linewidth = .25) +
    geom_text(data = subset(d, sig), aes(label = "***"),
              position = position_dodge(.78), hjust = -.2, size = 5) +
    coord_flip() +
    scale_fill_manual(values = condition_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, .12))) +
    annotate("text", x = 1, y = max(d$epsilon2, na.rm = TRUE),
             label = sprintf("%d of %d tests significant after BH correction",
                             nsig, nrow(d)),
             hjust = 1, vjust = -.2, size = 4.6, fontface = "italic",
             colour = "grey25") +
    labs(x = NULL, y = expression(bold(paste("Variance explained (", epsilon^2, ")")))) +
    plot_theme + (if (show_legend) NULL else guides(fill = "none"))
}

# PS is the single trait associated with genetic cluster, and only under LIN.
# Each trial keeps its own y-scale: the comparison being made is between
# clusters within a trial, never between trials.
ps_panel <- function(cond) {
  d <- dat %>% filter(Condition == cond, !is.na(KCluster))
  ps_lim <- quantile(d$PS, c(0, .99), na.rm = TRUE)
  s <- tests %>% filter(Grouping == "K.Cluster", Trait == "PS", Condition == cond)
  lab <- sprintf("epsilon^2 == %.3f ~ ~ italic(q) == %s",
                 s$epsilon2[1], format(signif(s$q_BH[1], 3), scientific = TRUE))
  ggplot(d, aes(KCluster, PS, fill = KCluster)) +
    geom_boxplot(width = .62, outlier.shape = NA, colour = "black", linewidth = .3) +
    geom_jitter(width = .14, size = .55, alpha = .35, colour = "grey20") +
    scale_fill_manual(values = k_colors, guide = "none") +
    coord_cartesian(ylim = c(ps_lim[1], ps_lim[2] * 1.18)) +
    annotate("text", x = 6.4, y = ps_lim[2] * 1.14, label = lab, parse = TRUE,
             hjust = 1, size = 4.6, colour = "grey25") +
    labs(x = sprintf("Genetic cluster, %s", cond), y = "PS (% TIC)") +
    plot_theme
}

# The Condition legend sits inside panel A, where the theme puts it; panel B
# would only repeat it.
fig1 <- (eps2_panel("RaceGroup", TRUE) | eps2_panel("K.Cluster", FALSE)) /
        (ps_panel("CTL") | ps_panel("LIN")) +
  plot_layout(heights = c(1.25, 1)) +
  plot_annotation(tag_levels = "A", theme = TAG_THEME)

save_fig(fig1, "Figure1_Population_Structure.png", width = 18, height = 14)

# ---- Supplementary Figure S3 -------------------------------------------------
# PCA of the 13 class percentages within each trial, on log10 values, scaled so
# that no single dominant class (MGDG, PC) sets the axes on its own.
trial_pca <- function(cond) {
  sub <- dat %>% filter(Condition == cond)
  m <- as.matrix(sub[, CLASS_ORDER, drop = FALSE])
  m <- log10(pmax(m, 0) + 1e-4)
  keep <- apply(m, 2, sd, na.rm = TRUE) > 0
  p <- prcomp(m[, keep, drop = FALSE], center = TRUE, scale. = TRUE)
  ve <- p$sdev^2 / sum(p$sdev^2)
  list(scores = sub %>% mutate(PC1 = p$x[, 1], PC2 = p$x[, 2]), ve = ve)
}
pca <- list(CTL = trial_pca("CTL"), LIN = trial_pca("LIN"))

# How much of PC1 and PC2 does each grouping explain? Same statistic as Figure 1,
# so "no visible clustering" is backed by a number rather than by eye.
pc_tests <- bind_rows(lapply(names(pca), function(cond) {
  s <- pca[[cond]]$scores
  bind_rows(lapply(c("RaceGroup", "KCluster"), function(g) {
    lab <- if (g == "KCluster") "K.Cluster" else g
    bind_rows(lapply(c("PC1", "PC2"), function(ax) {
      r <- kw_eps2(s[[ax]], s[[g]])
      if (is.null(r)) return(NULL)
      bind_cols(tibble(Condition = cond, Grouping = lab, Axis = ax), r)
    }))
  }))
})) %>%
  group_by(Condition, Grouping) %>%
  mutate(q_BH = p.adjust(p, method = "BH")) %>% ungroup()

save_table(pc_tests, "SuppTable_S25b_Population_Structure_PC_Tests.csv")

pca_panel <- function(cond, group_col, palette, show_legend = TRUE) {
  s <- pca[[cond]]$scores %>% filter(!is.na(.data[[group_col]]))
  ve <- pca[[cond]]$ve
  grouping_label <- if (group_col == "KCluster") "K.Cluster" else group_col
  e <- pc_tests %>% filter(Condition == cond, Grouping == grouping_label)
  # plotmath needs the axis names quoted and joined with ~, otherwise
  # PC1 and epsilon sit next to each other as two bare symbols and parse() fails.
  lab <- sprintf('list("PC1" ~ epsilon^2 == %.3f, "PC2" ~ epsilon^2 == %.3f)',
                 e$epsilon2[e$Axis == "PC1"], e$epsilon2[e$Axis == "PC2"])
  ggplot(s, aes(PC1, PC2, colour = .data[[group_col]])) +
    geom_point(size = 2.6, alpha = .8) +
    scale_colour_manual(values = palette, drop = FALSE) +
    annotate("text", x = max(s$PC1), y = min(s$PC2), label = lab, parse = TRUE,
             hjust = 1, vjust = 0, size = 4.6, colour = "grey25") +
    labs(x = sprintf("PC1 (%.1f%%), %s", 100 * ve[1], cond),
         y = sprintf("PC2 (%.1f%%)", 100 * ve[2])) +
    plot_theme +
    (if (show_legend) guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1)))
     else guides(colour = "none"))
}

figs3 <- (pca_panel("CTL", "RaceGroup", race_colors, TRUE) |
          pca_panel("LIN", "RaceGroup", race_colors, FALSE)) /
         (pca_panel("CTL", "KCluster", k_colors, TRUE) |
          pca_panel("LIN", "KCluster", k_colors, FALSE)) +
  plot_annotation(tag_levels = "A", theme = TAG_THEME)

save_fig(figs3, "SuppFig_S3_Class_PCA_Structure.png", width = 17, height = 14, subdir = "supp")

# ---- console summary for the manuscript --------------------------------------
cat("\n-- Kruskal-Wallis on class composition ------------------------------\n")
for (g in c("RaceGroup", "K.Cluster")) {
  for (cond in c("CTL", "LIN")) {
    d <- tests %>% filter(Grouping == g, Condition == cond)
    cat(sprintf("  %-10s %s  n=%d  traits=%d  max eps2=%.4f  min q=%.4g  sig=%d\n",
                g, cond, max(d$n), nrow(d), max(d$epsilon2), min(d$q_BH), sum(d$sig)))
  }
}
cat("\n-- significant tests -------------------------------------------------\n")
print(tests %>% filter(sig) %>%
        dplyr::select(Condition, Grouping, Trait, n, epsilon2, p, q_BH))
cat("\n-- PS by genetic cluster --------------------------------------------\n")
print(dat %>% filter(!is.na(KCluster)) %>% group_by(Condition, KCluster) %>%
        summarise(n = dplyr::n(), median_PS = median(PS), .groups = "drop"))
cat("\n-- grouping vs PC1/PC2 ----------------------------------------------\n")
print(pc_tests %>% dplyr::select(Condition, Grouping, Axis, n, epsilon2, q_BH))
