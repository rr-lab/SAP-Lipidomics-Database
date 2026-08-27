# ==============================================================================
# Figure 1 -- Lipid-class variation is largely independent of botanical race and
#             of genome-wide genetic cluster, except for PS.
#
#   A,B  variance explained by each grouping, for all 14 traits, per trial
#   C,D  PS abundance by genetic cluster in each trial
#
# Everything is computed WITHIN a trial. No CTL-LIN contrast, no compositional
# (CLR/ALR) transform. The class table, the race and cluster join, and the
# Kruskal-Wallis helper come from _common.R so this figure and its PCA
# supplement cannot drift apart.
#
# plot_theme is used unmodified. The one deviation is that panel B's legend is
# switched off, since it would only repeat panel A's.
#
# Inputs
#   data/SPATS_fitted/non_normalized_intensities/Final_subset_{control,lowinput}_*.csv
#   data/SAP_geoloc.csv
#
# Outputs
#   fig/main/Figure1_Population_Structure.png
#   table/supp/SuppTable_S25_Population_Structure_Lipid_Tests.csv
# ==============================================================================
source("scripts/new_new_script/_common.R")
suppressPackageStartupMessages({ library(tibble); library(tidyr) })

dat    <- population_table()
TRAITS <- c("TotalLipid", CLASS_ORDER)
trait_column <- function(tr) if (tr == "TotalLipid") "TotalLipid_log10" else tr

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
  bind_rows(rows) %>% group_by(Condition, Grouping) %>%
    mutate(q_BH = p.adjust(p, method = "BH")) %>% ungroup()
}

# "K.Cluster" is the label the published supplementary table already uses.
tests <- bind_rows(run_tests(dat, "RaceGroup", TRAITS, "RaceGroup"),
                   run_tests(dat, "KCluster",  TRAITS, "K.Cluster")) %>%
  mutate(Condition = factor(Condition, c("CTL", "LIN")),
         sig = !is.na(q_BH) & q_BH < 0.05)

save_table(tests, "SuppTable_S25_Population_Structure_Lipid_Tests.csv")

ord <- tests %>% group_by(Trait) %>%
  summarise(m = max(epsilon2, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(m)) %>% pull(Trait)
tests <- tests %>% mutate(Trait = factor(Trait, levels = rev(ord)))

# ---- A,B ---------------------------------------------------------------------
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

# ---- C,D ---------------------------------------------------------------------
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

fig1 <- (eps2_panel("RaceGroup", TRUE) | eps2_panel("K.Cluster", FALSE)) /
        (ps_panel("CTL") | ps_panel("LIN")) +
  plot_layout(heights = c(1.25, 1)) +
  plot_annotation(tag_levels = "A", theme = TAG_THEME)

save_fig(fig1, "Figure1_Population_Structure.png", width = 18, height = 14)

# ---- console summary ---------------------------------------------------------
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
