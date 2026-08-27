# ==============================================================================
# Supplementary Figure S3 -- Class-composition PCA of each field trial,
#                            coloured by botanical race and by genetic cluster.
#
#   A,B  CTL and LIN coloured by botanical race
#   C,D  the same two ordinations coloured by genetic cluster (K = 6)
#
# The point of the figure is that there is nothing to see: every grouping lands
# in one overlapping cloud. Each panel carries the variance in PC1 and PC2 that
# its grouping explains, by the same Kruskal-Wallis statistic Figure 1 uses, so
# "no clustering" is backed by a number rather than by eye.
#
# A PS-gradient overlay was tried and dropped. In CTL, PS tracks overall class
# composition (Spearman rho 0.56 on PC1), which is a separate observation from
# the population-structure question this figure answers.
#
# plot_theme is used unmodified. The one deviation is that the repeated legend
# is switched off in the right-hand panel of each row.
#
# Inputs
#   data/SPATS_fitted/non_normalized_intensities/Final_subset_{control,lowinput}_*.csv
#   data/SAP_geoloc.csv
#
# Outputs
#   fig/supp/SuppFig_S3_Class_PCA_Structure.png
#   table/supp/SuppTable_S25b_Population_Structure_PC_Tests.csv
# ==============================================================================
source("scripts/new_new_script/_common.R")
suppressPackageStartupMessages({ library(tibble); library(tidyr) })

dat <- population_table()
pca <- list(CTL = trial_class_pca(dat, "CTL"),
            LIN = trial_class_pca(dat, "LIN"))

# ---- how much of PC1 and PC2 does each grouping explain? ---------------------
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

# ---- panels ------------------------------------------------------------------
pca_panel <- function(cond, group_col, palette, show_legend = TRUE) {
  s  <- pca[[cond]]$scores %>% filter(!is.na(.data[[group_col]]))
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

cat("\n-- grouping vs PC1/PC2 ----------------------------------------------\n")
print(pc_tests %>% dplyr::select(Condition, Grouping, Axis, n, epsilon2, q_BH))
for (cond in names(pca)) {
  ve <- pca[[cond]]$ve
  cat(sprintf("  %s  PC1 %.1f%%  PC2 %.1f%%\n", cond, 100 * ve[1], 100 * ve[2]))
}
