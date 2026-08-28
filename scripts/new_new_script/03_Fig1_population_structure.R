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
#   table/supp/SuppTable_S25a_Group_Composition.csv
#     -- accessions per race group and per genetic cluster, per trial
#   table/supp/SuppTable_S25_Population_Structure_Lipid_Tests.csv
#     -- one row per trait x grouping x trial, carrying n, H, p, epsilon^2, the
#        BH q, which group is highest and lowest, and the median for every group
# ==============================================================================
source("scripts/new_new_script/_common.R")
suppressPackageStartupMessages({ library(tibble); library(tidyr) })

MIN_PER_GROUP <- 5          # the threshold kw_eps2() applies
dat    <- population_table()
TRAITS <- c("TotalLipid", CLASS_ORDER)
trait_column <- function(tr) if (tr == "TotalLipid") "TotalLipid_log10" else tr

# The per-group medians travel with the test, so the table alone shows which
# group is high and which is low. That is what makes the per-race boxplot grids
# unnecessary: the same numbers, without 28 panels of nothing.
group_medians <- function(values, groups, min_per_group = 5) {
  ok <- !is.na(values) & !is.na(groups)
  values <- values[ok]; groups <- droplevels(factor(groups[ok]))
  big <- names(which(table(groups) >= min_per_group))
  ok <- groups %in% big
  values <- values[ok]; groups <- droplevels(factor(groups[ok]))
  tapply(values, groups, median)
}

run_tests <- function(df, group_col, traits, label = group_col) {
  rows <- list()
  for (cond in levels(df$Condition)) {
    sub <- df %>% filter(Condition == cond, !is.na(.data[[group_col]]))
    for (tr in traits) {
      v <- sub[[trait_column(tr)]]
      g <- sub[[group_col]]
      r <- kw_eps2(v, g)
      if (is.null(r)) next
      med <- group_medians(v, g)
      row <- bind_cols(tibble(Condition = cond, Grouping = label, Trait = tr), r) %>%
        mutate(highest = names(med)[which.max(med)],
               lowest  = names(med)[which.min(med)])
      for (nm in names(med)) row[[paste0("median_", nm)]] <- unname(med[[nm]])
      rows[[length(rows) + 1L]] <- row
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

# ---- group composition -------------------------------------------------------
# How many accessions sit in each group, per trial. Needed to read everything
# else here: Kafir has only seven accessions in the whole passport table, so a
# race-level result rests on very different sample sizes from group to group.
# The raw passport designations folded into each race group are listed too,
# since "Mixed" is eleven different hyphenated labels rather than a category.
race_map <- dat %>%
  filter(!is.na(RaceGroup)) %>%
  distinct(LineRaw, Original_Race, RaceGroup) %>%
  group_by(RaceGroup) %>%
  summarise(Folded_labels = paste(sort(unique(Original_Race)), collapse = "; "),
            .groups = "drop")

composition <- bind_rows(
  dat %>% filter(!is.na(RaceGroup)) %>%
    count(Condition, Group = as.character(RaceGroup)) %>%
    mutate(Grouping = "RaceGroup"),
  dat %>% filter(!is.na(KCluster)) %>%
    count(Condition, Group = as.character(KCluster)) %>%
    mutate(Grouping = "K.Cluster")) %>%
  pivot_wider(names_from = Condition, values_from = n, names_prefix = "n_",
              values_fill = 0) %>%
  left_join(race_map, by = c("Group" = "RaceGroup")) %>%
  mutate(Folded_labels = ifelse(is.na(Folded_labels), Group, Folded_labels),
         tested_CTL = n_CTL >= MIN_PER_GROUP,
         tested_LIN = n_LIN >= MIN_PER_GROUP) %>%
  dplyr::select(Grouping, Group, n_CTL, n_LIN, tested_CTL, tested_LIN, Folded_labels) %>%
  arrange(Grouping, Group)

save_table(composition, "SuppTable_S25a_Group_Composition.csv")

cat("\n-- group composition -------------------------------------------------\n")
print(as.data.frame(composition %>% dplyr::select(-Folded_labels)))
cat(sprintf("\n  accessions with no usable race assignment: CTL %d, LIN %d\n",
            sum(dat$Condition == "CTL" & is.na(dat$RaceGroup)),
            sum(dat$Condition == "LIN" & is.na(dat$RaceGroup))))

save_table(tests %>% dplyr::select(Condition, Grouping, Trait, n, k_groups, H, p,
                                   epsilon2, q_BH, highest, lowest,
                                   dplyr::starts_with("median_")),
           "SuppTable_S25_Population_Structure_Lipid_Tests.csv")

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
