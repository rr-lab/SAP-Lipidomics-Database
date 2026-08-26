# ==============================================================================
# Figure 1 -- Lipid class variation is not explained by botanical race or by
#             genome-wide genetic cluster.
#
# Two panels, one per grouping, each showing variance explained (epsilon^2) for
# total lipid signal and the 13 lipid classes, separately for each trial.
# Both tests are run WITHIN a trial, so nothing here is a cross-trial contrast.
#
# The former panels A and B (top-variance species per trial) were removed: the
# "seven of ten shared" observation is a cross-trial comparison, and the
# high-variance ranking itself was judged not to carry biology.
#
# Input : table/race/race_structure_lipid_tests.csv
# Output: fig/main/Figure1_Race_Structure.png
# ==============================================================================
source("scripts/new_script/_common.R")

race_csv <- Sys.getenv("RACE_CSV", file.path(REPO, "table/race/race_structure_lipid_tests.csv"))
if (!file.exists(race_csv)) stop("Race test table not found: ", race_csv)

race <- vroom(race_csv, show_col_types = FALSE) %>%
  mutate(Condition = recode(Condition, Control = "CTL", LowInput = "LIN",
                            CTL = "CTL", LIN = "LIN"),
         Condition = factor(Condition, c("CTL", "LIN")),
         sig = !is.na(q_BH) & q_BH < 0.05)

ord <- race %>% group_by(Trait) %>%
  summarise(m = max(epsilon2, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(m)) %>% pull(Trait)
race <- race %>% mutate(Trait = factor(Trait, levels = rev(ord)))

grp_panel <- function(g) {
  d <- race %>% filter(Grouping == g)
  nsig <- sum(d$sig, na.rm = TRUE)
  ggplot(d, aes(Trait, epsilon2, fill = Condition)) +
    geom_col(position = position_dodge(.78), width = .7,
             colour = "black", linewidth = .25) +
    geom_text(data = subset(d, sig), aes(label = "***"),
              position = position_dodge(.78), hjust = -.2, size = 3.2) +
    coord_flip() +
    scale_fill_manual(values = condition_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, .12))) +
    annotate("text", x = 1, y = max(d$epsilon2, na.rm = TRUE),
             label = sprintf("%d trait%s significant after BH correction",
                             nsig, ifelse(nsig == 1, "", "s")),
             hjust = 1, vjust = -.2, size = 2.9, fontface = "italic", colour = "grey25") +
    labs(x = NULL, y = expression(bold(paste("Variance explained (", epsilon^2, ")")))) +
    plot_theme + theme(panel.grid.major.y = element_blank(), legend.position = "right")
}

groupings <- unique(race$Grouping)
pA <- grp_panel(groupings[1])
pB <- grp_panel(groupings[2])

fig <- (pA | pB) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") & theme(legend.position = "right")

save_fig(fig, "Figure1_Race_Structure.png", width = 14, height = 5.5)
cat("groupings:", paste(groupings, collapse = " | "), "\n")
for (g in groupings) {
  d <- race %>% filter(Grouping == g)
  cat(sprintf("  %-12s max eps2 = %.4f | min q = %.4g | significant = %d\n",
              g, max(d$epsilon2, na.rm = TRUE), min(d$q_BH, na.rm = TRUE), sum(d$sig)))
}
