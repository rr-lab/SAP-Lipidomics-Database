# ==============================================================================
# Figure 1 -- Cross-trial analytical comparability (CTL 2019 vs LIN 2022)
#
# The two trials were run three years apart, so every CTL-vs-LIN contrast in
# this study rests on the platform having behaved comparably. This figure makes
# that premise explicit.
#   A  QC-RSD distributions, both trials, before and after SERRF
#   B  the same comparison as summary metrics
#   C  named lipid species detected in each trial
#
# Inputs  : SoLD_paper/data/SEERF_normalized_intensities/{A,B}_QC-RSDs.csv
#           SoLD_paper/data/SPATS_fitted/non_normalized_intensities/Final_subset_*.csv
# Output  : fig/main/Figure1_Analytical_Comparability.png
# ==============================================================================
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(patchwork)
})

data_root <- Sys.getenv("SOLD_DATA", "/Users/nirwantandukar/Documents/Github/SoLD_paper/data")
out_file  <- Sys.getenv("FIG1_OUT",
  "/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database/fig/main/Figure1_Analytical_Comparability.png")
serrf_dir <- file.path(data_root, "SEERF_normalized_intensities")
spats_dir <- file.path(data_root, "SPATS_fitted/non_normalized_intensities")
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

# project condition palette (viridis endpoints, as used in Fig 1 and Supp S1-S6)
condition_colors <- c(CTL = "#440154FF", LIN = "#FDE725FF")

plot_theme <- theme_minimal(base_size = 13) +
  theme(axis.text = element_text(colour = "black", size = 9),
        axis.title = element_text(face = "bold", size = 11),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        plot.tag = element_text(face = "bold", size = 16))

# ---- A: QC-RSD distributions, both trials, before vs after SERRF -------------
rsd <- bind_rows(
  read_csv(file.path(serrf_dir, "A_QC-RSDs.csv"), show_col_types = FALSE) %>% mutate(Condition = "CTL"),
  read_csv(file.path(serrf_dir, "B_QC-RSDs.csv"), show_col_types = FALSE) %>% mutate(Condition = "LIN")
) %>%
  select(Condition, Before = QC_none, After = QC_SERRF) %>%
  pivot_longer(c(Before, After), names_to = "Stage", values_to = "RSD") %>%
  filter(is.finite(RSD), RSD > 0) %>%
  mutate(Stage = factor(Stage, c("Before", "After")), RSD = RSD * 100,
         Condition = factor(Condition, c("CTL", "LIN")))

pA <- ggplot(rsd, aes(RSD, colour = Condition, fill = Condition)) +
  geom_density(alpha = .40, linewidth = .55) +
  geom_vline(xintercept = 30, linetype = "dashed", colour = "grey35", linewidth = .4) +
  facet_wrap(~Stage, ncol = 1, scales = "free_y") +
  scale_x_log10(breaks = c(.1, 1, 10, 30, 100), labels = c("0.1", "1", "10", "30", "100")) +
  scale_colour_manual(values = c(CTL = "#440154FF", LIN = "#B8A800")) +
  scale_fill_manual(values = condition_colors) +
  labs(x = "QC relative standard deviation (%)", y = "Density") +
  plot_theme + theme(legend.position = "top",
                     strip.text = element_text(face = "bold", size = 9, hjust = 0))

# ---- B: analytical performance summary --------------------------------------
summ <- rsd %>% group_by(Condition, Stage) %>%
  summarise(`Median QC-RSD (%)` = median(RSD),
            `Features < 30% RSD (%)` = 100 * mean(RSD < 30), .groups = "drop") %>%
  pivot_longer(-c(Condition, Stage), names_to = "Metric", values_to = "Value")

pB <- ggplot(summ, aes(Stage, Value, fill = Condition)) +
  geom_col(position = position_dodge(.72), width = .64, colour = "black", linewidth = .25) +
  geom_text(aes(label = sprintf("%.1f", Value)), position = position_dodge(.72),
            vjust = -.4, size = 2.9) +
  facet_wrap(~Metric, scales = "free_y") +
  scale_fill_manual(values = condition_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, .18))) +
  labs(x = NULL, y = NULL) +
  plot_theme + theme(legend.position = "top",
                     strip.text = element_text(face = "bold", size = 9),
                     panel.grid.major.x = element_blank())

# ---- C: named lipid species detected in each trial --------------------------
sp <- function(f) { h <- names(read_csv(f, n_max = 0, show_col_types = FALSE))[-(1:4)]; h[grepl("\\(", h)] }
ctl <- sp(file.path(spats_dir, "Final_subset_control_all_lipids_fitted_phenotype_non_normalized.csv"))
lin <- sp(file.path(spats_dir, "Final_subset_lowinput_all_lipids_fitted_phenotype_non_normalized.csv"))

det <- tibble(Set = factor(c("CTL only", "Shared", "LIN only"), c("CTL only", "Shared", "LIN only")),
              n = c(length(setdiff(ctl, lin)), length(intersect(ctl, lin)), length(setdiff(lin, ctl))))

pC <- ggplot(det, aes(Set, n, fill = Set)) +
  geom_col(width = .6, colour = "black", linewidth = .25) +
  geom_text(aes(label = n), vjust = -.4, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("CTL only" = condition_colors[["CTL"]],
                               "Shared"   = "grey65",
                               "LIN only" = condition_colors[["LIN"]]), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, .20))) +
  labs(x = NULL, y = "Lipid species") +
  plot_theme + theme(panel.grid.major.x = element_blank())

fig <- (pA | pB) / (pC | plot_spacer()) +
  plot_layout(heights = c(1.35, 1), widths = c(1, 1)) +
  plot_annotation(tag_levels = "A")

ggsave(out_file, fig, width = 13.5, height = 8.4, dpi = 300, bg = "white")

cat(sprintf("species: union %d, shared %d, CTL-only %d, LIN-only %d\n",
            length(union(ctl, lin)), length(intersect(ctl, lin)),
            length(setdiff(ctl, lin)), length(setdiff(lin, ctl))))
print(as.data.frame(summ), digits = 4)
