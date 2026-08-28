# ==============================================================================
# Supplementary Figure S4 -- Lipid classes in a reduced chemical space.
#
# Each class sits at its abundance-weighted mean total carbon number and its
# abundance-weighted mean double-bond count, with an arrow running from its CTL
# position to its LIN position. One panel, colour by class, shape by trial.
#
# The arrow makes this a CTL-LIN contrast, and the two trials differ in year,
# field and acquisition batch, so the caption calls the shifts descriptive.
#
# Weights are the %TIC composition from _common.R, the same numbers Figure 2
# plots, so a class cannot sit in one place here and another there.
#
# plot_theme is used unmodified except for the legend, which sits outside on the
# right: thirteen classes plus two shapes do not fit inside a panel.
#
# Inputs
#   data/SPATS_fitted/non_normalized_intensities/Final_subset_{control,lowinput}_*.csv
#
# Outputs
#   fig/supp/SuppFig_S6_Chemical_Space.png   (prints as S4; filename kept
#                                             so the tex include path does not move)
#   table/supp/SuppTable_S5F_Chemical_Space.csv
# ==============================================================================
source("scripts/new_new_script/_common.R")
suppressPackageStartupMessages({ library(tibble) })

# Every n:m pair in a species name is summed, so MGDG(18:3/18:3) contributes 36
# carbons and 6 double bonds, and TG(16:0/18:2/18:3) contributes 52 and 5.
sum_pairs <- function(x, part) {
  x <- normalize_lipid_name(x)
  vapply(regmatches(x, gregexpr("[0-9]+:[0-9]+", x)), function(v) {
    if (!length(v)) return(NA_real_)
    sum(as.numeric(if (part == "C") sub(":.*", "", v) else sub(".*:", "", v)))
  }, numeric(1))
}

chem <- bind_rows(pct_tic(CTL_CSV, "CTL"), pct_tic(LIN_CSV, "LIN")) %>%
  mutate(Condition  = factor(Condition, c("CTL", "LIN")),
         FocusClass = lipid_class(Feature)) %>%
  filter(FocusClass %in% CLASS_ORDER) %>%
  mutate(total_c = sum_pairs(Feature, "C"), total_db = sum_pairs(Feature, "DB")) %>%
  filter(!is.na(total_c), pct > 0) %>%
  group_by(Condition, FocusClass) %>%
  summarise(n_species  = dplyr::n(),
            WeightedC  = weighted.mean(total_c, pct),
            WeightedDB = weighted.mean(total_db, pct), .groups = "drop") %>%
  mutate(FocusClass = factor(FocusClass, levels = CLASS_ORDER))

save_table(chem, "SuppTable_S5F_Chemical_Space.csv")

# geom_path draws in data order within each group, and summarise returns CTL
# before LIN, so every arrow points from the control to the low-input position.
figs6 <- ggplot(chem, aes(WeightedC, WeightedDB, colour = FocusClass, shape = Condition)) +
  geom_path(aes(group = FocusClass),
            arrow = arrow(type = "closed", length = unit(0.12, "inches")),
            linewidth = 0.8, alpha = 0.7) +
  geom_point(size = 4.5, alpha = 0.9) +
  scale_colour_manual(values = class_colors, name = "Class") +
  scale_shape_manual(values = c(CTL = 16, LIN = 17), name = "Condition") +
  labs(x = "Weighted Mean Total Carbons", y = "Weighted Mean Double Bonds") +
  plot_theme +
  theme(legend.position   = "right",
        legend.background = element_blank(),
        legend.title      = element_text(size = 14, face = "bold"))

save_fig(figs6, "SuppFig_S6_Chemical_Space.png", width = 12, height = 8, subdir = "supp")

cat("\n-- chemical space --\n")
print(as.data.frame(chem %>% mutate(across(where(is.numeric), ~round(.x, 2)))))
