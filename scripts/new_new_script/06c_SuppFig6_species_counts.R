# ==============================================================================
# Supplementary Figure S6 -- lipid species counts and overlap between trials.
#
#   Rscript scripts/new_new_script/06c_SuppFig6_species_counts.R
#
#   A  totals: species per trial, plus common and trial-specific
#   B  species counts by lipid class
#   C  species counts by lipid superclass
#
# Everything is read from the frozen species set written by
# 06b_SuppTableS6_species_inventory.R, so this figure and Supplementary Tables
# S6a-S6c can never disagree. Run 06b first.
#
# The counting rule and its caveats live in 06b's header. In short, a feature
# counts as a species if its annotated name carries lipid shorthand; class
# labels use the lyso normalisation, so the three CTL PC(x/0:0) species are
# counted as LPC and PC is 36 in CTL rather than 39.
#
# plot_theme unmodified except for the legend, which moves to the bottom on B
# and C (a legend inside the panel sits on top of the bars) and is dropped on A.
# No titles, per house style; panels are tagged A/B/C.
#
# Inputs
#   data/final_species_set/species_inventory.csv
#
# Outputs
#   fig/supp/SuppFig_S4_Lipid_Species_Counts.png  (prints as S6; filename kept
#                                                  so the tex include path holds)
# ==============================================================================
source("scripts/new_new_script/_common.R")
suppressPackageStartupMessages({ library(tidyr); library(forcats) })

inv_csv <- file.path(DATA_ROOT, "final_species_set/species_inventory.csv")
stopifnot(file.exists(inv_csv))
inv <- vroom(inv_csv, show_col_types = FALSE)

# ---- A: totals ---------------------------------------------------------------
tot <- tibble::tibble(
  Category = factor(c("CTL total", "LIN total", "Common", "CTL only", "LIN only"),
                    levels = c("CTL total", "LIN total", "Common", "CTL only", "LIN only")),
  Count = c(sum(inv$In_CTL), sum(inv$In_LIN), sum(inv$Status == "Common"),
            sum(inv$Status == "CTL only"), sum(inv$Status == "LIN only")),
  Fill  = c("CTL", "LIN", "Common", "CTL", "LIN"))

pA <- ggplot(tot, aes(Category, Count, fill = Fill)) +
  geom_col(width = .7, colour = "grey30", linewidth = .3) +
  geom_text(aes(label = Count), vjust = -0.4, size = 5.5, fontface = "bold") +
  scale_fill_manual(values = c(condition_colors, Common = "grey70"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, .12))) +
  labs(x = NULL, y = "Lipid species") +
  plot_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# ---- B and C: counts by class and by superclass -------------------------------
counts_by <- function(col) {
  inv %>%
    group_by(Group = .data[[col]]) %>%
    summarise(CTL = sum(In_CTL), LIN = sum(In_LIN), .groups = "drop") %>%
    mutate(Group = fct_reorder(Group, CTL + LIN)) %>%
    pivot_longer(c(CTL, LIN), names_to = "Condition", values_to = "Count") %>%
    mutate(Condition = factor(Condition, c("CTL", "LIN")))
}

bar_panel <- function(d) {
  ggplot(d, aes(Count, Group, fill = Condition)) +
    geom_col(position = position_dodge(width = .75), width = .7,
             colour = "grey30", linewidth = .3) +
    scale_fill_manual(values = condition_colors) +
    scale_x_continuous(expand = expansion(mult = c(0, .08))) +
    labs(x = "Lipid species", y = NULL) +
    plot_theme +
    theme(legend.position   = "bottom",
          legend.direction  = "horizontal",
          legend.background = element_blank())
}

pB <- bar_panel(counts_by("Class"))
pC <- bar_panel(counts_by("SuperClass"))

fig <- (pA / (pB | pC)) +
  plot_layout(heights = c(0.8, 1.4)) +
  plot_annotation(tag_levels = "A") & TAG_THEME

save_fig(fig, "SuppFig_S4_Lipid_Species_Counts.png",
         width = 14, height = 12, subdir = "supp")

cat("\n-- totals --\n");     print(as.data.frame(tot[, c("Category", "Count")]))
cat("\n-- by class --\n");   print(as.data.frame(counts_by("Class") %>%
  pivot_wider(names_from = Condition, values_from = Count) %>% arrange(desc(CTL + LIN))))
