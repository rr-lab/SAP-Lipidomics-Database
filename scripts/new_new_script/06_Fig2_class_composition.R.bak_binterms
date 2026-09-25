# ==============================================================================
# Figure 2 -- Lipid class composition of the SAP leaf lipidome, and the LION
#             ontology terms that separate the two trials.
#
#   A  every annotated superclass, on a log scale
#   B  the 13 focal lipid classes, stacked
#   C  LION enrichment
#
# Composition is %TIC, the share of total annotated signal, computed per sample
# and then averaged. Panels A and B use the same denominator (all annotated
# features), which is why the 13 classes in B stop near 93% rather than 100%.
#
# Panel C is the one CTL-LIN contrast in this figure and it is labelled as such
# in the caption.
#
# plot_theme is used unmodified except where its inside-panel legend would land
# on the data. Panel A moves the legend to the empty bottom-right corner, and
# panels B and C put theirs outside on the right, because a 13-row class table
# and a size key do not fit inside a panel.
#
# Inputs
#   data/SPATS_fitted/non_normalized_intensities/Final_subset_{control,lowinput}_*.csv
#   data/metadata/final_lipid_classes.csv
#
#   NOT data/lipid_class/final_lipid_classes.csv. Those two files are different.
#   data/lipid_class/ is the pre-2026-09-16 annotation, 325 rows, with the old
#   home-grown class vocabulary (Ether lipid, Betaine lipid, Terpenoid, Prenol).
#   data/metadata/ is the current one, 316 rows, carrying the eight LIPID MAPS
#   categories and a Category_Code column. This script defaulted to the old file
#   until 2026-09-17, so running it without LIPID_CLASS_CSV set silently
#   reproduced the superseded grouping while the manuscript quoted the new one.
#   table/Linex2/LION-enrichment.csv
#
# Outputs
#   fig/main/Figure2_Class_Composition.png
#   table/supp/SuppTable_S5D_Class_Composition_pctTIC.csv
#   table/supp/SuppTable_S5E_LION_Enrichment.csv
# SuperClass RENAMED TO Category, 2026-09-17. The values in this column have
# been the eight LIPID MAPS categories since the annotation moved over on
# 2026-09-16, and the manuscript calls them categories throughout -- it contains
# no occurrence of the word superclass, because LIPID MAPS has no such tier. The
# column, the sheet and the file kept the old name, so the supplement shipped a
# heading the paper never uses and implied a hierarchy that does not exist.
# ==============================================================================
source("scripts/new_new_script/_common.R")
suppressPackageStartupMessages({ library(tidyr); library(forcats); library(tibble) })

class_csv <- Sys.getenv("LIPID_CLASS_CSV",
  file.path(DATA_ROOT, "metadata/final_lipid_classes.csv"))
lion_csv  <- Sys.getenv("LION_CSV", file.path(REPO, "table/Linex2/LION-enrichment.csv"))
stopifnot(file.exists(class_csv), file.exists(lion_csv))

STACK_ORDER <- c("MGDG", "DGDG", "SQDG",                        # galactolipids
                 "PC", "PE", "PG", "PA", "PS", "LPC", "LPE",    # glycerophospholipids
                 "TG", "DG", "MG")                              # neutral glycerolipids

ann <- vroom(class_csv, show_col_types = FALSE) %>%
  transmute(key = tolower(normalize_lipid_name(Lipids)), Category = Class) %>%
  distinct(key, .keep_all = TRUE)

comp <- bind_rows(pct_tic(CTL_CSV, "CTL"), pct_tic(LIN_CSV, "LIN")) %>%
  mutate(FocusClass = lipid_class(Feature),
         key = tolower(normalize_lipid_name(Feature))) %>%
  left_join(ann, by = "key") %>%
  mutate(Category = ifelse(is.na(Category), "Unclassified", Category),
         Condition  = factor(Condition, c("CTL", "LIN")))

# ---- A: every annotated superclass -------------------------------------------
# Values span four orders of magnitude, so a log axis is needed, and on a log
# axis a bar's length is no longer proportional to its value. Hence points.
supers <- comp %>%
  group_by(Condition, Category) %>%
  summarise(pct = sum(pct), .groups = "drop") %>%
  filter(pct > 0, !is.na(pct)) %>%
  mutate(Category = fct_reorder(Category, pct, .fun = max, .desc = FALSE))

pA <- ggplot(supers, aes(pct, Category)) +
  geom_line(aes(group = Category), colour = "grey65", linewidth = .5) +
  geom_point(aes(fill = Condition), shape = 21, size = 4.5,
             colour = "black", stroke = .5) +
  geom_text(aes(label = ifelse(pct >= 0.01, sprintf("%.2f", pct), "<0.01"),
                vjust = ifelse(Condition == "CTL", -1.25, 2.15)),
            size = 4.2, colour = "grey25") +
  scale_fill_manual(values = condition_colors) +
  scale_x_log10(breaks = c(0.001, 0.1, 10), labels = c("0.001", "0.1", "10"),
                expand = expansion(mult = c(.12, .12))) +
  annotation_logticks(sides = "b", outside = TRUE,
                      short = unit(.05, "cm"), mid = unit(.1, "cm"), long = unit(.15, "cm")) +
  coord_cartesian(clip = "off") +
  labs(x = "Mean %TIC (log scale)", y = NULL) +
  plot_theme + theme(plot.margin = margin(15, 15, 26, 15),
                     legend.position = c(0.98, 0.04),
                     legend.justification = c("right", "bottom"))

# ---- B: the 13 focal classes -------------------------------------------------
zoom <- comp %>%
  filter(FocusClass %in% STACK_ORDER) %>%
  group_by(Condition, FocusClass) %>%
  summarise(pct = sum(pct), .groups = "drop") %>%
  mutate(FocusClass = factor(FocusClass, levels = STACK_ORDER))

zoom_lab <- zoom %>%
  pivot_wider(names_from = Condition, values_from = pct, values_fill = 0) %>%
  mutate(lab = sprintf("%-5s %5.2f %5.2f", FocusClass, CTL, LIN)) %>%
  arrange(factor(FocusClass, levels = STACK_ORDER))
zoom_labels <- setNames(zoom_lab$lab, as.character(zoom_lab$FocusClass))

pB <- ggplot(zoom, aes(Condition, pct, fill = FocusClass)) +
  geom_col(width = .62, colour = "black", linewidth = .25) +
  geom_text(data = subset(zoom, pct >= 2.5), aes(label = sprintf("%.1f", pct)),
            position = position_stack(vjust = .5), size = 4.8, colour = "white") +
  scale_fill_manual(values = class_colors, labels = zoom_labels,
                    name = sprintf("%-5s %5s %5s", "", "CTL", "LIN")) +
  scale_y_continuous(expand = expansion(mult = c(0, .03))) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = NULL, y = "Mean %TIC") +
  plot_theme +
  theme(legend.position   = "right",
        legend.background = element_blank(),
        legend.title      = element_text(family = "mono", size = 14, face = "bold"),
        legend.text       = element_text(family = "mono", size = 13),
        legend.key.height = unit(20, "pt"))

# ---- C: LION enrichment ------------------------------------------------------
# The LION run (job3, Sep 2025, post-deduplication input, two-tailed K-S) was
# submitted with "condition of interest" = CTL and "control condition" = LIN, so
# LION's Regulated column is referenced to CTL: UP means higher in CTL, that is
# LOWER in LIN. The mapping below inverts it so the figure reads in LIN terms
# like the rest of the paper.
#
# The convention is checked against the composition, not assumed. In the input
# matrix TG is 2.41 %TIC under CTL and 3.67 under LIN, and the whole
# glycerophospholipid pool is 30.23 under CTL and 28.46 under LIN. So the
# triacylglycerol term must carry LION's "higher in LIN" label and the
# glycerophospholipid term must carry "lower in LIN". If a future LION run is
# submitted the other way round, the stopifnot below fails and the ifelse needs
# flipping back.
#
# The ten strongest terms in each direction are shown; the full list of
# significant terms goes to the supplementary table.
lion <- vroom(lion_csv, show_col_types = FALSE) %>%
  rename(Term = `Term ID`, Description = Discription,
         p_value = `p-value`, q_value = `FDR q-value`) %>%
  filter(q_value < 0.05) %>%
  mutate(logQ = -log10(q_value),
         Direction = factor(ifelse(Regulated == "UP", "Lower in LIN", "Higher in LIN"),
                            c("Higher in LIN", "Lower in LIN")))

stopifnot(lion$Direction[lion$Description == "triacylglycerols [GL0301]"] ==
            "Higher in LIN")

save_table(lion %>% arrange(Direction, q_value) %>%
             dplyr::select(Term, Description, Direction, Annotated, ES, p_value, q_value),
           "SuppTable_S5E_LION_Enrichment.csv")

lion_top <- lion %>% group_by(Direction) %>% slice_min(q_value, n = 10) %>% ungroup() %>%
  mutate(Description = fct_reorder(Description, logQ))

pC <- ggplot(lion_top, aes(logQ, Description)) +
  geom_segment(aes(x = 0, xend = logQ, yend = Description),
               colour = "grey80", linewidth = .4) +
  geom_point(aes(size = Annotated, fill = Direction), shape = 21,
             colour = "black", stroke = .45) +
  scale_fill_manual(values = c("Higher in LIN" = "#FDE725FF", "Lower in LIN" = "#440154FF")) +
  scale_size_continuous(name = "Annotated lipids", range = c(3, 10)) +
  facet_grid(Direction ~ ., scales = "free_y", space = "free_y") +
  labs(x = expression(bold(-log[10]~(FDR~italic(q)))), y = NULL) +
  guides(fill = "none") +
  plot_theme +
  theme(legend.position   = "right",
        legend.background = element_blank(),
        legend.title      = element_text(size = 14, face = "bold"),
        axis.text.y       = element_text(size = 14),
        strip.text.y      = element_text(face = "bold", size = 15))

fig2 <- ((pA | pB) + plot_layout(widths = c(1.15, 1))) / pC +
  plot_layout(heights = c(1, 1.15)) +
  plot_annotation(tag_levels = "A", theme = TAG_THEME)
save_fig(fig2, "Figure2_Class_Composition.png", width = 19, height = 17)

# ---- composition table -------------------------------------------------------
class_tab <- comp %>% group_by(Condition, FocusClass) %>%
  summarise(pct = sum(pct), .groups = "drop")
save_table(bind_rows(
  supers    %>% transmute(Level = "Category",  Group = as.character(Category), Condition, pct_TIC = pct),
  class_tab %>% transmute(Level = "Lipid class", Group = as.character(FocusClass), Condition, pct_TIC = pct)
) %>% arrange(Level, Group, Condition), "SuppTable_S5D_Class_Composition_pctTIC.csv")

cat("\n-- lipid classes, mean %TIC --\n")
print(as.data.frame(class_tab %>% pivot_wider(names_from = Condition, values_from = pct) %>%
                    arrange(desc(CTL)) %>% mutate(across(where(is.numeric), ~round(.x, 2)))))
cat("\n-- LIPID MAPS categories, mean %TIC --\n")
print(as.data.frame(supers %>% pivot_wider(names_from = Condition, values_from = pct) %>%
                    arrange(desc(CTL)) %>% mutate(across(where(is.numeric), ~round(.x, 3)))))
cat("\n-- LION terms at q < 0.05 --\n")
print(as.data.frame(lion %>% count(Direction)))
