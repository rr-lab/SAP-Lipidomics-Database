# ==============================================================================
# Figure 2        -- Lipid class composition of the SAP leaf lipidome, and the
#                    LION ontology terms that separate the two trials.
#                      A  every annotated superclass, on a log scale
#                      B  the 13 focal lipid classes, stacked
#                      C  LION enrichment
# Supp Figure S6  -- Each class in a reduced chemical space defined by its
#                    abundance-weighted mean carbon number and double-bond
#                    count, drawn once per trial.
#
# Composition is %TIC, the share of total annotated signal, computed per sample
# and then averaged, so a few high-signal samples cannot dominate. Panels A and
# B use the same denominator (all annotated features), which is why the 13
# classes in B stop near 93% rather than reaching 100%.
#
# Panel C is the one CTL-LIN contrast in this figure and it is labelled as such.
# The chemical-space panels are deliberately drawn without arrows between the
# trials: each point describes one class in one trial and nothing is subtracted.
#
# Inputs
#   data/SPATS_fitted/non_normalized_intensities/Final_subset_{control,lowinput}_*.csv
#   data/lipid_class/final_lipid_classes.csv
#   table/Linex2/LION-enrichment.csv
#
# Outputs
#   fig/main/Figure2_Class_Composition.png
#   fig/supp/SuppFig_S6_Chemical_Space.png
#   table/supp/SuppTable_S5D_Class_Composition_pctTIC.csv
#   table/supp/SuppTable_S5E_LION_Enrichment.csv
#   table/supp/SuppTable_S5F_Chemical_Space.csv
# ==============================================================================
source("scripts/new_new_script/_common.R")
suppressPackageStartupMessages({ library(tidyr); library(forcats); library(tibble) })

class_csv <- Sys.getenv("LIPID_CLASS_CSV",
  file.path(DATA_ROOT, "lipid_class/final_lipid_classes.csv"))
lion_csv  <- Sys.getenv("LION_CSV", file.path(REPO, "table/Linex2/LION-enrichment.csv"))
stopifnot(file.exists(class_csv), file.exists(lion_csv))

CLASS_ORDER <- c("MGDG", "DGDG", "SQDG",                        # galactolipids
                 "PC", "PE", "PG", "PA", "PS", "LPC", "LPE",    # glycerophospholipids
                 "TG", "DG", "MG")                              # neutral glycerolipids

# ---- annotation --------------------------------------------------------------
ann <- vroom(class_csv, show_col_types = FALSE) %>%
  transmute(key = tolower(normalize_lipid_name(Lipids)), SuperClass = Class) %>%
  distinct(key, .keep_all = TRUE)

# ---- %TIC per trial ----------------------------------------------------------
pct_tic <- function(path, label) {
  d <- read_trial(path)
  feats <- names(d)[-1]
  m <- as.matrix(d[, feats]); storage.mode(m) <- "numeric"
  m[!is.finite(m)] <- 0
  share <- sweep(m, 1, pmax(rowSums(m), 1e-12), "/") * 100
  tibble(Feature = feats, pct = colMeans(share, na.rm = TRUE), Condition = label)
}

comp <- bind_rows(pct_tic(CTL_CSV, "CTL"), pct_tic(LIN_CSV, "LIN")) %>%
  mutate(FocusClass = lipid_class(Feature),
         key = tolower(normalize_lipid_name(Feature))) %>%
  left_join(ann, by = "key") %>%
  mutate(SuperClass = ifelse(is.na(SuperClass), "Unclassified", SuperClass),
         Condition  = factor(Condition, c("CTL", "LIN")))

# ---- A: every annotated superclass -------------------------------------------
# Values span four orders of magnitude, so a log axis is needed, and on a log
# axis a bar's length is no longer proportional to its value. Hence points.
supers <- comp %>%
  group_by(Condition, SuperClass) %>%
  summarise(pct = sum(pct), .groups = "drop") %>%
  filter(pct > 0, !is.na(pct)) %>%
  mutate(SuperClass = fct_reorder(SuperClass, pct, .fun = max, .desc = FALSE))

pA <- ggplot(supers, aes(pct, SuperClass)) +
  geom_line(aes(group = SuperClass), colour = "grey65", linewidth = .5) +
  geom_point(aes(fill = Condition), shape = 21, size = 3.4,
             colour = "black", stroke = .4) +
  geom_text(aes(label = ifelse(pct >= 0.01, sprintf("%.2f", pct), "<0.01"),
                vjust = ifelse(Condition == "CTL", -1.25, 2.15)),
            size = 2.5, colour = "grey25") +
  scale_fill_manual(values = condition_colors) +
  scale_x_log10(breaks = c(0.001, 0.1, 10), labels = c("0.001", "0.1", "10"),
                expand = expansion(mult = c(.12, .12))) +
  annotation_logticks(sides = "b", outside = TRUE,
                      short = unit(.05, "cm"), mid = unit(.1, "cm"), long = unit(.15, "cm")) +
  coord_cartesian(clip = "off") +
  labs(x = "Mean %TIC (log scale)", y = NULL) +
  plot_theme + theme(legend.position = "right",
                     panel.grid.major.y = element_line(colour = "grey92"),
                     panel.grid.major.x = element_line(colour = "grey92"),
                     plot.margin = margin(6, 10, 14, 6))

# ---- B: the 13 focal classes -------------------------------------------------
zoom <- comp %>%
  filter(FocusClass %in% CLASS_ORDER) %>%
  group_by(Condition, FocusClass) %>%
  summarise(pct = sum(pct), .groups = "drop") %>%
  mutate(FocusClass = factor(FocusClass, levels = CLASS_ORDER))

zoom_lab <- zoom %>%
  pivot_wider(names_from = Condition, values_from = pct, values_fill = 0) %>%
  mutate(lab = sprintf("%-5s %5.2f %5.2f", FocusClass, CTL, LIN)) %>%
  arrange(factor(FocusClass, levels = CLASS_ORDER))
zoom_labels <- setNames(zoom_lab$lab, as.character(zoom_lab$FocusClass))

pB <- ggplot(zoom, aes(Condition, pct, fill = FocusClass)) +
  geom_col(width = .62, colour = "black", linewidth = .25) +
  geom_text(data = subset(zoom, pct >= 2.5), aes(label = sprintf("%.1f", pct)),
            position = position_stack(vjust = .5), size = 3, colour = "white") +
  scale_fill_manual(values = class_colors, labels = zoom_labels,
                    name = sprintf("%-5s %5s %5s", "", "CTL", "LIN")) +
  scale_y_continuous(expand = expansion(mult = c(0, .03))) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = NULL, y = "Mean %TIC") +
  plot_theme +
  theme(legend.position    = "right",
        legend.title       = element_text(family = "mono", size = 9, face = "bold"),
        legend.text        = element_text(family = "mono", size = 8.5),
        legend.key.height  = unit(12, "pt"),
        panel.grid.major.x = element_blank())

# ---- C: LION enrichment ------------------------------------------------------
# LION compares the two trials directly, so DOWN means depleted in LIN relative
# to CTL. The ten strongest terms in each direction are shown; the full list of
# significant terms is written to the supplementary table.
lion <- vroom(lion_csv, show_col_types = FALSE) %>%
  rename(Term = `Term ID`, Description = Discription,
         p_value = `p-value`, q_value = `FDR q-value`) %>%
  filter(q_value < 0.05) %>%
  mutate(logQ = -log10(q_value),
         Direction = factor(ifelse(Regulated == "UP", "Higher in LIN", "Lower in LIN"),
                            c("Higher in LIN", "Lower in LIN")))

save_table(lion %>% arrange(Direction, q_value) %>%
             dplyr::select(Term, Description, Direction, Annotated, ES, p_value, q_value),
           "SuppTable_S5E_LION_Enrichment.csv")

lion_top <- lion %>% group_by(Direction) %>% slice_min(q_value, n = 10) %>% ungroup() %>%
  mutate(Description = fct_reorder(Description, logQ))

pC <- ggplot(lion_top, aes(logQ, Description)) +
  geom_segment(aes(x = 0, xend = logQ, yend = Description),
               colour = "grey80", linewidth = .4) +
  geom_point(aes(size = Annotated, fill = Direction), shape = 21,
             colour = "black", stroke = .35) +
  scale_fill_manual(values = c("Higher in LIN" = "#FDE725FF", "Lower in LIN" = "#440154FF")) +
  scale_size_continuous(name = "Annotated lipids", range = c(2, 7)) +
  facet_grid(Direction ~ ., scales = "free_y", space = "free_y") +
  labs(x = expression(bold(-log[10]~(FDR~italic(q)))), y = NULL) +
  guides(fill = "none") +
  plot_theme +
  theme(legend.position    = "right",
        legend.title       = element_text(size = 9, face = "bold"),
        axis.text.y        = element_text(size = 8),
        strip.text.y       = element_text(face = "bold", size = 9),
        panel.grid.major.y = element_line(colour = "grey94"),
        panel.grid.major.x = element_line(colour = "grey92"))

fig2 <- ((pA | pB) + plot_layout(widths = c(1.15, 1))) / pC +
  plot_layout(heights = c(1, 1.15)) +
  plot_annotation(tag_levels = "A")
save_fig(fig2, "Figure2_Class_Composition.png", width = 15, height = 13)

# ---- composition table -------------------------------------------------------
class_tab <- comp %>% group_by(Condition, FocusClass) %>%
  summarise(pct = sum(pct), .groups = "drop")
save_table(bind_rows(
  supers    %>% transmute(Level = "Superclass",  Group = as.character(SuperClass), Condition, pct_TIC = pct),
  class_tab %>% transmute(Level = "Lipid class", Group = as.character(FocusClass), Condition, pct_TIC = pct)
) %>% arrange(Level, Group, Condition), "SuppTable_S5D_Class_Composition_pctTIC.csv")

# ---- Supp Figure S6: reduced chemical space ----------------------------------
# Every n:m pair in a species name is summed, so MGDG(18:3/18:3) contributes 36
# carbons and 6 double bonds. Class means are weighted by %TIC, so the abundant
# species inside a class set its position.
sum_pairs <- function(x, part) {
  x <- normalize_lipid_name(x)
  vapply(regmatches(x, gregexpr("[0-9]+:[0-9]+", x)), function(v) {
    if (!length(v)) return(NA_real_)
    sum(as.numeric(if (part == "C") sub(":.*", "", v) else sub(".*:", "", v)))
  }, numeric(1))
}

chem <- comp %>%
  filter(FocusClass %in% CLASS_ORDER) %>%
  mutate(total_c = sum_pairs(Feature, "C"), total_db = sum_pairs(Feature, "DB")) %>%
  filter(!is.na(total_c), pct > 0) %>%
  group_by(Condition, FocusClass) %>%
  summarise(n_species = dplyr::n(),
            WeightedC  = weighted.mean(total_c, pct),
            WeightedDB = weighted.mean(total_db, pct), .groups = "drop") %>%
  mutate(FocusClass = factor(FocusClass, levels = CLASS_ORDER))

save_table(chem, "SuppTable_S5F_Chemical_Space.csv")

figs6 <- ggplot(chem, aes(WeightedC, WeightedDB, fill = FocusClass)) +
  geom_point(shape = 21, size = 4, colour = "black", stroke = .4, alpha = .95) +
  ggrepel::geom_text_repel(aes(label = FocusClass), size = 2.9, colour = "grey20",
                           min.segment.length = 0, segment.colour = "grey70",
                           max.overlaps = 30, seed = 1) +
  scale_fill_manual(values = class_colors, guide = "none") +
  facet_wrap(~ Condition, nrow = 1) +
  labs(x = "Abundance-weighted mean total carbons",
       y = "Abundance-weighted mean double bonds") +
  plot_theme + theme(panel.grid.major.x = element_line(colour = "grey92"),
                     strip.text = element_text(face = "bold", size = 11))

save_fig(figs6, "SuppFig_S6_Chemical_Space.png", width = 11, height = 5.5, subdir = "supp")

# ---- console summary ---------------------------------------------------------
cat("\n-- lipid classes, mean %TIC --\n")
print(as.data.frame(class_tab %>% pivot_wider(names_from = Condition, values_from = pct) %>%
                    arrange(desc(CTL)) %>% mutate(across(where(is.numeric), ~round(.x, 2)))))
cat("\n-- superclasses, mean %TIC --\n")
print(as.data.frame(supers %>% pivot_wider(names_from = Condition, values_from = pct) %>%
                    arrange(desc(CTL)) %>% mutate(across(where(is.numeric), ~round(.x, 3)))))
cat("\n-- LION terms at q < 0.05 --\n")
print(as.data.frame(lion %>% count(Direction)))
cat("\n-- chemical space --\n")
print(as.data.frame(chem %>% mutate(across(where(is.numeric), ~round(.x, 2)))))
