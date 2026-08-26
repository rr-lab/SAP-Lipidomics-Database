# ==============================================================================
# Figure 2 -- Lipid class composition of the SAP leaf lipidome in each trial.
#
# Composition is reported as %TIC, the share of total annotated signal, for each
# trial SEPARATELY. Nothing here is a cross-trial contrast: the two trials are
# described side by side and no difference is computed or tested, because CTL
# (2019) and LIN (2022) differ in year, planting date, field and acquisition
# batch, and a difference between them cannot be attributed to management.
#
#   A  every annotated superclass, including the minor ones (terpenoids,
#      sterols, prenols, fatty acids and the rest), on a log scale so classes
#      below 1% are legible rather than invisible inside a stacked bar
#   B  a zoom on the two dominant families, showing the individual classes
#      within the galactolipids and the glycerophospholipids
#
# This replaces the former lipidome landscape figure (CLR contrast, chemical space, LION), all
# three panels of which were LIN - CTL contrasts.
#
# Output: fig/main/Figure2_Class_Composition.png
#         table/supp/SuppTable_S5D_Class_Composition_pctTIC.csv
# ==============================================================================
source("scripts/pipeline/_common.R")
suppressPackageStartupMessages({ library(tidyr); library(forcats) })

class_csv <- Sys.getenv("LIPID_CLASS_CSV",
  file.path(DATA_ROOT, "lipid_class/final_lipid_classes.csv"))
out_png <- Sys.getenv("FIG2_OUT", file.path(FIG_MAIN, "Figure2_Class_Composition.png"))

# ---- annotation -------------------------------------------------------------
ann <- vroom(class_csv, show_col_types = FALSE) %>%
  transmute(key = tolower(normalize_lipid_name(Lipids)), SuperClass = Class) %>%
  distinct(key, .keep_all = TRUE)

# ---- %TIC per trial ---------------------------------------------------------
# TIC is the total annotated signal in a sample. Each feature's share is taken
# per sample first, then averaged across samples, so a few high-signal samples
# cannot dominate the mean composition.
pct_tic <- function(path, label) {
  d <- read_trial(path)
  feats <- names(d)[-1]
  m <- as.matrix(d[, feats])
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

# ---- A: every annotated superclass ------------------------------------------
# The overview. Values span four orders of magnitude, so a log axis is needed --
# and on a log axis a bar's length is no longer proportional to its value, which
# is why these are points rather than bars.
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
  scale_x_log10(breaks = c(0.001, 0.1, 10),
                labels = c("0.001", "0.1", "10"),
                expand = expansion(mult = c(.12, .12))) +
  annotation_logticks(sides = "b", outside = TRUE,
                      short = unit(.05,"cm"), mid = unit(.1,"cm"), long = unit(.15,"cm")) +
  coord_cartesian(clip = "off") +
  labs(x = "Mean %TIC (log scale)", y = NULL) +
  plot_theme + theme(legend.position = "right",
                     panel.grid.major.y = element_line(colour = "grey92"),
                     panel.grid.major.x = element_line(colour = "grey92"),
                     plot.margin = margin(6, 10, 14, 6))

# ---- B: zoom on the two dominant families -----------------------------------
# The individual classes inside the galactolipids and the glycerophospholipids,
# which together are roughly three quarters of the annotated signal. Everything
# outside these two families is in panel A and is deliberately not repeated
# here, so the bars stop short of 100%.
ZOOM_CLASSES <- c("MGDG", "DGDG", "SQDG",                       # galactolipids
                  "PC", "PE", "PG", "PA", "PS", "LPC", "LPE")   # glycerophospholipids

zoom <- comp %>%
  filter(FocusClass %in% ZOOM_CLASSES) %>%
  group_by(Condition, FocusClass) %>%
  summarise(pct = sum(pct), .groups = "drop") %>%
  mutate(FocusClass = factor(FocusClass, levels = ZOOM_CLASSES))

zoom_lab <- zoom %>%
  tidyr::pivot_wider(names_from = Condition, values_from = pct, values_fill = 0) %>%
  mutate(lab = sprintf("%-6s %5.2f %5.2f", FocusClass, CTL, LIN)) %>%
  arrange(factor(FocusClass, levels = ZOOM_CLASSES))
zoom_labels <- setNames(zoom_lab$lab, as.character(zoom_lab$FocusClass))

pB <- ggplot(zoom, aes(Condition, pct, fill = FocusClass)) +
  geom_col(width = .62, colour = "black", linewidth = .25) +
  geom_text(data = subset(zoom, pct >= 2.5),
            aes(label = sprintf("%.1f", pct)),
            position = position_stack(vjust = .5), size = 3, colour = "white") +
  scale_fill_manual(values = class_colors, labels = zoom_labels,
                    name = sprintf("%-6s %5s %5s", "", "CTL", "LIN")) +
  scale_y_continuous(expand = expansion(mult = c(0, .03))) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = NULL, y = "Mean %TIC") +
  plot_theme +
  theme(legend.position    = "right",
        legend.title       = element_text(family = "mono", size = 9, face = "bold"),
        legend.text        = element_text(family = "mono", size = 8.5),
        legend.key.height  = unit(12, "pt"),
        panel.grid.major.x = element_blank())

fig <- (pA | pB) + plot_layout(widths = c(1.15, 1)) +
  plot_annotation(tag_levels = "A")
save_fig(fig, basename(out_png), width = 15, height = 7)

# ---- the numbers behind the figure ------------------------------------------
class_tab <- comp %>% group_by(Condition, FocusClass) %>%
  summarise(pct = sum(pct), .groups = "drop")
tab <- bind_rows(
  supers    %>% transmute(Level = "Superclass",  Group = as.character(SuperClass), Condition, pct_TIC = pct),
  class_tab %>% transmute(Level = "Lipid class", Group = as.character(FocusClass), Condition, pct_TIC = pct)
) %>% arrange(Level, Group, Condition)
save_table(tab, "SuppTable_S5D_Class_Composition_pctTIC.csv")

cat("\n-- individual classes --\n")
print(as.data.frame(class_tab %>% pivot_wider(names_from = Condition, values_from = pct) %>%
                    arrange(desc(CTL)) %>% mutate(across(where(is.numeric), ~round(.x, 2)))))
cat("\n-- superclasses, %TIC (panel A) --\n")
print(as.data.frame(supers %>% pivot_wider(names_from = Condition, values_from = pct) %>%
                    arrange(desc(CTL)) %>% mutate(across(where(is.numeric), ~round(.x, 3)))))
