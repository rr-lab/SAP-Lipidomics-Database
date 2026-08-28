# ==============================================================================
# Supplementary Figure S8 -- CTL and LIN GWAS identify largely non-overlapping
#                            candidate loci.
#
#   A  candidate-set overlap at gene level and at 250 kb locus level
#   B  fold enrichment of the overlap against the resolution it is counted at
#   C  candidate genes per lipid class, split CTL-only / shared / LIN-only
#   D  candidate-gene counts against the independent 250 kb loci they occupy
#
# Replaces the Python version (_legacy_pipeline/41_figure_overlap.py), which
# carried an overall title and a title on every panel. Panel identity comes from
# the tag letter and the caption, as everywhere else in the paper.
#
# Panel A draws its circles from first principles rather than pulling in a Venn
# package, so the script has no dependency beyond what _common.R already loads.
# The circles are fixed-size and are not area-proportional; the counts are
# printed, which is what the panel is for.
#
# Inputs (all under table/overlap/)
#   gwas_overlap_overall.csv          gene-level overlap per trait layer
#   gwas_overlap_locus_level.csv      the same after collapsing to 100/250/500 kb
#   gwas_overlap_by_class.csv         per lipid class, individual-lipid layer
#   gwas_gene_to_locus_inflation.csv  genes and loci per condition and layer
#
# Output
#   fig/supp/Figure7_CTL_LIN_Overlap.png   (prints as Supp Fig S8)
# ==============================================================================
source("scripts/new_new_script/_common.R")
suppressPackageStartupMessages({ library(tidyr); library(forcats); library(tibble) })

ov_dir   <- Sys.getenv("OVERLAP_DIR", file.path(REPO, "table/overlap"))
out_name <- "Figure7_CTL_LIN_Overlap.png"
stopifnot(dir.exists(ov_dir))

ctl_col <- unname(condition_colors["CTL"])
lin_col <- unname(condition_colors["LIN"])
shared_col <- "#7A7A7A"

stars <- function(p) ifelse(p < .001, "***", ifelse(p < .01, "**",
                     ifelse(p < .05, "*", "n.s.")))

# ---- A: overlap at two resolutions -------------------------------------------
overall <- vroom(file.path(ov_dir, "gwas_overlap_overall.csv"), show_col_types = FALSE)
locus   <- vroom(file.path(ov_dir, "gwas_overlap_locus_level.csv"), show_col_types = FALSE)

g <- overall %>% filter(Layer == "All layers")
l <- locus   %>% filter(Layer == "All layers", window_kb == 250)

venn_counts <- tribble(
  ~panel,                ~ctl_only,      ~shared,      ~lin_only,      ~fold,           ~p,
  "Gene level",          g$CTL_only[1],  g$n_shared[1], g$LIN_only[1], g$fold_enrichment[1], g$p_hypergeom[1],
  "250 kb locus level",  l$n_CTL[1] - l$n_shared[1], l$n_shared[1],
                         l$n_LIN[1] - l$n_shared[1], l$fold_enrichment[1], l$p_hypergeom[1]) %>%
  mutate(panel = factor(panel, c("Gene level", "250 kb locus level")))

circle <- function(x0, r = 1, n = 240) {
  t <- seq(0, 2 * pi, length.out = n)
  data.frame(x = x0 + r * cos(t), y = r * sin(t))
}
circles <- bind_rows(
  lapply(levels(venn_counts$panel), function(pn)
    bind_rows(cbind(circle(-0.42), set = "CTL", panel = pn),
              cbind(circle( 0.42), set = "LIN", panel = pn)))) %>%
  mutate(panel = factor(panel, levels(venn_counts$panel)))

lab_a <- venn_counts %>%
  transmute(panel,
            lab = sprintf("%.2f%s chance,  %s", fold, "×", stars(p)))

pA <- ggplot() +
  geom_polygon(data = circles, aes(x, y, group = interaction(panel, set), fill = set),
               colour = "black", linewidth = .4, alpha = .55) +
  geom_text(data = venn_counts, aes(-0.92, 0, label = ctl_only), size = 5, fontface = "bold", colour = "white") +
  geom_text(data = venn_counts, aes( 0.00, 0, label = shared),   size = 5, fontface = "bold") +
  geom_text(data = venn_counts, aes( 0.92, 0, label = lin_only), size = 5, fontface = "bold") +
  geom_text(data = venn_counts, aes(-0.72, 1.16, label = "CTL"), size = 4.6, fontface = "bold", colour = ctl_col) +
  geom_text(data = venn_counts, aes( 0.72, 1.16, label = "LIN"), size = 4.6, fontface = "bold", colour = lin_col) +
  geom_text(data = lab_a, aes(0, -1.42, label = lab), size = 4.2, colour = "grey25") +
  facet_wrap(~ panel, nrow = 1) +
  scale_fill_manual(values = c(CTL = ctl_col, LIN = lin_col), guide = "none") +
  coord_fixed(xlim = c(-1.7, 1.7), ylim = c(-1.6, 1.4)) +
  labs(x = NULL, y = NULL) +
  plot_theme +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank(), panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 13))

# ---- B: fold enrichment against the resolution it is counted at --------------
res <- bind_rows(
  overall %>% filter(Layer != "All layers") %>%
    transmute(Layer, res = "Gene", fold = fold_enrichment, p = p_hypergeom),
  locus %>% filter(Layer != "All layers") %>%
    transmute(Layer, res = paste0(window_kb, " kb"), fold = fold_enrichment, p = p_hypergeom)) %>%
  mutate(res = factor(res, c("Gene", "100 kb", "250 kb", "500 kb")),
         Layer = factor(Layer, c("Individual lipids", "Class sums / ratios")))

pB <- ggplot(res, aes(res, fold, fill = Layer)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey45", linewidth = .5) +
  geom_col(position = position_dodge(.8), width = .7, colour = "black", linewidth = .25) +
  geom_text(aes(label = stars(p)), position = position_dodge(.8),
            vjust = -0.35, size = 4.2, fontface = "bold") +
  scale_fill_manual(values = c("Individual lipids" = "#2166AC",
                               "Class sums / ratios" = "#D55E00")) +
  scale_y_continuous(expand = expansion(mult = c(0, .16))) +
  labs(x = "Resolution at which overlap is counted",
       y = "Fold enrichment over chance") +
  plot_theme + theme(panel.grid.major.x = element_blank())

# ---- C: candidate genes per lipid class --------------------------------------
by_class <- vroom(file.path(ov_dir, "gwas_overlap_by_class.csv"), show_col_types = FALSE) %>%
  filter(Layer == "individual") %>%
  mutate(Class = fct_reorder(Class, n_CTL + n_LIN))

cls_long <- by_class %>%
  transmute(Class, `CTL only` = CTL_only, Shared = n_shared, `LIN only` = LIN_only) %>%
  pivot_longer(-Class, names_to = "part", values_to = "n") %>%
  mutate(part = factor(part, c("CTL only", "Shared", "LIN only")))

cls_lab <- by_class %>%
  transmute(Class, tot = CTL_only + n_shared + LIN_only,
            lab = sprintf("%d shared, %.1f%s %s", n_shared, fold_enrichment, "×", stars(q_BH)))

pC <- ggplot(cls_long, aes(n, Class, fill = part)) +
  geom_col(width = .68, colour = "black", linewidth = .25) +
  geom_text(data = cls_lab, aes(tot, Class, label = lab), inherit.aes = FALSE,
            hjust = -0.06, size = 3.9, colour = "grey25") +
  scale_fill_manual(values = c("CTL only" = ctl_col, "Shared" = shared_col, "LIN only" = lin_col)) +
  scale_x_continuous(expand = expansion(mult = c(0, .34))) +
  labs(x = "LD-mapped candidate genes (individual-lipid GWAS)", y = NULL) +
  plot_theme + theme(panel.grid.major.y = element_blank())

# ---- D: genes against the independent loci they occupy -----------------------
infl <- vroom(file.path(ov_dir, "gwas_gene_to_locus_inflation.csv"), show_col_types = FALSE) %>%
  filter(Layer != "all layers") %>%
  mutate(grp = paste(Condition, sub("class sums/ratios", "sums/ratios", Layer)),
         grp = factor(grp, unique(grp)))

infl_long <- infl %>%
  transmute(grp, `Candidate genes` = n_genes, `Independent 250 kb loci` = n_loci) %>%
  pivot_longer(-grp, names_to = "what", values_to = "n") %>%
  mutate(what = factor(what, c("Candidate genes", "Independent 250 kb loci")))

pD <- ggplot(infl_long, aes(grp, n, fill = what)) +
  geom_col(position = position_dodge(.8), width = .7, colour = "black", linewidth = .25) +
  geom_text(aes(label = n), position = position_dodge(.8),
            vjust = -0.35, size = 3.9, fontface = "bold") +
  geom_text(data = infl, aes(grp, n_genes, label = sprintf("%.1f genes/locus", genes_per_locus)),
            inherit.aes = FALSE, vjust = -2.1, size = 3.7, colour = "grey25") +
  scale_fill_manual(values = c("Candidate genes" = "grey72",
                               "Independent 250 kb loci" = "#1B9E77")) +
  scale_y_log10(expand = expansion(mult = c(0, .22))) +
  labs(x = NULL, y = "Count (log scale)") +
  plot_theme + theme(panel.grid.major.x = element_blank())

fig <- (pA | pB) / (pC | pD) +
  plot_annotation(tag_levels = "A", theme = TAG_THEME)

save_fig(fig, out_name, width = 17, height = 12, subdir = "supp")

cat("\n-- overlap, all layers ------------------------------------------------\n")
cat(sprintf("  gene level      CTL-only %d  shared %d  LIN-only %d   %.2fx  p = %.3g\n",
            g$CTL_only[1], g$n_shared[1], g$LIN_only[1], g$fold_enrichment[1], g$p_hypergeom[1]))
cat(sprintf("  250 kb loci     CTL-only %d  shared %d  LIN-only %d   %.2fx  p = %.3g\n",
            l$n_CTL[1]-l$n_shared[1], l$n_shared[1], l$n_LIN[1]-l$n_shared[1],
            l$fold_enrichment[1], l$p_hypergeom[1]))
cat("\n-- classes sharing no genes -------------------------------------------\n")
print(as.data.frame(by_class %>% filter(n_shared == 0) %>%
        dplyr::select(Class, n_CTL, n_LIN, n_shared)))
