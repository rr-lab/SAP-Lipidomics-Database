#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Chapter 2 GO enrichment outputs, built from the LD-aware permutation results.
#
#   Fig_GO_BP_main.png    main text  -- GO-BP, both trait layers, terms carried
#                                       by >= 3 independent 250 kb intervals
#   SuppFig_GO_BP_all.png supplement -- GO-BP, every significant term
#   SuppFig_GO_MF_all.png supplement -- GO-MF, every significant term
#   Table_GO_enrichment_all.tsv / .tex  -- every term, both ontologies
#
# No titles or subtitles are drawn: that text belongs in the figure caption.
# Panel tags A/B/C are kept because captions refer to them.
#
# Usage:  Rscript 36_GO_figures_and_table.R [repo_root]
#         (run from anywhere inside the repo and the root is found for you)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(patchwork); library(scales)
})

find_root <- function(start = getwd()) {
  d <- normalizePath(start, mustWork = FALSE)
  for (i in 1:8) {
    if (dir.exists(file.path(d, "data/LD_mapped/go_enrichment"))) return(d)
    p <- dirname(d); if (identical(p, d)) break; d <- p
  }
  NA_character_
}
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[1] else find_root()
if (is.na(root) || !dir.exists(file.path(root, "data/LD_mapped/go_enrichment")))
  stop("Could not locate the repo. setwd() inside SAP-Lipidomics-Database, or pass the root.\n",
       "  cwd: ", getwd(), call. = FALSE)
message("repo root: ", root)

in_dir  <- file.path(root, "data/LD_mapped/go_enrichment")
fig_dir <- file.path(root, "fig/go_ldaware");  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
tab_dir <- file.path(root, "table/go_enrichment"); dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

plot_theme <- theme_minimal(base_size = 14) +
  theme(axis.text    = element_text(colour = 'black', size = 11),
        axis.title   = element_text(face = 'bold', size = 14),
        axis.line    = element_line(colour = 'black', linewidth = .5),
        panel.grid   = element_blank(),
        legend.title = element_text(size = 12),
        legend.text  = element_text(size = 11),
        plot.margin  = margin(12, 100, 12, 12),
        plot.tag     = element_text(face = 'bold', size = 16),
        plot.tag.position = c(0, 1))

COND_COL <- c(CTL = '#2a78d6', LIN = '#eb6834')   # CVD-validated categorical pair
MIN_INTERVALS <- 3

# ---- load -----------------------------------------------------------------
# term_bg_count lives in the primary enrichment table, the permutation columns
# in the LD-aware one; join so the table can report both.
read_layer <- function(onto, layer) {
  ld  <- file.path(in_dir, sprintf('GO_%s_enrichment_by_lipid_class_%s_LDaware.tsv', onto, layer))
  raw <- file.path(in_dir, sprintf('GO_%s_enrichment_by_lipid_class_%s_LD.tsv',      onto, layer))
  if (!file.exists(ld)) { warning('missing: ', ld, call. = FALSE); return(NULL) }
  d <- read_tsv(ld, show_col_types = FALSE)
  if (file.exists(raw)) {
    bg <- read_tsv(raw, show_col_types = FALSE) |>
      select(condition, lipid_class, go_term, term_bg_count) |> distinct()
    d <- left_join(d, bg, by = c('condition', 'lipid_class', 'go_term'))
  } else d$term_bg_count <- NA_integer_
  d |> mutate(ontology = onto,
              layer = if (layer == 'individual') 'Individual lipid species'
                      else 'Class sums and ratios')
}

load_onto <- function(onto)
  bind_rows(read_layer(onto, 'individual'), read_layer(onto, 'sumratio')) |>
    mutate(intervals = blocks_carrying_term,
           solid     = intervals >= MIN_INTERVALS,
           stroke_col = COND_COL[condition],
           fill_col   = ifelse(solid, COND_COL[condition], 'white'))

# ---- figure ---------------------------------------------------------------
# rows are ordered within each layer facet; y is a plain integer so all three
# panels share one scale, and the facet split is by trait layer.
order_rows <- function(d, wrap)
  d |>
    mutate(layer = factor(layer, levels = c('Individual lipid species', 'Class sums and ratios'))) |>
    arrange(layer, desc(condition), desc(lipid_class), fold_enrichment) |>
    group_by(layer) |>
    mutate(y_id = row_number()) |>
    ungroup() |>
    mutate(row_key = factor(paste0(as.integer(layer), '_', y_id),
                            levels = paste0(as.integer(layer), '_', y_id)),
           label   = paste0(lipid_class, "  ", str_trunc(go_term, wrap)))

pt <- function(base, left = 8, right = 8)
  plot_theme +
  theme(plot.margin = margin(10, right, 10, left),
        plot.title  = element_text(face = 'bold', size = base + .5, hjust = 0, margin = margin(b = 8)),
        axis.text   = element_text(colour = 'black', size = base - 1),
        axis.title.x = element_text(face = 'bold', size = base + 1, margin = margin(t = 9)),
        axis.title.y = element_blank(),
        strip.text.y.left = element_text(angle = 90, face = 'bold', size = base, margin = margin(r = 6)),
        strip.placement = 'outside',
        panel.spacing.y = unit(14, 'pt'))

facet_layers <- function(p, show_strip) {
  p <- p + facet_grid(layer ~ ., scales = 'free_y', space = 'free_y',
                      switch = 'y', labeller = label_wrap_gen(18))
  if (show_strip) p else p + theme(strip.text.y.left = element_blank())
}

panelA <- function(d, onto, base) {
  ont <- if (onto == 'BP') 'biological-process' else 'molecular-function'
  p <- ggplot(d, aes(fold_enrichment, row_key)) +
    geom_vline(xintercept = 1, linetype = 'dashed', colour = 'grey45', linewidth = .45) +
    geom_segment(aes(x = 1, xend = fold_enrichment, y = row_key, yend = row_key, colour = stroke_col),
                 linewidth = .45, alpha = .45) +
    geom_point(aes(size = term_test_count, colour = stroke_col, fill = fill_col),
               shape = 21, stroke = .9) +
    scale_colour_identity() + scale_fill_identity() +
    scale_size_continuous(range = c(2.0, 8), guide = 'none') +
    scale_x_log10(breaks = c(1, 2, 5, 10, 20, 40), labels = c('1','2','5','10','20','40')) +
    scale_y_discrete(labels = setNames(d$label, as.character(d$row_key)),
                     expand = expansion(add = .6)) +
    labs(x = 'Fold enrichment (log scale)', y = NULL,
         title = paste0('A   GO ', ont, ' enrichment')) +
    pt(base, left = 12)
  facet_layers(p, TRUE)
}

panelB <- function(d, base) {
  p <- ggplot(d, aes(intervals, row_key)) +
    geom_vline(xintercept = MIN_INTERVALS, linetype = 'dashed', colour = 'grey45', linewidth = .45) +
    geom_segment(aes(x = 0, xend = intervals, y = row_key, yend = row_key, colour = stroke_col),
                 linewidth = .45, alpha = .45) +
    geom_point(aes(colour = stroke_col, fill = fill_col), shape = 21, size = base * .26, stroke = .9) +
    scale_colour_identity() + scale_fill_identity() +
    scale_x_continuous(limits = c(0, max(d$intervals) + 1.5), expand = expansion(mult = c(.02, .02))) +
    scale_y_discrete(labels = NULL, expand = expansion(add = .6)) +
    labs(x = 'Independent 250 kb intervals', y = NULL, title = 'B   Distributed support') +
    pt(base)
  facet_layers(p, FALSE)
}

panelC <- function(d, base) {
  d <- d |> mutate(gene_q = -log10(p_adj_bh), ld_q = -log10(q_perm_bh))
  p <- ggplot(d) +
    geom_vline(xintercept = -log10(.05), linetype = 'dashed', colour = '#c1272d', linewidth = .55) +
    geom_segment(aes(x = gene_q, xend = ld_q, y = row_key, yend = row_key), colour = 'grey70', linewidth = .45) +
    geom_point(aes(gene_q, row_key), shape = 21, fill = 'white', colour = 'grey45',
               size = base * .22, stroke = .7) +
    geom_point(aes(ld_q, row_key, colour = stroke_col), shape = 18, size = base * .34) +
    scale_colour_identity() +
    scale_y_discrete(labels = NULL, expand = expansion(add = .6)) +
    labs(x = expression(-log[10]*(q)), y = NULL, title = 'C   LD-aware correction') +
    pt(base, right = 16)
  facet_layers(p, FALSE)
}

legend_strip <- function(base, any_hollow = TRUE) {
  k <- tibble(x = c(1, 2.0, 3.2, 6.6, 8.1), y = 1,
              lab = c('CTL', 'LIN', sprintf('fewer than %d independent intervals', MIN_INTERVALS),
                      'gene-level q', 'LD-aware q'))
  if (!any_hollow) k <- k[-3, ] |> mutate(x = c(1, 2.0, 3.4, 4.9))
  ggplot(k, aes(x, y)) +
    annotate('point', x = 1,   y = 1, shape = 21, size = 4, fill = COND_COL['CTL'], colour = COND_COL['CTL']) +
    annotate('point', x = 2.0, y = 1, shape = 21, size = 4, fill = COND_COL['LIN'], colour = COND_COL['LIN']) +
    {if (any_hollow) annotate('point', x = 3.2, y = 1, shape = 21, size = 4,
                              fill = 'white', colour = 'grey25') else NULL} +
    annotate('point', x = if (any_hollow) 6.6 else 3.4, y = 1, shape = 21, size = 3,
             fill = 'white', colour = 'grey45') +
    annotate('point', x = if (any_hollow) 8.1 else 4.9, y = 1, shape = 18, size = 4,
             colour = 'grey35') +
    geom_text(aes(label = lab), hjust = 0, nudge_x = .10, size = base * .27, colour = 'grey20') +
    scale_x_continuous(limits = c(.85, if (any_hollow) 10.2 else 6.6)) +
    theme_void() + theme(plot.margin = margin(2, 12, 8, 12))
}

build_fig <- function(d, onto, outfile, pitch, base, wrap, width) {
  if (!nrow(d)) { warning('nothing to plot for ', outfile, call. = FALSE); return(invisible(NULL)) }
  d <- order_rows(d, wrap)
  fig <- ((panelA(d, onto, base) | panelB(d, base) | panelC(d, base)) +
            plot_layout(widths = c(2.50, 1.10, 1.30))) /
         legend_strip(base, any(!d$solid)) +
         plot_layout(heights = unit(c(1, 0.42), c('null', 'in')))
  h <- pitch * nrow(d) + 2.6
  ggsave(outfile, fig, width = width, height = h, dpi = 350, bg = 'white', limitsize = FALSE)
  message(sprintf('%-52s %2d terms  %4.1f x %4.1f in', basename(outfile), nrow(d), width, h))
}

BP <- load_onto('BP')
MF <- load_onto('MF')

# main text: only the terms that keep >= 3 independent intervals
build_fig(filter(BP, intervals >= MIN_INTERVALS), 'BP',
          file.path(fig_dir, 'Fig_GO_BP_main.png'),
          pitch = .34, base = 12, wrap = 58, width = 15.5)

# supplement: everything, tighter pitch so a full page stays legible
build_fig(BP, 'BP', file.path(fig_dir, 'SuppFig_GO_BP_all.png'),
          pitch = .21, base = 9.5, wrap = 54, width = 14.0)
build_fig(MF, 'MF', file.path(fig_dir, 'SuppFig_GO_MF_all.png'),
          pitch = .21, base = 9.5, wrap = 54, width = 14.0)

# ---- table ----------------------------------------------------------------
tab <- bind_rows(BP, MF) |>
  transmute(Ontology = ontology, Layer = layer, Condition = condition,
            Class = lipid_class, Term = go_term, GO = go_id,
            Genes = term_test_count, Background = term_bg_count,
            Intervals = intervals,
            Fold = round(as.numeric(fold_enrichment), 1),
            q_gene = signif(as.numeric(p_adj_bh), 2),
            q_LD   = signif(as.numeric(q_perm_bh), 2),
            Overlap_genes = overlap_genes) |>
  arrange(Ontology, Layer, Condition, Class, desc(Intervals), q_LD)

write_tsv(tab, file.path(tab_dir, 'Table_GO_enrichment_all.tsv'))

esc <- function(x) x |> str_replace_all('\\\\', '\\\\textbackslash{}') |>
  str_replace_all('([&%$#_{}])', '\\\\\\1')
body <- tab |>
  select(Ontology, Layer, Condition, Class, Term, Genes, Background, Intervals, Fold, q_gene, q_LD) |>
  mutate(across(where(is.character), esc),
         Layer = recode(Layer, 'Individual lipid species' = 'Individual',
                        'Class sums and ratios' = 'Sum/ratio')) |>
  apply(1, function(r) paste0('  ', paste(trimws(r), collapse = ' & '), ' \\\\')) |>
  paste(collapse = '\n')

tex <- paste0(
'% Auto-generated by 36_GO_figures_and_table.R -- do not edit by hand.\n',
'% Requires: \\usepackage{longtable,booktabs}\n',
'\\begin{longtable}{llllp{5.2cm}rrrrrr}\n',
'\\caption[GO enrichment of GWAS candidate sets]{\\textbf{Gene Ontology enrichment of\n',
'LD-mapped GWAS candidate sets, all significant terms.} \\emph{Genes} is the number of\n',
'candidate genes carrying the term, \\emph{Background} the number of sorghum genes\n',
'carrying it genome-wide, and \\emph{Intervals} the number of independent 250\\,kb\n',
'genomic intervals those candidate genes occupy. $q_{\\mathrm{gene}}$ is the\n',
'Benjamini--Hochberg value from the standard gene-sampling test and\n',
'$q_{\\mathrm{LD}}$ the value from the LD-aware permutation null. Terms with fewer\n',
'than three independent intervals are single loci and should be read as such.}\\\\\n',
'\\label{tab:go_enrichment_all}\\\\\n',
'\\toprule\n',
'Ont. & Layer & Cond. & Class & Term & Genes & Bg & Int. & Fold & $q_{\\mathrm{gene}}$ & $q_{\\mathrm{LD}}$ \\\\\n',
'\\midrule\n\\endfirsthead\n',
'\\toprule\n',
'Ont. & Layer & Cond. & Class & Term & Genes & Bg & Int. & Fold & $q_{\\mathrm{gene}}$ & $q_{\\mathrm{LD}}$ \\\\\n',
'\\midrule\n\\endhead\n',
'\\midrule \\multicolumn{11}{r}{\\emph{continued on next page}}\\\\ \\endfoot\n',
'\\bottomrule\n\\endlastfoot\n',
body, '\n\\end{longtable}\n')
writeLines(tex, file.path(tab_dir, 'Table_GO_enrichment_all.tex'))

message(sprintf('\ntable: %d rows -> %s{.tsv,.tex}', nrow(tab),
                file.path(tab_dir, 'Table_GO_enrichment_all')))
message(sprintf('main figure keeps %d of %d GO-BP terms (>= %d intervals)',
                sum(BP$intervals >= MIN_INTERVALS), nrow(BP), MIN_INTERVALS))
