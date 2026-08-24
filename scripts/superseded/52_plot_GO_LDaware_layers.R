#!/usr/bin/env Rscript
# LD-aware GO enrichment figures (GO-BP and GO-MF), one per trait layer.
#
# Three panels sharing the GO-term axis:
#   A  fold enrichment (log x); point size = candidate genes carrying the term;
#      hollow point = fewer than three independent 250 kb intervals
#   B  number of independent 250 kb intervals carrying the term
#   C  gene-level q (open grey circle) -> LD-aware permutation q (filled diamond)
#
# Reads the *_LDaware.tsv written by 29/33; writes fig/go_ldaware/*.png.
# Usage:  Rscript 35_plot_GO_LDaware.R [repo_root]

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(ggplot2)
  library(patchwork); library(scales)
})

# Repo root: an explicit argument wins; otherwise walk up from the working
# directory looking for data/LD_mapped/go_enrichment. This makes the script work
# both from Rscript at the repo root and from an interactive console started
# anywhere inside the repo.
find_root <- function(start = getwd()) {
  d <- normalizePath(start, mustWork = FALSE)
  for (i in 1:8) {
    if (dir.exists(file.path(d, "data/LD_mapped/go_enrichment"))) return(d)
    parent <- dirname(d); if (identical(parent, d)) break; d <- parent
  }
  NA_character_
}

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[1] else find_root()
if (is.na(root) || !dir.exists(file.path(root, "data/LD_mapped/go_enrichment")))
  stop("Could not locate the repo. Either setwd() to somewhere inside\n",
       "  SAP-Lipidomics-Database, or pass the repo root explicitly:\n",
       "    Rscript 35_plot_GO_LDaware.R /path/to/SAP-Lipidomics-Database\n",
       "  Current working directory: ", getwd(), call. = FALSE)
message("repo root: ", root)

in_dir <- file.path(root, "data/LD_mapped/go_enrichment")
out_dir<- file.path(root, "fig/go_ldaware"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

# panel-local override: plot_theme's 100pt right margin is meant for a single
# wide plot; inside a 3-panel row it strands the panel C title off-canvas.
panel_theme <- function(left = 8, right = 8)
  plot_theme + theme(plot.margin = margin(12, right, 12, left),
                     plot.title  = element_text(face = 'bold', size = 11.5, hjust = 0,
                                                margin = margin(b = 10)))

COND_COL <- c(CTL = '#2a78d6', LIN = '#eb6834')   # CVD-validated categorical pair

prep <- function(f) {
  read_tsv(f, show_col_types = FALSE) |>
    mutate(blocks = blocks_carrying_term,
           solid  = blocks >= 3,
           stroke_col = COND_COL[condition],
           # shape 21 draws fill; hollow marks the <3-interval terms
           fill_col   = ifelse(solid, COND_COL[condition], 'white'),
           label  = paste0(lipid_class, "  ", str_trunc(go_term, 60))) |>
    arrange(desc(condition), desc(lipid_class), fold_enrichment) |>
    mutate(y_id = row_number())
}

divider <- function(d) {
  b <- which(diff(as.integer(factor(d$condition))) != 0)
  if (length(b)) b[1] + .5 else NA_real_
}

panelA <- function(d, onto) {
  ont <- if (onto == "BP") "biological-process" else "molecular-function"
  ggplot(d, aes(fold_enrichment, y_id)) +
    geom_vline(xintercept = 1, linetype = 'dashed', colour = 'grey45', linewidth = .45) +
    geom_segment(aes(x = 1, xend = fold_enrichment, y = y_id, yend = y_id, colour = stroke_col),
                 linewidth = .5, alpha = .45) +
    geom_point(aes(size = term_test_count, colour = stroke_col, fill = fill_col),
               shape = 21, stroke = 1.0) +
    scale_colour_identity() + scale_fill_identity() +
    scale_size_continuous(range = c(2.4, 9), guide = 'none') +
    scale_x_log10(breaks = c(1, 2, 5, 10, 20, 40), labels = c('1','2','5','10','20','40')) +
    scale_y_continuous(breaks = d$y_id, labels = d$label,
                       expand = expansion(mult = c(.015, .015))) +
    labs(x = 'Fold enrichment (log scale)', y = NULL,
         title = paste('A   Enrichment of GO', ont, 'terms')) +
    panel_theme(left = 12) +
    theme(axis.text.y = element_text(size = 10))
}

panelB <- function(d) {
  ggplot(d, aes(blocks, y_id)) +
    geom_vline(xintercept = 3, linetype = 'dashed', colour = 'grey45', linewidth = .45) +
    geom_segment(aes(x = 0, xend = blocks, y = y_id, yend = y_id, colour = stroke_col),
                 linewidth = .5, alpha = .45) +
    geom_point(aes(colour = stroke_col, fill = fill_col), shape = 21, size = 2.9, stroke = 1.0) +
    scale_colour_identity() + scale_fill_identity() +
    scale_x_continuous(limits = c(0, max(d$blocks) + 1.5), expand = expansion(mult = c(.02, .02))) +
    scale_y_continuous(breaks = d$y_id, labels = NULL, expand = expansion(mult = c(.015, .015))) +
    labs(x = 'Independent 250 kb intervals', y = NULL, title = 'B   Distributed support') +
    panel_theme()
}

panelC <- function(d) {
  d <- d |> mutate(gene_q = -log10(p_adj_bh), ld_q = -log10(q_perm_bh))
  ggplot(d) +
    geom_vline(xintercept = -log10(.05), linetype = 'dashed', colour = '#c1272d', linewidth = .6) +
    geom_segment(aes(x = gene_q, xend = ld_q, y = y_id, yend = y_id),
                 colour = 'grey70', linewidth = .5) +
    geom_point(aes(gene_q, y_id), shape = 21, fill = 'white', colour = 'grey45',
               size = 2.4, stroke = .8) +
    geom_point(aes(ld_q, y_id, colour = stroke_col), shape = 18, size = 3.8) +
    scale_colour_identity() +
    scale_y_continuous(breaks = d$y_id, labels = NULL, expand = expansion(mult = c(.015, .015))) +
    labs(x = expression(-log[10]*(q)), y = NULL, title = 'C   Gene-level vs LD-aware q') +
    panel_theme(right = 16)
}

# scale_*_identity drops the legend, so the key is drawn as a standalone strip
legend_strip <- function() {
  k <- tibble(x = c(1, 2.0, 3.2, 6.4, 7.9), y = 1,
              lab = c('CTL', 'LIN', 'fewer than 3 independent intervals',
                      'gene-level q', 'LD-aware q'))
  ggplot(k, aes(x, y)) +
    annotate('point', x = 1, y = 1, shape = 21, size = 4, fill = COND_COL['CTL'], colour = COND_COL['CTL']) +
    annotate('point', x = 2.0, y = 1, shape = 21, size = 4, fill = COND_COL['LIN'], colour = COND_COL['LIN']) +
    annotate('point', x = 3.2, y = 1, shape = 21, size = 4, fill = 'white', colour = 'grey25') +
    annotate('point', x = 6.4, y = 1, shape = 21, size = 3, fill = 'white', colour = 'grey45') +
    annotate('point', x = 7.9, y = 1, shape = 18, size = 4, colour = 'grey35') +
    geom_text(aes(label = lab), hjust = 0, nudge_x = .10, size = 3.7, colour = 'grey20') +
    scale_x_continuous(limits = c(.85, 9.9)) +
    theme_void() + theme(plot.margin = margin(0, 12, 6, 12))
}

add_divider <- function(p, yb)
  if (is.na(yb)) p else p + geom_hline(yintercept = yb, colour = 'grey72', linewidth = .35)

build <- function(onto, layer, title, outfile) {
  f <- file.path(in_dir, sprintf('GO_%s_enrichment_by_lipid_class_%s_LDaware.tsv', onto, layer))
  if (!file.exists(f)) { warning('missing input, skipped: ', f, call. = FALSE); return(invisible(NULL)) }
  d <- prep(f); yb <- divider(d)
  row <- add_divider(panelA(d, onto), yb) | add_divider(panelB(d), yb) | add_divider(panelC(d), yb)
  fig <- (row + plot_layout(widths = c(2.55, 1.12, 1.33))) / legend_strip() +
    plot_layout(heights = c(1, .022)) +
    plot_annotation(
      title = title,
      subtitle = sprintf(paste('Point area in panel A is the number of candidate genes carrying the term (%d–%d).',
                               'Hollow points rest on fewer than three independent 250 kb intervals.'),
                         min(d$term_test_count), max(d$term_test_count)),
      theme = theme(plot.title    = element_text(face = 'bold', size = 16),
                    plot.subtitle = element_text(size = 11, colour = 'grey30',
                                                 margin = margin(b = 6)),
                    plot.margin   = margin(12, 12, 6, 12)))
  ggsave(outfile, fig, width = 15.5, height = max(6, .34 * nrow(d) + 3.2),
         dpi = 350, bg = 'white', limitsize = FALSE)
  message(sprintf('%s  (%d terms)', outfile, nrow(d)))
}

for (onto in c('BP', 'MF')) {
  lbl <- if (onto == 'BP') 'GO-BP' else 'GO-MF'
  build(onto, 'individual',
        paste(lbl, 'enrichment of individual-lipid GWAS candidates, by lipid class'),
        file.path(out_dir, sprintf('SuppFig_GO%s_LDaware_Individual.png', onto)))
  build(onto, 'sumratio',
        paste(lbl, 'enrichment of class-sum and class-ratio GWAS candidates, by lipid class'),
        file.path(out_dir, sprintf('SuppFig_GO%s_LDaware_SumRatio.png', onto)))
}
