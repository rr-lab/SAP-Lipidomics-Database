#!/usr/bin/env Rscript
# Heritability figure.
#   A  class-sum h2 with profile-likelihood 95% CI, CTL vs LIN (13 focal classes)
#   B  per-species paired h2, CTL vs LIN (152 species shared by both trials)
#
# Marks: LIN's fill (#FDE725) sits at 1.23:1 against white, so every point is a
# filled shape 21 with a thin black stroke and the intervals are drawn in grey.
# Identity comes from the fill; legibility from the stroke.

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(readr)
  library(stringr); library(ggplot2); library(patchwork)
})

root   <- '/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database'
tabdir <- file.path(root, 'table', 'new_table')
figdir <- file.path(root, 'fig', 'new_figures')
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

plot_theme <- theme_minimal(base_size = 13) +
  theme(axis.text = element_text(colour = 'black', size = 9),
        axis.title = element_text(face = 'bold', size = 11),
        axis.line = element_line(colour = 'black'),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank())

condition_colors <- c(CTL = '#440154FF', LIN = '#FDE725FF')

## ---- A: class sums with CI ------------------------------------------------
cls <- read_csv(file.path(tabdir, 'SuppTable_Genomic_Heritability_ProfileCI.csv'),
                show_col_types = FALSE)
ord <- cls |> filter(trial == 'LIN') |> arrange(h2) |> pull(Class)
cls <- cls |> mutate(Class = factor(Class, levels = ord),
                     trial = factor(trial, levels = c('CTL','LIN')))

pA <- ggplot(cls, aes(h2, Class, group = trial)) +
  geom_linerange(aes(xmin = lo, xmax = hi),
                 position = position_dodge(width = 0.65),
                 colour = 'grey55', linewidth = 0.5) +
  geom_point(aes(fill = trial), position = position_dodge(width = 0.65),
             shape = 21, size = 2.9, stroke = 0.4, colour = 'black') +
  scale_fill_manual(values = condition_colors) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0.01, 0.03))) +
  labs(x = expression('Genomic heritability ('*h^2*') of class sum'), y = NULL) +
  plot_theme + theme(legend.position = 'bottom')

## ---- B: per-species paired ------------------------------------------------
sp <- read_csv(file.path(tabdir, 'SuppTable_Heritability_PerSpecies_Paired.csv'),
               show_col_types = FALSE) |>
  select(trial, Species, h2) |>
  pivot_wider(names_from = trial, values_from = h2) |>
  filter(is.finite(CTL), is.finite(LIN))

n_up <- sum(sp$LIN > sp$CTL); n_dn <- sum(sp$LIN < sp$CTL)
lim  <- c(0, max(c(sp$CTL, sp$LIN), na.rm = TRUE) * 1.04)

pB <- ggplot(sp, aes(CTL, LIN)) +
  geom_abline(slope = 1, intercept = 0, colour = 'grey60', linetype = 2, linewidth = 0.4) +
  geom_point(shape = 21, fill = 'grey75', colour = 'black',
             stroke = 0.25, size = 2.0, alpha = 0.85) +
  annotate('text', x = lim[1] + 0.04 * diff(lim), y = lim[2] * 0.96,
           label = paste0('higher under LIN\n', n_up, ' species'),
           hjust = 0, vjust = 1, size = 3.1, colour = 'grey20') +
  annotate('text', x = lim[2] * 0.97, y = lim[1] + 0.02 * diff(lim),
           label = paste0('higher under CTL\n', n_dn, ' species'),
           hjust = 1, vjust = 0, size = 3.1, colour = 'grey20') +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = expression('CTL  '*h^2), y = expression('LIN  '*h^2)) +
  plot_theme + theme(panel.grid.major.y = element_line(colour = 'grey92'))

fig <- (pA | pB) +
  plot_layout(widths = c(1.05, 1)) +
  plot_annotation(tag_levels = 'A',
                  theme = theme(plot.tag = element_text(face = 'bold', size = 15)))

ggsave(file.path(figdir, 'Fig_Heritability_CTL_LIN.png'), fig,
       width = 12, height = 5.6, dpi = 350, bg = 'white')
ggsave(file.path(figdir, 'Fig_Heritability_CTL_LIN.pdf'), fig,
       width = 12, height = 5.6, bg = 'white')
cat('species plotted:', nrow(sp), ' LIN higher:', n_up, ' CTL higher:', n_dn, '\n')
message('wrote fig/new_figures/Fig_Heritability_CTL_LIN.{png,pdf}')
