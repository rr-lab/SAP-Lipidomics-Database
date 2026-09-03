#!/usr/bin/env Rscript
# Genomic heritability of the 13 focal lipid class sums, CTL vs LIN.
#
# Two things this adds over 22c:
#   (1) class sums are built from the species DETECTED IN BOTH TRIALS, so the
#       trait compared across trials is the same measurement in each;
#   (2) bootstrap percentile intervals over accessions, because mixed.solve
#       returns no standard error and h2 near zero is otherwise unreadable.
#
# Inputs in data/kinship/ : sap_grm.rel , sap_grm.rel.id   (plink2 --make-rel square)
# Phenotypes             : data/raw/Final_subset_{control,lowinput}_*.csv

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(readr)
  library(stringr); library(rrBLUP); library(ggplot2)
})

root   <- '/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database'
kdir   <- file.path(root, 'data', 'kinship')
tabdir <- file.path(root, 'table', 'new_table')
figdir <- file.path(root, 'fig', 'new_figures')
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

NBOOT <- as.integer(Sys.getenv('NBOOT', '300'))
set.seed(1)

plot_theme <- theme_minimal(base_size = 13) +
  theme(axis.text = element_text(colour = 'black', size = 9),
        axis.title = element_text(face = 'bold', size = 11),
        axis.line = element_line(colour = 'black'),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank())

condition_colors <- c(CTL = '#440154FF', LIN = '#FDE725FF')
focal <- c('DG','DGDG','LPC','LPE','MG','MGDG','PA','PC','PE','PG','PS','SQDG','TG')

## ---- GRM ------------------------------------------------------------------
load_k <- function() {
  K   <- as.matrix(fread(file.path(kdir, 'sap_grm.rel'), header = FALSE, data.table = FALSE))
  idf <- fread(file.path(kdir, 'sap_grm.rel.id'), header = FALSE, data.table = FALSE)
  idf <- idf[!startsWith(trimws(as.character(idf[[1]])), '#'), , drop = FALSE]
  ids <- as.character(idf[[ncol(idf)]])
  stopifnot(nrow(K) == length(ids))
  K <- (K + t(K)) / 2; dimnames(K) <- list(ids, ids)
  message(sprintf('GRM: %d genotypes | diag mean %.3f | offdiag mean %.4f',
                  nrow(K), mean(diag(K)), mean(K[upper.tri(K)])))
  list(K = K, ids = ids)
}

h2_of <- function(y, K) {
  q <- tryCatch(mixed.solve(y, K = K), error = function(e) NULL)
  if (is.null(q)) return(NA_real_)
  h <- q$Vu / (q$Vu + q$Ve)
  if (!is.finite(h)) NA_real_ else if (h < 1e-6) 0 else h
}

# bootstrap over accessions: resample lines, subset K to the resampled set.
# duplicated rows make K singular, so we jitter the diagonal slightly.
boot_ci <- function(y, ids, KIN, nboot = NBOOT) {
  i  <- match(ids, KIN$ids); ok <- !is.na(i) & is.finite(y)
  y0 <- y[ok]; i0 <- i[ok]; n <- length(y0)
  if (n < 30) return(c(NA_real_, NA_real_))
  out <- numeric(nboot)
  for (b in seq_len(nboot)) {
    s  <- sample.int(n, n, replace = TRUE)
    Kb <- KIN$K[i0[s], i0[s], drop = FALSE]
    diag(Kb) <- diag(Kb) + 1e-6
    out[b] <- h2_of(y0[s], Kb)
  }
  as.numeric(quantile(out, c(0.025, 0.975), na.rm = TRUE))
}

## ---- phenotypes -----------------------------------------------------------
read_spats <- function(cond) {
  f <- file.path(root, 'data','SPATS_fitted','non_normalized_intensities',
        sprintf('Final_subset_%s_all_lipids_fitted_phenotype_non_normalized.csv',
                if (cond == 'CTL') 'control' else 'lowinput'))
  x   <- fread(f, data.table = FALSE, check.names = FALSE)
  lip <- setdiff(names(x), c('LineRaw','PlotID','row','col'))
  lip <- lip[grepl('\\(', lip)]
  m <- as.matrix(x[, lip, drop = FALSE]); m[!is.finite(m)] <- NA_real_
  list(line = as.character(x$LineRaw), mat = m, lipids = lip)
}

CTLd <- read_spats('CTL'); LINd <- read_spats('LIN'); KIN <- load_k()
shared <- intersect(CTLd$lipids, LINd$lipids)
message(sprintf('species: CTL %d, LIN %d, shared %d',
                length(CTLd$lipids), length(LINd$lipids), length(shared)))

cls <- function(v) str_match(v, '^([A-Za-z]+)\\(')[, 2]

estimate <- function(d, cond, species_set, label) {
  use <- d$lipids %in% species_set
  k   <- cls(d$lipids)
  bind_rows(lapply(focal, function(cl) {
    idx <- which(use & !is.na(k) & k == cl)
    if (!length(idx)) return(NULL)
    y  <- log10(rowSums(d$mat[, idx, drop = FALSE]) + 1)
    i  <- match(d$line, KIN$ids); ok <- !is.na(i) & is.finite(y)
    h  <- h2_of(y[ok], KIN$K[i[ok], i[ok], drop = FALSE])
    ci <- boot_ci(y, d$line, KIN)
    tibble(trial = cond, species_set = label, Class = cl, n_species = length(idx),
           n = sum(ok), h2 = h, lo = ci[1], hi = ci[2])
  }))
}

all_res <- bind_rows(
  estimate(CTLd, 'CTL', CTLd$lipids, 'all detected'),
  estimate(LINd, 'LIN', LINd$lipids, 'all detected'),
  estimate(CTLd, 'CTL', shared,      'shared only'),
  estimate(LINd, 'LIN', shared,      'shared only')
)

write_csv(all_res, file.path(tabdir, 'SuppTable_Genomic_Heritability_Focal_Bootstrap.csv'))

sh <- all_res |> filter(species_set == 'shared only')
cat('\n========== h2, class sums over SHARED species only (', NBOOT, ' bootstraps) ==========\n', sep = '')
print(as.data.frame(sh |> select(trial, Class, n_species, n, h2, lo, hi) |>
                    arrange(trial, desc(h2))), row.names = FALSE, digits = 3)

wide <- sh |> select(trial, Class, h2, lo, hi) |>
  pivot_wider(names_from = trial, values_from = c(h2, lo, hi)) |>
  mutate(nonzero_CTL = lo_CTL > 0, nonzero_LIN = lo_LIN > 0,
         LIN_gt_CTL  = lo_LIN > hi_CTL)
cat('\n---- per class, shared species ----\n')
print(as.data.frame(wide |> select(Class, h2_CTL, lo_CTL, hi_CTL, h2_LIN, lo_LIN, hi_LIN,
                                   nonzero_CTL, nonzero_LIN, LIN_gt_CTL)),
      row.names = FALSE, digits = 3)
cat('\nclasses whose CI excludes zero   CTL:', sum(wide$nonzero_CTL, na.rm = TRUE),
    ' LIN:', sum(wide$nonzero_LIN, na.rm = TRUE), 'of', nrow(wide), '\n')
cat('classes where LIN CI is entirely above CTL CI:',
    sum(wide$LIN_gt_CTL, na.rm = TRUE), 'of', nrow(wide), '\n')
cat('median h2   CTL:', round(median(sh$h2[sh$trial=='CTL'], na.rm=TRUE), 3),
    '  LIN:',           round(median(sh$h2[sh$trial=='LIN'], na.rm=TRUE), 3), '\n')

p <- ggplot(sh, aes(h2, reorder(Class, h2), colour = trial)) +
  geom_linerange(aes(xmin = lo, xmax = hi),
                 position = position_dodge(width = 0.6), linewidth = 0.6) +
  geom_point(position = position_dodge(width = 0.6), size = 2.4) +
  scale_colour_manual(values = condition_colors) +
  labs(x = expression('Genomic heritability ('*h^2*')'), y = NULL) +
  plot_theme + theme(legend.position = 'bottom')
ggsave(file.path(figdir, 'SuppFig_Heritability_SharedSpecies_CI.png'), p,
       width = 8, height = 6, dpi = 350, bg = 'white')
message('done')
