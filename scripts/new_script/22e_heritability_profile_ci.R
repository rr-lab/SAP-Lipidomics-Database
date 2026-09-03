#!/usr/bin/env Rscript
# Genomic heritability of the 13 focal lipid class sums, CTL vs LIN, on the
# species detected in BOTH trials, with profile-likelihood confidence intervals.
#
# Why not bootstrap: resampling accessions with replacement creates duplicate
# individuals whose kinship with each other is 1 and whose phenotypes are
# identical, so the genetic term fits them exactly, Ve collapses and h2 -> 1 in
# every replicate. Profile likelihood avoids resampling entirely.
#
# Model  y = 1*mu + u + e ,  var(u) = sg2 * K ,  var(e) = se2 * I
# Reparameterised by h2 = sg2/(sg2+se2); REML profiled over a grid of h2, and
# the 95% CI is the set with 2*(llmax - ll) <= qchisq(0.95, 1).
#
# Inputs in data/kinship/ : sap_grm.rel , sap_grm.rel.id

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(readr)
  library(stringr); library(ggplot2)
})

root   <- '/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database'
kdir   <- file.path(root, 'data', 'kinship')
tabdir <- file.path(root, 'table', 'new_table')
figdir <- file.path(root, 'fig', 'new_figures')
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

plot_theme <- theme_minimal(base_size = 13) +
  theme(axis.text = element_text(colour = 'black', size = 9),
        axis.title = element_text(face = 'bold', size = 11),
        axis.line = element_line(colour = 'black'),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank())

condition_colors <- c(CTL = '#440154FF', LIN = '#FDE725FF')
focal <- c('DG','DGDG','LPC','LPE','MG','MGDG','PA','PC','PE','PG','PS','SQDG','TG')
GRID  <- seq(0, 0.999, by = 0.001)

load_k <- function() {
  K   <- as.matrix(fread(file.path(kdir, 'sap_grm.rel'), header = FALSE, data.table = FALSE))
  idf <- fread(file.path(kdir, 'sap_grm.rel.id'), header = FALSE, data.table = FALSE)
  idf <- idf[!startsWith(trimws(as.character(idf[[1]])), '#'), , drop = FALSE]
  ids <- as.character(idf[[ncol(idf)]]); stopifnot(nrow(K) == length(ids))
  K <- (K + t(K)) / 2; dimnames(K) <- list(ids, ids)
  message(sprintf('GRM: %d genotypes | diag mean %.3f | offdiag mean %.4f',
                  nrow(K), mean(diag(K)), mean(K[upper.tri(K)])))
  list(K = K, ids = ids)
}

# REML profile over h2 using one eigendecomposition of K
prof_h2 <- function(y, K) {
  n <- length(y); p <- 1
  e <- eigen(K, symmetric = TRUE)
  d <- pmax(e$values, 0)                       # clamp tiny negative eigenvalues
  yt <- as.numeric(crossprod(e$vectors, y))
  xt <- as.numeric(crossprod(e$vectors, rep(1, n)))
  ll <- vapply(GRID, function(h) {
    v <- h * d + (1 - h)
    if (any(v <= 1e-12)) return(-Inf)
    w <- 1 / v
    XtWX <- sum(w * xt * xt)
    if (!is.finite(XtWX) || XtWX <= 0) return(-Inf)
    b  <- sum(w * xt * yt) / XtWX
    r  <- yt - xt * b
    s2 <- sum(w * r * r) / (n - p)
    if (!is.finite(s2) || s2 <= 0) return(-Inf)
    -0.5 * (sum(log(v)) + (n - p) * log(s2) + log(XtWX))
  }, numeric(1))
  i    <- which.max(ll)
  keep <- which(ll >= ll[i] - qchisq(0.95, 1) / 2)
  c(h2 = GRID[i], lo = GRID[min(keep)], hi = GRID[max(keep)],
    ll_flat = as.numeric(ll[i] - ll[1]))       # gain over h2 = 0
}

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

estimate <- function(d, cond) {
  k <- cls(d$lipids); use <- d$lipids %in% shared
  bind_rows(lapply(focal, function(cl) {
    idx <- which(use & !is.na(k) & k == cl)
    if (!length(idx)) return(NULL)
    y <- log10(rowSums(d$mat[, idx, drop = FALSE]) + 1)
    i <- match(d$line, KIN$ids); ok <- !is.na(i) & is.finite(y)
    if (sum(ok) < 30) return(NULL)
    r <- prof_h2(y[ok], KIN$K[i[ok], i[ok], drop = FALSE])
    tibble(trial = cond, Class = cl, n_species = length(idx), n = sum(ok),
           h2 = r[['h2']], lo = r[['lo']], hi = r[['hi']], ll_gain = r[['ll_flat']])
  }))
}

res <- bind_rows(estimate(CTLd, 'CTL'), estimate(LINd, 'LIN')) |>
  mutate(nonzero = lo > 0)
write_csv(res, file.path(tabdir, 'SuppTable_Genomic_Heritability_ProfileCI.csv'))

cat('\n===== h2 with profile-likelihood 95% CI, shared species only =====\n')
print(as.data.frame(res |> select(trial, Class, n_species, n, h2, lo, hi, nonzero) |>
                    arrange(trial, desc(h2))), row.names = FALSE, digits = 3)

w <- res |> select(trial, Class, h2, lo, hi) |>
  pivot_wider(names_from = trial, values_from = c(h2, lo, hi)) |>
  mutate(LIN_gt_CTL = lo_LIN > hi_CTL, diff = h2_LIN - h2_CTL)
cat('\n---- paired by class ----\n')
print(as.data.frame(w |> select(Class, h2_CTL, lo_CTL, hi_CTL,
                                h2_LIN, lo_LIN, hi_LIN, diff, LIN_gt_CTL)),
      row.names = FALSE, digits = 3)

cat('\nCI excludes zero    CTL:', sum(res$nonzero[res$trial=='CTL']),
    ' LIN:', sum(res$nonzero[res$trial=='LIN']), 'of', nrow(w), '\n')
cat('LIN CI entirely above CTL CI:', sum(w$LIN_gt_CTL, na.rm = TRUE), 'of', nrow(w), '\n')
cat('classes with LIN > CTL point estimate:', sum(w$diff > 0, na.rm = TRUE), 'of', nrow(w), '\n')
cat('median h2   CTL:', round(median(res$h2[res$trial=='CTL']), 3),
    '  LIN:',           round(median(res$h2[res$trial=='LIN']), 3), '\n')
st <- suppressWarnings(wilcox.test(w$h2_LIN, w$h2_CTL, paired = TRUE))
cat(sprintf('paired Wilcoxon across the %d classes: V = %s, p = %.4g\n',
            nrow(w), st$statistic, st$p.value))

p <- ggplot(res, aes(h2, reorder(Class, h2), colour = trial)) +
  geom_linerange(aes(xmin = lo, xmax = hi),
                 position = position_dodge(width = 0.6), linewidth = 0.6) +
  geom_point(position = position_dodge(width = 0.6), size = 2.4) +
  scale_colour_manual(values = condition_colors) +
  labs(x = expression('Genomic heritability ('*h^2*')  with profile-likelihood 95% CI'),
       y = NULL) +
  plot_theme + theme(legend.position = 'bottom')
ggsave(file.path(figdir, 'SuppFig_Heritability_ProfileCI.png'), p,
       width = 8, height = 6, dpi = 350, bg = 'white')
message('done')
