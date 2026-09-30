#!/usr/bin/env Rscript
# Per-species genomic heritability, CTL vs LIN, as a paired design.
#
# Fixes the two weaknesses of the class-sum analysis (22e):
#   - 152 shared species instead of 13 correlated class sums, so the paired
#     test has many more, and far less dependent, units;
#   - the SAME genotypes in both trials, so precision is identical and the
#     comparison is genuinely paired.
#
# h2 by REML profile over a grid; 95% CI from 2*(llmax - ll) <= qchisq(.95,1).

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
        panel.grid.minor = element_blank(), legend.title = element_blank())
condition_colors <- c(CTL = '#440154FF', LIN = '#FDE725FF')
GRID <- seq(0, 0.999, by = 0.001)

load_k <- function() {
  K   <- as.matrix(fread(file.path(kdir,'sap_grm.rel'), header=FALSE, data.table=FALSE))
  idf <- fread(file.path(kdir,'sap_grm.rel.id'), header=FALSE, data.table=FALSE)
  idf <- idf[!startsWith(trimws(as.character(idf[[1]])), '#'), , drop=FALSE]
  ids <- as.character(idf[[ncol(idf)]]); stopifnot(nrow(K)==length(ids))
  K <- (K+t(K))/2; dimnames(K) <- list(ids, ids); list(K=K, ids=ids)
}

prof_h2 <- function(y, K) {
  n <- length(y); p <- 1
  e <- eigen(K, symmetric = TRUE); d <- pmax(e$values, 0)
  yt <- as.numeric(crossprod(e$vectors, y))
  xt <- as.numeric(crossprod(e$vectors, rep(1, n)))
  ll <- vapply(GRID, function(h) {
    v <- h*d + (1-h); if (any(v <= 1e-12)) return(-Inf)
    w <- 1/v; XtWX <- sum(w*xt*xt); if (!is.finite(XtWX) || XtWX <= 0) return(-Inf)
    b <- sum(w*xt*yt)/XtWX; r <- yt - xt*b
    s2 <- sum(w*r*r)/(n-p); if (!is.finite(s2) || s2 <= 0) return(-Inf)
    -0.5*(sum(log(v)) + (n-p)*log(s2) + log(XtWX))
  }, numeric(1))
  i <- which.max(ll); keep <- which(ll >= ll[i] - qchisq(0.95,1)/2)
  c(GRID[i], GRID[min(keep)], GRID[max(keep)])
}

read_spats <- function(cond) {
  f <- file.path(root,'data','raw',
        sprintf('Final_subset_%s_all_lipids_fitted_phenotype_non_normalized.csv',
                if (cond=='CTL') 'control' else 'lowinput'))
  x <- fread(f, data.table=FALSE, check.names=FALSE)
  lip <- setdiff(names(x), c('LineRaw','PlotID','row','col')); lip <- lip[grepl('\\(', lip)]
  m <- as.matrix(x[, lip, drop=FALSE]); m[!is.finite(m)] <- NA_real_
  rownames(m) <- as.character(x$LineRaw)
  list(mat = m, lipids = lip)
}

CTLd <- read_spats('CTL'); LINd <- read_spats('LIN'); KIN <- load_k()
species <- intersect(CTLd$lipids, LINd$lipids)
lines   <- Reduce(intersect, list(rownames(CTLd$mat), rownames(LINd$mat), KIN$ids))
message(sprintf('shared species: %d | shared genotypes with GRM: %d',
                length(species), length(lines)))

Ksub <- KIN$K[lines, lines, drop = FALSE]
run <- function(d, cond) {
  M <- d$mat[lines, species, drop = FALSE]
  bind_rows(lapply(seq_along(species), function(j) {
    y <- log10(M[, j] + 1); ok <- is.finite(y)
    if (sum(ok) < 100) return(NULL)
    r <- prof_h2(y[ok], Ksub[ok, ok, drop = FALSE])
    tibble(trial = cond, Species = species[j],
           Class = str_match(species[j], '^([A-Za-z]+)\\(')[,2],
           n = sum(ok), h2 = r[1], lo = r[2], hi = r[3])
  }))
}

res <- bind_rows(run(CTLd,'CTL'), run(LINd,'LIN')) |> mutate(nonzero = lo > 0)
write_csv(res, file.path(tabdir,'SuppTable_Heritability_PerSpecies_Paired.csv'))

w <- res |> select(trial, Species, Class, h2) |>
  pivot_wider(names_from = trial, values_from = h2) |>
  filter(is.finite(CTL), is.finite(LIN)) |> mutate(diff = LIN - CTL)

cat('\n================ per-species heritability, paired ================\n')
cat('species compared      :', nrow(w), '\n')
cat('genotypes per trial   :', length(lines), '\n')
cat('median h2   CTL:', round(median(w$CTL),3), '  LIN:', round(median(w$LIN),3), '\n')
cat('mean   h2   CTL:', round(mean(w$CTL),3),   '  LIN:', round(mean(w$LIN),3), '\n')
cat('CI excludes zero  CTL:', sum(res$nonzero[res$trial=='CTL']),
    ' LIN:', sum(res$nonzero[res$trial=='LIN']), 'of', nrow(w), '\n')
cat('LIN > CTL         :', sum(w$diff > 0), ' CTL > LIN:', sum(w$diff < 0),
    ' ties:', sum(w$diff == 0), '\n')
nb <- sum(w$diff != 0)
bt <- binom.test(sum(w$diff > 0), nb, 0.5)
wt <- suppressWarnings(wilcox.test(w$LIN, w$CTL, paired = TRUE))
cat(sprintf('sign test (%d non-ties): p = %.3g\n', nb, bt$p.value))
cat(sprintf('paired Wilcoxon        : V = %s, p = %.3g\n', wt$statistic, wt$p.value))

cat('\n---- by class (median h2) ----\n')
byc <- w |> group_by(Class) |>
  summarise(n_species = n(), CTL = median(CTL), LIN = median(LIN),
            LIN_higher = sum(diff > 0), .groups='drop') |> arrange(desc(n_species))
print(as.data.frame(byc), row.names = FALSE, digits = 3)

p1 <- ggplot(w, aes(CTL, LIN)) +
  geom_abline(slope=1, intercept=0, colour='grey60', linetype=2) +
  geom_point(alpha=.65, size=1.8, colour='#440154FF') +
  coord_equal() +
  labs(x=expression('CTL  '*h^2), y=expression('LIN  '*h^2)) + plot_theme
p2 <- res |> ggplot(aes(h2, fill = trial)) +
  geom_histogram(bins=30, position='identity', alpha=.6, colour='black', linewidth=.2) +
  scale_fill_manual(values=condition_colors) +
  labs(x=expression('per-species '*h^2), y='species') +
  plot_theme + theme(legend.position='bottom')
ggsave(file.path(figdir,'SuppFig_Heritability_PerSpecies_Scatter.png'), p1,
       width=6, height=6, dpi=350, bg='white')
ggsave(file.path(figdir,'SuppFig_Heritability_PerSpecies_Hist.png'), p2,
       width=7, height=5, dpi=350, bg='white')
message('done')
