#!/usr/bin/env Rscript
# Genomic heritability of the 13 focal lipid class sums, computed from our own
# GRM (LD-pruned markers from the SAP VCF, built by grm/build_grm.R on the
# server) and the SpATS-adjusted phenotypes.
#
# Inputs expected in data/kinship/ :
#   sap_grm.rel      square GRM from  plink2 --make-rel square
#   sap_grm.rel.id   FID/IID, one line per row of the matrix
# Phenotypes read from data/raw/Final_subset_{control,lowinput}_*.csv
# One GRM serves both trials; only the phenotype vector changes.

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

plot_theme <- theme_minimal(base_size = 13) +
  theme(axis.text = element_text(colour = 'black', size = 9),
        axis.title = element_text(face = 'bold', size = 11),
        axis.line = element_line(colour = 'black'),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank())

condition_colors <- c(CTL = '#440154FF', LIN = '#FDE725FF')
focal <- c('DG','DGDG','LPC','LPE','MG','MGDG','PA','PC','PE','PG','PS','SQDG','TG')

## ---- our own GRM ----------------------------------------------------------
load_k <- function() {
  fk <- file.path(kdir, 'sap_grm.rel'); fi <- file.path(kdir, 'sap_grm.rel.id')
  stopifnot(file.exists(fk), file.exists(fi))
  K   <- as.matrix(fread(fk, header = FALSE, data.table = FALSE))
  idf <- fread(fi, header = FALSE, data.table = FALSE)
  idf <- idf[!startsWith(trimws(as.character(idf[[1]])), '#'), , drop = FALSE]  # drop '#FID IID'
  ids <- as.character(idf[[ncol(idf)]])          # IID is the last column
  if (nrow(K) != length(ids))
    stop('GRM is ', nrow(K), ' rows but .rel.id has ', length(ids), ' lines')
  K <- (K + t(K)) / 2
  dimnames(K) <- list(ids, ids)
  message(sprintf('GRM: %d genotypes | diag mean %.3f | offdiag mean %.4f',
                  nrow(K), mean(diag(K)), mean(K[upper.tri(K)])))
  list(K = K, ids = ids)
}

fit_h2 <- function(y, ids, kin) {
  i  <- match(ids, kin$ids)
  ok <- !is.na(i) & is.finite(y)
  y  <- y[ok]; i <- i[ok]
  if (length(y) < 30 || sd(y) == 0)
    return(tibble(n = length(y), Vu = NA_real_, Ve = NA_real_, h2 = NA_real_,
                  status = 'insufficient_variation'))
  q <- tryCatch(mixed.solve(y, K = kin$K[i, i, drop = FALSE]), error = function(e) e)
  if (inherits(q, 'error'))
    return(tibble(n = length(y), Vu = NA_real_, Ve = NA_real_, h2 = NA_real_,
                  status = conditionMessage(q)))
  tibble(n = length(y), Vu = q$Vu, Ve = q$Ve, h2 = q$Vu / (q$Vu + q$Ve), status = 'ok')
}

## ---- SpATS-adjusted phenotypes -------------------------------------------
read_spats <- function(cond) {
  f <- file.path(root, 'data','SPATS_fitted','non_normalized_intensities',
        sprintf('Final_subset_%s_all_lipids_fitted_phenotype_non_normalized.csv',
                if (cond == 'CTL') 'control' else 'lowinput'))
  x   <- fread(f, data.table = FALSE, check.names = FALSE)
  lip <- setdiff(names(x), c('LineRaw','PlotID','row','col'))
  lip <- lip[grepl('\\(', lip)]
  list(line = as.character(x$LineRaw), mat = as.matrix(x[, lip, drop = FALSE]), lipids = lip)
}

estimate <- function(cond) {
  d <- read_spats(cond)
  cls <- str_match(d$lipids, '^([A-Za-z]+)\\(')[, 2]
  mat <- d$mat; mat[!is.finite(mat)] <- NA_real_
  message(sprintf('%s: %d lines, %d species, %d matched to GRM',
                  cond, length(d$line), length(d$lipids), sum(d$line %in% KIN$ids)))
  bind_rows(lapply(focal, function(k) {
    idx <- which(!is.na(cls) & cls == k)
    if (!length(idx)) return(NULL)
    # no na.rm: a missing species makes the line's class sum NA rather than
    # silently smaller, which would push technical dropout into Ve
    y <- rowSums(mat[, idx, drop = FALSE])
    fit_h2(log10(y + 1), d$line, KIN) |>
      mutate(trial = cond, Class = k, n_species = length(idx), .before = 1)
  }))
}

KIN <- load_k()
res <- bind_rows(estimate('CTL'), estimate('LIN')) |>
  mutate(h2 = ifelse(h2 < 1e-6, 0, h2))

## ---- compare with GEMMA PVE ----------------------------------------------
pf <- file.path(kdir, 'pve_all_traits.tsv')
if (file.exists(pf)) {
  pve <- read_tsv(pf, show_col_types = FALSE) |>
    mutate(Class = str_replace(trait, '^Sum_', '') |> str_replace('_.*$', '')) |>
    filter(Class %in% focal, trial %in% c('CTL','LIN')) |>
    group_by(trial, Class) |>
    summarise(pve_gemma = mean(pve, na.rm = TRUE),
              se_gemma  = mean(se_pve, na.rm = TRUE),
              n_gemma_traits = n(), .groups = 'drop')
  res <- left_join(res, pve, by = c('trial','Class')) |>
    mutate(sig = ifelse(!is.na(pve_gemma) & pve_gemma > 2 * se_gemma, 'yes', 'no'))
} else {
  message('note: ', pf, ' not found - skipping GEMMA comparison')
}

write_csv(res, file.path(tabdir, 'SuppTable_Genomic_Heritability_FocalClasses.csv'))

cat('\n=============== genomic heritability, 13 focal classes ===============\n')
print(as.data.frame(res |> arrange(trial, desc(h2))), row.names = FALSE, digits = 3)
cat('\nmedian h2   CTL:', round(median(res$h2[res$trial=='CTL'], na.rm=TRUE), 4),
    '  LIN:',            round(median(res$h2[res$trial=='LIN'], na.rm=TRUE), 4), '\n')
cat('classes at zero   CTL:', sum(res$h2[res$trial=='CTL'] == 0, na.rm=TRUE),
    ' LIN:',                  sum(res$h2[res$trial=='LIN'] == 0, na.rm=TRUE), '\n')

p <- ggplot(res |> filter(is.finite(h2)),
            aes(h2, reorder(Class, h2), fill = trial)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7,
           colour = 'black', linewidth = 0.2) +
  scale_fill_manual(values = condition_colors) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = expression('Genomic heritability ('*h^2*')'), y = NULL) +
  plot_theme + theme(legend.position = 'bottom')
ggsave(file.path(figdir, 'SuppFig_Heritability_FocalClasses.png'), p,
       width = 8, height = 6, dpi = 350, bg = 'white')
message('done')
