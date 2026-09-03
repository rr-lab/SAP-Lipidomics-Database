#!/usr/bin/env Rscript
# Genomic heritability of lipid class sums computed on the SpATS-adjusted
# phenotypes actually used for GWAS, compared against the SERRF-based estimates
# from 22_SuppTable16_genomic_heritability.R.
#
# Same model as script 22: rrBLUP::mixed.solve with the LD-thinned GRM,
# h2 = Vu / (Vu + Ve). Only the input phenotype changes.

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(readr)
  library(stringr); library(rrBLUP); library(ggplot2); library(jsonlite); library(reticulate)
})

out    <- '/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database'
tabdir <- file.path(out, 'table', 'new_table')
figdir <- file.path(out, 'fig', 'new_figures')
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

plot_theme <- theme_minimal(base_size = 13) +
  theme(axis.text = element_text(colour = 'black', size = 9),
        axis.title = element_text(face = 'bold', size = 11),
        axis.line = element_line(colour = 'black'),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank())

# 13 focal classes reported in the paper (Supplementary Table S5D)
focal <- c('DG','DGDG','LPC','LPE','MG','MGDG','PA','PC','PE','PG','PS','SQDG','TG')

## ---- GRM (identical to script 22) -----------------------------------------
load_grm <- function(gd = file.path(out, 'data', 'genomic_grm_serrf')) {
  stopifnot(dir.exists(gd))
  ids <- read_tsv(file.path(gd, 'samples.tsv'), show_col_types = FALSE)$sample_id
  f_tsv <- file.path(gd, 'K_weighted.tsv'); f_npy <- file.path(gd, 'K_weighted.npy')
  if (file.exists(f_tsv)) {
    K <- as.matrix(fread(f_tsv, header = FALSE, data.table = FALSE))
  } else if (file.exists(f_npy)) {
    K <- reticulate::import('numpy', convert = TRUE)$load(f_npy)   # builder writes .npy
    storage.mode(K) <- 'numeric'
  } else stop('no K_weighted.tsv or K_weighted.npy in ', gd)
  # Confirm the GRM is unweighted: a p-value-weighted GRM makes h2 circular.
  sj <- file.path(gd, 'grm_summary.json')
  nm <- NA_integer_
  if (file.exists(sj)) {
    j <- fromJSON(sj)
    nm <- tryCatch(as.integer(j$counts$snps_kept), error = function(e) NA_integer_)
    cfg <- tryCatch(j$config, error = function(e) NULL)
    if (!is.null(cfg)) {
      a <- suppressWarnings(as.numeric(cfg$alpha))
      gw <- cfg$gwas_file
      if ((!is.na(a) && a > 0) && !is.null(gw) && !is.na(gw) && nzchar(as.character(gw)))
        warning('GRM appears p-value WEIGHTED (alpha=', a, ', gwas_file=', gw,
                '). Heritability from this GRM is circular. Rebuild with no --gwas-file.',
                immediate. = TRUE)
    }
  }
  list(K = (K + t(K)) / 2, ids = ids, nmarkers = nm)
}

fit_h2 <- function(y, ids, grm) {
  i  <- match(ids, grm$ids)
  ok <- !is.na(i) & is.finite(y)
  y  <- y[ok]; i <- i[ok]
  if (length(y) < 30 || sd(y) == 0)
    return(tibble(n_genotypes = length(y), Vu = NA_real_, Ve = NA_real_,
                  h2 = NA_real_, model_status = 'insufficient_variation'))
  q <- tryCatch(mixed.solve(y, K = grm$K[i, i, drop = FALSE]), error = function(e) e)
  if (inherits(q, 'error'))
    return(tibble(n_genotypes = length(y), Vu = NA_real_, Ve = NA_real_,
                  h2 = NA_real_, model_status = conditionMessage(q)))
  tibble(n_genotypes = length(y), Vu = q$Vu, Ve = q$Ve,
         h2 = q$Vu / (q$Vu + q$Ve), model_status = 'ok')
}

## ---- SpATS-adjusted phenotypes --------------------------------------------
# Layout: LineRaw, PlotID, row, col, <lipid columns>
read_spats <- function(cond) {
  f <- file.path(out, 'data','SPATS_fitted','non_normalized_intensities',
                 sprintf('Final_subset_%s_all_lipids_fitted_phenotype_non_normalized.csv',
                         if (cond == 'CTL') 'control' else 'lowinput'))
  x <- fread(f, data.table = FALSE, check.names = FALSE)
  lip <- setdiff(names(x), c('LineRaw', 'PlotID', 'row', 'col'))
  lip <- lip[grepl('\\(', lip)]           # species-level names only
  list(line = as.character(x$LineRaw), mat = as.matrix(x[, lip, drop = FALSE]), lipids = lip)
}

class_of <- function(lipid) str_match(lipid, '^([A-Za-z]+)\\(')[, 2]

estimate <- function(cond, grm) {
  d   <- read_spats(cond)
  cls <- class_of(d$lipids)
  keep <- !is.na(cls) & cls %in% focal
  message(sprintf('%s: %d samples, %d species, %d in focal classes',
                  cond, length(d$line), length(d$lipids), sum(keep)))
  mat <- d$mat; mat[!is.finite(mat)] <- NA_real_
  res <- lapply(focal, function(k) {
    idx <- which(keep & cls == k)
    if (length(idx) == 0) return(NULL)
    # NA-aware class sum: a species missing for a line makes that line's sum NA
    # rather than silently smaller (na.rm = TRUE would push technical dropout into Ve)
    y <- rowSums(mat[, idx, drop = FALSE])
    fit_h2(log10(y + 1), d$line, grm) |>
      mutate(Condition = cond, Class = k, n_species = length(idx), .before = 1)
  })
  bind_rows(res)
}

grm <- load_grm()
message('GRM: ', length(grm$ids), ' genotypes, ', grm$nmarkers, ' markers')

spats <- bind_rows(estimate('CTL', grm), estimate('LIN', grm)) |>
  mutate(Input = 'SpATS-adjusted')

## ---- merge with the existing SERRF estimates ------------------------------
serrf_file <- file.path(tabdir, 'SuppTable_Genomic_Heritability_Lipid_Class_Sums.csv')
cmp <- spats |> select(Condition, Class, h2_spats = h2, n_genotypes)
if (file.exists(serrf_file)) {
  serrf <- read_csv(serrf_file, show_col_types = FALSE) |>
    filter(Class %in% focal) |> select(Condition, Class, h2_serrf = h2)
  cmp <- full_join(cmp, serrf, by = c('Condition', 'Class'))
} else {
  cmp$h2_serrf <- NA_real_
}
cmp <- cmp |>
  mutate(across(c(h2_spats, h2_serrf), ~ ifelse(.x < 1e-6, 0, .x)),   # boundary -> 0
         delta = h2_spats - h2_serrf) |>
  arrange(Condition, desc(h2_spats))

write_csv(spats, file.path(tabdir, 'SuppTable_Genomic_Heritability_ClassSums_SpATS.csv'))
write_csv(cmp,   file.path(tabdir, 'SuppTable_Genomic_Heritability_SpATS_vs_SERRF.csv'))

cat('\n================ h2, SpATS-adjusted vs SERRF ================\n')
print(as.data.frame(cmp), row.names = FALSE, digits = 3)
cat('\nmedian h2  SpATS  CTL:', round(median(cmp$h2_spats[cmp$Condition=='CTL'], na.rm=TRUE), 4),
    ' LIN:',                   round(median(cmp$h2_spats[cmp$Condition=='LIN'], na.rm=TRUE), 4), '\n')
cat('median h2  SERRF  CTL:',  round(median(cmp$h2_serrf[cmp$Condition=='CTL'], na.rm=TRUE), 4),
    ' LIN:',                   round(median(cmp$h2_serrf[cmp$Condition=='LIN'], na.rm=TRUE), 4), '\n')
cat('classes with h2 == 0 (SpATS)  CTL:', sum(cmp$h2_spats[cmp$Condition=='CTL'] == 0, na.rm=TRUE),
    'of', sum(cmp$Condition=='CTL'),
    ' LIN:', sum(cmp$h2_spats[cmp$Condition=='LIN'] == 0, na.rm=TRUE),
    'of', sum(cmp$Condition=='LIN'), '\n')

p <- cmp |>
  pivot_longer(c(h2_spats, h2_serrf), names_to = 'Input', values_to = 'h2') |>
  mutate(Input = recode(Input, h2_spats = 'SpATS-adjusted', h2_serrf = 'SERRF only')) |>
  filter(is.finite(h2)) |>
  ggplot(aes(h2, reorder(Class, h2), fill = Input)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7,
           colour = 'black', linewidth = 0.2) +
  facet_wrap(~ Condition) +
  scale_fill_manual(values = c('SpATS-adjusted' = '#440154FF', 'SERRF only' = '#FDE725FF')) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = expression('Genomic heritability ('*h^2*')'), y = NULL) +
  plot_theme + theme(legend.position = 'bottom')

ggsave(file.path(figdir, 'SuppFig_Heritability_SpATS_vs_SERRF.png'), p,
       width = 10, height = 6, dpi = 350, bg = 'white')
message('done')
