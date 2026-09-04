# ==============================================================================
# 42_heritability_structure_conditioned.R
#
#   Rscript scripts/new_script/42_heritability_structure_conditioned.R
#
# Genomic heritability of each lipid-class sum is estimated twice: once with an
# intercept only, exactly as in 22e, and once with the leading eigenvectors of
# the genomic relationship matrix as fixed covariates.
#
# Why. Heritability from a GRM and variance explained by genetic cluster are not
# independent quantities. A GRM encodes relatedness, and population structure is
# relatedness, so a trait that differs between clusters will show a raised h2
# whether or not any variance segregates within clusters. Conditioning on the
# leading eigenvectors removes the between-cluster component and leaves the part
# that does. A trait whose interval still excludes zero after conditioning is
# carrying additive variance beyond structure.
#
# Everything about the plain fit matches 22e: the same shared-species set, the
# same log10(sum + 1) phenotype, the same profile-likelihood interval. Only the
# fixed-effect design changes between the two fits, so the pair is comparable.
#
# Input   data/SPATS_fitted/non_normalized_intensities/Final_subset_{control,lowinput}_*.csv
#         data/kinship/sap_grm.rel, sap_grm.rel.id
# Output  table/new_table/SuppTable_Heritability_StructureConditioned.csv
#
# Run from the repository root.
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(readr); library(stringr)
})

root   <- Sys.getenv('SOLD_DB_REPO', '.')
kdir   <- file.path(root, 'data', 'kinship')
tabdir <- file.path(root, 'table', 'new_table')
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)

N_PC  <- 3                      # matches the GWAS, which fitted three genotype PCs
GRID  <- seq(0, 1, by = 0.002)
focal <- c('PC','PE','PG','PA','PS','LPC','LPE','DG','MG','TG','MGDG','DGDG','SQDG')

load_k <- function() {
  K   <- as.matrix(fread(file.path(kdir, 'sap_grm.rel'), header = FALSE, data.table = FALSE))
  idf <- fread(file.path(kdir, 'sap_grm.rel.id'), header = FALSE, data.table = FALSE, sep = '\t')
  idf <- idf[!grepl('^#', idf[[1]]), , drop = FALSE]      # plink2 writes a #FID IID header
  ids <- as.character(idf[[2]])
  stopifnot(nrow(K) == length(ids))
  dimnames(K) <- list(ids, ids)
  list(K = K, ids = ids)
}

# profile-likelihood h2 for an arbitrary fixed-effect design X
prof_h2 <- function(y, K, X) {
  n <- length(y); p <- ncol(X)
  e <- eigen(K, symmetric = TRUE); d <- pmax(e$values, 0)
  yt <- crossprod(e$vectors, y); Xt <- crossprod(e$vectors, X)
  ll <- vapply(GRID, function(h) {
    v <- h * d + (1 - h); if (any(v <= 1e-12)) return(-Inf)
    w <- 1 / v
    XtWX <- crossprod(Xt * w, Xt)
    b <- tryCatch(solve(XtWX, crossprod(Xt * w, yt)), error = function(e) NULL)
    if (is.null(b)) return(-Inf)
    r <- yt - Xt %*% b; s2 <- sum(w * r * r) / (n - p)
    if (!is.finite(s2) || s2 <= 0) return(-Inf)
    -0.5 * (sum(log(v)) + (n - p) * log(s2) +
            as.numeric(determinant(XtWX, logarithm = TRUE)$modulus))
  }, numeric(1))
  i <- which.max(ll); keep <- which(ll >= ll[i] - qchisq(0.95, 1) / 2)
  c(h2 = GRID[i], lo = GRID[min(keep)], hi = GRID[max(keep)])
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
cls <- function(v) str_match(v, '^([A-Za-z]+)\\(')[, 2]

estimate <- function(d, cond) {
  k <- cls(d$lipids); use <- d$lipids %in% shared
  bind_rows(lapply(focal, function(cl) {
    idx <- which(use & !is.na(k) & k == cl)
    if (!length(idx)) return(NULL)
    y <- log10(rowSums(d$mat[, idx, drop = FALSE]) + 1)     # as in 22e
    i <- match(d$line, KIN$ids); ok <- !is.na(i) & is.finite(y)
    if (sum(ok) < 30) return(NULL)
    Ksub <- KIN$K[i[ok], i[ok], drop = FALSE]; yy <- y[ok]
    ev   <- eigen(Ksub, symmetric = TRUE)$vectors[, seq_len(N_PC), drop = FALSE]
    a <- prof_h2(yy, Ksub, matrix(1, length(yy), 1))
    b <- prof_h2(yy, Ksub, cbind(1, ev))
    tibble(trial = cond, Class = cl, n_species = length(idx), n = sum(ok),
           h2 = a[['h2']],    lo = a[['lo']],    hi = a[['hi']],
           h2_pc = b[['h2']], lo_pc = b[['lo']], hi_pc = b[['hi']])
  }))
}

res <- bind_rows(estimate(CTLd, 'CTL'), estimate(LINd, 'LIN')) |>
  mutate(nonzero = lo > 0, nonzero_pc = lo_pc > 0,
         h2_lost_to_structure = round(h2 - h2_pc, 3),
         single_species = n_species == 1)

write_csv(res, file.path(tabdir, 'SuppTable_Heritability_StructureConditioned.csv'))

cat('\nGenomic heritability of class sums, before and after conditioning on',
    N_PC, 'GRM eigenvectors\n\n')
print(as.data.frame(res |>
  arrange(trial, desc(h2)) |>
  transmute(trial, Class, n_sp = n_species,
            h2 = round(h2, 3), CI = sprintf('[%.2f,%.2f]', lo, hi), sig = nonzero,
            h2_PC = round(h2_pc, 3), CI_PC = sprintf('[%.2f,%.2f]', lo_pc, hi_pc),
            sig_PC = nonzero_pc, single_species)), row.names = FALSE)
cat('\nTraits still non-zero after conditioning:\n')
print(as.data.frame(res |> filter(nonzero_pc) |>
  transmute(trial, Class, n_species, h2 = round(h2,3), h2_pc = round(h2_pc,3))), row.names = FALSE)
