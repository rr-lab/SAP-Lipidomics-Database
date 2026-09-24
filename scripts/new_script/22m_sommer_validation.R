#!/usr/bin/env Rscript
# ==============================================================================
# sommer validation of the SoLD class-sum genomic heritability
#
# Run on the machine that has sommer installed.
#   setwd('/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database')
#   source('scripts/new_script/22m_sommer_validation.R')
#
# Same model fitted three ways, printed side by side.
#   ours    profile likelihood over h2, the pipeline estimator
#   sommer  mmer(y ~ 1, random = ~vsr(Genotype, Gu = K), rcov = ~units),
#           h2 taken as V_G / (V_G + V_E)
#   lme4    the same model with K replaced by its symmetric square root, which
#           is the construction sommer applies internally
#
# Phenotype is log10(summed SpATS fitted intensity of the class + 1), on every
# accession of the trial that is in the GRM and every species of that class
# detected in that trial. Same construction as the per-species estimates.
#
# Writes table/new_table/h2_sommer_validation.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(readr); library(stringr)
})
have_sommer <- requireNamespace('sommer', quietly = TRUE)
have_lme4   <- requireNamespace('lme4',   quietly = TRUE)
if (have_sommer) message('sommer ', as.character(packageVersion('sommer')))
if (!have_sommer) message('sommer not found, that column will be skipped')

ROOT <- Sys.getenv('SOLD_REPO', getwd())
kdir <- file.path(ROOT, 'data', 'kinship')
sdir <- file.path(ROOT, 'data', 'SPATS_fitted', 'non_normalized_intensities')
odir <- file.path(ROOT, 'table', 'new_table')
dir.create(odir, recursive = TRUE, showWarnings = FALSE)

CLASSES <- c('DGDG','MGDG','SQDG','GalCer','LPC','LPE','Cer','AEG','FA',
             'TG','DG','MG','PC','PE','PG','PA','PS')
GRID <- seq(0, 0.999, by = 0.001)

# ---- inputs ------------------------------------------------------------------
K   <- as.matrix(fread(file.path(kdir,'sap_grm.rel'), header=FALSE, data.table=FALSE))
idf <- fread(file.path(kdir,'sap_grm.rel.id'), header=FALSE, data.table=FALSE, sep='\t')
idf <- idf[!startsWith(trimws(as.character(idf[[1]])),'#'), , drop=FALSE]
ids <- as.character(idf[[ncol(idf)]])
stopifnot(nrow(K) == length(ids))
K <- (K + t(K))/2; dimnames(K) <- list(ids, ids)
message(sprintf('GRM %d genotypes, diag mean %.3f', nrow(K), mean(diag(K))))

read_spats <- function(cond) {
  f <- file.path(sdir, sprintf(
        'Final_subset_%s_all_lipids_fitted_phenotype_non_normalized.csv',
        if (cond == 'CTL') 'control' else 'lowinput'))
  x   <- fread(f, data.table = FALSE, check.names = FALSE)
  lip <- setdiff(names(x), c('LineRaw','PlotID','row','col'))
  lip <- lip[grepl('\\(', lip)]
  m <- as.matrix(x[, lip, drop = FALSE]); m[!is.finite(m)] <- NA_real_
  list(line = as.character(x$LineRaw), mat = m, lipids = lip)
}
cls <- function(v) str_match(v, '^([A-Za-z]+)\\(')[, 2]

# ---- the three estimators ----------------------------------------------------
prof_h2 <- function(y, Ksub) {
  n <- length(y); p <- 1
  e <- eigen(Ksub, symmetric = TRUE); d <- pmax(e$values, 0)
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
  c(h2 = GRID[i], lo = GRID[min(keep)], hi = GRID[max(keep)])
}

sommer_h2 <- function(y, gid, Ksub) {
  if (!have_sommer) return(c(h2 = NA, vg = NA, ve = NA))
  d <- data.frame(y = y, Genotype = factor(gid, levels = rownames(Ksub)))
  A <- Ksub
  fit <- NULL
  # sommer 4.4 renamed mmer to mmes and vsr to vsm(ism(.)), so try both.
  if ('mmer' %in% getNamespaceExports('sommer')) {
    fit <- try(sommer::mmer(y ~ 1,
                 random = stats::as.formula('~ sommer::vsr(Genotype, Gu = A)'),
                 rcov = ~ units, data = d, verbose = FALSE), silent = TRUE)
  }
  if (is.null(fit) || inherits(fit, 'try-error')) {
    fit <- try(sommer::mmes(y ~ 1,
                 random = stats::as.formula('~ sommer::vsm(sommer::ism(Genotype), Gu = A)'),
                 rcov = ~ units, data = d, verbose = FALSE), silent = TRUE)
  }
  if (is.null(fit) || inherits(fit, 'try-error')) return(c(h2 = NA, vg = NA, ve = NA))
  v <- as.numeric(unlist(if (!is.null(fit$sigma)) fit$sigma else fit$theta))
  if (length(v) < 2) return(c(h2 = NA, vg = NA, ve = NA))
  vg <- v[1]; ve <- v[length(v)]
  c(h2 = vg/(vg+ve), vg = vg, ve = ve)
}

lme4_h2 <- function(y, Ksub) {
  if (!have_lme4) return(c(h2 = NA))
  n  <- nrow(Ksub)
  ev <- eigen(Ksub, symmetric = TRUE)
  L  <- ev$vectors %*% (sqrt(pmax(ev$values, 0)) * t(ev$vectors))
  ct <- lme4::lmerControl(check.nobs.vs.nlev='ignore', check.nobs.vs.nRE='ignore',
                          check.nlev.gtreq.5='ignore')
  d  <- data.frame(y = y, id = factor(sprintf('g%04d', seq_len(n))))
  lf <- lme4::lFormula(y ~ 1 + (1|id), data = d, REML = TRUE, control = ct)
  lf$reTrms$Zt <- Matrix::Matrix(crossprod(L, as.matrix(lf$reTrms$Zt)), sparse = TRUE)
  devf <- do.call(lme4::mkLmerDevfun, lf)
  opt  <- lme4::optimizeLmer(devf, control = list(maxfun = 1e5))
  fit  <- lme4::mkMerMod(environment(devf), opt, lf$reTrms, fr = lf$fr)
  vc   <- as.data.frame(lme4::VarCorr(fit))
  sg2  <- vc$vcov[vc$grp == 'id']; se2 <- vc$vcov[vc$grp == 'Residual']
  c(h2 = sg2/(sg2+se2))
}

int <- function(v) qnorm((rank(v) - 0.5)/length(v))

# ---- run ---------------------------------------------------------------------
res <- bind_rows(lapply(c('CTL','LIN'), function(cond) {
  d <- read_spats(cond); k <- cls(d$lipids)
  keep <- d$line %in% ids
  gid  <- d$line[keep]
  Ks   <- K[gid, gid, drop = FALSE]
  bind_rows(lapply(CLASSES, function(cl) {
    idx <- which(!is.na(k) & k == cl); if (!length(idx)) return(NULL)
    yraw <- rowSums(d$mat[keep, idx, drop = FALSE])
    y    <- log10(yraw + 1)
    a  <- prof_h2(y, Ks)
    sm <- sommer_h2(y, gid, Ks)
    l4 <- lme4_h2(y, Ks)
    message(sprintf('  %s %-6s  ours %.3f  sommer %.3f  lme4 %.3f',
                    cond, cl, a[['h2']], sm[['h2']], l4[['h2']]))
    tibble(Condition = cond, Class = cl, N_species = length(idx), n = sum(keep),
           h2_ours = a[['h2']], CI_lo = a[['lo']], CI_hi = a[['hi']],
           h2_sommer = sm[['h2']], Vg_sommer = sm[['vg']], Ve_sommer = sm[['ve']],
           h2_lme4 = l4[['h2']],
           h2_untransformed = prof_h2(yraw, Ks)[['h2']],
           h2_inverse_normal = prof_h2(int(yraw), Ks)[['h2']])
  }))
}))

res <- res %>% mutate(d_sommer = h2_sommer - h2_ours, d_lme4 = h2_lme4 - h2_ours)
write_csv(res, file.path(odir, 'h2_sommer_validation.csv'))

cat('\n===== class-sum h2, log10(summed fitted intensity + 1) =====\n')
print(as.data.frame(res %>% select(Condition, Class, N_species, n, h2_ours, CI_lo, CI_hi,
                                   h2_sommer, h2_lme4, d_sommer, d_lme4)),
      row.names = FALSE, digits = 4)
ok <- is.finite(res$d_sommer)
if (any(ok)) {
  cat(sprintf('\n  max |ours - sommer|   %.5f   over %d traits\n',
              max(abs(res$d_sommer[ok])), sum(ok)))
  cat(sprintf('  correlation           %.8f\n', cor(res$h2_ours[ok], res$h2_sommer[ok])))
}
ok4 <- is.finite(res$d_lme4)
if (any(ok4)) cat(sprintf('  max |ours - lme4|     %.5f\n', max(abs(res$d_lme4[ok4]))))

cat('\n===== scale sensitivity =====\n')
print(as.data.frame(res %>% select(Condition, Class, h2_ours, h2_untransformed,
                                   h2_inverse_normal)), row.names = FALSE, digits = 3)
cat('\nwritten ', file.path(odir, 'h2_sommer_validation.csv'), '\n')
