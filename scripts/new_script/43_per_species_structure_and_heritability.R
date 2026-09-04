# ==============================================================================
# 43_per_species_structure_and_heritability.R
#
#   Rscript scripts/new_script/43_per_species_structure_and_heritability.R
#
# The class-level analyses hide a problem: PS, PA and LPE each contain a single
# molecular species, so their "class sum" is that species and a class-level
# statement about them is really a statement about one measurement. This script
# drops to the species level and runs both analyses there, so every claim is
# made on a stated number of species rather than on a class label.
#
# The structure test runs on every species detected in a trial, because that
# question does not require the trials to be paired. The heritability columns
# are reported for the same species, and an in_both flag marks the 146 species
# shared by the two trials, which are the ones the paired CTL/LIN comparison in
# 22f is built on.
#
# For each lipid species, separately within each trial:
#
#   * variance explained by botanical race and by marker-based genetic cluster
#     (Kruskal-Wallis epsilon^2, BH-corrected within trial x grouping), the
#     species-level counterpart of the class-level tests behind Figure 1;
#
#   * genomic heritability with a 95% profile-likelihood interval, fitted with
#     an intercept only and again with the leading eigenvectors of the genomic
#     relationship matrix as fixed effects, so that variance between genetic
#     clusters is separated from variance segregating within them.
#
# The per-class summary then counts, rather than asserts: how many species in a
# class carry non-zero heritability, how many survive conditioning on structure,
# and how many show a significant ancestry effect.
#
# Input   data/SPATS_fitted/non_normalized_intensities/Final_subset_{control,lowinput}_*.csv
#         data/kinship/sap_grm.rel, sap_grm.rel.id
#         data/SAP_geoloc.csv
# Output  table/new_table/SuppTable_PerSpecies_Structure_and_Heritability.csv
#         table/new_table/SuppTable_PerSpecies_Summary_byClass.csv
#
# Run from the repository root.
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(readr); library(stringr); library(tibble)
})

root   <- Sys.getenv('SOLD_DB_REPO', '.')
kdir   <- file.path(root, 'data', 'kinship')
tabdir <- file.path(root, 'table', 'new_table')
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)

N_PC <- 3
GRID <- seq(0, 1, by = 0.002)
PURE <- c('Bicolor','Caudatum','Durra','Guinea','Kafir')

load_k <- function() {
  K   <- as.matrix(fread(file.path(kdir,'sap_grm.rel'), header = FALSE, data.table = FALSE))
  idf <- fread(file.path(kdir,'sap_grm.rel.id'), header = FALSE, data.table = FALSE, sep = '\t')
  idf <- idf[!grepl('^#', idf[[1]]), , drop = FALSE]
  ids <- as.character(idf[[2]]); dimnames(K) <- list(ids, ids)
  list(K = K, ids = ids)
}

read_spats <- function(cond) {
  f <- file.path(root,'data','SPATS_fitted','non_normalized_intensities',
        sprintf('Final_subset_%s_all_lipids_fitted_phenotype_non_normalized.csv',
                if (cond == 'CTL') 'control' else 'lowinput'))
  x   <- fread(f, data.table = FALSE, check.names = FALSE)
  lip <- setdiff(names(x), c('LineRaw','PlotID','row','col')); lip <- lip[grepl('\\(', lip)]
  m <- as.matrix(x[, lip, drop = FALSE]); m[!is.finite(m)] <- NA_real_
  list(line = as.character(x$LineRaw), mat = m, lipids = lip)
}

# epsilon^2 from Kruskal-Wallis: (H - k + 1) / (n - k)
eps2_kw <- function(y, g) {
  ok <- !is.na(y) & !is.na(g); y <- y[ok]; g <- droplevels(factor(g[ok]))
  k <- nlevels(g); n <- length(y)
  if (k < 2 || n < 30) return(c(eps2 = NA_real_, p = NA_real_, k = k, n = n))
  kt <- kruskal.test(y ~ g)
  c(eps2 = max(0, (unname(kt$statistic) - k + 1) / (n - k)),
    p = kt$p.value, k = k, n = n)
}

# profile-likelihood h2 reusing one eigendecomposition for every species
prof_h2 <- function(yt, Xt, d, p) {
  n <- length(yt)
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
  c(GRID[i], GRID[min(keep)], GRID[max(keep)])
}

geo <- fread(file.path(root,'data','SAP_geoloc.csv'))
setnames(geo, 1, 'Taxa')
geo <- unique(geo[, .(Taxa = as.character(Taxa), K.Cluster, Original_Race)], by = 'Taxa')
geo[, RaceGroup := fifelse(is.na(Original_Race) | Original_Race %in% c('','NA') |
                           grepl('verticilliflorum', Original_Race, ignore.case = TRUE),
                           NA_character_,
                           fifelse(Original_Race %in% PURE, Original_Race, 'Mixed'))]

CTLd <- read_spats('CTL'); LINd <- read_spats('LIN'); KIN <- load_k()
shared <- intersect(CTLd$lipids, LINd$lipids)
message(sprintf('species: CTL %d, LIN %d, shared %d',
                length(CTLd$lipids), length(LINd$lipids), length(shared)))
cls <- function(v) str_match(v, '^([A-Za-z]+)\\(')[, 2]

run_trial <- function(d, cond) {
  keep <- d$line %in% KIN$ids
  line <- d$line[keep]; M <- d$mat[keep, d$lipids, drop = FALSE]
  i <- match(line, KIN$ids); Ksub <- KIN$K[i, i, drop = FALSE]
  e  <- eigen(Ksub, symmetric = TRUE); dd <- pmax(e$values, 0); V <- e$vectors
  X0 <- crossprod(V, matrix(1, nrow(M), 1))
  X1 <- crossprod(V, cbind(1, V[, seq_len(N_PC), drop = FALSE]))
  g  <- geo[match(line, Taxa)]
  message(sprintf('  %s: %d genotypes with a GRM entry', cond, nrow(M)))

  rbindlist(lapply(d$lipids, function(sp) {
    y <- M[, sp]; ok <- is.finite(y)
    if (sum(ok) < 30) return(NULL)
    yl <- log10(y + 1)
    r  <- eps2_kw(yl, g$RaceGroup)
    kc <- eps2_kw(yl, g$K.Cluster)
    yt <- crossprod(V, yl)
    a  <- prof_h2(yt, X0, dd, 1)
    b  <- prof_h2(yt, X1, dd, 1 + N_PC)
    data.table(trial = cond, Species = sp, Class = cls(sp),
               in_both = sp %in% shared, n = sum(ok),
               eps2_race = r[['eps2']], p_race = r[['p']],
               eps2_cluster = kc[['eps2']], p_cluster = kc[['p']],
               h2 = a[1], lo = a[2], hi = a[3],
               h2_pc = b[1], lo_pc = b[2], hi_pc = b[3])
  }))
}

res <- rbind(run_trial(CTLd, 'CTL'), run_trial(LINd, 'LIN'))
res[, `:=`(q_race    = p.adjust(p_race, 'BH'),
           q_cluster = p.adjust(p_cluster, 'BH')), by = trial]
res[, `:=`(nonzero = lo > 0, nonzero_pc = lo_pc > 0,
           sig_race = q_race < 0.05, sig_cluster = q_cluster < 0.05)]
write_csv(res, file.path(tabdir, 'SuppTable_PerSpecies_Structure_and_Heritability.csv'))

summ <- res[, .(n_species = .N,
                median_h2 = round(median(h2), 3),
                n_h2_nonzero = sum(nonzero),
                n_survives_structure = sum(nonzero_pc),
                n_sig_cluster = sum(sig_cluster, na.rm = TRUE),
                n_sig_race = sum(sig_race, na.rm = TRUE)),
            by = .(trial, Class)][order(trial, -n_h2_nonzero, -median_h2)]
write_csv(summ, file.path(tabdir, 'SuppTable_PerSpecies_Summary_byClass.csv'))

cat('\n== per-class counts, species level ==\n'); print(as.data.frame(summ), row.names = FALSE)
cat('\n== totals ==\n')
print(as.data.frame(res[, .(species = .N, h2_nonzero = sum(nonzero),
        survives_structure = sum(nonzero_pc),
        sig_cluster = sum(sig_cluster, na.rm = TRUE),
        sig_race = sum(sig_race, na.rm = TRUE)), by = trial]), row.names = FALSE)
