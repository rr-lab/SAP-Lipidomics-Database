# ==============================================================================
# 44_heritability_leave_one_race_out.R
#
#   Rscript scripts/new_script/44_heritability_leave_one_race_out.R
#
# Robustness of the heritability contrast to ancestry composition.
#
# 05_SuppTableS5G_composition_stability.R already asks this of lipid-class
# composition: drop one botanical race group, or 20 per cent of the panel, and
# see whether class shares move. The heritability result needs the same
# treatment for a different reason. Profile-likelihood intervals express
# sampling error at a fixed panel, but they do not say whether the count of
# heritable species depends on one ancestry group being present.
#
# Two checks, per trial:
#   * leave-one-race-out: drop all accessions of one botanical race group and
#     recount how many species have a 95% interval excluding zero;
#   * 80 per cent subsamples: draw random subsets of genotypes and do the same,
#     which also shows how much of any drop is loss of sample size alone.
#
# Only the intercept-only fit is run here. The structure-conditioned fit is in
# script 42 and answers a different question.
#
# Input   data/SPATS_fitted/non_normalized_intensities/Final_subset_{control,lowinput}_*.csv
#         data/kinship/sap_grm.rel, sap_grm.rel.id ; data/SAP_geoloc.csv
# Output  table/new_table/SuppTable_Heritability_AncestryRobustness.csv
#
# Run from the repository root.
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(readr); library(stringr)
})

root   <- Sys.getenv('SOLD_DB_REPO', '.')
kdir   <- file.path(root,'data','kinship'); tabdir <- file.path(root,'table','new_table')
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)
GRID <- seq(0, 1, by = 0.004)          # coarser than 22e; we count, not report, h2
N_SUB <- 20; SUB_FRAC <- 0.8
PURE <- c('Bicolor','Caudatum','Durra','Guinea','Kafir')
set.seed(1)

load_k <- function() {
  K <- as.matrix(fread(file.path(kdir,'sap_grm.rel'), header=FALSE, data.table=FALSE))
  idf <- fread(file.path(kdir,'sap_grm.rel.id'), header=FALSE, data.table=FALSE, sep='\t')
  idf <- idf[!grepl('^#', idf[[1]]), , drop=FALSE]
  ids <- as.character(idf[[2]]); dimnames(K) <- list(ids,ids); list(K=K, ids=ids)
}
read_spats <- function(cond) {
  f <- file.path(root,'data','SPATS_fitted','non_normalized_intensities',
        sprintf('Final_subset_%s_all_lipids_fitted_phenotype_non_normalized.csv',
                if (cond=='CTL') 'control' else 'lowinput'))
  x <- fread(f, data.table=FALSE, check.names=FALSE)
  lip <- setdiff(names(x), c('LineRaw','PlotID','row','col')); lip <- lip[grepl('\\(', lip)]
  m <- as.matrix(x[, lip, drop=FALSE]); m[!is.finite(m)] <- NA_real_
  list(line=as.character(x$LineRaw), mat=m, lipids=lip)
}
# returns TRUE when the 95% profile interval excludes zero
nonzero_h2 <- function(yt, xt, d) {
  n <- length(yt)
  ll <- vapply(GRID, function(h){ v <- h*d+(1-h); if(any(v<=1e-12)) return(-Inf)
    w <- 1/v; XtWX <- sum(w*xt*xt); if(!is.finite(XtWX)||XtWX<=0) return(-Inf)
    b <- sum(w*xt*yt)/XtWX; r <- yt-xt*b; s2 <- sum(w*r*r)/(n-1)
    if(!is.finite(s2)||s2<=0) return(-Inf)
    -0.5*(sum(log(v)) + (n-1)*log(s2) + log(XtWX)) }, numeric(1))
  i <- which.max(ll); keep <- which(ll >= ll[i]-qchisq(0.95,1)/2)
  GRID[min(keep)] > 0
}
count_nonzero <- function(M, Ksub) {
  e <- eigen(Ksub, symmetric=TRUE); d <- pmax(e$values,0); V <- e$vectors
  xt <- as.numeric(crossprod(V, rep(1, nrow(M))))
  sum(vapply(seq_len(ncol(M)), function(j){
    y <- M[,j]; if(!all(is.finite(y))) return(FALSE)
    nonzero_h2(as.numeric(crossprod(V, log10(y+1))), xt, d) }, logical(1)))
}

geo <- fread(file.path(root,'data','SAP_geoloc.csv')); setnames(geo,1,'Taxa')
geo <- unique(geo[, .(Taxa=as.character(Taxa), Original_Race)], by='Taxa')
geo[, RaceGroup := fifelse(is.na(Original_Race)|Original_Race %in% c('','NA')|
      grepl('verticilliflorum',Original_Race,ignore.case=TRUE), NA_character_,
      fifelse(Original_Race %in% PURE, Original_Race, 'Mixed'))]

KIN <- load_k(); out <- list()
for (cond in c('CTL','LIN')) {
  d <- read_spats(cond)
  keep <- d$line %in% KIN$ids
  line <- d$line[keep]; M <- d$mat[keep, , drop=FALSE]
  i <- match(line, KIN$ids); K <- KIN$K[i,i,drop=FALSE]
  rg <- geo$RaceGroup[match(line, geo$Taxa)]
  full <- count_nonzero(M, K)
  out[[length(out)+1]] <- data.table(trial=cond, analysis='full panel', dropped=NA_character_,
                                     n_genotypes=nrow(M), n_species=ncol(M), n_h2_nonzero=full)
  message(sprintf('%s full panel: %d of %d species', cond, full, ncol(M)))
  for (g in sort(unique(na.omit(rg)))) {
    s <- which(rg != g | is.na(rg))
    out[[length(out)+1]] <- data.table(trial=cond, analysis='leave-one-race-out', dropped=g,
      n_genotypes=length(s), n_species=ncol(M), n_h2_nonzero=count_nonzero(M[s,,drop=FALSE], K[s,s,drop=FALSE]))
  }
  for (b in seq_len(N_SUB)) {
    s <- sort(sample(nrow(M), floor(SUB_FRAC*nrow(M))))
    out[[length(out)+1]] <- data.table(trial=cond, analysis='80% subsample', dropped=sprintf('rep%02d',b),
      n_genotypes=length(s), n_species=ncol(M), n_h2_nonzero=count_nonzero(M[s,,drop=FALSE], K[s,s,drop=FALSE]))
  }
}
res <- rbindlist(out)
write_csv(res, file.path(tabdir,'SuppTable_Heritability_AncestryRobustness.csv'))
cat('\n== species with non-zero heritability, by panel ==\n')
print(as.data.frame(res[analysis!='80% subsample']), row.names=FALSE)
cat('\n== 80% subsamples ==\n')
print(as.data.frame(res[analysis=='80% subsample',
  .(reps=.N, median=median(n_h2_nonzero), min=min(n_h2_nonzero), max=max(n_h2_nonzero)), by=trial]), row.names=FALSE)
