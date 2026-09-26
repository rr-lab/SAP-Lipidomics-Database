# ==============================================================================
# Supplementary Table S25 -- per-species paired h2, rebuilt on 164 species.
#
# WHY. 22q built the shared-species set with intersect() on RAW feature names,
# so CTL's PC(18:2/0:0) and LIN's LPC(18:2) were treated as different species and
# the intersection came to 163. They are the same lipid. normalize_lipid_name()
# resolves the two spellings, giving 164 -- the number the species inventory
# (S6a) reports as Common, and the number that makes 164 + 50 + 52 = 266 add up,
# and the number the LION input already used.
#
# Only LPC(18:2) is added. Nothing else moves, because the rename affects only
# the lyso spellings. Features whose names collapse onto one species are summed
# before the fit so nothing is counted twice.
#
# Everything else is 22q unchanged: log10(fitted + 1), the profile-likelihood h2
# with V = sigma^2 (h K + (1-h) I) over the genotypes phenotyped in both trials,
# and the CI as the set of h2 within chi2_1/2 = 1.92 of the maximum restricted
# log-likelihood.
# ==============================================================================
source('/home/claude/h2/check/_engine.R')
KIN <- load_k('sap_grm_september')

rd <- function(cond){
  f <- file.path(sdir, sprintf('Final_subset_%s_all_lipids_fitted_phenotype_non_normalized.csv',
        if (cond=='CTL') 'control' else 'lowinput'))
  x <- fread(f, data.table=FALSE, check.names=FALSE)
  lip <- setdiff(names(x), c('LineRaw','PlotID','row','col'))
  m <- as.matrix(x[, lip, drop=FALSE]); m[!is.finite(m)] <- NA_real_
  # collapse features that normalise to the same species
  nm <- normalize_lipid_name(lip)
  m <- t(rowsum(t(m), nm, na.rm = FALSE))
  list(line=as.character(x$LineRaw), mat=m, lipids=colnames(m))
}
C <- rd('CTL'); L <- rd('LIN')
shared_sp <- sort(intersect(C$lipids, L$lipids))
shared_ge <- intersect(intersect(C$line, L$line), KIN$ids)
message(sprintf('shared species %d | shared genotypes %d', length(shared_sp), length(shared_ge)))

one <- function(d, cond){
  keep <- d$line %in% shared_ge
  i <- match(d$line[keep], KIN$ids)
  e <- eigen(KIN$K[i,i,drop=FALSE], symmetric=TRUE)
  EG <- list(d=pmax(e$values,0), U=e$vectors, xt=as.numeric(crossprod(e$vectors, rep(1,sum(keep)))))
  bind_rows(lapply(shared_sp, function(sp){
    y <- log10(d$mat[keep, sp] + 1)
    if (!all(is.finite(y))) return(NULL)
    r <- prof_h2_eig(y, EG)
    tibble(trial=cond, Species=sp, Class=cls(sp), n=sum(keep),
           h2=r[['h2']], lo=r[['lo']], hi=r[['hi']], nonzero=r[['lo']]>0)
  }))
}
r <- bind_rows(one(C,'CTL'), one(L,'LIN'))
r$Class[is.na(r$Class)] <- ''
write_csv(r, '/home/claude/h2/check/SuppTable_S25_Heritability_PerSpecies_Paired.csv')

w <- r %>% select(trial,Species,h2) %>% pivot_wider(names_from=trial, values_from=h2)
cat(sprintf('\nspecies        %d\n', nrow(w)))
cat(sprintf('higher LIN     %d\n', sum(w$LIN > w$CTL)))
cat(sprintf('higher CTL     %d\n', sum(w$LIN < w$CTL)))
cat(sprintf('equal          %d\n', sum(w$LIN == w$CTL)))
cat(sprintf('median h2      CTL %.3f   LIN %.3f\n', median(w$CTL), median(w$LIN)))
cat(sprintf('CI excludes 0  CTL %d   LIN %d\n',
            sum(r$nonzero[r$trial=='CTL']), sum(r$nonzero[r$trial=='LIN'])))
cat('\nLPC(18:2):\n'); print(as.data.frame(r[r$Species=='LPC(18:2)',]))
