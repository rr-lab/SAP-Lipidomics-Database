# ==============================================================================
# Supplementary Table S5G -- Is the class composition of a trial a property of
#                            the whole panel, or of a subset of accessions?
#
# Two robustness checks, both computed ENTIRELY WITHIN a trial. Nothing is
# subtracted between CTL and LIN at any point, so neither result depends on the
# two trials being comparable.
#
#   1. Subsample stability. Draw 80% of the accessions without replacement,
#      recompute the mean %TIC of every class, and repeat 1000 times. Report the
#      full-panel value and the 2.5th and 97.5th percentiles of the resampled
#      means. A rank-order-preserved statistic was tried and dropped: it returns
#      100% in both trials, so like the old jackknife it cannot discriminate.
#
#   2. Leave-one-race-out. Recompute the same means with each botanical race
#      group removed in turn, and report the smallest and largest value any
#      omission produces. This one is deterministic, six omissions per trial.
#
# This replaces the leave-one-sample jackknife, which deleted one of 757
# observations from a median and so returned a perfect score for every feature,
# including features with no significant effect at all.
#
# Inputs
#   data/SPATS_fitted/non_normalized_intensities/Final_subset_{control,lowinput}_*.csv
#   data/SAP_geoloc.csv
#
# Output
#   table/supp/SuppTable_S5G_Composition_Stability.csv
# ==============================================================================
source("scripts/new_new_script/_common.R")
suppressPackageStartupMessages({ library(tibble); library(tidyr) })

B      <- as.integer(Sys.getenv("N_RESAMPLE", "1000"))
FRAC   <- 0.80
set.seed(1)

dat <- population_table()

# Composition of a set of accessions: the mean of their per-accession %TIC.
class_means <- function(d) vapply(CLASS_ORDER, function(cl) mean(d[[cl]]), numeric(1))

stability_for <- function(cond) {
  d <- dat %>% filter(Condition == cond)
  n <- nrow(d)
  k <- floor(FRAC * n)

  full <- class_means(d)

  # ---- 1. subsample ----------------------------------------------------------
  draws <- matrix(NA_real_, nrow = B, ncol = length(CLASS_ORDER),
                  dimnames = list(NULL, CLASS_ORDER))
  for (b in seq_len(B)) draws[b, ] <- class_means(d[sample.int(n, k), , drop = FALSE])
  lo <- apply(draws, 2, quantile, .025)
  hi <- apply(draws, 2, quantile, .975)

  # ---- 2. leave one race group out -------------------------------------------
  races <- setdiff(unique(as.character(d$RaceGroup)), NA_character_)
  loo <- vapply(races, function(r) class_means(d %>% filter(is.na(RaceGroup) | RaceGroup != r)),
                numeric(length(CLASS_ORDER)))

  tibble(
    Condition        = cond,
    Class            = CLASS_ORDER,
    n_accessions     = n,
    full_panel_pctTIC = unname(full),
    resample_lo95    = unname(lo),
    resample_hi95    = unname(hi),
    resample_halfwidth = unname((hi - lo) / 2),
    loo_race_min     = unname(apply(loo, 1, min)),
    loo_race_max     = unname(apply(loo, 1, max)),
    loo_race_maxshift = unname(pmax(abs(apply(loo, 1, max) - full),
                                    abs(apply(loo, 1, min) - full)))
  )
}

out <- bind_rows(lapply(c("CTL", "LIN"), stability_for))
save_table(out, "SuppTable_S5G_Composition_Stability.csv")

# ---- console summary ---------------------------------------------------------
cat(sprintf("\n%d resamples at %.0f%% of accessions, seed 1\n", B, 100 * FRAC))
for (cond in c("CTL", "LIN")) {
  d <- out %>% filter(Condition == cond)
  cat(sprintf("\n-- %s (n = %d accessions) --\n", cond, d$n_accessions[1]))
  cat(sprintf("   widest 95%% resample half-width: %.3f pp (%s)\n",
              max(d$resample_halfwidth), d$Class[which.max(d$resample_halfwidth)]))
  cat(sprintf("   largest shift from dropping any one race group: %.3f pp (%s)\n",
              max(d$loo_race_maxshift), d$Class[which.max(d$loo_race_maxshift)]))
  cat("\n")
  print(as.data.frame(d %>% dplyr::select(Class, full_panel_pctTIC, resample_lo95,
                                          resample_hi95, loo_race_min, loo_race_max) %>%
                        mutate(across(where(is.numeric), ~round(.x, 3)))))
}
