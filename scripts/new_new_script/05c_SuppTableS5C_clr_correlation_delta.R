# ==============================================================================
# Supplementary Table S5C -- within-trial CLR correlations between lipid
# classes, and how they change from CTL to LIN.
#
#   Rscript scripts/new_new_script/05c_SuppTableS5C_SuppFigS6_clr_correlations.R
#
# WHY THIS EXISTS. Both were produced by scripts/_legacy_pipeline/
# 28_class_logratio_stats.R, a 50 KB monolith that in the same run also rebuilds
# Figure2_OPLS_DA and writes SuppTable_S7_OPLS_VIP_ratios.csv -- a table that was
# removed from the supplement in September 2026. Re-running it to refresh S5C
# put the removed table back and overwrote a main figure. This is that section
# lifted out, the same way 20_SuppTable1to3_ratio_species_stats.R was lifted.
#
# The statistics are unchanged from the original. Two things are not.
#
#   Lyso names are normalised before the class sums are taken. CTL writes a
#   lyso species as PC(18:2/0:0) and LIN writes it as LPC(18:2); without this
#   the CTL spelling is summed into the diacyl PC pool. The shipped S5C predates
#   the fix and predates the 2026-09-03 species deduplication, so every
#   correlation involving LPC or PC in it is off.
#
# Input   data/SPATS_fitted/non_normalized_intensities/Final_subset_{control,lowinput}_*.csv
# Output  table/supp/SuppTable_S5C_Class_CLR_Correlation_Delta.csv
# ==============================================================================
source("scripts/new_new_script/_common.R")
suppressPackageStartupMessages({ library(tidyr) })

SPATS_DIR <- file.path(DATA_ROOT, "SPATS_fitted/non_normalized_intensities")
TRIALS <- c(CTL = file.path(SPATS_DIR, "Final_subset_control_all_lipids_fitted_phenotype_non_normalized.csv"),
            LIN = file.path(SPATS_DIR, "Final_subset_lowinput_all_lipids_fitted_phenotype_non_normalized.csv"))

FOCAL <- c("MGDG", "DGDG", "SQDG", "PC", "PE", "PG", "PA", "PS",
           "LPC", "LPE", "TG", "DG", "MG")

normalize_lyso <- function(x) sub("^([A-Z]+)\\(([^/]+)/0:0\\)$", "L\\1(\\2)", x)

clr_matrix <- function(path) {
  d <- vroom(path, show_col_types = FALSE, progress = FALSE)
  m <- as.matrix(d[, -1]); colnames(m) <- normalize_lyso(names(d)[-1])
  s <- sapply(FOCAL, function(cl) {
    sp <- grep(paste0("^", cl, "\\("), colnames(m), value = TRUE)
    if (!length(sp)) return(rep(NA_real_, nrow(m)))
    rowSums(m[, sp, drop = FALSE], na.rm = TRUE)
  })
  s[s <= 0] <- NA
  L <- log(s)
  L[complete.cases(L), , drop = FALSE] - rowMeans(L[complete.cases(L), , drop = FALSE])
}

cor_and_p <- function(X) {
  k <- ncol(X); r <- matrix(NA_real_, k, k, dimnames = list(colnames(X), colnames(X)))
  p <- r
  for (i in seq_len(k)) for (j in seq_len(k)) {
    ct <- suppressWarnings(cor.test(X[, i], X[, j]))
    r[i, j] <- ct$estimate; p[i, j] <- ct$p.value
  }
  list(r = r, p = p, n = nrow(X))
}

CTL <- cor_and_p(clr_matrix(TRIALS[["CTL"]]))
LIN <- cor_and_p(clr_matrix(TRIALS[["LIN"]]))
stopifnot(identical(rownames(CTL$r), rownames(LIN$r)))

# ---- S5C ---------------------------------------------------------------------
pairs <- t(combn(rownames(CTL$r), 2))
s5c <- data.frame(
  Class_A = pairs[, 1], Class_B = pairs[, 2],
  r_CTL = CTL$r[pairs], r_LIN = LIN$r[pairs],
  p_CTL = CTL$p[pairs], p_LIN = LIN$p[pairs],
  stringsAsFactors = FALSE)
s5c$delta_r       <- s5c$r_LIN - s5c$r_CTL
s5c$sign_reversal <- sign(s5c$r_CTL) != sign(s5c$r_LIN)
s5c$q_CTL <- p.adjust(s5c$p_CTL, "BH")
s5c$q_LIN <- p.adjust(s5c$p_LIN, "BH")
s5c <- s5c[order(-abs(s5c$delta_r)), ]
save_table(s5c, "SuppTable_S5C_Class_CLR_Correlation_Delta.csv")

# NO FIGURE IS DRAWN HERE. 08_SuppFig7_class_correlations.R owns
# fig/supp/SuppFig_S6_CLR_Correlations.png and draws it as two panels, CTL and
# LIN, with no delta-r panel. That third panel is a CTL-LIN difference in
# compositional space, which the manuscript deliberately stopped reporting, and
# the caption on that figure now says in as many words that the two matrices are
# not differenced. An earlier draft of this script redrew all three panels and
# would have put the removed comparison back into the supplement.
#
# The delta_r column below stays in the table, because the sign reversals it
# records are cited in the Results and the Discussion. A number in a table a
# reader looks up is not the same claim as a panel in a figure.

# A reversal only counts where the correlation is significant in BOTH trials; a
# pair that is noise in one of them has not reversed, it was never there.
n_real <- sum(s5c$sign_reversal & s5c$q_CTL < .05 & s5c$q_LIN < .05)
message(sprintf("S5C  %d class pairs, %d sign reversals, %d of them significant in both trials",
                nrow(s5c), sum(s5c$sign_reversal), n_real))
cat("\n-- largest changes in correlation --\n")
print(head(within(s5c, { r_CTL <- round(r_CTL, 3); r_LIN <- round(r_LIN, 3)
                         delta_r <- round(delta_r, 3)
                         p_CTL <- NULL; p_LIN <- NULL; q_CTL <- NULL; q_LIN <- NULL }), 10),
      row.names = FALSE)
