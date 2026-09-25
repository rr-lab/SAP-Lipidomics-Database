# ==============================================================================
# Class composition contrasts -- Supplementary Tables S5A, S5B and S1.
#
#   Rscript scripts/new_new_script/75_class_composition_contrasts.R
#
# One script, one definition of a lipid class, so the CLR contrast, the ALR
# contrast, the ratio statistics, Figure 2 and the composition table cannot
# disagree with each other again.
#
# WHY THIS EXISTS. Three of these tables were written by two legacy monoliths
# and had drifted apart.
#
#   S5A/S5B came from _legacy_pipeline/22_lipidome_class_composition.R, last run
#   2026-08-24. It reads the class off a feature name with "^(CLASS)(?=\\()", so
#   PC(18:2/0:0) and PC(22:0/0:0) -- which are lyso species -- were counted as
#   PC. The file therefore reports LPC at +0.264 on the CLR scale, a 1.30-fold
#   INCREASE, while the current naming gives -0.138, a 0.87-fold decrease. The
#   manuscript argues that LPC does not accumulate, so the shipped table said
#   the opposite of the paper.
#
#   S1 came from new_script/20_SuppTable1to3_ratio_species_stats.R, which
#   computes a class as mean(log10 relative abundance) over that class's
#   species. That is a per-species geometric mean, not a class total, and the
#   two can differ in SIGN. On that scale MG rises and TG/MG falls; on class
#   totals MG falls from 1.30% to 1.05% of TIC and TG/MG rises. The manuscript
#   quoted the geometric-mean ratios in one sentence and class-total values in
#   the next.
#
# WHAT THIS SCRIPT DOES INSTEAD. Every quantity below is built from the same
# object: the per-sample composition of the 13 focal classes, closed within
# those 13, with normalize_lipid_name applied first so a lyso species is a lyso
# species everywhere.
#
#   S5A  CLR contrast, mean(LIN) - mean(CTL), 500-bootstrap CI
#   S5B  ALR contrast against TG, same estimator
#   S1   class-total ratio statistics, the same 33 ratios as before
#
# Two artefact features are dropped first, as the legacy script did, for the
# reasons recorded there. Phytosphingosine goes from 0.035% of TIC under CTL to
# 16.9% under LIN, which no leaf lipidome supports, and SM(d18:1/17:0) is the
# only SM species in a plant, an odd-chain internal-standard chemotype. Neither
# produced a GWAS candidate.
#
# Inputs
#   data/SPATS_fitted/non_normalized_intensities/Final_subset_{control,lowinput}_*.csv
# Outputs
#   table/supp/SuppTable_S5A_Class_CLR_Contrast.csv
#   table/supp/SuppTable_S5B_Class_ALR_Contrast.csv
#   table/supp/SuppTable_S1_Ratio_Statistics.csv
# ==============================================================================
source("scripts/new_new_script/_common.R")
suppressPackageStartupMessages({ library(tibble) })

FOCAL    <- c("MGDG","DGDG","SQDG","PC","PE","PG","PA","PS","LPC","LPE","TG","DG","MG")
ARTEFACT <- c("Phytosphingosine", "SM(d18:1/17:0)")
N_BOOT   <- 500

# ---- per-sample composition over the 13 focal classes ------------------------
composition <- function(path) {
  x <- read_trial(path)
  m <- as.matrix(x[, -1, drop = FALSE]); storage.mode(m) <- "numeric"
  m[!is.finite(m)] <- 0
  drop <- colnames(m) %in% ARTEFACT
  if (any(drop)) message("  dropping artefact feature(s): ",
                         paste(colnames(m)[drop], collapse = ", "))
  m <- m[, !drop, drop = FALSE]
  colnames(m) <- normalize_lipid_name(colnames(m))
  cl <- sub("\\(.*$", "", colnames(m))
  s <- sapply(FOCAL, function(k) {
    j <- which(cl == k); if (length(j)) rowSums(m[, j, drop = FALSE]) else rep(0, nrow(m))
  })
  rownames(s) <- as.character(x[[1]])
  s / rowSums(s)                       # closed within the 13 focal classes
}

# Zeros are replaced by half the smallest positive part in that sample and the
# row re-closed, so the log ratios below are defined. Same rule as the legacy
# script's close_replace_rows.
close_replace <- function(mat, delta_frac = 0.5) {
  t(apply(mat, 1, function(x) {
    x[!is.finite(x) | x < 0] <- 0
    if (sum(x) <= 0) return(x)
    if (any(x == 0)) x[x == 0] <- min(x[x > 0]) * delta_frac
    x / sum(x)
  }))
}

A <- close_replace(composition(CTL_CSV))
B <- close_replace(composition(LIN_CSV))
message("CTL samples: ", nrow(A), "   LIN samples: ", nrow(B))

# ---- bootstrap contrast ------------------------------------------------------
boot_contrast <- function(a, b, seed) {
  est <- colMeans(b) - colMeans(a)
  set.seed(seed)
  bo <- t(replicate(N_BOOT,
    colMeans(b[sample(nrow(b), replace = TRUE), , drop = FALSE]) -
    colMeans(a[sample(nrow(a), replace = TRUE), , drop = FALSE])))
  tibble(Class = colnames(a), Effect = est,
         CI_Low  = apply(bo, 2, quantile, 0.025),
         CI_High = apply(bo, 2, quantile, 0.975),
         AbsEffect = abs(est), FoldChange = exp(est))
}

clr <- function(x) log(x) - rowMeans(log(x))
s5a <- boot_contrast(clr(A), clr(B), 1101)
s5a <- s5a[order(-s5a$AbsEffect), ]
write.csv(s5a, file.path(TAB_SUPP, "SuppTable_S5A_Class_CLR_Contrast.csv"), row.names = FALSE)
cat("\n-- S5A, CLR contrast (natural log; FoldChange = exp(Effect)) --\n")
print(as.data.frame(s5a %>% mutate(across(where(is.numeric), ~round(.x, 4)))))

REF <- "TG"
alr <- function(x) log(sweep(x[, setdiff(FOCAL, REF), drop = FALSE], 1, x[, REF], "/"))
s5b <- boot_contrast(alr(A), alr(B), 1201)
s5b <- s5b[order(-s5b$AbsEffect), ]
write.csv(s5b, file.path(TAB_SUPP, "SuppTable_S5B_Class_ALR_Contrast.csv"), row.names = FALSE)
message("Saved: SuppTable_S5B_Class_ALR_Contrast.csv  (reference class ", REF, ")")

# ---- S1, ratio statistics on class totals ------------------------------------
# The same 33 ratios the previous S1 carried, so nothing the manuscript cites
# disappears; only the definition of a class changes, from a per-species
# geometric mean to the class total.
RATIOS <- c(
  "MG/SQDG","MG/MGDG","PS/SQDG","DGDG/MG","PG/SQDG","DG/MG","PE/SQDG",
  "LPC/PS","DGDG/PS","DG/PS","PC/PS","LPE/SQDG","LPC/MG","PE/PS",
  "MGDG/PG","PG/PS","LPC/LPE","MGDG/PE","MG/PG","SQDG/TG","PA/PS",
  "LPE/MGDG","PC/SQDG",
  "TG/PE","TG/MGDG","TG/MG","TG/DG","PS/PC","PS/MGDG","PS/DG",
  "PG/MGDG","PG/DGDG","PE/PA")

s1 <- do.call(rbind, lapply(RATIOS, function(r) {
  p <- strsplit(r, "/", fixed = TRUE)[[1]]
  a <- log10(A[, p[1]] / A[, p[2]]); b <- log10(B[, p[1]] / B[, p[2]])
  w <- wilcox.test(b, a)
  eff <- median(b) - median(a)
  data.frame(Ratio = r, n_C = length(a), n_LI = length(b),
             median_C = round(median(a), 4), median_LI = round(median(b), 4),
             effect_log10 = round(eff, 4), effect_fc = round(10^eff, 2),
             direction = if (eff > 0) "LI > C" else "LI < C",
             p_wilcox = w$p.value)
}))
s1$p_adj_BH <- p.adjust(s1$p_wilcox, "BH")
s1$significance <- cut(s1$p_adj_BH, c(-Inf, 1e-4, 1e-3, 0.01, 0.05, Inf),
                       labels = c("****", "***", "**", "*", "ns"))
s1 <- s1[order(-abs(s1$effect_log10)), ]
write.csv(s1, file.path(TAB_SUPP, "SuppTable_S1_Ratio_Statistics.csv"), row.names = FALSE)
cat("\n-- S1, ratio statistics on class totals --\n")
print(s1[, c("Ratio","median_C","median_LI","effect_log10","effect_fc","direction","significance")],
      row.names = FALSE)

cat("\n-- class composition, % of the 13 focal classes --\n")
print(data.frame(Class = FOCAL,
                 CTL = round(100 * colMeans(A), 4),
                 LIN = round(100 * colMeans(B), 4)), row.names = FALSE)
