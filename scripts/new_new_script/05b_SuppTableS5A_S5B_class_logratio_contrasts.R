# ==============================================================================
# Supplementary Tables S5A and S5B -- the CTL/LIN class contrast on the centred
# and additive log-ratio scales.
#
#   Rscript scripts/new_new_script/05b_SuppTableS5A_S5B_class_logratio_contrasts.R
#
# WHY THIS EXISTS. Both tables shipped in the supplement for months with no
# generating script anywhere in the repository. They were reproduced here on
# 2026-09-17 by fitting the definition below against the shipped numbers, which
# it recovers for ten of the thirteen classes to within 0.005 on the effect and
# on both interval bounds. The three it does not recover -- PA, PS and TG -- are
# the classes whose species membership changed in the 2026-09-03 deduplication,
# which the shipped tables predate. That is the point of rebuilding them.
#
# THE DEFINITION, stated once.
#
#   Class sums are taken per sample over the thirteen focal classes, from the
#   SpATS-fitted non-normalised intensities. A sample contributing a zero or
#   negative sum in any class is dropped from that class.
#
#   CLR   clr_i = log(x_i) - mean_j log(x_j), natural log, the mean running over
#         the same thirteen classes. The geometric mean is the reference, so the
#         thirteen values sum to zero and no class is privileged.
#
#   ALR   alr_i = log(x_i / x_TG) = clr_i - clr_TG. TG is the reference, so TG
#         has no row of its own and the table carries twelve.
#
#   Effect is mean(LIN) - mean(CTL) on whichever scale, and the interval is the
#   95% Welch two-sample interval for that difference. Welch rather than pooled
#   because the two trials differ in sample size (394 against 363) and in
#   variance.
#
# READ THE TWO TABLES TOGETHER, NOT SEPARATELY. A CLR effect is a statement
# about a class relative to the average class, and an ALR effect is a statement
# about it relative to TG. A class can rise on one and fall on the other without
# either being wrong, and LPC does exactly that. Neither is the %TIC pool share
# in Supplementary Table S5D, which is a third quantity again.
#
# Input   data/SPATS_fitted/non_normalized_intensities/Final_subset_{control,lowinput}_*.csv
# Output  table/supp/SuppTable_S5A_Class_CLR_Contrast.csv
#         table/supp/SuppTable_S5B_Class_ALR_Contrast.csv
# ==============================================================================
suppressPackageStartupMessages({ library(vroom) })

REPO      <- Sys.getenv("SOLD_REPO", ".")
DATA_ROOT <- Sys.getenv("SOLD_DATA", file.path(REPO, "data"))
SPATS     <- file.path(DATA_ROOT, "SPATS_fitted/non_normalized_intensities")
OUT       <- Sys.getenv("SUPP_OUT", file.path(REPO, "table/supp"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

TRIALS <- c(CTL = file.path(SPATS, "Final_subset_control_all_lipids_fitted_phenotype_non_normalized.csv"),
            LIN = file.path(SPATS, "Final_subset_lowinput_all_lipids_fitted_phenotype_non_normalized.csv"))

# The thirteen focal classes, in the order Figure 2B stacks them.
FOCAL   <- c("MGDG", "DGDG", "SQDG", "PC", "PE", "PG", "PA", "PS",
             "LPC", "LPE", "TG", "DG", "MG")
ALR_REF <- "TG"

# CTL writes a lyso species as PC(18:2/0:0); LIN writes it as LPC(18:2). Both
# spellings must fold to one before any class sum is taken, or CTL's lyso
# species land in the diacyl class and LIN's do not.
normalize_lipid_name <- function(x) sub("^([A-Z]+)\\(([^/]+)/0:0\\)$", "L\\1(\\2)", x)

class_sums <- function(path) {
  d <- vroom(path, show_col_types = FALSE, progress = FALSE)
  m <- as.matrix(d[, -1])
  colnames(m) <- normalize_lipid_name(names(d)[-1])
  s <- sapply(FOCAL, function(cl) {
    sp <- grep(paste0("^", cl, "\\("), colnames(m), value = TRUE)
    if (!length(sp)) return(rep(NA_real_, nrow(m)))
    rowSums(m[, sp, drop = FALSE], na.rm = TRUE)
  })
  s[s <= 0] <- NA
  s
}

log_ratio <- function(s, scale) {
  L <- log(s)
  if (scale == "clr") L - rowMeans(L, na.rm = TRUE) else L - L[, ALR_REF]
}

contrast <- function(scale) {
  a <- log_ratio(class_sums(TRIALS[["CTL"]]), scale)
  b <- log_ratio(class_sums(TRIALS[["LIN"]]), scale)
  cls <- if (scale == "alr") setdiff(FOCAL, ALR_REF) else FOCAL
  out <- do.call(rbind, lapply(cls, function(cl) {
    tt <- t.test(b[, cl], a[, cl])          # Welch, LIN - CTL
    eff <- unname(tt$estimate[1] - tt$estimate[2])
    data.frame(Class = cl, Effect = eff,
               CI_Low = tt$conf.int[1], CI_High = tt$conf.int[2],
               AbsEffect = abs(eff), stringsAsFactors = FALSE)
  }))
  out[order(-out$AbsEffect), ]
}

s5a <- contrast("clr")
s5b <- contrast("alr")

write.csv(s5a[, c("Class", "Effect", "CI_Low", "CI_High", "AbsEffect")],
          file.path(OUT, "SuppTable_S5A_Class_CLR_Contrast.csv"), row.names = FALSE)
# S5B keeps the shipped column order, which leads with Feature and repeats the
# class name, because the manuscript cites it by those headings.
write.csv(data.frame(Feature = s5b$Class, Effect = s5b$Effect,
                     CI_Low = s5b$CI_Low, CI_High = s5b$CI_High,
                     Class = s5b$Class, AbsEffect = s5b$AbsEffect),
          file.path(OUT, "SuppTable_S5B_Class_ALR_Contrast.csv"), row.names = FALSE)

message(sprintf("S5A  %d classes on the CLR scale", nrow(s5a)))
message(sprintf("S5B  %d classes on the ALR scale, reference %s", nrow(s5b), ALR_REF))
cat("\n-- CLR, ordered by absolute effect --\n")
print(within(s5a, { Effect <- round(Effect, 4); CI_Low <- round(CI_Low, 4)
                    CI_High <- round(CI_High, 4); AbsEffect <- NULL }), row.names = FALSE)
cat("\n-- ALR against ", ALR_REF, " --\n", sep = "")
print(within(s5b, { Effect <- round(Effect, 4); CI_Low <- round(CI_Low, 4)
                    CI_High <- round(CI_High, 4); AbsEffect <- NULL }), row.names = FALSE)
