# ==============================================================================
# 12c  Class sums and non-redundant class log-ratios for the sum/ratio GWAS
#
# Replaces  Final_{control,lowinput}_BLUPs_class_sums_and_ratios.csv
#
# Why this rebuild exists. The shipped ratio phenotypes were
#
#     r_AB = log10(A + shift_A) - log10( max(B + shift_B, floor_B) )
#
# with shift_c = -min(c) + 1e-6 and floor_c the 1st percentile. Two defects
# follow from that definition.
#
#   (1) Redundancy. All 306 ordered pairs were built, but log(A/B) = -log(B/A),
#       and a GWAS on y and on -y is the same test. 153 of the 306 were the
#       mirror of another. Non-redundant count is 18 sums + C(18,2) = 171.
#
#   (2) A floor/shift artifact. shift_c = -min(c) + 1e-6 places the lowest
#       genotype of each class at 1e-6 while the class median sits near 1e7,
#       i.e. 10.4 to 13.9 log10 units below the rest of the panel. That
#       genotype spikes when its class is the numerator and is clamped back by
#       floor_c when its class is the denominator, so the two directions are
#       NOT mirrors and give different candidate genes. 136 of 306 ratios carry
#       such a spike.
#
# Fix. BLUPs are deviations from the class mean, so ~62% of values are negative
# and cannot be logged directly. Rather than shift by -min, add the class
# intercept back: the mean class sum of the SpATS fitted intensities. That puts
# each genotype on the original intensity scale, keeps the spatial correction,
# and is strictly positive with margin, so log10(A/B) needs no shift and no
# floor and is exactly antisymmetric. Only 153 ratios are then built.
#
# Input   data/SPATS_fitted/BLUP_GWAS_phenotype/Final_{control,lowinput}_all_lipids_BLUPs.csv
#         data/summed_lipid_intensities/{A,B}_final_summed_lipids.csv   (intercept source)
# Output  data/SPATS_fitted/BLUP_GWAS_phenotype/Final_{control,lowinput}_BLUPs_class_sums_and_logratios.csv
#         data/SPATS_fitted/BLUP_GWAS_phenotype/Final_all_lipids_BLUPs_class_logratios.meta.csv
#         data/SPATS_fitted/BLUP_GWAS_phenotype/Final_all_lipids_BLUPs_class_logratios.qc.csv
# ==============================================================================
suppressPackageStartupMessages({ library(vroom) })

REPO      <- Sys.getenv("SOLD_REPO", ".")
DATA_ROOT <- Sys.getenv("SOLD_DATA", file.path(REPO, "data"))
BLUP_DIR  <- file.path(DATA_ROOT, "SPATS_fitted/BLUP_GWAS_phenotype")
SUM_DIR   <- file.path(DATA_ROOT, "summed_lipid_intensities")

TRIALS <- list(
  control  = list(blup = file.path(BLUP_DIR, "Final_control_all_lipids_BLUPs.csv"),
                  int  = file.path(SUM_DIR,  "A_final_summed_lipids.csv")),
  lowinput = list(blup = file.path(BLUP_DIR, "Final_lowinput_all_lipids_BLUPs.csv"),
                  int  = file.path(SUM_DIR,  "B_final_summed_lipids.csv")))

# The summed-intensity matrix is stored with compounds as rows and genotypes as
# columns, but it also carries annotation and instrument columns (MZErrorPPM and
# friends). Selecting "everything that is not annotation" pulls those in and
# skews the mean, so the sample columns are taken as the intersection with the
# genotypes actually present in the BLUP matrix.

# Pipeline class order. Numerator of a kept ratio is whichever class comes first
# here, so the retained direction is deterministic and documented.
CLASSES <- c("DGDG","MGDG","SQDG","GalCer","LPC","LPE","Cer","AEG","SM","FA",
             "TG","DG","MG","PC","PE","PG","PA","PS")

# A class is usable for ratios only if BLUP + intercept is strictly positive for
# every genotype. A class that fails is NOT rescued with a shift, because a shift
# is what produced the artifact this script exists to remove. It is dropped and
# reported. With MATCH_TRIALS the drop is applied to both trials so CTL and LIN
# carry the same trait set and recurrence counts stay comparable.
MATCH_TRIALS <- TRUE
ID_COLS <- c("LineRaw","Sample","Line","PlotID","row","col")

# CTL writes lyso species as CLASS(x:y/0:0), LIN as LCLASS(x:y). Same molecule.
# The 3-position guard keeps DG(18:0/18:2/0:0) a diacylglycerol.
normalize_lipid_name <- function(x) {
  x <- sub("^(PC|PE|PG|PS|PA)\\(([^/()]+)/0:0\\)$", "L\\1(\\2)", x)
  x <- sub("^(LPC|LPE)\\(([^/()]+)/0:0\\)$", "\\1(\\2)", x)
  x
}
lipid_class <- function(x) sub("\\(.*$", "", normalize_lipid_name(x))

# Class sums per genotype from a genotype x species matrix.
class_sums <- function(path) {
  d   <- vroom(path, show_col_types = FALSE, progress = FALSE)
  key <- names(d)[1]
  lip <- setdiff(names(d), ID_COLS)
  k   <- lipid_class(lip)
  line <- d[[key]]
  m <- sapply(CLASSES, function(cl) {
    sp <- lip[k == cl]
    if (!length(sp)) stop("no species in class ", cl, " for ", basename(path))
    rowSums(as.matrix(d[, sp, drop = FALSE]), na.rm = TRUE)
  })
  ag <- aggregate(as.data.frame(m), list(LineRaw = line), mean, na.rm = TRUE)
  list(line = ag$LineRaw,
       m    = as.matrix(ag[, CLASSES, drop = FALSE]),
       sp_by_class = split(lip, factor(k, levels = CLASSES))[CLASSES])
}

# Class intercept, i.e. the panel mean class abundance on the input intensity
# scale. Taken from the intensity matrix the SpATS model was fitted to, which is
# stored with compounds as rows. Spot checks against the fitted matrix agree to
# within 0.7%, so the two are on the same scale; this file is used because it
# also carries species the artifact-filtered fitted matrix dropped.
class_intercepts <- function(path, sp_by_class, genotypes) {
  d <- vroom(path, show_col_types = FALSE, progress = FALSE)
  samp <- intersect(names(d), genotypes)
  if (length(samp) < 0.9 * length(genotypes))
    stop("intensity matrix ", basename(path), " covers only ", length(samp),
         " of ", length(genotypes), " genotypes")
  mu <- setNames(
    suppressWarnings(apply(as.matrix(d[, samp, drop = FALSE]), 1,
                           function(v) mean(as.numeric(v), na.rm = TRUE))),
    d$Compound_Name)
  v <- sapply(CLASSES, function(cl) {
    sp <- sp_by_class[[cl]]
    miss <- setdiff(sp, names(mu))
    if (length(miss))
      stop("class ", cl, ": no intensity row for ", paste(miss, collapse = ", "))
    sum(mu[sp], na.rm = TRUE)
  })
  attr(v, "n_samples") <- length(samp)
  v
}

build_trial <- function(tag, paths) {
  b <- class_sums(paths$blup)                              # BLUPs, centred near zero
  intercept <- class_intercepts(paths$int, b$sp_by_class, b$line)  # intensity scale

  # Recentred class abundance. Strictly positive, spatially corrected, and on
  # the original intensity scale. This is the quantity the ratios are taken of.
  A <- sweep(b$m, 2, intercept[CLASSES], "+")

  bad <- CLASSES[apply(!is.finite(A) | A <= 0, 2, any)]
  if (length(bad))
    message("  !! ", tag, ": not strictly positive in ", paste(bad, collapse = ", "),
            " (", paste(sprintf("%d/%d genotypes", colSums(A[, bad, drop = FALSE] <= 0),
                                nrow(A)), collapse = ", "), ")")

  qc <- data.frame(
    Trial = tag, Class = CLASSES,
    N_species = sapply(b$sp_by_class, length)[CLASSES],
    N_genotypes_for_intercept = attr(intercept, "n_samples"),
    Class_intercept = unname(intercept[CLASSES]),
    BLUP_min = apply(b$m, 2, min), BLUP_max = apply(b$m, 2, max),
    Recentred_min = apply(A, 2, min), Recentred_median = apply(A, 2, median),
    N_genotypes_nonpositive = colSums(A <= 0),
    Usable_for_ratios = !(CLASSES %in% bad),
    row.names = NULL)
  qc$Log10_spread_min_to_median <- ifelse(qc$Usable_for_ratios,
    log10(apply(A, 2, median)) - log10(pmax(apply(A, 2, min), .Machine$double.xmin)), NA_real_)

  list(line = b$line, sums = b$m, A = A, qc = qc, bad = bad, n_geno = length(b$line))
}

# Assemble the phenotype matrix once the usable class set is known.
assemble <- function(r, keep) {
  logA <- log10(r$A[, keep, drop = FALSE])
  out  <- data.frame(Line = r$line, r$sums[, keep, drop = FALSE], check.names = FALSE)
  pairs <- t(combn(keep, 2))
  for (i in seq_len(nrow(pairs))) {
    A1 <- pairs[i, 1]; B1 <- pairs[i, 2]
    out[[paste0("Sum_", A1, "_over_", B1, "_log10ratio")]] <- logA[, A1] - logA[, B1]
  }
  list(out = out, pairs = pairs, logA = logA)
}

res <- lapply(names(TRIALS), function(tg) build_trial(tg, TRIALS[[tg]]))
names(res) <- names(TRIALS)

dropped <- if (MATCH_TRIALS) unique(unlist(lapply(res, `[[`, "bad"))) else NULL
KEEP <- setdiff(CLASSES, dropped)
if (length(dropped))
  message("\nDropped from the ratio set in BOTH trials: ", paste(dropped, collapse = ", "),
          "\n  (a shift would hide this, so the class is excluded instead)")

asm <- lapply(names(res), function(tg) assemble(res[[tg]],
  if (MATCH_TRIALS) KEEP else setdiff(CLASSES, res[[tg]]$bad)))
names(asm) <- names(res)

# ---- checks ------------------------------------------------------------------
# Antisymmetry is now exact by construction. Assert it rather than assume it.
for (tg in names(asm)) {
  a <- asm[[tg]]; r <- res[[tg]]
  worst <- max(apply(a$pairs, 1, function(p)
    max(abs((a$logA[, p[1]] - a$logA[, p[2]]) + (a$logA[, p[2]] - a$logA[, p[1]])))))
  if (worst > 1e-12) stop(tg, ": ratios are not antisymmetric, worst = ", worst)
  nk <- ncol(a$logA); nt <- ncol(a$out) - 1
  stopifnot(nt == nk + nrow(a$pairs), nrow(a$pairs) == choose(nk, 2))
  message(sprintf("%-9s %d genotypes   %d traits (%d sums + %d ratios)   antisymmetry residual %.1e",
                  tg, r$n_geno, nt, nk, nrow(a$pairs), worst))
}

# ---- write -------------------------------------------------------------------
for (tg in names(asm)) {
  p <- file.path(BLUP_DIR, sprintf("Final_%s_BLUPs_class_sums_and_logratios.csv", tg))
  write.csv(asm[[tg]]$out, p, row.names = FALSE)
  message("Saved: ", p, "  (", nrow(asm[[tg]]$out), " rows x ", ncol(asm[[tg]]$out), " cols)")
}

qc <- do.call(rbind, lapply(res, `[[`, "qc"))
qc$Kept_in_output <- qc$Class %in% KEEP
write.csv(qc, file.path(BLUP_DIR, "Final_all_lipids_BLUPs_class_logratios.qc.csv"), row.names = FALSE)

pr <- asm[[1]]$pairs
meta <- data.frame(
  Trait = c(KEEP, paste0("Sum_", pr[,1], "_over_", pr[,2], "_log10ratio")),
  Type  = c(rep("class_sum", length(KEEP)), rep("class_log10ratio", nrow(pr))),
  Definition = c(rep("sum of member-species BLUPs, unchanged from the previous file", length(KEEP)),
                 paste0("log10(", pr[,1], " + class intercept) - log10(", pr[,2], " + class intercept)")))
write.csv(meta, file.path(BLUP_DIR, "Final_all_lipids_BLUPs_class_logratios.meta.csv"), row.names = FALSE)
message("Saved: meta and qc.  Traits per trial = ", nrow(meta))

# ---- split into 10-phenotype chunks -----------------------------------------
# The GWAS is driven one chunk at a time, so each part carries the Line key plus
# at most CHUNK trait columns. Traits keep the order of the combined file, so
# part k holds traits ((k-1)*CHUNK + 1) .. (k*CHUNK).
CHUNK <- 10
for (tg in names(asm)) {
  d  <- asm[[tg]]$out
  tr <- setdiff(names(d), "Line")
  idx <- split(seq_along(tr), ceiling(seq_along(tr) / CHUNK))
  for (k in seq_along(idx)) {
    p <- file.path(BLUP_DIR,
      sprintf("Final_%s_BLUPs_class_sums_and_logratios_%d.csv", tg, k))
    write.csv(d[, c("Line", tr[idx[[k]]]), drop = FALSE], p, row.names = FALSE)
  }
  message(sprintf("%-9s split into %d parts of <=%d traits (%d traits total)",
                  tg, length(idx), CHUNK, length(tr)))
}

# Index so a result file can be traced back to the part it came from.
part_index <- do.call(rbind, lapply(names(asm), function(tg) {
  tr <- setdiff(names(asm[[tg]]$out), "Line")
  data.frame(Trial = tg, Part = ceiling(seq_along(tr) / CHUNK),
             File = sprintf("Final_%s_BLUPs_class_sums_and_logratios_%d.csv",
                            tg, ceiling(seq_along(tr) / CHUNK)),
             Trait = tr, row.names = NULL)
}))
write.csv(part_index, file.path(BLUP_DIR, "Final_all_lipids_BLUPs_class_logratios.parts.csv"),
          row.names = FALSE)
message("Saved: parts index (", nrow(part_index), " rows)")
