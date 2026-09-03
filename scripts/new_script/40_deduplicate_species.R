# ==============================================================================
# 40_deduplicate_species.R
#
# Remove redundant lipid species from the SpATS-fitted phenotype matrices.
#
# Sixteen pairs of "species" in the two matrices correlate at r = 1.0000 with a
# constant ratio between them. They are not two lipids. They are one MS feature
# that received two library annotations and was then carried through the
# pipeline twice. Evidence, from data/summed_lipid_intensities/{A,B}_final_summed_lipids.csv:
#
#   DG(18:1/18:1) and MG(18:1) share scan 11369 (CTL) / 10144 (LIN) at
#   m/z 339.289. That mass is a monoacylglycerol; DG(18:1/18:1) would be near
#   m/z 621. The DG row is flagged ambiguous with three candidate names, one of
#   which is "Monoolein". The DG annotation is wrong and the feature is an MG.
#
#   13-Keto-9Z,11E-octadecadienoic acid and 13S-HOTrE share scan 3647 at
#   m/z 277.216 and have the same formula. They cannot be distinguished.
#
#   The TG pairs are isobaric or near-isobaric acyl assignments of one feature.
#
#   LPC(18:3)/PC(18:3/0:0), LPE(16:0)/PE(16:0/0:0) and SPB 18:0;2OH/SPHINGANINE
#   are the same compound written under two nomenclatures.
#
#   Lutein/Zeaxanthin and alpha/beta-Carotene are isomer pairs that this
#   chromatography does not resolve.
#
# One member of each pair is kept. Values are left unchanged; only the redundant
# column is removed. The retained member is chosen on MS evidence where the
# evidence discriminates (a real scan, a non-ambiguous library match, the higher
# MQScore) and on nomenclature otherwise.
#
# Magnitude. Within each pair one member is an exact 2x or 3x multiple of the
# other. Tracing the pipeline shows why. data/raw_lipid_intensities/A_cleaned_lipids.csv
# holds 617 rows under 501 distinct names, because 116 names carry more than one
# row: alternative library annotations of the same lipid, collapsed onto a common
# name and then summed. LPC(18:3) is the sum of three such annotations, one of
# which is PC(18:3/0:0); LPE(16:0) is the sum of two, one of which is
# PE(16:0/0:0); SPB 18:0;2OH is the sum of two, one of which is SPHINGANINE. The
# merged name therefore carries the intended total and the other member of the
# pair is a single un-merged variant that survived alongside it.
#
# The retained column is given the merged total, so these species stay on the
# same footing as the other 113 merged names already in the data set. Where the
# retained name is not the merged member, its values are replaced by the merged
# member's, which is exact because the two are proportional.
#
# One caveat, recorded rather than acted on. Of the 116 merged groups, 81 combine
# rows that are exact scalar multiples of one another and 35 combine genuinely
# different measurements. Summing is unambiguously right for the 35. For the 81
# the summed total counts one measurement more than once, so the absolute
# per cent-of-TIC of those species is high by the merge factor. Correcting only
# the pairs found here would make them inconsistent with the other merged names,
# which is why the merged total is kept throughout. It affects composition
# percentages only; GWAS, heritability, correlations and every log-ratio are
# scale-invariant and unaffected.
#
# Deduplication is per trial. A pair that is redundant under CTL is not always
# redundant under LIN, so only the column that is actually redundant in a given
# matrix is removed from it.
#
# Both phenotype families are treated, because both are read downstream. The
# SpATS-fitted matrices feed every figure and table; the BLUP matrices under
# BLUP_GWAS_phenotype/ are what the GWAS and the class sum/ratio builder
# (SoLD_paper/scripts/4.1_sum_ratios_lipids.R) actually read. Deduplicating one
# and not the other would leave the GWAS running on the duplicated set.
#
# Input   data/SPATS_fitted/non_normalized_intensities/Final_subset_{control,lowinput}_*.csv
#         data/SPATS_fitted/BLUP_GWAS_phenotype/Final_{control,lowinput}_all_lipids_BLUPs.csv
#         data/final_species_set/species_inventory.csv
# Output  those four phenotype files, deduplicated, with the originals archived
#         alongside under _pre_dedup/
#         data/final_species_set/species_inventory.csv, rebuilt from them
#         decision, count and class tables in table/supp/
#
# Run from the repository root.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table) })

REPO  <- Sys.getenv("SOLD_DB_REPO", ".")
SPATS <- file.path(REPO, "data/SPATS_fitted/non_normalized_intensities")
BLUP  <- file.path(REPO, "data/SPATS_fitted/BLUP_GWAS_phenotype")
OUT   <- file.path(REPO, "table/supp")
META  <- c("LineRaw", "PlotID", "row", "col")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

FILES <- c(CTL = "Final_subset_control_all_lipids_fitted_phenotype_non_normalized.csv",
           LIN = "Final_subset_lowinput_all_lipids_fitted_phenotype_non_normalized.csv")

# every phenotype matrix that carries the duplicated species
targets <- data.table(
  family = c("spats", "spats", "blup", "blup"),
  trial  = c("CTL", "LIN", "CTL", "LIN"),
  dir    = c(SPATS, SPATS, BLUP, BLUP),
  file   = c(FILES[["CTL"]], FILES[["LIN"]],
             "Final_control_all_lipids_BLUPs.csv",
             "Final_lowinput_all_lipids_BLUPs.csv"))
for (d in unique(targets$dir)) dir.create(file.path(d, "_pre_dedup"), recursive = TRUE, showWarnings = FALSE)
ARCH <- file.path(SPATS, "_pre_dedup")   # the inventory archive lives with the SpATS files

# ---- decisions ---------------------------------------------------------------
# trial "BOTH" applies wherever the column is present.
decisions <- fread(text = '
trial|keep|drop|rename_keep_to|reason
BOTH|MG(18:1)|DG(18:1/18:1)||one feature at m/z 339.289 is a monoacylglycerol; the DG call is the ambiguous one and lists Monoolein among its candidates
BOTH|LPC(18:3)|PC(18:3/0:0)||same compound; PC(18:3/0:0) is lyso-PC written in full-position notation
BOTH|LPE(16:0)|PE(16:0/0:0)||same compound; PE(16:0/0:0) is lyso-PE written in full-position notation
BOTH|SPB 18:0;2OH|SPHINGANINE||same compound; a dihydroxy C18:0 sphingoid base is sphinganine
BOTH|beta-Carotene|alpha-Carotene||beta-Carotene has a real scan and a non-ambiguous match; the alpha call has neither
BOTH|beta-Carotene|Alpha-Carotene||as above, LIN spelling of the same redundant call
CTL|Lutein|Zeaxanthin||unresolved isomers of one C40H56O2 feature under CTL, where the two correlate at r = 1; the LIN run resolves them (r = 0.48) and both are kept there
LIN|13S-HOTrE|13-Keto-9Z,11E-octadecadienoic acid||same scan, same m/z, same formula; HOTrE has the higher MQScore
CTL|TG(12:0/12:0/14:0)|TG(10:0/12:0/16:0)||isobaric TG 38:0; the retained call is non-ambiguous, the dropped one has three candidates
CTL|TG(14:0/14:0/16:1)|TG(12:0/14:0/18:1)||isobaric TG 44:1; the retained call is non-ambiguous
CTL|TG(12:0/12:0/18:2)|TG(8:0/16:1/18:1)||isobaric TG 42:2; an 8:0 acyl chain is implausible in sorghum leaf
CTL|TG(12:0/16:0/16:0)|TG(14:0/14:0/16:0)||isobaric TG 44:0; retained call also survives deduplication in LIN
LIN|TG(10:0/14:0/16:0)|TG(12:0/12:0/18:1)||one feature carried under two acyl assignments; the retained call has fewer candidates
LIN|TG(12:0/16:0/16:0)|TG(14:0/14:0/18:0)||one feature carried under two acyl assignments; retained call matches the CTL decision
LIN|TG(18:1/18:2/18:3)|TG(18:0/18:2/20:1)||one feature carried under two acyl assignments; TG 54:6 is far more plausible in sorghum leaf than TG 56:3
', sep = "|")

R_DUP <- 0.9999   # correlation above which two columns are one measurement
log <- list(); counts <- list(); skipped <- list()

for (ti in seq_len(nrow(targets))) {
  tr  <- targets$trial[ti]; fam <- targets$family[ti]
  f   <- file.path(targets$dir[ti], targets$file[ti])
  if (!file.exists(f)) { warning("missing: ", f, call. = FALSE); next }
  d   <- fread(f)
  arc <- file.path(targets$dir[ti], "_pre_dedup", targets$file[ti])
  if (!file.exists(arc)) file.copy(f, arc)

  dec <- decisions[trial %in% c("BOTH", tr)]
  n_before <- length(setdiff(names(d), META))

  for (i in seq_len(nrow(dec))) {
    k <- dec$keep[i]; dr <- dec$drop[i]
    rn <- dec$rename_keep_to[i]; if (is.na(rn)) rn <- ""
    if (!(dr %in% names(d))) next
    if (!(k %in% names(d))) next
    r <- suppressWarnings(cor(d[[k]], d[[dr]], use = "pairwise.complete.obs"))
    # guard: only ever remove a column that really is the same measurement.
    if (is.na(r) || r <= R_DUP) {
      warning(sprintf("%s/%s: %s and %s correlate at r = %.4f, below the duplicate threshold; kept both",
                      fam, tr, k, dr, r), call. = FALSE)
      skipped[[length(skipped) + 1]] <- data.table(family = fam, trial = tr,
                                                   kept = k, not_dropped = dr, r = round(r, 6))
      next
    }
    # keep the merged total under the retained name. The pair is proportional, so
    # taking the larger member's vector is exact, not an approximation.
    # BLUPs are centred and go negative, so column sums cannot say which member
    # carries the merged total. Spread does, and is exact here because the two
    # columns are proportional.
    sk <- stats::sd(d[[k]],  na.rm = TRUE)
    sr <- stats::sd(d[[dr]], na.rm = TRUE)
    factor_applied <- 1
    if (is.finite(sk) && is.finite(sr) && sk > 0 && sr > sk * (1 + 1e-9)) {
      factor_applied <- sr / sk
      set(d, j = k, value = d[[dr]])
    }
    d[, (dr) := NULL]
    if (nzchar(rn) && rn != k && k %in% names(d)) setnames(d, k, rn)
    log[[length(log) + 1]] <- data.table(family = fam, trial = tr,
                                         kept = if (nzchar(rn)) rn else k,
                                         dropped = dr, r = round(r, 6),
                                         merged_total_from = if (factor_applied > 1) dr else k,
                                         factor_applied = round(factor_applied, 4),
                                         reason = dec$reason[i])
  }

  n_after <- length(setdiff(names(d), META))
  counts[[ti]] <- data.table(family = fam, trial = tr, file = targets$file[ti],
                             species_before = n_before, species_after = n_after,
                             removed = n_before - n_after)
  fwrite(d, f)
}

log    <- rbindlist(log)
counts <- rbindlist(counts)
fwrite(log,    file.path(OUT, "SuppTable_SpeciesDeduplication_Decisions.csv"))
fwrite(counts, file.path(OUT, "SuppTable_SpeciesDeduplication_Counts.csv"))

cat("\n-- species removed --\n")
print(log[, .(family, trial, kept, dropped, r, merged_total_from, factor_applied)])
cat("\n-- species counts --\n");  print(counts)
if (length(skipped)) {
  skipped <- rbindlist(skipped)
  fwrite(skipped, file.path(OUT, "SuppTable_SpeciesDeduplication_NotRemoved.csv"))
  cat("\n-- listed as duplicates but NOT removed (correlation below threshold) --\n"); print(skipped)
}
# ---- rebuild the species inventory from the deduplicated matrices ------------
INV <- file.path(REPO, "data/final_species_set/species_inventory.csv")
if (file.exists(INV)) {
  inv <- fread(INV)
  if (!file.exists(file.path(ARCH, "species_inventory.csv")))
    file.copy(INV, file.path(ARCH, "species_inventory.csv"))

  present <- lapply(FILES, function(f) {
    x <- fread(file.path(SPATS, f)); n <- setdiff(names(x), META)
    n[vapply(x[, ..n], is.numeric, logical(1))]
  })

  before <- copy(inv)[, `:=`(In_CTL = Status %in% c("Common", "CTL only"),
                             In_LIN = Status %in% c("Common", "LIN only"))]

  inv[, `:=`(In_CTL = Species %in% present$CTL, In_LIN = Species %in% present$LIN)]
  inv <- inv[In_CTL | In_LIN]
  inv[, Status := fifelse(In_CTL & In_LIN, "Common",
                   fifelse(In_CTL, "CTL only", "LIN only"))]
  fwrite(inv, INV)

  cls <- rbindlist(lapply(sort(unique(before$Class)), function(cl) data.table(
    Class = cl,
    CTL_before = before[Class == cl & In_CTL == TRUE, .N],
    CTL_after  = inv   [Class == cl & In_CTL == TRUE, .N],
    LIN_before = before[Class == cl & In_LIN == TRUE, .N],
    LIN_after  = inv   [Class == cl & In_LIN == TRUE, .N])))
  fwrite(cls, file.path(OUT, "SuppTable_SpeciesDeduplication_ClassCounts.csv"))

  cat("\n-- species inventory --\n")
  cat(sprintf("  total     %3d -> %3d\n", nrow(before), nrow(inv)))
  for (s_ in c("Common", "CTL only", "LIN only"))
    cat(sprintf("  %-8s  %3d -> %3d\n", s_, sum(before$Status == s_), sum(inv$Status == s_)))
  cat("\n-- species per class --\n"); print(cls[CTL_before != CTL_after | LIN_before != LIN_after])
}

cat("\nOriginals archived in", ARCH, "\n")
