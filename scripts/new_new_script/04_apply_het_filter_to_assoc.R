# =============================================================================
# Check, and optionally apply, the high-heterozygosity SNP removal against the
# sum/ratio GWAS results.
#
# Why this is a check first and an edit second. The new run's files are named
# SAP_bialleles_MAF_0.05, with no _filtered_80perc_hets_ token, so it is not
# clear from the filenames whether the filtered or the unfiltered genotypes
# were used. This script answers that from the data. If none of the removed
# SNPs appear in the .assoc.txt files, the GWAS already used the filtered
# genotypes and nothing needs doing. If they do appear, the run used the
# unfiltered genotypes and the SNPs must be dropped post hoc, which is what
# 0_remove_bad_snps.R was written for.
#
# It also reports how many of the removed SNPs are currently Bonferroni
# significant. That is the direct test of whether high-heterozygosity markers
# are producing the inflated signals.
#
# Matching is by SNP id and, independently, by chromosome plus position, since
# the HapMap rs# and the GEMMA rs column do not have to use the same naming.
#
# Usage   Rscript 04_apply_het_filter_to_assoc.R              # check only
#         APPLY=1 Rscript 04_apply_het_filter_to_assoc.R      # write filtered copies
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })

ROOT    <- Sys.getenv("ROOT", "/rsstu/users/r/rrellan/DOE_CAREER/SAP/results/spats_corrected/BLUP/sum_ratio_lipids_final")
REMOVED <- Sys.getenv("REMOVED", file.path(ROOT, "genotype_qc/removed_high_het_SNPs_80perc.tsv"))
APPLY   <- nzchar(Sys.getenv("APPLY", ""))
TRIALS  <- c(CTL = "control/raw_gwas", LIN = "lowinput/raw_gwas")

if (!file.exists(REMOVED)) stop("removal list not found, run 03_het_filter_audit.sh first: ", REMOVED)
bad <- fread(REMOVED, showProgress = FALSE)
bad[, key := paste(chr, pos, sep = ":")]
message("removed-SNP list: ", nrow(bad), " SNPs over ", uniqueN(bad$chr), " chromosomes")

report <- rbindlist(lapply(names(TRIALS), function(tg) {
  dir <- file.path(ROOT, TRIALS[[tg]])
  fs  <- list.files(dir, pattern = "\\.assoc\\.txt$", full.names = TRUE)
  if (!length(fs)) { message("no .assoc.txt in ", dir); return(NULL) }
  message("\n== ", tg, "  ", length(fs), " files")

  # one file is enough to decide which genotype set was used, since every trait
  # in a trial was run on the same markers
  probe <- fread(fs[1], select = c("chr","rs","ps","p_wald"), showProgress = FALSE)
  probe[, key := paste(chr, ps, sep = ":")]
  hit_rs  <- sum(probe$rs  %in% bad$rs)
  hit_pos <- sum(probe$key %in% bad$key)
  message("  markers in assoc      : ", nrow(probe))
  message("  matching removed by id: ", hit_rs)
  message("  matching removed by chr:pos: ", hit_pos)

  use_pos <- hit_pos >= hit_rs
  n_hit   <- max(hit_rs, hit_pos)
  if (n_hit == 0L) {
    message("  -> the filtered genotypes were used, nothing to remove")
  } else {
    message("  -> the UNFILTERED genotypes were used, ", n_hit, " high-het markers present")
    message("     matching on ", if (use_pos) "chr:pos" else "SNP id")
  }

  # Bonferroni denominator is the number of markers GEMMA actually TESTED, which
  # is well below the HapMap total because MAF and missingness are re-evaluated
  # within each trial's genotype subset. Never take it from the HapMap count.
  bonf <- 0.05 / nrow(probe)
  n_after <- nrow(probe) - n_hit
  message("  Bonferroni now         : ", signif(bonf, 3), "  over ", nrow(probe), " markers")
  message("  Bonferroni after strip : ", signif(0.05 / n_after, 3), "  over ", n_after, " markers")
  rbindlist(lapply(fs, function(f) {
    d <- fread(f, showProgress = FALSE)
    d[, key := paste(chr, ps, sep = ":")]
    drop <- if (use_pos) d$key %in% bad$key else d$rs %in% bad$rs
    row <- data.table(
      Trial = tg, Trait = sub("_mod_sub_.*$", "", basename(f)),
      N_markers = nrow(d), N_high_het = sum(drop),
      N_sig_before   = sum(d$p_wald <= bonf, na.rm = TRUE),
      N_sig_high_het = sum(drop & d$p_wald <= bonf, na.rm = TRUE),
      Min_p_high_het = if (any(drop)) min(d$p_wald[drop], na.rm = TRUE) else NA_real_)
    if (APPLY && any(drop)) {
      out <- file.path(ROOT, dirname(TRIALS[[tg]]), "raw_gwas_hetfiltered")
      dir.create(out, recursive = TRUE, showWarnings = FALSE)
      d[, key := NULL]
      fwrite(d[!drop], file.path(out, basename(f)), sep = "\t", quote = FALSE, na = "NA")
    }
    row
  }))
}))

if (is.null(report) || !nrow(report)) stop("no assoc files found under ", ROOT)

cat("\n-- per trial --\n")
print(report[, .(traits = .N,
                 markers = unique(N_markers),
                 high_het_markers = unique(N_high_het),
                 sig_total = sum(N_sig_before),
                 sig_that_are_high_het = sum(N_sig_high_het),
                 pct_of_sig = round(100 * sum(N_sig_high_het) / pmax(sum(N_sig_before), 1), 2)),
              by = Trial])

cat("\n-- traits where high-het markers contribute the most significant hits --\n")
print(head(report[order(-N_sig_high_het),
                  .(Trial, Trait, N_high_het, N_sig_before, N_sig_high_het, Min_p_high_het)], 15),
      row.names = FALSE)

out <- file.path(ROOT, "genotype_qc/het_filter_check_by_trait.csv")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
fwrite(report, out)
cat("\nSaved:", out, "\n")

if (!APPLY && any(report$N_high_het > 0))
  cat("\nHigh-het markers are present. Re-run with  APPLY=1 Rscript 04_apply_het_filter_to_assoc.R\n",
      "to write filtered copies into <trial>/raw_gwas_hetfiltered, leaving the originals untouched.\n")
if (APPLY)
  cat("\nFiltered copies written. Use the reduced marker count as the new Bonferroni denominator.\n")
