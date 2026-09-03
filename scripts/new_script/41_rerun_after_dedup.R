# ==============================================================================
# 41_rerun_after_dedup.R
#
# Rerun every figure and table that reads the SpATS-fitted phenotype matrices,
# after 40_deduplicate_species.R has removed the redundant species from them.
#
#   Rscript scripts/new_script/41_rerun_after_dedup.R
#
# Run from the repository root. Each step is a script that already exists; this
# only fixes the order and reports which ones succeeded, so nothing is silently
# left on stale data.
#
# NOT rerun here, and why.
#
#   GWAS of any kind. No candidate gene is lost by the deduplication, and every
#   class-sum and class-ratio phenotype correlates with its pre-deduplication
#   version at r >= 0.999 in both trials (196 of 324 under CTL and 256 of 324
#   under LIN are affine-identical, so their p-values are unchanged exactly).
#   Rerunning the GWAS is therefore not warranted. Steps 1 and 2 below relabel
#   the phenotype names in the existing candidate tables instead.
#
#   GO enrichment (10, 11) and the overlap figure (13) read GWAS output only,
#   and are unaffected until the sum/ratio GWAS is rerun.
#
# Set SOLD_RERUN_STOP_ON_ERROR=1 to halt at the first failure instead of
# carrying on and reporting at the end.
# ==============================================================================

STOP_ON_ERROR <- nzchar(Sys.getenv("SOLD_RERUN_STOP_ON_ERROR"))

steps <- list(
  # --- candidate tables: relabel duplicate phenotype names ------------------
  list(s = "scripts/new_script/40b_relabel_candidate_phenotypes.R",
       w = "master candidate table, phenotype names relabelled"),
  list(s = "scripts/new_script/38_export_supp_tables_S7toS10.R",
       w = "Supp Tables S7-S10 rebuilt from it"),
  list(s = "scripts/new_new_script/06e_MainTable_top_recurrent_loci.R",
       w = "main table of most recurrent loci"),

  # --- lipidome composition and structure, all read the phenotype matrices ---
  list(s = "scripts/new_new_script/06b_SuppTableS6_species_inventory.R",
       w = "species inventory + Supp Tables S6a-c   (run first, others use the frozen set)"),
  list(s = "scripts/new_new_script/03_Fig1_population_structure.R",
       w = "Figure 1 + Supp Table S25"),
  list(s = "scripts/new_new_script/04_SuppFig3_class_pca.R",
       w = "Supp Fig S3 + Supp Table S25b"),
  list(s = "scripts/new_new_script/05_SuppTableS5G_composition_stability.R",
       w = "Supp Table S5G"),
  list(s = "scripts/new_new_script/06_Fig2_class_composition.R",
       w = "Figure 2 + Supp Tables S5D, S5E"),
  list(s = "scripts/new_new_script/06a_Fig2C_LION_input.R",
       w = "Figure 2C LION input"),
  list(s = "scripts/new_new_script/06c_SuppFig6_species_counts.R",
       w = "Supp Fig S6 species counts"),
  list(s = "scripts/new_new_script/07_SuppFig4_chemical_space.R",
       w = "Supp Fig S4 + Supp Table S5F"),
  list(s = "scripts/new_new_script/08_SuppFig7_class_correlations.R",
       w = "Supp Fig S7"),
  list(s = "scripts/new_new_script/12_Fig5_linex.R",
       w = "Figure 5, LINEX and reaction balance"),

  # --- genomic heritability, now reading the deduplicated matrices -----------
  list(s = "scripts/new_script/22e_heritability_profile_ci.R",
       w = "class-sum heritability with profile-likelihood CIs"),
  list(s = "scripts/new_script/22f_heritability_species_paired.R",
       w = "per-species paired heritability, CTL vs LIN"),
  list(s = "scripts/new_script/22g_Fig_heritability.R",
       w = "heritability figure   (needs 22e and 22f)")
)

run_one <- function(path) {
  if (!file.exists(path)) return(list(ok = FALSE, msg = "script not found"))
  res <- tryCatch({
    system2("Rscript", path, stdout = TRUE, stderr = TRUE)
  }, error = function(e) structure(conditionMessage(e), status = 1L))
  st <- attr(res, "status")
  list(ok = is.null(st) || st == 0L, msg = paste(utils::tail(res, 12), collapse = "\n"))
}

cat(strrep("=", 78), "\n")
cat("Rerunning", length(steps), "scripts after species deduplication\n")
cat(strrep("=", 78), "\n")

out <- vector("list", length(steps))
for (i in seq_along(steps)) {
  s <- steps[[i]]
  cat(sprintf("\n[%2d/%2d] %s\n         %s\n", i, length(steps), basename(s$s), s$w))
  t0 <- Sys.time()
  r  <- run_one(s$s)
  el <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  cat(sprintf("         %s  (%.1fs)\n", if (r$ok) "ok" else "FAILED", el))
  if (!r$ok) {
    cat(paste0("         | ", strsplit(r$msg, "\n")[[1]], collapse = "\n"), "\n")
    if (STOP_ON_ERROR) stop("stopping at ", s$s, call. = FALSE)
  }
  out[[i]] <- data.frame(step = i, script = s$s, produces = s$w,
                         status = if (r$ok) "ok" else "FAILED", seconds = el)
}

out <- do.call(rbind, out)
dir.create("table/supp", recursive = TRUE, showWarnings = FALSE)
write.csv(out, "table/supp/rerun_after_dedup_log.csv", row.names = FALSE)

cat("\n", strrep("=", 78), "\n", sep = "")
print(out[, c("step", "status", "seconds", "produces")], row.names = FALSE)
bad <- out[out$status != "ok", ]
if (nrow(bad)) {
  cat("\n", nrow(bad), "script(s) failed:\n")
  for (p in bad$script) cat("  ", p, "\n")
} else {
  cat("\nAll steps completed.\n")
}
cat("\nNothing further is outstanding. The GWAS itself does not need rerunning;\n")
cat("see the header of this script for the numbers behind that.\n")
