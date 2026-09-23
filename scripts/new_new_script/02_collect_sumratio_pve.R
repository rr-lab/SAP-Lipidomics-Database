# =============================================================================
# Harvest the GEMMA null-model PVE (SNP heritability) for every sum/ratio trait.
#
# GEMMA writes the PVE into the per-trait .log.txt that sits beside the
# .assoc.txt, not into the .assoc.txt itself and not into the Summary folder.
# The lines are
#     ## pve estimate in the null model = 0.xxxx
#     ## se(pve) in the null model      = 0.xxxx
# PVE here is the proportion of phenotypic variance explained by the centred
# kinship matrix, after the fitted covariates, so it is a marker-based h^2 and
# is NOT the same quantity as a broad-sense H^2 from a replicated design.
#
# Run this on the cluster, against the ORIGINAL vcf2gwas tree, because
# 01_collect_sumratio_assoc.sh copies only the .assoc.txt files.
#
# Usage   Rscript 02_collect_sumratio_pve.R
#         SRC=/path/to/Linear_Mixed_Model Rscript 02_collect_sumratio_pve.R
# =============================================================================
SRC  <- Sys.getenv("SRC",  "/share/maize/ntanduk/SoLD/ratios/Output/Linear_Mixed_Model")
DEST <- Sys.getenv("DEST", "/rsstu/users/r/rrellan/DOE_CAREER/SAP/results/spats_corrected/BLUP/sum_ratio_lipids_final")
dir.create(DEST, recursive = TRUE, showWarnings = FALSE)

logs <- list.files(SRC, pattern = "_mod_sub_.*\\.log\\.txt$",
                   full.names = TRUE, recursive = TRUE)
# the run-level logs sit at the top of the tree and have no _mod_sub_, so the
# pattern above already excludes them; drop anything not beside an .assoc.txt
logs <- logs[file.exists(sub("\\.log\\.txt$", ".assoc.txt", logs))]
message("GEMMA per-trait logs found: ", length(logs))
if (!length(logs)) stop("no per-trait .log.txt found under ", SRC)

grab <- function(txt, key) {
  i <- grep(key, txt, fixed = TRUE)
  if (!length(i)) return(NA_real_)
  suppressWarnings(as.numeric(trimws(sub("^.*=", "", txt[i[1]]))))
}

rows <- lapply(logs, function(f) {
  txt <- readLines(f, warn = FALSE)
  b   <- basename(f)
  trial <- if (grepl("Final_control_",  b)) "CTL"
      else if (grepl("Final_lowinput_", b)) "LIN" else NA_character_
  data.frame(
    Trial      = trial,
    Trait      = sub("_mod_sub_.*$", "", b),
    PVE        = grab(txt, "pve estimate in the null model"),
    SE_PVE     = grab(txt, "se(pve) in the null model"),
    VG         = grab(txt, "vg estimate in the null model"),
    VE         = grab(txt, "ve estimate in the null model"),
    N_analyzed = grab(txt, "number of analyzed individuals"),
    N_SNPs     = grab(txt, "number of analyzed SNPs"),
    LogFile    = b,
    stringsAsFactors = FALSE)
})
pve <- do.call(rbind, rows)
pve <- pve[order(pve$Trial, -pve$PVE), ]

# ---- checks -----------------------------------------------------------------
cat("\n-- rows per trial --\n");  print(table(pve$Trial, useNA = "ifany"))
cat("\nPVE missing (GEMMA did not report it):", sum(is.na(pve$PVE)), "\n")
if (any(is.na(pve$Trial))) cat("Trial not resolved for", sum(is.na(pve$Trial)), "files\n")
d <- pve[duplicated(pve[, c("Trial","Trait")]), ]
if (nrow(d)) { cat("\n!! duplicated Trial/Trait rows:\n"); print(d[, c("Trial","Trait","LogFile")]) }

cat("\n-- PVE summary by trial --\n")
print(do.call(rbind, lapply(split(pve, pve$Trial), function(x) data.frame(
  n = nrow(x), median = median(x$PVE, na.rm = TRUE), mean = mean(x$PVE, na.rm = TRUE),
  min = min(x$PVE, na.rm = TRUE), max = max(x$PVE, na.rm = TRUE),
  n_gt_0.1 = sum(x$PVE > 0.10, na.rm = TRUE),
  n_gt_0.2 = sum(x$PVE > 0.20, na.rm = TRUE)))), digits = 3)

cat("\n-- 10 highest PVE per trial --\n")
for (tg in sort(unique(na.omit(pve$Trial))))
  print(head(pve[pve$Trial == tg, c("Trial","Trait","PVE","SE_PVE","N_analyzed")], 10),
        row.names = FALSE, digits = 3)

out <- file.path(DEST, "sumratio_PVE_by_trait.csv")
write.csv(pve, out, row.names = FALSE)
cat("\nSaved:", out, " (", nrow(pve), "rows )\n")
