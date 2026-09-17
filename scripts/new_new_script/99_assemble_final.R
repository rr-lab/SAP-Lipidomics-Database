# ==============================================================================
# Assemble final/ -- the one folder the manuscript compiles against.
#
#   Rscript scripts/new_new_script/99_assemble_final.R
#
# Run this last, after the figure and table scripts. It copies the fifteen files
# main.tex includes, the thirty-nine supplementary tables and the ten workbooks
# into final/, and refuses to finish if any of the fifteen is missing.
#
# WHY A SEPARATE FOLDER. fig/ and table/supp/ at the repository root accumulate
# intermediate and superseded output -- figures under two numbering schemes,
# tables under their pre-consolidation numbers, files from runs that used the old
# lipid class annotation. Several of those had names close enough to the current
# ones to be picked up by mistake, and at least one was. final/ contains what the
# manuscript reads and nothing else, so there is no second copy to choose wrongly.
#
# THE FILE LIST IS DELIBERATELY HARD-CODED. It is the list main.tex includes. If
# a figure is renamed in the manuscript it must be renamed here too, and the
# check below will fail loudly until it is, which is the point.
# ==============================================================================
REPO  <- Sys.getenv("SOLD_REPO", ".")
FINAL <- file.path(REPO, "final")

MAIN <- c("Figure1_Population_Structure.png", "Figure2_Class_Composition.png",
          "Figure3_GWAS_Manhattan.png", "Figure5_GO_BP_LDaware.png",
          "Figure5_LINEX.png", "Figure7_Shiny_App.png")

SUPP <- c("Figure7_CTL_LIN_Overlap.png", "SuppFig_S1_Workflow.png",
          "SuppFig_S2_QC_RunOrder_SERRF_PCA_SpATS.png",
          "SuppFig_S3_Class_PCA_Structure.png", "SuppFig_S4_Heritability.png",
          "SuppFig_S5_Lipid_Species_Counts.png", "SuppFig_S6_CLR_Correlations.png",
          "SuppFig_S7_PCA_Lipids.png", "SuppFig_S8_Chemical_Space.png")

for (d in c("fig/main", "fig/supp", "table/workbooks"))
  dir.create(file.path(FINAL, d), recursive = TRUE, showWarnings = FALSE)

copy_in <- function(files, from, to) {
  src <- file.path(REPO, from, files)
  missing <- files[!file.exists(src)]
  ok <- file.copy(src[file.exists(src)], file.path(FINAL, to), overwrite = TRUE)
  list(n = sum(ok), missing = missing)
}

m <- copy_in(MAIN, "fig/main", "fig/main")
s <- copy_in(SUPP, "fig/supp", "fig/supp")

# ONLY THE TEN WORKBOOKS GO IN. The thirty-nine CSV and TSV files under
# table/supp are the sources those workbooks are packed from, and they still
# carry the pre-consolidation numbers -- S5A, S6c, S24b, S31. The manuscript
# cites Supplementary Tables S1 to S10 and names a sheet inside them; it never
# cites a source file. Copying both sets into final/ put two numbering schemes
# side by side in the folder that is supposed to end that confusion.
wb <- list.files(file.path(REPO, "table/supp/workbooks"),
                 pattern = "\\.(xlsx|csv)$", full.names = TRUE)
invisible(file.copy(wb, file.path(FINAL, "table/workbooks"), overwrite = TRUE))

message(sprintf("final/fig/main         %2d of %d", m$n, length(MAIN)))
message(sprintf("final/fig/supp         %2d of %d", s$n, length(SUPP)))
message(sprintf("final/table/workbooks  %2d files", length(wb)))

gone <- c(m$missing, s$missing)
if (length(gone)) {
  message("\nMISSING, and main.tex includes every one of these:")
  for (g in gone) message("  ", g)
  stop(length(gone), " file(s) the manuscript needs are not in fig/. ",
       "Run the script that writes them, or see final/MANIFEST.md for the two ",
       "that need the GWAS result files and cannot be built without them.")
}
message("\nfinal/ is complete -- every file main.tex includes is present.")
