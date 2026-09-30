# ==============================================================================
# Figure 2C -- LION ranking-mode input
#
#   source("scripts/new_new_script/_common.R")   # via the header below
#   Rscript scripts/new_new_script/06a_Fig2C_LION_input.R
#
# Writes table/Linex2/LION_input.csv in the two-header-row layout the LION web
# tool expects in Ranking mode (lipidontology.com -> Ranking mode -> (i) process
# input -> file input). LION computes the CTL vs LIN ranking itself; this script
# only supplies the matrix.
#
# Layout, matching LION's own example file
#   row 1   empty cell, then the group label of every sample column
#   row 2   empty cell, then a sample id
#   rows 3+ lipid name, then that lipid's value in every sample
#
# Choices made here, all of which belong in Methods
#   * values are per-sample %TIC, not raw intensity. The two trials are separate
#     acquisition batches, so raw totals differ for instrument reasons; %TIC puts
#     both on the compositional scale the rest of the paper reports.
#   * only species detected in both trials are included (171 of 243). A species
#     seen in one trial is below the detection threshold in the other, not known
#     to be absent, and feeding LION a near-infinite fold change for those would
#     drive the ranking on a detection artefact.
#   * features whose names normalise to the same species are summed, so the
#     lyso species are not double counted (see normalize_lipid_name in _common.R).
#
# Known limitation, worth stating in the paper: LION's ontology is built on
# LIPID MAPS and does not carry the plant galactolipids or the sulfolipid, so
# MGDG, DGDG and SQDG (22 of the 171 species) map to no term and are absent from
# every LION result. The classes carrying the strongest compositional signal are
# therefore outside the enrichment universe.
# ==============================================================================
source("scripts/new_new_script/_common.R")

OUT <- file.path(REPO, "table/Linex2/LION_input.csv")
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)

# per-sample %TIC, with same-named features summed
trial_matrix <- function(path) {
  x <- read_trial(path)
  m <- as.matrix(x[, -1, drop = FALSE]); storage.mode(m) <- "numeric"
  m[!is.finite(m)] <- 0
  s <- sweep(m, 1, pmax(rowSums(m), 1e-12), "/") * 100
  s <- t(rowsum(t(s), normalize_lipid_name(colnames(m))))
  rownames(s) <- as.character(x[[1]])
  s
}

A <- trial_matrix(CTL_CSV)
B <- trial_matrix(LIN_CSV)

keep <- intersect(colnames(A), colnames(B))
A <- A[, keep, drop = FALSE]
B <- B[, keep, drop = FALSE]

vals   <- cbind(t(A), t(B))
groups <- c(rep("CTL", nrow(A)), rep("LIN", nrow(B)))
ids    <- paste0("#", seq_len(ncol(vals)))

out <- rbind(c("", groups),
             c("", ids),
             cbind(rownames(vals),
                   format(vals, scientific = TRUE, digits = 8, trim = TRUE)))

write.table(out, OUT, sep = ",", row.names = FALSE, col.names = FALSE,
            quote = FALSE, na = "")

message("Saved: ", OUT, "  (", length(keep), " species x ",
        ncol(vals), " samples; CTL ", nrow(A), ", LIN ", nrow(B), ")")

# what LION will and will not be able to map
cls <- lipid_class(keep)
unmappable <- c("MGDG", "DGDG", "SQDG")
cat("\n-- species by class in the input --\n")
print(sort(table(cls), decreasing = TRUE))
cat("\nOutside LION's ontology (", sum(cls %in% unmappable), " species): ",
    paste(unmappable, collapse = ", "), "\n", sep = "")
