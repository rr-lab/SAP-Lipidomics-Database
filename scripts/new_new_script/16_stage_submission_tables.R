# ==============================================================================
# Stage the supplementary tables for submission and regenerate the numbering map.
#
#   Rscript scripts/new_new_script/16_stage_submission_tables.R
#
# 15_SuppTables_build_flat.R writes table/supp/tables/ under descriptive names.
# The journal wants one file per table under its own number, so this copies the
# same files to final_submission/supplementary_tables/ as SN_Table.xlsx and to
# final/table/supp_tables/ as the working set, then writes SUPP_TABLE_NUMBERING.md
# with the citation count and first-cited line taken from main.tex itself.
#
# Nothing is rebuilt here. If a table is wrong, fix its generator and re-run
# 15_SuppTables_build_flat.R first.
#
# Input   table/supp/tables/SuppTable_S*.xlsx and SuppTable_index.csv
#         main.tex, for the citation map
# Output  final/table/supp_tables/
#         final_submission/supplementary_tables/
#         SUPP_TABLE_NUMBERING.md
# ==============================================================================
REPO <- Sys.getenv("SOLD_REPO", ".")
SRC  <- file.path(REPO, "table/supp/tables")
MAIN <- file.path(REPO, Sys.getenv("SOLD_MAIN", "main.tex"))
WORK <- file.path(REPO, "final/table/supp_tables")
SUB  <- file.path(REPO, "final_submission/supplementary_tables")

idx <- read.csv(file.path(SRC, "SuppTable_index.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(idx) > 0, all(file.exists(file.path(SRC, idx$file))))

for (d in c(WORK, SUB)) {
  unlink(list.files(d, full.names = TRUE))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---- citation map, read out of main.tex -------------------------------------
tex <- readLines(MAIN, warn = FALSE)
cite_lines <- function(n) {
  # a citation is Supplementary Table~Sn or Supplementary Tables~Sa, Sb and Sn,
  # and Sn must not be the prefix of a longer number
  pat <- sprintf("Supplementary Tables?~(S[0-9]+(, |, and | and ))*S%d(?![0-9])", n)
  which(grepl(pat, tex, perl = TRUE))
}
hits <- lapply(as.integer(sub("^S", "", idx$table)), cite_lines)

idx$line  <- vapply(hits, function(h) if (length(h)) h[1] else NA_integer_, 1L)
idx$cites <- lengths(hits)

# ---- copy --------------------------------------------------------------------
for (i in seq_len(nrow(idx))) {
  from <- file.path(SRC, idx$file[i])
  stopifnot(file.copy(from, file.path(WORK, idx$file[i]), overwrite = TRUE))
  stopifnot(file.copy(from, file.path(SUB, sprintf("%s_Table.xlsx", idx$table[i])),
                      overwrite = TRUE))
}
write.csv(idx, file.path(SUB, "S_Table_index.csv"), row.names = FALSE)
write.csv(idx, file.path(WORK, "SuppTable_index.csv"), row.names = FALSE)

# ---- the numbering document --------------------------------------------------
n <- nrow(idx)
md <- c(
  sprintf("# Supplementary tables, S1 to S%d", n), "",
  sprintf("One table per file, built by `scripts/new_new_script/15_SuppTables_build_flat.R` and staged by `scripts/new_new_script/16_stage_submission_tables.R` on %s.",
          format(Sys.Date())), "",
  "Files live in `table/supp/tables/` under descriptive names, in `final/table/supp_tables/` as the working set,",
  "and in `final_submission/supplementary_tables/` under the journal's `SN_Table.xlsx` names.", "",
  "`was` is the number the table carried in the 41-table scheme that preceded the ten workbooks.",
  "`line` is the first line of `main.tex` that cites it, and `cites` how many lines cite it.", "",
  "| # | Title | was | rows | cols | line | cites |", "|---|---|---|---|---|---|---|",
  sprintf("| %s | %s | %s | %d | %d | %s | %d |", idx$table, gsub("_", " ", idx$title),
          idx$was, idx$rows, idx$cols,
          ifelse(is.na(idx$line), "not cited", idx$line), idx$cites),
  "",
  sprintf("%d tables, %s data rows in total.", n, format(sum(idx$rows), big.mark = ",")))
writeLines(md, file.path(REPO, "SUPP_TABLE_NUMBERING.md"))

message(sprintf("staged %d tables to final/table/supp_tables/ and final_submission/supplementary_tables/", n))
if (any(idx$cites == 0))
  message("NOT CITED in main.tex: ", paste(idx$table[idx$cites == 0], collapse = ", "))

# The journal numbers supporting items in order of first citation, so the table
# whose first mention is earliest must be S1. Checking it here means a later edit
# that moves a section, as moving Materials and methods above Results did on
# 2026-10-01, cannot leave the numbering silently out of order.
ord_ok <- order(idx$line, na.last = TRUE)
if (!identical(ord_ok, seq_len(nrow(idx)))) {
  message("\nNUMBERING IS NOT IN CITATION ORDER. By first mention it should run")
  message("  ", paste(idx$table[ord_ok], collapse = " "))
  stop("renumber the TABLES list in 15_SuppTables_build_flat.R, then rebuild")
}
message("numbering is in order of first citation, S1 to S", nrow(idx))
print(idx[, c("table", "title", "rows", "line", "cites")], row.names = FALSE)
