# ==============================================================================
# Supplementary Tables S1-S10 -- pack the per-analysis CSV/TSV files into ten
# numbered workbooks, one sheet per file.
#
#   Rscript scripts/new_new_script/15_SuppTables_build_workbooks.R
#
# WHY THIS EXISTS. The supplement had grown to 41 separately numbered tables,
# most of them tiny -- six rows here, three rows there -- which is a lot of
# numbers for about ten topics. Grouping them by topic gives ten workbooks whose
# sheets carry the detail, so nothing is lost and the numbering is readable.
# Letter suffixes (S5A, S6a, S24b) are gone; every table is a plain number.
#
# THE ORDER IS FIRST MENTION IN main.tex, not the order the analyses were run.
# S1 is the population-structure group because it is cited first, at line 134;
# S10 is the frozen species set because it is not cited until the Methods. If a
# citation moves in the manuscript, the number that table should carry moves
# with it, and this file is where that gets re-decided.
#
# Sheet names are capped at 31 characters by the xlsx format, so they are short
# labels rather than the source filenames. The mapping from sheet back to source
# file is the GROUPS list below and is also written to the index sheet.
#
# Input   table/supp/SuppTable_S*.csv, table/supp/SuppTable_S*.tsv  (41 files)
# Output  table/supp/workbooks/SuppTable_S1..S10_*.xlsx
#         table/supp/workbooks/SuppTable_index.csv
# ==============================================================================
suppressPackageStartupMessages({ library(openxlsx); library(vroom) })

REPO <- Sys.getenv("SOLD_REPO", ".")
SRC  <- file.path(REPO, "table/supp")
OUT  <- file.path(SRC, "workbooks")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# new number -> title, and the sheets in the order they should appear.
# Each sheet is  "Sheet name" = "old table number".
GROUPS <- list(
  list(n = 1, title = "Population_structure_and_ancestry", sheets = c(
    "Lipid tests"                  = "S24",
    "Group composition"            = "S24a",
    "PC tests"                     = "S24b")),
  list(n = 2, title = "Genomic_heritability", sheets = c(
    "Per-species paired"           = "S25",
    "Class sums"                   = "S26",
    "Structure-conditioned"        = "S27",
    "Per-species structure and h2" = "S28",
    "Ancestry robustness"          = "S29")),
  list(n = 3, title = "Species_inventory_and_stability", sheets = c(
    "Species summary"              = "S6a",
    "By class"                     = "S6b",
    "By superclass"                = "S6c",
    "Composition stability"        = "S5G")),
  list(n = 4, title = "Class_composition_and_contrasts", sheets = c(
    "Composition pctTIC"           = "S5D",
    "CLR contrast"                 = "S5A",
    "ALR contrast"                 = "S5B",
    "CLR correlation delta"        = "S5C",
    "Ratio statistics"             = "S1",
    "Top-variance lipids"          = "S4")),
  list(n = 5, title = "Lipid_ontology_and_chemical_space", sheets = c(
    "LION enrichment"              = "S5E",
    "Chemical space"               = "S5F")),
  list(n = 6, title = "GWAS_candidate_genes", sheets = c(
    "CTL individual"               = "S7",
    "CTL sum-ratio"                = "S8",
    "LIN individual"               = "S9",
    "LIN sum-ratio"                = "S10")),
  list(n = 7, title = "GO_enrichment", sheets = c(
    "Loci collapsed"               = "S16",
    "BP all terms"                 = "S17",
    "MF all terms"                 = "S18",
    "Genes in enriched terms"      = "S19")),
  list(n = 8, title = "CTL_LIN_candidate_overlap", sheets = c(
    "Gene level"                   = "S20",
    "Locus level"                  = "S21",
    "Shared individual"            = "S22",
    "Shared sum-ratio"             = "S23",
    "By lipid class"               = "S30",
    "Shared ranked"                = "S31")),
  list(n = 9, title = "LINEX_reaction_mapping", sheets = c(
    "Reactions"                    = "S11",
    "Reaction balance"             = "S12",
    "GWAS gene support"            = "S13",
    "Branch summary"               = "S14")),
  list(n = 10, title = "Frozen_lipid_species_set", sheets = c(
    "Final lipid classes"          = "S15"))
)

# ---- locate each source file by its old number -------------------------------
avail <- list.files(SRC, pattern = "^SuppTable_S[0-9]+[a-zA-Z]?_.*\\.(csv|tsv)$")
avail <- avail[!grepl("^SuppTable_S[0-9]+to", avail)]
key   <- sub("^SuppTable_(S[0-9]+[a-zA-Z]?)_.*$", "\\1", avail)
stopifnot(!anyDuplicated(key))
path_of <- setNames(file.path(SRC, avail), key)

# DROPPED FROM THE SUPPLEMENT, 2026-09-16.
#
# S2 and S3 were the per-species CTL/LIN Wilcoxon contrast and its per-class
# summary, each carrying a leave-one-genotype-out jackknife stability column.
# They are out for three reasons. The manuscript reports no per-species
# between-trial test -- every trial contrast it makes is at class level -- so
# they backed a claim the paper does not make. The sentence that cited them
# described them as a per-species DETECTION result, which is not what they
# measure. And they predate the species deduplication, testing 152 species on
# 394 CTL and 363 LIN genotypes against the current 146 on 389 and 357.
#
# The files move to table/archive/ rather than being deleted, so the analysis
# can be recovered if a reviewer asks for it. Listing them here keeps the
# unassigned-file guard below honest instead of silently widening it.
DROPPED <- c("S2", "S3")

wanted <- unlist(lapply(GROUPS, function(g) unname(g$sheets)))
stopifnot(length(wanted) == length(unique(wanted)))
missing <- setdiff(wanted, names(path_of))
orphan  <- setdiff(names(path_of), c(wanted, DROPPED))
if (length(missing)) stop("no file for old table(s) ", paste(missing, collapse = ", "))
if (length(orphan))  stop("file(s) not assigned to any group ", paste(orphan, collapse = ", "))
message(sprintf("%d source files, all assigned", length(wanted)))

# ---- build ------------------------------------------------------------------
index <- list()
for (g in GROUPS) {
  wb <- createWorkbook()
  for (i in seq_along(g$sheets)) {
    nm  <- names(g$sheets)[i]; old <- unname(g$sheets)[i]
    stopifnot(nchar(nm) <= 31)
    f   <- path_of[[old]]
    dat <- vroom(f, delim = if (grepl("\\.tsv$", f)) "\t" else ",",
                 show_col_types = FALSE, progress = FALSE)
    addWorksheet(wb, nm)
    writeData(wb, nm, dat)
    freezePane(wb, nm, firstRow = TRUE)
    index[[length(index) + 1]] <- data.frame(
      table = sprintf("S%d", g$n), title = g$title, sheet = nm,
      was = old, rows = nrow(dat), cols = ncol(dat), source_file = basename(f))
  }
  f_out <- file.path(OUT, sprintf("SuppTable_S%d_%s.xlsx", g$n, g$title))
  saveWorkbook(wb, f_out, overwrite = TRUE)
  message(sprintf("S%-3d %-36s %d sheets", g$n, g$title, length(g$sheets)))
}
idx <- do.call(rbind, index)
write.csv(idx, file.path(OUT, "SuppTable_index.csv"), row.names = FALSE)
message(sprintf("\n%d workbooks, %d sheets, %d data rows",
                length(GROUPS), nrow(idx), sum(idx$rows)))
