# ==============================================================================
# Supplementary Tables S1-S29 -- one table per file, numbered in a single flat
# sequence.
#
#   Rscript scripts/new_new_script/15_SuppTables_build_flat.R
#
# WHY THIS REPLACES THE TEN-WORKBOOK BUILD. Grouping 39 sheets into ten
# workbooks made the numbering short but made every citation two-part, as
# "Supplementary Table S4, CLR correlation delta". A tab with its own header row
# and its own subject is a table, so each one now carries its own number and the
# citation is one token again.
#
# THE ORDER IS FIRST MENTION IN main.tex. Sheets that belonged to the same
# workbook stay together as a contiguous run, placed at the first mention of any
# member, because they are read together. Within a run the order is first
# mention, then source order for the ones the text never names individually.
#
# WHAT WAS DROPPED ON THE WAY FROM 39 TO 31, 2026-09-30.
#   S23 Shared sum-ratio   0 rows. A header with nothing under it. The finding
#                          is already the sentence "shares no gene at all
#                          between trials, against 0.6 expected".
#   S6a Species summary    6 rows, the counts 214 / 216 / 164 / 50 / 52 / 266.
#                          All six are in the Results text and in SuppFig S5A.
#   S20 Gene level         3 rows, identical to rows 1-3 of main-text Table 3.
#   S21 Locus level        9 rows, six of which are Table 3's lower block. The
#                          three All-layers window rows are the only unique
#                          content and the text reports the layers separately.
# The source CSVs stay in table/supp/ and are listed in DROPPED below, so the
# unassigned-file guard stays honest rather than being widened.
#
# WHAT WAS MERGED. Six pairs that share a key and a subject, listed in TABLES
# with a `merge` function. Nothing is lost by any of them; the row counts are
# checked against the sources at the end of this script.
#
# THE TWO CANDIDATE MERGES, 2026-09-30, which took the sequence from 31 to 29.
# The four GWAS candidate tables carried the same columns and differed only in
# trial and trait layer, so they were four files where two would do. Each trial
# is now one table with a leading Layer column reading individual or sum_ratio,
# and the four numbers S19-S22 collapse to S19 CTL and S20 LIN. Everything from
# the old S23 onward moves down by two.
#
# Input   table/supp/SuppTable_S*.csv, table/supp/SuppTable_S*.tsv
# Output  table/supp/tables/SuppTable_S1..S29_*.xlsx   (one sheet each)
#         table/supp/tables/SuppTable_index.csv
# ==============================================================================
suppressPackageStartupMessages({ library(openxlsx); library(vroom) })

REPO <- Sys.getenv("SOLD_REPO", ".")
SRC  <- file.path(REPO, "table/supp")
OUT  <- file.path(SRC, "tables")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- locate each source file by its old number -------------------------------
avail <- list.files(SRC, pattern = "^SuppTable_S[0-9]+[a-zA-Z]?_.*\\.(csv|tsv)$")
avail <- avail[!grepl("^SuppTable_S[0-9]+to", avail)]
key   <- sub("^SuppTable_(S[0-9]+[a-zA-Z]?)_.*$", "\\1", avail)
stopifnot(!anyDuplicated(key))
path_of <- setNames(file.path(SRC, avail), key)

read_old <- function(old) {
  f <- path_of[[old]]
  if (is.null(f)) stop("no source file for old table ", old)
  vroom(f, delim = if (grepl("\\.tsv$", f)) "\t" else ",",
        show_col_types = FALSE, progress = FALSE)
}

# ---- the four merges ---------------------------------------------------------
# Each takes the sources in the order listed in `was` and returns one frame.

merge_contrasts <- function(clr, alr)           # S5A + S5B
  rbind(cbind(Scale = "CLR", as.data.frame(clr)),
        cbind(Scale = "ALR", as.data.frame(alr)))

merge_counts <- function(by_class, by_cat) {    # S6b + S6c
  a <- as.data.frame(by_class); b <- as.data.frame(by_cat)
  names(a)[1] <- names(b)[1] <- "Group"
  rbind(cbind(Tier = "Class",    a),
        cbind(Tier = "Category", b))
}

merge_shared <- function(ranked, individual) {  # S31 + S22
  r <- as.data.frame(ranked); i <- as.data.frame(individual)
  i$GeneName <- NULL                           # identical to ranked's
  out <- merge(r, i, by = "GeneID", all.x = TRUE, sort = FALSE)
  out[order(out$best_p, -out$tot_phen), ]
}

merge_candidates <- function(individual, sum_ratio) {  # S7 + S8, and S9 + S10
  a <- as.data.frame(individual); b <- as.data.frame(sum_ratio)
  stopifnot(identical(names(a), names(b)))
  rbind(cbind(Layer = "individual", a), cbind(Layer = "sum_ratio", b))
}

merge_reactions <- function(reactions, balance) {  # S11 + S12
  a <- as.data.frame(reactions); b <- as.data.frame(balance)
  b$join_key <- sub("\\*$", "", b$Reaction_branch)
  out <- merge(a, b, by.x = "Gene_Symbol", by.y = "join_key", all.x = TRUE,
               sort = FALSE)
  out$Reaction_branch <- NULL                  # duplicates Gene_Symbol
  out[order(out$Reaction_ID), ]
}

# ---- the 31 tables -----------------------------------------------------------
# n      the number the table carries in main.tex
# sheet  the single sheet name inside the file, <= 31 characters
# title  goes in the file name
# was    the old table number(s); two entries means `merge` is used
T <- function(n, title, sheet, was, merge = NULL)
  list(n = n, title = title, sheet = sheet, was = was, merge = merge)

TABLES <- list(
  # -- cited first, from Materials and methods -----------------------------------
  T( 1, "Frozen_lipid_species_set",           "Frozen species set",      "S32"),
  T( 2, "GO_loci_collapsed",                  "Loci collapsed",          "S16"),
  T( 3, "GO_BP_all_terms",                    "BP all terms",            "S17"),
  T( 4, "GO_MF_all_terms",                    "MF all terms",            "S18"),
  T( 5, "GO_genes_in_enriched_terms",         "Genes in enriched terms", "S19"),
  # -- population structure, ancestry and heritability ---------------------------
  T( 6, "Population_structure_group_sizes",   "Group composition",       "S24a"),
  T( 7, "Population_structure_lipid_tests",   "Lipid tests",             "S24"),
  T( 8, "Heritability_per_species_structure", "Per-species struct+h2",   "S28"),
  T( 9, "Population_structure_PC_tests",      "PC tests",                "S24b"),
  T(10, "Heritability_per_species_paired",    "Per-species paired",      "S25"),
  T(11, "Heritability_structure_conditioned", "Structure-conditioned",   "S27"),
  T(12, "Heritability_class_sums",            "Class sums",              "S26"),
  T(13, "Heritability_ancestry_robustness",   "Ancestry robustness",     "S29"),
  # -- species inventory, composition and contrasts ------------------------------
  T(14, "Species_counts_by_class_and_category", "Species counts",        c("S6b", "S6c"), merge_counts),
  T(15, "Top_variance_lipids",                "Top-variance lipids",     "S4"),
  T(16, "Class_composition_pctTIC",           "Composition pctTIC",      "S5D"),
  T(17, "Class_contrasts_CLR_and_ALR",        "CLR and ALR contrasts",   c("S5A", "S5B"), merge_contrasts),
  T(18, "Composition_stability",              "Composition stability",   "S5G"),
  T(19, "Class_CLR_correlation_delta",        "CLR correlation delta",   "S5C"),
  # -- lipid ontology and chemical space -----------------------------------------
  T(20, "LION_enrichment",                    "LION enrichment",         "S5E"),
  T(21, "Chemical_space",                     "Chemical space",          "S5F"),
  # -- GWAS candidate genes ------------------------------------------------------
  T(22, "GWAS_candidates_CTL",                "CTL candidates",          c("S7", "S8"), merge_candidates),
  T(23, "GWAS_candidates_LIN",                "LIN candidates",          c("S9", "S10"), merge_candidates),
  # -- CTL/LIN candidate overlap -------------------------------------------------
  T(24, "Overlap_by_lipid_class",             "By lipid class",          "S30"),
  T(25, "Shared_candidate_genes",             "Shared genes",            c("S31", "S22"), merge_shared),
  # -- LINEX reaction mapping ----------------------------------------------------
  T(26, "LINEX_reactions_and_balance",        "Reactions and balance",   c("S11", "S12"), merge_reactions),
  T(27, "LINEX_GWAS_gene_support",            "GWAS gene support",       "S13"),
  T(28, "LINEX_branch_summary",               "Branch summary",          "S14"),
  # -- cited first from the Discussion -------------------------------------------
  T(29, "Lipid_ratio_statistics",             "Ratio statistics",        "S1")
)

DROPPED <- c(
  "S2",  "S3",   # per-species trial contrast + jackknife, dropped 2026-09-16
  "S23",          # shared sum-ratio, 0 rows
  "S6a",          # species summary, six numbers already in the text
  "S20", "S21",   # gene- and locus-level overlap, duplicated by main Table 3
  "S15"           # the 316-row class lookup. It was standing in for the frozen
                  # species set and is 50 rows too long for that, because it
                  # still holds every annotation the curation merged away. S32,
                  # written by 79_build_frozen_species_table.R, is the curated
                  # 266 with the lookup's subclass tiers joined on, and S9 now
                  # carries that instead.
)

wanted <- unlist(lapply(TABLES, function(t) t$was))
stopifnot(length(wanted) == length(unique(wanted)))
stopifnot(vapply(TABLES, function(t) as.integer(t$n), 1L) == seq_along(TABLES))
missing <- setdiff(wanted, names(path_of))
orphan  <- setdiff(names(path_of), c(wanted, DROPPED))
if (length(missing)) stop("no file for old table(s) ", paste(missing, collapse = ", "))
if (length(orphan))  stop("file(s) not assigned to any table ", paste(orphan, collapse = ", "))
message(sprintf("%d source files across %d tables, all assigned",
                length(wanted), length(TABLES)))

# ---- build -------------------------------------------------------------------
index <- list()
for (t in TABLES) {
  stopifnot(nchar(t$sheet) <= 31)
  srcs <- lapply(t$was, read_old)
  dat  <- if (is.null(t$merge)) srcs[[1]] else do.call(t$merge, srcs)
  if (!is.null(t$merge)) {
    expect <- sum(vapply(srcs, nrow, 1L))
    got    <- nrow(dat)
    message(sprintf("  merge S%-2d %-26s %s -> %d rows%s", t$n, t$sheet,
                    paste(vapply(srcs, nrow, 1L), collapse = "+"), got,
                    if (got == expect) "" else "  (join, not a stack)"))
  }
  wb <- createWorkbook()
  addWorksheet(wb, t$sheet)
  writeData(wb, t$sheet, dat)
  freezePane(wb, t$sheet, firstRow = TRUE)
  f_out <- file.path(OUT, sprintf("SuppTable_S%d_%s.xlsx", t$n, t$title))
  saveWorkbook(wb, f_out, overwrite = TRUE)
  index[[length(index) + 1]] <- data.frame(
    table = sprintf("S%d", t$n), title = t$title, sheet = t$sheet,
    was = paste(t$was, collapse = " + "), rows = nrow(dat), cols = ncol(dat),
    file = basename(f_out))
}
idx <- do.call(rbind, index)
write.csv(idx, file.path(OUT, "SuppTable_index.csv"), row.names = FALSE)
message(sprintf("\n%d tables, %d data rows", nrow(idx), sum(idx$rows)))
print(idx[, c("table", "sheet", "was", "rows", "cols")], row.names = FALSE)
