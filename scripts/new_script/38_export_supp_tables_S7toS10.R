# ==============================================================================
# 38_export_supp_tables_S7toS10.R
#
# Rebuild Supplementary Tables S7-S10 (the GWAS candidate-gene tables) directly
# from the LD-mapping master table, so that every number quoted in the manuscript
# comes from the same run as the tables the reader downloads.
#
#   S7  CTL, individual lipid traits
#   S8  CTL, lipid-class sums and ratios
#   S9  LIN, individual lipid traits
#   S10 LIN, lipid-class sums and ratios
#
# The previously shipped S7-S10 were produced by an earlier GWAS run over a
# larger phenotype set and at a suggestive significance cutoff, so their gene
# counts and phenotype-recurrence maxima did not match the manuscript text.
# This script replaces them.
#
# Phenotype relabelling. 40_deduplicate_species.R removed species that were
# duplicate annotations of a retained species. An association reported against a
# removed name is the same association against the name that was kept, since the
# two columns were the same measurement. Those phenotype names are therefore
# relabelled to their retained partner and de-duplicated within each gene, rather
# than dropped. Best_P_Value is unaffected: the relabelled phenotype carries the
# same p-value it always did. N_Phenotypes is recounted after relabelling.
#
# Input   data/LD_mapped/candidate_tables/ALL_LD_candidate_genes_master.tsv
#         table/supp/SuppTable_SpeciesDeduplication_Decisions.csv  (optional)
# Output  table/supp/SuppTable_S7..S10_*.tsv          (canonical, feeds the workbook)
#         new_new_tables/supp/SuppTable_S7..S10_*.tsv (mirror)
#
# Run from the repository root.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table) })
`%||%` <- function(a, b) if (is.null(a)) b else a

# "; " (semicolon SPACE) is the Phenotypes delimiter. Splitting on ";" alone
# would tear the species "SPB 18:0;2OH" in half. See 40b for the full note.
SEP <- "; "

REPO      <- Sys.getenv("SOLD_DB_REPO", ".")
MASTER    <- file.path(REPO, "data/LD_mapped/candidate_tables/ALL_LD_candidate_genes_master.tsv")
OUT_MAIN  <- file.path(REPO, "table/supp")
OUT_MIRROR<- file.path(REPO, "new_new_tables/supp")
ARCHIVE   <- file.path(REPO, "table/supp/_superseded_JUN05_run")

BONF <- 8.19e-9   # 0.05 / 6,105,000 markers; -log10 = 8.09

for (d in c(OUT_MAIN, OUT_MIRROR, ARCHIVE)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

master <- fread(MASTER, sep = "\t", quote = "", colClasses = "character")
# read as character so p-values and r2 keep the exact text of the master table;
# numeric copies are made only for sorting and for the summary below.

# ---- column hygiene ----------------------------------------------------------
# Family_Sufamily is a typo carried through the annotation pipeline.
setnames(master, "Family_Sufamily", "Family_SubFamily", skip_absent = TRUE)
setnames(master, "ProteinClass",    "Protein_class",    skip_absent = TRUE)

KEEP <- c("GeneID", "GeneName", "Family_SubFamily", "Protein_class",
          "GO_MF", "GO_BP", "GO_CC",
          "Chromosome", "Gene_Start", "Gene_End",
          "Best_P_Value", "N_sig_SNPs", "Max_r2",
          "N_Phenotypes", "Phenotypes")
stopifnot(all(KEEP %in% names(master)))

# ---- the four tables ---------------------------------------------------------
spec <- data.table(
  table_id  = c("S7", "S8", "S9", "S10"),
  condition = c("CTL", "CTL", "LIN", "LIN"),
  layer     = c("individual", "sumratio", "individual", "sumratio"),
  fname     = c("SuppTable_S7_CTL_GWAS_candidate_genes_for_individual_lipid_traits.tsv",
                "SuppTable_S8_CTL_GWAS_candidate_genes_for_lipid_sum_ratio_traits.tsv",
                "SuppTable_S9_LIN_GWAS_candidate_genes_for_individual_lipid_traits.tsv",
                "SuppTable_S10_LIN_GWAS_candidate_genes_for_lipid_sum_ratio_traits.tsv"),
  mirror    = c("SuppTable_S7_CTL_GWAS_candidates_ind.tsv",
                "SuppTable_S8_CTL_GWAS_candidates_sumratio.tsv",
                "SuppTable_S9_LIN_GWAS_candidates_ind.tsv",
                "SuppTable_S10_LIN_GWAS_candidate_sumratio.tsv")
)

# ---- phenotype relabelling map, from the deduplication decision table --------
DEC <- file.path(OUT_MAIN, "SuppTable_SpeciesDeduplication_Decisions.csv")
relabel <- list()
if (file.exists(DEC)) {
  dd <- unique(fread(DEC)[, .(trial, kept, dropped)])
  # candidate tables write acyl separators as "_", the species tables as "/"
  u <- function(x) gsub("/", "_", x, fixed = TRUE)
  for (tr in unique(dd$trial))
    relabel[[tr]] <- setNames(u(dd[trial == tr]$kept), u(dd[trial == tr]$dropped))
}

apply_relabel <- function(tab, tr) {
  m <- relabel[[tr]]
  if (is.null(m) || !length(m)) return(tab)
  fixed <- vapply(tab$Phenotypes, function(ph) {
    v <- trimws(strsplit(ph, SEP, fixed = TRUE)[[1]]); v <- v[nzchar(v)]
    hit <- v %in% names(m); if (any(hit)) v[hit] <- m[v[hit]]
    paste(unique(v), collapse = SEP)
  }, character(1), USE.NAMES = FALSE)
  n <- vapply(fixed, function(ph) length(strsplit(ph, SEP, fixed = TRUE)[[1]]),
              integer(1), USE.NAMES = FALSE)
  changed <- sum(fixed != tab$Phenotypes)
  tab[, Phenotypes := fixed][, N_Phenotypes := as.character(n)]
  attr(tab, "relabelled") <- changed
  tab
}

summ <- vector("list", nrow(spec))
relabel_log <- integer(nrow(spec))

for (i in seq_len(nrow(spec))) {
  s   <- spec[i]
  tab <- master[condition == s$condition & layer == s$layer, ..KEEP]

  # sort strongest signal first, then by genome position
  tab[, `:=`(..p = as.numeric(Best_P_Value), ..c = as.integer(Chromosome), ..s = as.numeric(Gene_Start))]
  setorder(tab, ..p, ..c, ..s)
  tab[, c("..p", "..c", "..s") := NULL]

  tab <- apply_relabel(tab, s$condition)
  relabel_log[i] <- as.integer(attr(tab, "relabelled") %||% 0L)

  # archive whatever is currently shipped under that name before overwriting
  for (old in c(file.path(OUT_MAIN, s$fname), file.path(OUT_MIRROR, s$mirror))) {
    if (file.exists(old)) file.rename(
      old, file.path(ARCHIVE, paste0(basename(dirname(old)), "__", basename(old))))
  }

  fwrite(tab, file.path(OUT_MAIN,   s$fname),  sep = "\t", quote = FALSE)
  fwrite(tab, file.path(OUT_MIRROR, s$mirror), sep = "\t", quote = FALSE)

  phen <- unique(trimws(unlist(strsplit(tab$Phenotypes, SEP, fixed = TRUE))))
  phen <- phen[nzchar(phen)]

  summ[[i]] <- data.table(
    Table              = s$table_id,
    Trial              = s$condition,
    Layer              = s$layer,
    Genes              = nrow(tab),
    Phenotypes_tested  = length(phen),
    Max_phenotypes_per_gene = max(as.integer(tab$N_Phenotypes)),
    Min_P              = min(as.numeric(tab$Best_P_Value)),
    Rows_above_Bonferroni   = sum(as.numeric(tab$Best_P_Value) > BONF),
    Genes_with_relabelled_phenotypes = relabel_log[i]
  )
}

summ <- rbindlist(summ)
fwrite(summ, file.path(OUT_MAIN, "SuppTable_S7toS10_summary.tsv"), sep = "\t", quote = FALSE)

cat("\nRegenerated S7-S10 from", basename(MASTER), "\n\n")
print(summ)
cat("\nSuperseded files moved to", ARCHIVE, "\n")
