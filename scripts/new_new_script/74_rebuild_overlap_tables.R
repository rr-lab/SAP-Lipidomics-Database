# ==============================================================================
# The CTL/LIN candidate overlap -- rebuild the master table and Supplementary
# Tables S20 and S21.
#
#   Rscript scripts/new_new_script/74_rebuild_overlap_tables.R
#
# Why this exists. Every overlap output on disk descends from
# data/LD_mapped/candidate_tables/ALL_LD_candidate_genes_master.tsv, written
# 2026-09-16, which predates the high-heterozygosity marker filter and the class
# sum and ratio rebuild. It carries 1,100 / 115 / 4,319 / 812 candidate genes
# where the current annotation tables give 1,062 / 54 / 4,370 / 385. S20 and S21
# were never written by an R script at all -- they came from the first-pass
# Python and have no generator in the repository.
#
# This script rebuilds the master from data/gene_annotation_final at
# Max_r2 >= 0.4, the same source and cutoff as Supplementary Tables S7-S10, and
# then writes S20 and S21 from it. Running
# 14_SuppTableS30S31_overlap_by_class_and_shared_ranked.R afterwards regenerates
# S22, S23, S30 and S31 against the rebuilt master with its class vocabulary
# untouched.
#
# The gene universe is every gene with coordinates, 34,027, which is what the
# manuscript quotes. Overlap is a right-tailed hypergeometric test. Windows are
# consecutive non-overlapping bins on gene start, and n_windows_genome counts the
# bins holding at least one gene model, so the null is over windows that could
# have been hit rather than over the whole genome length.
#
# Inputs
#   data/gene_annotation_final/{CTL,LIN}_gene_annotation.tsv
#   data/gene_annotation_final/GO_terms_all.txt
#   data/gene_annotation_final/sorghum_all_genes_with_coords.tsv
# Outputs
#   data/LD_mapped/candidate_tables/ALL_LD_candidate_genes_master.tsv
#   table/supp/SuppTable_S20_Overlap_gene_level.csv
#   table/supp/SuppTable_S21_Overlap_locus_level.csv
# ==============================================================================
suppressPackageStartupMessages(library(data.table))

REPO   <- Sys.getenv("SOLD_REPO", ".")
ANN    <- file.path(REPO, "data/gene_annotation_final")
MASTER <- file.path(REPO, "data/LD_mapped/candidate_tables/ALL_LD_candidate_genes_master.tsv")
SUPP   <- file.path(REPO, "table/supp")
R2_CUT <- 0.4
WINDOWS <- c(100, 250, 500)
dir.create(dirname(MASTER), recursive = TRUE, showWarnings = FALSE)
dir.create(SUPP, recursive = TRUE, showWarnings = FALSE)

A <- rbindlist(lapply(c("CTL_gene_annotation.tsv", "LIN_gene_annotation.tsv"),
                      function(f) fread(file.path(ANN, f), sep = "\t", colClasses = "character")))
A[, `:=`(Best_P = as.numeric(Best_P), Max_r2 = as.numeric(Max_r2),
         N_sig_SNPs = as.integer(N_sig_SNPs),
         Gene_Start = as.integer(Gene_Start), Gene_End = as.integer(Gene_End))]
A <- A[Max_r2 >= R2_CUT]
A[, layer := ifelse(grepl("individual", Layer), "individual", "sumratio")]

func <- unique(fread(file.path(ANN, "GO_terms_all.txt"), sep = "\t", colClasses = "character"),
               by = "GeneID")
setnames(func, c("Family_Subfamily", "Protein_Class"), c("Family_Sufamily", "ProteinClass"),
         skip_absent = TRUE)
coords <- unique(fread(file.path(ANN, "sorghum_all_genes_with_coords.tsv")), by = "GeneID")
coords[, Chr := as.character(Chr)]
N_GENES <- nrow(coords)
message("gene universe: ", N_GENES)

# ---- master ------------------------------------------------------------------
mst <- A[, .(Best_P_Value = min(Best_P), N_sig_SNPs = max(N_sig_SNPs),
             Max_r2 = max(Max_r2), N_Phenotypes = uniqueN(Trait),
             Phenotypes = paste(sort(unique(Trait)), collapse = "; "),
             Gene_Chr = Gene_Chr[1], Gene_Start = Gene_Start[1], Gene_End = Gene_End[1]),
         by = .(condition = Condition, layer, GeneID)]
mst <- merge(mst, func, by = "GeneID", all.x = TRUE)
mst <- merge(mst, coords[, .(GeneID, Chr, Start, End)], by = "GeneID", all.x = TRUE)
mst[, `:=`(Chromosome = fifelse(is.na(Chr), Gene_Chr, Chr),
           Gene_Start = fifelse(is.na(Start), Gene_Start, Start),
           Gene_End   = fifelse(is.na(End),   Gene_End,   End))]
for (cl in c("GeneName", "Family_Sufamily", "ProteinClass", "GO_MF", "GO_BP", "GO_CC")) {
  if (!cl %in% names(mst)) mst[[cl]] <- NA_character_
  mst[[cl]][is.na(mst[[cl]])] <- ""
}
COLS <- c("condition", "layer", "GeneID", "GeneName", "Family_Sufamily", "ProteinClass",
          "GO_MF", "GO_BP", "Chromosome", "Gene_Start", "Gene_End", "Best_P_Value",
          "N_sig_SNPs", "Max_r2", "N_Phenotypes", "Phenotypes", "GO_CC")
setorder(mst, condition, layer, Best_P_Value, GeneID)
fwrite(mst[, ..COLS], MASTER, sep = "\t", quote = FALSE)
message("Saved: ", MASTER, "  (", nrow(mst), " rows)")
print(mst[, .N, by = .(condition, layer)])

# ---- overlap helpers ---------------------------------------------------------
hyper <- function(a, b, s, N) phyper(s - 1, a, N - a, b, lower.tail = FALSE)
ovl <- function(lab, a, b, N, extra = NULL) {
  s <- length(intersect(a, b)); e <- length(a) * length(b) / N
  cbind(extra, data.table(Layer = lab, n_CTL = length(a), n_LIN = length(b), n_shared = s,
        expected_shared = e, fold_enrichment = s / e,
        jaccard = s / length(union(a, b)), p_hypergeom = hyper(length(a), length(b), s, N)))
}
LAYERS <- list("All layers" = NULL, "Individual lipids" = "individual",
               "Class sums / ratios" = "sumratio")
genes <- function(cond, lay) {
  d <- A[Condition == cond]; if (!is.null(lay)) d <- d[layer == lay]; unique(d$GeneID)
}

# ---- S20, gene level ---------------------------------------------------------
s20 <- rbindlist(lapply(names(LAYERS), function(nm) {
  a <- genes("CTL", LAYERS[[nm]]); b <- genes("LIN", LAYERS[[nm]])
  r <- ovl(nm, a, b, N_GENES)
  r[, `:=`(CTL_only = length(setdiff(a, b)), LIN_only = length(setdiff(b, a)))]
  setcolorder(r, c("Layer", "n_CTL", "n_LIN", "n_shared", "CTL_only", "LIN_only"))[]
}))
fwrite(s20, file.path(SUPP, "SuppTable_S20_Overlap_gene_level.csv"))
cat("\n-- S20, gene level --\n"); print(s20)

# ---- S21, window level -------------------------------------------------------
bins <- function(gs, kb) unique(coords[GeneID %in% gs, paste(Chr, Start %/% (kb * 1000L))])
s21 <- rbindlist(lapply(WINDOWS, function(kb) {
  NW <- uniqueN(coords[, paste(Chr, Start %/% (kb * 1000L))])
  rbindlist(lapply(names(LAYERS), function(nm) {
    a <- bins(genes("CTL", LAYERS[[nm]]), kb); b <- bins(genes("LIN", LAYERS[[nm]]), kb)
    ovl(nm, a, b, NW, data.table(window_kb = kb, n_windows_genome = NW))
  }))
}))
setcolorder(s21, c("window_kb", "Layer", "n_windows_genome"))
fwrite(s21, file.path(SUPP, "SuppTable_S21_Overlap_locus_level.csv"))
cat("\n-- S21, window level --\n"); print(s21)

# ---- gene-to-locus inflation, panel D of Supp Fig S8 -------------------------
# Was last written by the first-pass Python in May 2026 and carries the old
# candidate counts, which is why panel D disagreed with panels A to C.
INFL_DIR <- file.path(REPO, "table/overlap")
dir.create(INFL_DIR, recursive = TRUE, showWarnings = FALSE)
infl <- rbindlist(lapply(c("CTL", "LIN"), function(cond)
  rbindlist(lapply(list(c("individual lipids", "individual"),
                        c("class sums/ratios", "sumratio"),
                        c("all layers", NA)), function(z) {
    gs <- genes(cond, if (is.na(z[2])) NULL else z[2])
    lo <- uniqueN(coords[GeneID %in% gs, paste(Chr, Start %/% 250000L)])
    data.table(Condition = cond, Layer = z[1], n_genes = length(gs), n_loci = lo,
               genes_per_locus = length(gs) / lo)
  }))))
fwrite(infl, file.path(INFL_DIR, "gwas_gene_to_locus_inflation.csv"))
cat("\n-- gene-to-locus inflation --\n"); print(infl)

# ---- the two claims the paragraph makes that are not in either table ----------
sh <- intersect(genes("CTL", NULL), genes("LIN", NULL))
pairs_both <- A[GeneID %in% sh, .(k = uniqueN(Condition)), by = .(GeneID, Trait)][k == 2]
traits_both <- intersect(unique(A[Condition == "CTL"]$Trait), unique(A[Condition == "LIN"]$Trait))
cat("\nshared genes:", length(sh),
    "| gene-trait pairs recovered in both trials:", nrow(pairs_both),
    "| traits significant in both trials:", length(traits_both), "\n")
