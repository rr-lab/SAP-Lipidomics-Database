# ==============================================================================
# Supplementary Tables S20 and S21 -- the CTL/LIN candidate-gene overlap at gene
# resolution and after collapsing to 100, 250 and 500 kb windows. Also writes the
# gene-to-locus inflation table that Supplementary Figure S8 panel D reads.
#
#   Rscript scripts/new_new_script/13a_SuppTableS20S21_overlap_gene_and_locus.R
#
# WHY THIS EXISTS. Both tables came from the first-pass Python
# (scripts/chapter2_addons/gwas_overlap.py and locus_overlap.py) and had no R
# equivalent, so they could not be refreshed when the candidate master changed.
# They went stale twice over -- through the 2026-09-03 species deduplication and
# again when the phantom GWAS phenotypes were removed on 2026-09-16 -- while
# Supplementary Figure S8 kept drawing from them.
#
# This reproduces every cell of the Python output that the intervening changes do
# not touch, to fourteen significant digits, which is what licenses the cells they
# do. It is the same construction script 14 uses for the per-class breakdown, so
# S20, S21, S30 and S31 now come from one definition of the candidate set.
#
# THE UNIVERSE. Gene level tests against every gene with coordinates in
# genes.range (34,027), not the annotated subset. Locus level tests against the
# number of windows of that size containing at least one such gene. A gene is
# placed in a window by its start coordinate.
#
# THE LAYERS. "All layers" is the union of a condition's individual and sum/ratio
# candidates, so it is not the sum of the two rows below it -- a gene called in
# both layers is counted once. That is why LIN all-layers can equal LIN
# individual even when the individual count drops.
#
# Input   data/LD_mapped/candidate_tables/ALL_LD_candidate_genes_master.tsv
#         data/LD_mapped/genes_ranges/genes.range
# Output  table/supp/SuppTable_S20_Overlap_gene_level.csv
#         table/supp/SuppTable_S21_Overlap_locus_level.csv
#         table/overlap/gwas_gene_to_locus_inflation.csv
# ==============================================================================
suppressPackageStartupMessages({ library(vroom); library(dplyr) })

REPO      <- Sys.getenv("SOLD_REPO", ".")
DATA_ROOT <- Sys.getenv("SOLD_DATA", file.path(REPO, "data"))
MASTER    <- file.path(DATA_ROOT, "LD_mapped/candidate_tables/ALL_LD_candidate_genes_master.tsv")
GENES     <- file.path(DATA_ROOT, "LD_mapped/genes_ranges/genes.range")
SUPP      <- file.path(REPO, "table/supp")
OVERLAP   <- file.path(REPO, "table/overlap")
dir.create(SUPP, recursive = TRUE, showWarnings = FALSE)
dir.create(OVERLAP, recursive = TRUE, showWarnings = FALSE)

mst <- vroom(MASTER, delim = "\t", show_col_types = FALSE, progress = FALSE)
gr  <- vroom(GENES, delim = "\t", col_names = c("chr", "start", "end", "gene"),
             show_col_types = FALSE, progress = FALSE)

LAYERS <- list("All layers"          = c("individual", "sumratio"),
               "Individual lipids"   = "individual",
               "Class sums / ratios" = "sumratio")

genes_of <- function(cond, ly) unique(mst$GeneID[mst$condition == cond & mst$layer %in% ly])

# hyper: P(X >= k) for k shared out of nA and nB drawn from N. phyper's q is
# P(X > q), hence k - 1.
hyp <- function(nA, nB, k, N) list(
  expected = nA * nB / N, fold = k / (nA * nB / N),
  jaccard  = k / (nA + nB - k),
  p        = phyper(k - 1, nA, N - nA, nB, lower.tail = FALSE))

# ---- S20, gene level ---------------------------------------------------------
N <- nrow(gr)
s20 <- bind_rows(lapply(names(LAYERS), function(nm) {
  ctl <- genes_of("CTL", LAYERS[[nm]]); lin <- genes_of("LIN", LAYERS[[nm]])
  sh  <- intersect(ctl, lin); h <- hyp(length(ctl), length(lin), length(sh), N)
  data.frame(Layer = nm, n_CTL = length(ctl), n_LIN = length(lin),
             n_shared = length(sh), CTL_only = length(setdiff(ctl, lin)),
             LIN_only = length(setdiff(lin, ctl)), expected_shared = h$expected,
             fold_enrichment = h$fold, jaccard = h$jaccard, p_hypergeom = h$p)
}))

# ---- S21, locus level --------------------------------------------------------
window_of <- function(kb) setNames(paste0(gr$chr, ":", gr$start %/% (kb * 1000)), gr$gene)
s21 <- bind_rows(lapply(c(100, 250, 500), function(kb) {
  win <- window_of(kb); NW <- length(unique(win))
  bind_rows(lapply(names(LAYERS), function(nm) {
    ctl <- unique(win[genes_of("CTL", LAYERS[[nm]])]); ctl <- ctl[!is.na(ctl)]
    lin <- unique(win[genes_of("LIN", LAYERS[[nm]])]); lin <- lin[!is.na(lin)]
    sh  <- intersect(ctl, lin); h <- hyp(length(ctl), length(lin), length(sh), NW)
    data.frame(window_kb = kb, Layer = nm, n_windows_genome = NW,
               n_CTL = length(ctl), n_LIN = length(lin), n_shared = length(sh),
               expected_shared = h$expected, fold_enrichment = h$fold,
               jaccard = h$jaccard, p_hypergeom = h$p)
  }))
}))

# ---- gene-to-locus inflation, at 250 kb --------------------------------------
LAB  <- c("All layers" = "all layers", "Individual lipids" = "individual lipids",
          "Class sums / ratios" = "class sums/ratios")
win250 <- window_of(250)
infl <- bind_rows(lapply(c("CTL", "LIN"), function(cond)
  bind_rows(lapply(c("Individual lipids", "Class sums / ratios", "All layers"), function(nm) {
    g <- genes_of(cond, LAYERS[[nm]]); w <- unique(win250[g]); w <- w[!is.na(w)]
    data.frame(Condition = cond, Layer = unname(LAB[nm]), n_genes = length(g),
               n_loci = length(w), genes_per_locus = length(g) / length(w))
  }))))

write.csv(s20,  file.path(SUPP, "SuppTable_S20_Overlap_gene_level.csv"),   row.names = FALSE)
write.csv(s21,  file.path(SUPP, "SuppTable_S21_Overlap_locus_level.csv"),  row.names = FALSE)
write.csv(infl, file.path(OVERLAP, "gwas_gene_to_locus_inflation.csv"),    row.names = FALSE)

message(sprintf("Gene universe N = %d", N))
message(sprintf("S20  %d rows    S21  %d rows    inflation  %d rows",
                nrow(s20), nrow(s21), nrow(infl)))
cat("\n-- gene level --\n")
print(within(s20, { expected_shared <- round(expected_shared, 1)
                    fold_enrichment <- round(fold_enrichment, 3)
                    jaccard <- round(jaccard, 4); p_hypergeom <- signif(p_hypergeom, 3) }),
      row.names = FALSE)
cat("\n-- genes per 250 kb locus --\n")
print(within(infl, genes_per_locus <- round(genes_per_locus, 2)), row.names = FALSE)
