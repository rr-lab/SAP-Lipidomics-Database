# ==============================================================================
# GO enrichment -- collapse redundant terms, then rank by fold enrichment.
#
# Feeds Figure 4 and the GO supplementary tables. Three changes to how the
# results are selected:
#
#   1. Significance is q_LD alone. The LD-aware permutation with gene-density-
#      matched intervals is the test built for the "is this just one LD block"
#      worry, so an interval-count threshold on top of it tests the same thing
#      twice. The old intervals >= 3 cut also selected for large, vague GO terms:
#      terms passing it averaged 175 background genes at 10x enrichment, those
#      failing it averaged 24 genes at 23x. A small, specific term cannot span
#      three intervals because it does not contain enough genes.
#
#   2. Rows whose candidate-gene set is identical are one finding, not several.
#      GO nests parent and child terms, so a single locus can surface three
#      times: SORBI_3006G195000/195100/195200 appears as purine nucleobase
#      transport, purine nucleobase transmembrane transport and purine
#      nucleoside transmembrane transport. Such rows are merged, the most
#      specific term (smallest background) leads, and every merged term is kept
#      in a column so nothing is lost.
#
#   3. What remains is ranked by fold enrichment, which is meaningful only after
#      step 2, since before it the ranking compares synonyms with each other.
#
# Interval count is carried through as a column. It is a property of the
# evidence worth showing, not a filter.
#
# Input : table/go_enrichment/Table_GO_enrichment_all.tsv
# Output: table/go_enrichment/Table_GO_enrichment_collapsed.tsv
# ==============================================================================
source("scripts/new_new_script/_common.R")
suppressPackageStartupMessages({ library(tidyr); library(stringr); library(tibble) })

in_file  <- Sys.getenv("GO_TABLE",
  file.path(REPO, "table/go_enrichment/Table_GO_enrichment_all.tsv"))
out_file <- file.path(REPO, "table/go_enrichment/Table_GO_enrichment_collapsed.tsv")
Q_MAX    <- 0.05

stopifnot(file.exists(in_file))
go <- vroom(in_file, delim = "\t", show_col_types = FALSE)

# ---- 1. significance ---------------------------------------------------------
kept <- go %>% filter(!is.na(q_LD), q_LD < Q_MAX)

# ---- 2. collapse identical gene sets -----------------------------------------
# Terms sharing every candidate gene inside one stratum describe one locus.
# The label keeps the words the terms share and joins what differs, so
# "nitrate assimilation" + "nitrate transmembrane transport" becomes
# "nitrate assimilation and transmembrane transport".
merge_labels <- function(terms, background) {
  terms <- terms[order(background)]          # most specific first
  terms <- unique(terms)
  if (length(terms) == 1L) return(terms)
  w      <- strsplit(terms, " ", fixed = TRUE)
  n_min  <- min(lengths(w))
  shared <- 0L
  while (shared < n_min &&
         length(unique(vapply(w, function(x) x[shared + 1L], character(1)))) == 1L) {
    shared <- shared + 1L
  }
  stem <- if (shared) paste(w[[1]][seq_len(shared)], collapse = " ") else ""
  rest <- vapply(w, function(x)
    paste(x[seq.int(shared + 1L, length(x))], collapse = " "), character(1))
  if (length(rest) == 2L) {
    lab <- paste(rest, collapse = " and ")
  } else {
    lab <- sprintf("%s (+%d related terms)", rest[1], length(rest) - 1L)
  }
  str_squish(paste(stem, lab))
}

collapsed <- kept %>%
  group_by(Ontology, Layer, Condition, Class, Overlap_genes) %>%
  summarise(
    Term          = merge_labels(Term, Background),
    N_terms       = dplyr::n(),
    GO            = paste(sort(unique(GO)), collapse = ";"),
    Genes         = max(Genes),
    Background    = min(Background),
    Intervals     = max(Intervals),
    Fold          = max(Fold),
    q_gene        = min(q_gene),
    q_LD          = min(q_LD),
    .groups = "drop"
  )

# Term is rebuilt inside summarise, so the original names are collected separately.
merged_terms <- kept %>%
  group_by(Ontology, Layer, Condition, Class, Overlap_genes) %>%
  summarise(Merged_terms = paste(sort(unique(Term)), collapse = "; "), .groups = "drop")
collapsed <- collapsed %>%
  left_join(merged_terms, by = c("Ontology", "Layer", "Condition", "Class", "Overlap_genes"))

# ---- 3. rank by fold ---------------------------------------------------------
collapsed <- collapsed %>%
  arrange(Ontology, Layer, Condition, desc(Fold)) %>%
  dplyr::select(Ontology, Layer, Condition, Class, Term, N_terms, GO,
                Genes, Background, Intervals, Fold, q_gene, q_LD,
                Merged_terms, Overlap_genes)

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
write.table(collapsed, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
message("Saved: ", out_file, "  (", nrow(collapsed), " rows)")

# ---- what changed ------------------------------------------------------------
cat("\n-- selection ---------------------------------------------------------\n")
cat(sprintf("  all terms in input          %d\n", nrow(go)))
cat(sprintf("  q_LD < %.2f                  %d\n", Q_MAX, nrow(kept)))
cat(sprintf("  after collapsing duplicates %d   (%d rows merged away)\n",
            nrow(collapsed), nrow(kept) - nrow(collapsed)))
cat("\n  by ontology and layer:\n")
print(as.data.frame(collapsed %>% count(Ontology, Layer, Condition)))

cat("\n-- rows that were merged --------------------------------------------\n")
print(as.data.frame(collapsed %>% filter(N_terms > 1) %>%
        dplyr::select(Ontology, Condition, Class, Term, N_terms, Fold, Intervals)))

for (ont in c("BP", "MF")) {
  cat(sprintf("\n-- %s, ranked by fold enrichment ------------------------------\n", ont))
  print(as.data.frame(collapsed %>% filter(Ontology == ont) %>%
          mutate(Term = substr(Term, 1, 52)) %>%
          dplyr::select(Layer, Condition, Class, Term, Genes, Background,
                        Intervals, Fold, q_LD)))
}

# ==============================================================================
# Pass 2 -- collapse across lipid class, trait layer and ontology.
#
# Pass 1 merged terms sharing a gene set inside one stratum. The same gene set
# also recurs BETWEEN strata: SORBI_3006G195000/195100/195200 is one purine
# permease array, and it surfaces ten times across DG, SM, SQDG and FA, in both
# layers and both ontologies, described once as a process and once as an
# activity. Ten rows, one locus. Zinc ion binding does the same across seven
# classes under CTL.
#
# The unit of a finding here is therefore a gene set within a trial, not a row.
# Every stratum a gene set appears in is kept in the label, because that is
# informative in itself: a locus that surfaces in one lipid class only is
# class-specific biology, while one that surfaces across seven classes is a
# generic block that happens to sit near many candidates.
#
# Interval count travels through as a column for the supplementary table. It is
# deliberately NOT used to select or to rank: it scales with how large a GO term
# is, so filtering on it removes exactly the small, specific terms worth having.
# Significance remains q_LD, the LD-aware permutation, which is the test that
# actually addresses whether an enrichment is an artefact of linkage.
#
# Output: table/go_enrichment/Table_GO_enrichment_loci.tsv
# ==============================================================================
loci_file <- file.path(REPO, "table/go_enrichment/Table_GO_enrichment_loci.tsv")

loci <- collapsed %>%
  group_by(Condition, Overlap_genes) %>%
  arrange(Background, .by_group = TRUE) %>%
  summarise(
    Term        = dplyr::first(Term),          # most specific term in the group
    Classes     = paste(sort(unique(Class)), collapse = ", "),
    N_classes   = dplyr::n_distinct(Class),
    Ontologies  = paste(sort(unique(Ontology)), collapse = "/"),
    Layers      = paste(sort(unique(Layer)), collapse = "; "),
    N_strata    = dplyr::n(),
    Genes       = max(Genes),
    Background  = min(Background),
    Intervals   = max(Intervals),
    Fold_min    = min(Fold),
    Fold_max    = max(Fold),
    q_LD        = min(q_LD),
    All_terms   = paste(sort(unique(Merged_terms)), collapse = " | "),
    .groups = "drop"
  ) %>%
  arrange(Condition, desc(Fold_max))

write.table(loci, loci_file, sep = "\t", row.names = FALSE, quote = FALSE)
message("Saved: ", loci_file, "  (", nrow(loci), " rows)")

cat("\n======================================================================\n")
cat(sprintf("PASS 2  %d stratum rows  ->  %d distinct loci\n", nrow(collapsed), nrow(loci)))
cat("======================================================================\n")

cat("\n-- gene sets recurring across lipid classes (generic blocks) ----------\n")
print(as.data.frame(loci %>% filter(N_classes > 1) %>%
        mutate(Term = substr(Term, 1, 46)) %>%
        dplyr::select(Condition, Term, Classes, N_classes, N_strata,
                      Genes, Fold_min, Fold_max, q_LD)))

cat("\n-- class-specific loci, ranked by fold (the specific findings) --------\n")
print(as.data.frame(loci %>% filter(N_classes == 1) %>%
        mutate(Term = substr(Term, 1, 50)) %>%
        dplyr::select(Condition, Classes, Term, Ontologies, Genes,
                      Background, Fold_max, q_LD)), max = 1000)
