#!/usr/bin/env Rscript

# GO biological-process over-representation using Storey's q-value method.
# This is a sensitivity analysis parallel to go_enrichment_LD_by_condition_layer.py.
# It keeps the same candidate sets, GO background, term filters, and one-sided
# Fisher tests, changing only the multiple-testing correction.

suppressPackageStartupMessages({
  library(qvalue)
})

repo <- "/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database"
annot_file <- file.path(repo, "data/annotation/gene_annotation.txt")
master_file <- file.path(repo, "data/LD_mapped/candidate_tables/ALL_LD_candidate_genes_master.tsv")
out_dir <- file.path(repo, "data/LD_mapped/go_enrichment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

parse_go_bp <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character())
  terms <- strsplit(x, ";", fixed = TRUE)[[1]]
  terms <- trimws(terms)
  terms <- terms[nzchar(terms)]
  sub("\\(GO:[0-9]+\\)$", "", terms)
}

annotation <- read.delim(annot_file, sep = "\t", quote = "", comment.char = "",
                         check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c("GeneID", "GO_BP") %in% names(annotation)))

gene_terms <- lapply(annotation$GO_BP, parse_go_bp)
names(gene_terms) <- trimws(annotation$GeneID)
gene_terms <- gene_terms[lengths(gene_terms) > 0]
background <- names(gene_terms)
n_background <- length(background)

term_genes <- list()
for (gene in background) {
  for (term in gene_terms[[gene]]) {
    term_genes[[term]] <- unique(c(term_genes[[term]], gene))
  }
}

master <- read.delim(master_file, sep = "\t", quote = "", comment.char = "",
                     check.names = FALSE, stringsAsFactors = FALSE)
sets <- split(trimws(master$GeneID), paste(master$condition, master$layer, sep = "_"))
sets <- lapply(sets, unique)

enrich_one <- function(candidate_genes) {
  tested <- intersect(unique(candidate_genes), background)
  n_tested <- length(tested)
  rows <- list()

  for (term in names(term_genes)) {
    term_set <- term_genes[[term]]
    term_bg <- length(term_set)
    overlap <- intersect(tested, term_set)
    n_overlap <- length(overlap)
    if (term_bg < 5 || n_overlap < 3) next

    a <- n_overlap
    b <- n_tested - n_overlap
    c <- term_bg - n_overlap
    d <- n_background - term_bg - b
    p_value <- fisher.test(matrix(c(a, b, c, d), nrow = 2),
                           alternative = "greater")$p.value
    expected <- n_tested * term_bg / n_background
    rows[[length(rows) + 1L]] <- data.frame(
      go_term = term,
      term_bg_count = term_bg,
      term_test_count = n_overlap,
      expected = expected,
      fold_enrichment = n_overlap / expected,
      p_value = p_value,
      overlap_genes = paste(sort(overlap), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }

  if (!length(rows)) {
    return(list(n_tested = n_tested, results = data.frame()))
  }

  results <- do.call(rbind, rows)
  # Small term sets can make the package's pi0 estimator undefined. Use the
  # standard estimator when valid, and conservatively fix pi0 = 1 otherwise.
  qv <- tryCatch(
    qvalue::qvalue(results$p_value)$qvalues,
    error = function(e) qvalue::qvalue(results$p_value, pi0 = 1)$qvalues
  )
  results$q_value <- qv
  results <- results[order(results$q_value, -results$fold_enrichment), ]
  rownames(results) <- NULL
  list(n_tested = n_tested, results = results)
}

layer_order <- c("CTL_individual", "CTL_sumratio", "LIN_individual", "LIN_sumratio")
summary_rows <- list()

for (layer_name in layer_order) {
  result <- enrich_one(sets[[layer_name]])
  output <- result$results
  output_file <- file.path(out_dir, paste0(layer_name, "_GO_BP_enrichment_LD_qvalue.tsv"))

  if (nrow(output)) {
    write.table(output, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    n_significant <- sum(output$q_value < 0.05, na.rm = TRUE)
  } else {
    write.table(data.frame(), output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    n_significant <- 0L
  }

  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    layer = layer_name,
    candidate_genes = length(sets[[layer_name]]),
    genes_tested_with_GO_BP = result$n_tested,
    terms_tested = nrow(output),
    qvalue_below_0.05 = n_significant,
    stringsAsFactors = FALSE
  )

  message(sprintf("%s: %d candidate genes, %d GO-BP-testable genes, %d terms, %d q < 0.05",
                  layer_name, length(sets[[layer_name]]), result$n_tested,
                  nrow(output), n_significant))
}

summary_file <- file.path(out_dir, "GO_BP_enrichment_LD_qvalue_summary.tsv")
write.table(do.call(rbind, summary_rows), summary_file, sep = "\t",
            quote = FALSE, row.names = FALSE)
message("Outputs written to: ", out_dir)
