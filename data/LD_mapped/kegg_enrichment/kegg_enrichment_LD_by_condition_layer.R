#!/usr/bin/env Rscript

# KEGG pathway over-representation for the four LD-mapped GWAS layers.
#
# The local UniProt export provides SORBI to UniProt links. KEGGREST then
# retrieves the current UniProt to Sorghum KEGG mapping from KEGG online.

suppressPackageStartupMessages({
  library(KEGGREST)
})

repo <- "/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database"
master_file <- file.path(repo, "data/LD_mapped/candidate_tables/ALL_LD_candidate_genes_master.tsv")
mapping_file <- file.path(repo, "data/annotation/SORBI_to_KEGG.tsv")
uniprot_file <- file.path(repo, "data/LD_mapped/candidate_tables/uniprot.csv")
out_dir <- file.path(repo, "data/LD_mapped/kegg_enrichment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

master <- read.delim(master_file, sep = "\t", quote = "", comment.char = "",
                     check.names = FALSE, stringsAsFactors = FALSE)
sets <- split(trimws(master$GeneID), paste(master$condition, master$layer, sep = "_"))
sets <- lapply(sets, unique)

if (!file.exists(uniprot_file)) {
  stop("Missing UniProt mapping file: ", uniprot_file)
}

uniprot <- read.csv(uniprot_file, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(ncol(uniprot) >= 2)
uniprot_ids <- trimws(uniprot[[1]])
sorbi_ids <- trimws(uniprot[[2]])
uniprot_to_sorbi <- data.frame(
  uniprot_id = paste0("up:", uniprot_ids),
  SORBI_ID = sorbi_ids,
  stringsAsFactors = FALSE
)
uniprot_to_sorbi <- uniprot_to_sorbi[
  nzchar(uniprot_to_sorbi$uniprot_id) & nzchar(uniprot_to_sorbi$SORBI_ID), ]

message("Retrieving UniProt to Sorghum KEGG mappings from KEGG")
kegg_conversion <- keggConv("uniprot", "sbi")
mapping <- merge(
  uniprot_to_sorbi,
  data.frame(
    KEGG_GENE_ID = names(kegg_conversion),
    uniprot_id = unname(kegg_conversion),
    stringsAsFactors = FALSE
  ),
  by = "uniprot_id"
)
mapping <- unique(mapping[c("SORBI_ID", "KEGG_GENE_ID")])

# Use the mapped master-gene universe unless a broader tested-gene mapping is
# supplied. A full GWAS-tested background is preferable for publication.
background_sorbi <- unique(master$GeneID)
background_sorbi <- intersect(background_sorbi, mapping$SORBI_ID)
candidate_to_kegg <- function(genes) {
  unique(mapping$KEGG_GENE_ID[mapping$SORBI_ID %in% genes])
}

background_kegg <- candidate_to_kegg(background_sorbi)
if (length(background_kegg) < 10) {
  stop("Fewer than 10 mapped KEGG background genes were found. Check the mapping file.")
}

message("Querying KEGG pathway links for ", length(background_kegg), " mapped genes")
# Download the organism-wide link table once, then filter locally. This avoids
# KEGG's request-URI limit when querying thousands of genes individually.
all_pathway_links <- keggLink("pathway", "sbi")
pathway_links <- all_pathway_links[names(all_pathway_links) %in% background_kegg]
if (is.null(pathway_links) || !length(pathway_links)) {
  stop("KEGG returned no pathway links for the mapped genes.")
}

link_df <- data.frame(
  kegg_gene = names(pathway_links),
  pathway_id = unname(pathway_links),
  stringsAsFactors = FALSE
)
pathway_ids <- unique(link_df$pathway_id)
all_pathway_names <- keggList("pathway", "sbi")
pathway_names <- setNames(pathway_ids, pathway_ids)
pathway_names[sub("^path:", "", names(all_pathway_names))] <- unname(all_pathway_names)

pathway_genes <- split(link_df$kegg_gene, link_df$pathway_id)
n_background <- length(background_kegg)

enrich_one <- function(candidate_genes) {
  candidate_kegg <- candidate_to_kegg(candidate_genes)
  tested <- intersect(candidate_kegg, background_kegg)
  n_tested <- length(tested)
  rows <- list()

  for (pid in names(pathway_genes)) {
    pathway_set <- intersect(unique(pathway_genes[[pid]]), background_kegg)
    pathway_bg <- length(pathway_set)
    overlap <- intersect(tested, pathway_set)
    n_overlap <- length(overlap)
    if (pathway_bg < 5 || n_overlap < 3) next

    a <- n_overlap
    b <- n_tested - n_overlap
    c <- pathway_bg - n_overlap
    d <- n_background - pathway_bg - b
    p_value <- fisher.test(matrix(c(a, b, c, d), nrow = 2),
                           alternative = "greater")$p.value
    expected <- n_tested * pathway_bg / n_background
    sorbi_overlap <- mapping$SORBI_ID[match(overlap, mapping$KEGG_GENE_ID)]
    rows[[length(rows) + 1L]] <- data.frame(
      pathway_id = pid,
      pathway_name = unname(pathway_names[[sub("^path:", "", pid)]]),
      pathway_bg_count = pathway_bg,
      pathway_test_count = n_overlap,
      expected = expected,
      fold_enrichment = n_overlap / expected,
      p_value = p_value,
      overlap_genes = paste(sort(unique(sorbi_overlap)), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }

  if (!length(rows)) return(list(n_tested = n_tested, results = data.frame()))
  results <- do.call(rbind, rows)
  results$p_adj_bh <- p.adjust(results$p_value, method = "BH")
  results <- results[order(results$p_adj_bh, -results$fold_enrichment), ]
  rownames(results) <- NULL
  list(n_tested = n_tested, results = results)
}

layer_order <- c("CTL_individual", "CTL_sumratio", "LIN_individual", "LIN_sumratio")
summary_rows <- list()
for (layer_name in layer_order) {
  result <- enrich_one(sets[[layer_name]])
  output <- result$results
  output_file <- file.path(out_dir, paste0(layer_name, "_KEGG_enrichment.tsv"))
  write.table(output, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  n_sig <- if (nrow(output)) sum(output$p_adj_bh < 0.05) else 0L
  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    layer = layer_name,
    candidate_genes = length(sets[[layer_name]]),
    mapped_candidate_genes = length(candidate_to_kegg(sets[[layer_name]])),
    mapped_background_genes = length(background_kegg),
    pathways_tested = nrow(output),
    bh_below_0.05 = n_sig,
    stringsAsFactors = FALSE
  )
  message(sprintf("%s: %d candidate genes, %d mapped, %d pathways, %d BH-significant",
                  layer_name, length(sets[[layer_name]]),
                  length(candidate_to_kegg(sets[[layer_name]])), nrow(output), n_sig))
}

write.table(do.call(rbind, summary_rows),
            file.path(out_dir, "KEGG_enrichment_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Outputs written to: ", out_dir)
