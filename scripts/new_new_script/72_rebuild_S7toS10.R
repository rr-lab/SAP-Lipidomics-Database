# ==============================================================================
# Supplementary Tables S7-S10 -- GWAS candidate genes, rebuilt from the current
# annotation tables.
#
#   Rscript scripts/new_new_script/72_rebuild_S7toS10.R
#
# Why this exists. The S7-S10 files shipped in table/supp/ were built on 2025-08-09
# and 2025-09-03, before two things changed. The high-heterozygosity marker filter
# removed 83,640 markers, and the class-sum and ratio phenotypes were rebuilt. The
# old files therefore disagree with the manuscript on every count.
#
#              old file    manuscript
#   S7  CTL individual      1,100        1,062
#   S8  CTL sum/ratio         115           54
#   S9  LIN individual      4,319        4,370
#   S10 LIN sum/ratio         812          385
#
# The old S7 also still carries SORBI_3001G522600 at p = 4.5e-35, which is the
# signal traced to 23 removed markers sitting in ~300 bp at heterozygosity 0.81
# to 1.00. Its best p after the filter is 1.15e-05, well below significance.
#
# Source of truth is data/gene_annotation_final/{CTL,LIN}_gene_annotation.tsv,
# one row per gene per phenotype, which is what every number in the GWAS section
# is computed from. Rows are kept at Max_r2 >= 0.4, the LD cutoff the paper uses.
#
# Per-gene aggregation follows the convention of the files being replaced, checked
# against them gene by gene where the underlying data did not change
#   Best_P_Value   minimum Best_P over the gene's phenotypes
#   N_sig_SNPs     maximum N_sig_SNPs over the gene's phenotypes, not the sum,
#                  since the same SNP tags the gene for every phenotype it hits
#   Max_r2         maximum Max_r2 over the gene's phenotypes
#   N_Phenotypes   number of distinct phenotypes
#   Phenotypes     those phenotypes, alphabetical, "; " separated
#
# Functional annotation is joined from data/gene_annotation_final/GO_terms_all.txt
# and coordinates from sorghum_all_genes_with_coords.tsv. 83 of the 5,176 candidate
# genes carry no functional annotation in that file; their GeneName and GO columns
# are left empty rather than filled in, and the count is reported in the summary.
#
# Outputs (table/supp/)
#   SuppTable_S7_CTL_GWAS_candidate_genes_for_individual_lipid_traits.tsv
#   SuppTable_S8_CTL_GWAS_candidate_genes_for_lipid_sum_ratio_traits.tsv
#   SuppTable_S9_LIN_GWAS_candidate_genes_for_individual_lipid_traits.tsv
#   SuppTable_S10_LIN_GWAS_candidate_genes_for_lipid_sum_ratio_traits.tsv
#   SuppTable_S7toS10_summary.tsv
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)
})

REPO    <- Sys.getenv("SOLD_REPO", ".")
ANNDIR  <- file.path(REPO, "data/gene_annotation_final")
OUTDIR  <- file.path(REPO, "table/supp")
R2_CUT  <- 0.4
WIN     <- 250000
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

read_ann <- function(f) fread(file.path(ANNDIR, f), sep = "\t", colClasses = "character")

gwas <- rbindlist(list(read_ann("CTL_gene_annotation.tsv"),
                       read_ann("LIN_gene_annotation.tsv")))
gwas[, `:=`(Best_P = as.numeric(Best_P), Max_r2 = as.numeric(Max_r2),
            N_sig_SNPs = as.integer(N_sig_SNPs),
            Gene_Start = as.integer(Gene_Start), Gene_End = as.integer(Gene_End))]
gwas <- gwas[Max_r2 >= R2_CUT]

func <- fread(file.path(ANNDIR, "GO_terms_all.txt"), sep = "\t", colClasses = "character")
setnames(func, c("Family_Subfamily", "Protein_Class"), c("Family_SubFamily", "Protein_class"),
         skip_absent = TRUE)
func <- unique(func, by = "GeneID")

coords <- fread(file.path(ANNDIR, "sorghum_all_genes_with_coords.tsv"), sep = "\t",
                colClasses = list(character = "GeneID"))
coords <- unique(coords, by = "GeneID")

COLS <- c("GeneID", "GeneName", "Family_SubFamily", "Protein_class",
          "GO_MF", "GO_BP", "GO_CC", "Chromosome", "Gene_Start", "Gene_End",
          "Best_P_Value", "N_sig_SNPs", "Max_r2", "N_Phenotypes", "Phenotypes")

build <- function(cond, layer) {
  d <- gwas[Condition == cond & Layer == layer]
  if (!nrow(d)) return(NULL)
  g <- d[, .(Best_P_Value = min(Best_P),
             N_sig_SNPs   = max(N_sig_SNPs),
             Max_r2       = max(Max_r2),
             N_Phenotypes = uniqueN(Trait),
             Phenotypes   = paste(sort(unique(Trait)), collapse = "; "),
             Gene_Chr     = Gene_Chr[1],
             Gene_Start   = Gene_Start[1],
             Gene_End     = Gene_End[1]),
           by = GeneID]
  g <- merge(g, func, by = "GeneID", all.x = TRUE)
  g <- merge(g, coords[, .(GeneID, Chr, Start, End)], by = "GeneID", all.x = TRUE)
  # coordinates come from the genome annotation where available, from the GWAS
  # table otherwise, so no gene is written without a position
  g[, `:=`(Chromosome = fifelse(is.na(Chr), Gene_Chr, as.character(Chr)),
           Gene_Start = fifelse(is.na(Start), Gene_Start, Start),
           Gene_End   = fifelse(is.na(End),   Gene_End,   End))]
  for (cl in c("GeneName", "Family_SubFamily", "Protein_class", "GO_MF", "GO_BP", "GO_CC")) {
    if (!cl %in% names(g)) g[[cl]] <- NA_character_
    g[[cl]][is.na(g[[cl]])] <- ""
  }
  setorder(g, Best_P_Value, GeneID)
  g[, ..COLS]
}

SPEC <- list(
  list(tab = "S7",  cond = "CTL", layer = "individual_lipids_final",
       file = "SuppTable_S7_CTL_GWAS_candidate_genes_for_individual_lipid_traits.tsv"),
  list(tab = "S8",  cond = "CTL", layer = "sum_ratio_lipids_final",
       file = "SuppTable_S8_CTL_GWAS_candidate_genes_for_lipid_sum_ratio_traits.tsv"),
  list(tab = "S9",  cond = "LIN", layer = "individual_lipids_final",
       file = "SuppTable_S9_LIN_GWAS_candidate_genes_for_individual_lipid_traits.tsv"),
  list(tab = "S10", cond = "LIN", layer = "sum_ratio_lipids_final",
       file = "SuppTable_S10_LIN_GWAS_candidate_genes_for_lipid_sum_ratio_traits.tsv")
)

summ <- rbindlist(lapply(SPEC, function(s) {
  g <- build(s$cond, s$layer)
  fwrite(g, file.path(OUTDIR, s$file), sep = "\t", quote = FALSE)
  message("Saved: ", s$file, "  (", nrow(g), " genes)")
  d <- gwas[Condition == s$cond & Layer == s$layer]
  data.table(
    Table                    = s$tab,
    Trial                    = s$cond,
    Layer                    = sub("_lipids_final", "", s$layer),
    Candidate_genes          = nrow(g),
    Phenotypes_with_candidates = uniqueN(d$Trait),
    Genes_multi_phenotype    = sum(g$N_Phenotypes > 1),
    Median_phenotypes_per_gene = median(g$N_Phenotypes),
    Max_phenotypes_per_gene  = max(g$N_Phenotypes),
    Min_P                    = min(g$Best_P_Value),
    LD_cutoff_r2             = R2_CUT,
    Loci_250kb               = uniqueN(d[, paste(Gene_Chr, Gene_Start %/% WIN)]),
    Genes_without_functional_annotation = sum(g$GeneName == "")
  )
}))

# Loci pooled over both layers within a trial, which is what the manuscript quotes
pooled <- gwas[, .(Loci_250kb_both_layers = uniqueN(paste(Gene_Chr, Gene_Start %/% WIN)),
                   Candidate_genes_both_layers = uniqueN(GeneID)), by = Condition]
summ <- merge(summ, pooled, by.x = "Trial", by.y = "Condition", all.x = TRUE)
setorder(summ, Trial, Table)
fwrite(summ, file.path(OUTDIR, "SuppTable_S7toS10_summary.tsv"), sep = "\t", quote = FALSE)
message("Saved: SuppTable_S7toS10_summary.tsv")
print(as.data.frame(summ))
