# =============================================================================
# Supplementary GO tables S16-S19, built from the rebuilt enrichment.
#
# Script 10 writes the collapsed and loci tables; this turns them into the
# SuppTable_* files that 15_SuppTables_build_workbooks.R packages.
#
#   S16  loci collapsed          <- Table_GO_enrichment_loci.tsv
#   S17  BP all terms            <- Table_GO_enrichment_all.tsv, BP only
#   S18  MF all terms            <- Table_GO_enrichment_all.tsv, MF only
#   S19  genes in enriched terms <- one row per gene appearing in a term that
#                                   passed, with where it came from
#
# S17 and S18 in the previous supplement were byte-identical and each held the
# whole table, BP and MF together, so the sheet names did not describe the
# contents. They are split properly here.
#
# Input   table/go_enrichment/{Table_GO_enrichment_all,Table_GO_enrichment_loci}.tsv
#         data/gene_annotation_final/{CTL,LIN}_gene_annotation.tsv
#         data/gene_annotation_final/sorghum_GO_universe.txt
# Output  table/supp/SuppTable_S1[6-9]_*.tsv
#         table/go_enrichment/genes_in_enriched_GO_terms.tsv
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })

REPO  <- Sys.getenv("SOLD_REPO", ".")
GO    <- file.path(REPO, "table/go_enrichment")
ANN   <- Sys.getenv("GO_IN", file.path(REPO, "data/gene_annotation_final"))
SUPP  <- file.path(REPO, "table/supp")
Q_MAX <- 0.05
dir.create(SUPP, recursive = TRUE, showWarnings = FALSE)

all  <- fread(file.path(GO, "Table_GO_enrichment_all.tsv"))
loci <- fread(file.path(GO, "Table_GO_enrichment_loci.tsv"))

# ---- S16, S17, S18 -----------------------------------------------------------
fwrite(loci, file.path(SUPP, "SuppTable_S16_GO_loci_collapsed.tsv"), sep = "\t")
for (o in c("BP", "MF")) {
  x <- all[Ontology == o & q_LD < Q_MAX][order(Condition, Layer, Class, q_LD, -Fold)]
  n <- if (o == "BP") "SuppTable_S17_GO_BP_all_terms.tsv" else "SuppTable_S18_GO_MF_all_terms.tsv"
  fwrite(x, file.path(SUPP, n), sep = "\t")
  message(n, ": ", nrow(x), " rows")
}
message("SuppTable_S16_GO_loci_collapsed.tsv: ", nrow(loci), " rows")

# ---- S19, genes inside the terms that passed --------------------------------
sig <- all[q_LD < Q_MAX]
g <- sig[, .(GeneID = trimws(unlist(strsplit(Overlap_genes, ";", fixed = TRUE)))),
         by = .(Condition, Layer, Class, Ontology, Term, GO, Fold, Intervals, q_LD)]

ann <- rbind(fread(file.path(ANN, "CTL_gene_annotation.tsv")),
             fread(file.path(ANN, "LIN_gene_annotation.tsv")))
gw <- ann[, .(Best_GWAS_p = min(Best_P), Max_r2 = max(Max_r2),
              Max_recurrence = uniqueN(Trait), Chr = Gene_Chr[1], Start = Gene_Start[1]),
          by = .(Condition, GeneID)]

uni <- fread(file.path(ANN, "sorghum_GO_universe.txt"), quote = "")
setnames(uni, c("GeneID","GeneName","Family_Subfamily","Protein_Class","GO_MF","GO_BP","GO_CC"))

s19 <- g[, .(Conditions     = paste(sort(unique(Condition)), collapse = ";"),
             Ontologies     = paste(sort(unique(Ontology)),  collapse = ";"),
             Classes        = paste(sort(unique(Class)),     collapse = ";"),
             Max_fold       = max(Fold),
             Max_intervals  = max(Intervals),
             N_terms        = uniqueN(GO),
             Specific       = if (max(Intervals) >= 2) "yes" else "single locus",
             Terms          = paste(sort(unique(Term)), collapse = " | ")),
         by = GeneID]
s19 <- merge(s19, gw[, .(Best_GWAS_p = min(Best_GWAS_p), Max_r2 = max(Max_r2),
                         Max_recurrence = max(Max_recurrence), Chr = Chr[1], Start = Start[1]),
                     by = GeneID], by = "GeneID", all.x = TRUE)
s19 <- merge(s19, uni[, .(GeneID, GeneName, Family = Family_Subfamily, Protein_Class)],
             by = "GeneID", all.x = TRUE)

# Flag reproduces the previous table's coarse grouping, from the PANTHER class
flag <- function(pc, nm) {
  pc <- ifelse(is.na(pc), "", pc); nm <- ifelse(is.na(nm), "", nm)
  fifelse(grepl("lipid|lipase|acyltransferase|phospholip", paste(pc, nm), ignore.case = TRUE),
          "LIPID enzyme/transport",
  fifelse(grepl("transporter|transmembrane|channel|carrier", paste(pc, nm), ignore.case = TRUE),
          "transport",
  fifelse(grepl("transferase|kinase|hydrolase|oxidase|reductase|synthase|ligase|isomerase|peptidase",
                paste(pc, nm), ignore.case = TRUE), "metabolic enzyme",
  fifelse(grepl("repeat|domain-containing|zinc finger", paste(pc, nm), ignore.case = TRUE),
          "repeat-associated domain", ""))))
}
s19[, Flag := flag(Protein_Class, GeneName)]
setcolorder(s19, c("GeneID","GeneName","Flag","Conditions","Ontologies","Classes","Max_recurrence",
                   "Best_GWAS_p","Max_r2","Max_fold","Max_intervals","Specific","N_terms",
                   "Chr","Start","Family","Terms"))
s19[, Protein_Class := NULL]
setorder(s19, -Max_fold, Best_GWAS_p)
fwrite(s19, file.path(GO, "genes_in_enriched_GO_terms.tsv"), sep = "\t")
fwrite(s19, file.path(SUPP, "SuppTable_S19_genes_in_enriched_GO_terms.tsv"), sep = "\t")
message("S19: ", nrow(s19), " genes")

cat("\n-- S19 composition --\n")
print(s19[, .N, by = Flag][order(-N)])
print(s19[, .N, by = Specific])
cat("\n-- top genes by fold --\n")
print(head(s19[, .(GeneID, GeneName = substr(GeneName,1,40), Classes, Max_fold, Max_intervals, Best_GWAS_p)], 12),
      row.names = FALSE)
