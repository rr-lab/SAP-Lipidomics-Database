# ─────────────────────────────────────────────────────────────────────────────
# 78_GO_enriched_genes_two_column.R
#
# Two-column list of the candidate genes that fall inside GO terms passing the
# LD-aware permutation null, one column per trial. Genes enriched in both
# trials appear in both columns.
#
# Input : table/supp/SuppTable_S19_genes_in_enriched_GO_terms.tsv
# Output: table/supp/GO_enriched_genes_CTL_vs_LIN.csv
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages(library(data.table))

s19 <- fread("table/supp/SuppTable_S19_genes_in_enriched_GO_terms.tsv")

ctl <- sort(unique(s19[grepl("CTL", Conditions), GeneID]))
lin <- sort(unique(s19[grepl("LIN", Conditions), GeneID]))

message("CTL: ", length(ctl),
        "   LIN: ", length(lin),
        "   in both: ", length(intersect(ctl, lin)),
        "   unique genes: ", length(union(ctl, lin)))

n <- max(length(ctl), length(lin))
out <- data.table(
  CTL = c(ctl, rep("", n - length(ctl))),
  LIN = c(lin, rep("", n - length(lin)))
)

fwrite(out, "table/supp/GO_enriched_genes_CTL_vs_LIN.csv")
message("Saved: table/supp/GO_enriched_genes_CTL_vs_LIN.csv  (", n, " rows)")
