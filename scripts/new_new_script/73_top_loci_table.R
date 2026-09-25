# ==============================================================================
# Table 4 -- the most recurrent GWAS loci under CTL and LIN.
#
#   Rscript scripts/new_new_script/73_top_loci_table.R
#
# A locus is one signal, defined as every candidate gene sharing a trial, trait
# layer, chromosome and gene-level best p-value, so the genes tagged by one lead
# SNP occupy a single row. Within a signal
#   Interval   min Gene_Start to max Gene_End, in Mb
#   Genes      number of candidate genes the signal tags
#   Gene       first gene by position carrying a real functional annotation, used
#              as a label for the locus and not as a causal assignment. "Uncharacterized
#              protein" does not count as one, since it labels nothing; a signal made
#              only of those falls back to its first gene.
#   Traits     number of distinct phenotypes over all genes in the signal
#   P          the shared best p-value
# The five signals with the most traits are taken per trial and trait layer,
# ties broken by p.
#
# Built from data/gene_annotation_final at Max_r2 >= 0.4, the same source and
# cutoff as Supplementary Tables S7-S10.
#
# Outputs
#   table/supp/Table4_top_loci.tsv
#   table/main/Table4_top_loci.tex   (the tabular body only)
# ==============================================================================
suppressPackageStartupMessages(library(data.table))

REPO   <- Sys.getenv("SOLD_REPO", ".")
ANNDIR <- file.path(REPO, "data/gene_annotation_final")
R2_CUT <- 0.4
dir.create(file.path(REPO, "table/main"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(REPO, "table/supp"), recursive = TRUE, showWarnings = FALSE)

g <- rbindlist(lapply(c("CTL_gene_annotation.tsv", "LIN_gene_annotation.tsv"),
                      function(f) fread(file.path(ANNDIR, f), sep = "\t", colClasses = "character")))
g[, `:=`(Best_P = as.numeric(Best_P), Max_r2 = as.numeric(Max_r2),
         Gene_Start = as.integer(Gene_Start), Gene_End = as.integer(Gene_End))]
g <- g[Max_r2 >= R2_CUT]

func <- unique(fread(file.path(ANNDIR, "GO_terms_all.txt"), sep = "\t",
                     colClasses = "character")[, .(GeneID, GeneName)], by = "GeneID")

# gene level first, so one gene contributes one p-value to its signal
gene <- g[, .(Best_P = min(Best_P), Gene_Start = Gene_Start[1], Gene_End = Gene_End[1],
              Traits = list(unique(Trait))),
          by = .(Condition, Layer, Gene_Chr, GeneID)]
gene <- merge(gene, func, by = "GeneID", all.x = TRUE)
gene[is.na(GeneName), GeneName := ""]

loci <- gene[, {
  o <- order(Gene_Start)
  informative <- GeneName[o] != "" & GeneName[o] != "Uncharacterized protein"
  lab <- GeneID[o][informative][1]
  if (is.na(lab)) lab <- GeneID[o][1]
  .(Start = min(Gene_Start), End = max(Gene_End), NGenes = .N,
    Gene = lab, Function = GeneName[GeneID == lab][1],
    NTraits = uniqueN(unlist(Traits)))
}, by = .(Condition, Layer, Gene_Chr, Best_P)]

# The annotation file writes slashes in protein names as underscores
# (GDSL esterase_lipase, serine_threonine protein kinase). Put them back, and
# capitalise the first letter so the column reads consistently.
loci[, Function := sub("^(.)", "\\U\\1", gsub("_", "/", Function), perl = TRUE)]

setorder(loci, Condition, Layer, -NTraits, Best_P)
top <- loci[, head(.SD, 5), by = .(Condition, Layer)]

top[, Interval := sprintf("%s:%.2f--%.2f", Gene_Chr, Start/1e6, End/1e6)]
out <- top[, .(Condition, Layer = sub("_lipids_final", "", Layer),
               Interval, Genes = NGenes, Gene, Function, Traits = NTraits, P = Best_P)]
fwrite(out, file.path(REPO, "table/supp/Table4_top_loci.tsv"), sep = "\t", quote = FALSE)

# ---- LaTeX body --------------------------------------------------------------
sci <- function(p) {
  e <- floor(log10(p)); m <- p/10^e
  sprintf("$%.2f\\times10^{%d}$", m, e)
}
esc <- function(x) gsub("_", "\\\\_", x)
lines <- character()
for (cond in c("CTL", "LIN")) {
  lines <- c(lines, sprintf("\\multicolumn{7}{l}{\\textbf{%s}} \\\\", cond))
  for (lay in c("individual", "sum_ratio")) {
    d <- out[Condition == cond & Layer == lay]
    if (!nrow(d)) next
    lab <- if (lay == "individual") "Individual lipids" else "Sums / ratios"
    for (i in seq_len(nrow(d))) {
      lines <- c(lines, sprintf("%s & %s & %d & \\texttt{%s} & %s & %d & %s \\\\",
                                if (i == 1) lab else "", d$Interval[i], d$Genes[i],
                                esc(d$Gene[i]),
                                if (d$Function[i] == "") "Not annotated" else d$Function[i],
                                d$Traits[i], sci(d$P[i])))
    }
    if (lay == "individual") lines <- c(lines, "\\addlinespace")
  }
  if (cond == "CTL") lines <- c(lines, "\\midrule")
}
writeLines(lines, file.path(REPO, "table/main/Table4_top_loci.tex"))
message("Saved: table/supp/Table4_top_loci.tsv and table/main/Table4_top_loci.tex")
print(as.data.frame(out))
