# ---------------------------------------------------------------------------
# Main table: most recurrent GWAS loci under CTL and LIN
#
# Recurrence is counted per SIGNAL, not per gene. A single lead SNP tags every
# gene inside its LD window, so ranking raw genes just lists one locus many
# times over. A signal is therefore defined as one (trial, trait layer,
# chromosome, best p-value) combination, which groups all genes carried by the
# same lead SNP into one row. N_genes reports how many genes that signal tagged.
#
# The gene named on each row is the first annotated gene in the interval. It is
# a label for the locus, not a causal assignment: LD cannot resolve which gene
# in the window drives the association.
#
# Source  data/LD_mapped/candidate_tables/ALL_LD_candidate_genes_master.tsv
#         (the same candidate set the Results text counts: 1,100 / 115 CTL and
#         4,323 / 812 LIN genes)
# Outputs new_new_tables/main/MainTable_TopRecurrentLoci.csv
#         new_new_tables/main/MainTable_TopRecurrentLoci.tex
# ---------------------------------------------------------------------------

suppressMessages(library(data.table))

in_file  <- "data/LD_mapped/candidate_tables/ALL_LD_candidate_genes_master.tsv"
out_dir  <- "new_new_tables/main"
n_top    <- 5L
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

d <- fread(in_file)

# ---- putative function, with fallbacks ------------------------------------
is_blank <- function(x) is.na(x) | trimws(x) == "" | trimws(x) == "NA"
d[, Fn := trimws(GeneName)]
d[is_blank(Fn), Fn := trimws(Family_Sufamily)]
d[is_blank(Fn), Fn := trimws(ProteinClass)]
d[is_blank(Fn), Fn := "Uncharacterized protein"]
d[, informative := !(Fn %in% c("Uncharacterized protein", "Uncharacterised protein"))]

# shorten the few PANTHER-style names that are too long for a printed table
abbrev <- c(
  "Major facilitator superfamily (MFS) profile domain-containing protein" =
    "MFS transporter domain protein",
  "Cytoplasmic tRNA 2-thiolation protein 1 C-terminal domain-containing protein" =
    "tRNA 2-thiolation protein 1",
  "Aminotransferase-like plant mobile domain-containing protein" =
    "Aminotransferase-like domain protein",
  "SS18 N-terminal domain-containing protein" = "SS18 N-terminal domain protein"
)
d[Fn %in% names(abbrev), Fn := abbrev[Fn]]

# ---- collapse genes into signals -------------------------------------------
sig <- d[, {
  m    <- max(N_Phenotypes)
  cand <- .SD[N_Phenotypes == m][order(-informative, GeneID)]
  .(N_Phenotypes = m,
    N_genes      = .N,
    Max_r2       = max(Max_r2),
    Start_Mb     = min(Gene_Start) / 1e6,
    End_Mb       = max(Gene_End)   / 1e6,
    Gene         = cand$GeneID[1],
    Function     = cand$Fn[1])
}, by = .(condition, layer, Chromosome, Best_P_Value)]

cat("Signals recovered per trial and layer\n")
print(sig[, .N, by = .(condition, layer)])

# ---- top n per trial x layer -----------------------------------------------
top <- sig[order(condition, layer, -N_Phenotypes, Best_P_Value)][
  , head(.SD, n_top), by = .(condition, layer)]

top[, Layer := fifelse(layer == "individual", "Individual lipids", "Sums / ratios")]
top[, Interval := sprintf("%d:%.2f--%.2f", Chromosome, Start_Mb, End_Mb)]
setorder(top, condition, -N_Phenotypes, Best_P_Value)

out <- top[, .(Trial = condition, Layer, Interval_Mb = Interval,
               N_genes, Gene, Function,
               N_Phenotypes, Best_P_Value, Max_r2 = round(Max_r2, 2))]
fwrite(out, file.path(out_dir, "MainTable_TopRecurrentLoci.csv"))

# ---- LaTeX body -------------------------------------------------------------
sci <- function(p) {
  e <- floor(log10(p)); m <- p / 10^e
  sprintf("$%.2f\\times10^{%d}$", m, e)
}
esc <- function(x) gsub("_", "\\\\_", x)

lines <- character(0)
for (cond in c("CTL", "LIN")) {
  lines <- c(lines, sprintf("\\multicolumn{7}{l}{\\textbf{%s}} \\\\", cond))
  blk <- out[Trial == cond]
  for (lay in c("Individual lipids", "Sums / ratios")) {
    sub <- blk[Layer == lay]
    for (i in seq_len(nrow(sub))) {
      lines <- c(lines, sprintf(
        "%s & %s & %d & \\texttt{%s} & %s & %d & %s \\\\",
        if (i == 1) lay else "",
        sub$Interval_Mb[i], sub$N_genes[i], esc(sub$Gene[i]),
        sub$Function[i], sub$N_Phenotypes[i], sci(sub$Best_P_Value[i])))
    }
    if (lay == "Individual lipids") lines <- c(lines, "\\addlinespace")
  }
  if (cond == "CTL") lines <- c(lines, "\\midrule")
}
writeLines(lines, file.path(out_dir, "MainTable_TopRecurrentLoci.tex"))

cat("\nWrote", file.path(out_dir, "MainTable_TopRecurrentLoci.csv"), "and .tex\n")
print(out)
