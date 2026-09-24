# =============================================================================
# Rebuild the gene-counts tables from the annotation tables.
#
# The counts files shipped with a broken N_sig_SNPs. It was built by summing the
# per-trait marker counts, and a marker significant for 11 traits was counted 11
# times, which turned the 2 markers tagging the chromosome 3 block into 313.
#
# This rebuilds them from CTL/LIN_gene_annotation.tsv alone, so it needs no
# server files. The annotation table stores the marker count per trait x gene,
# which is genuinely distinct, and the lead marker, but not the marker
# identities, so the exact union across traits is not recoverable from it. The
# broken column is therefore replaced with three quantities that are, and each
# says something the sum did not:
#
#   Max_sig_SNPs_per_trait  the largest marker count for any single trait. This
#                           is a lower bound on the union and is the number the
#                           manuscript quotes, e.g. 52 for FA(22:1) at the
#                           chromosome 3 block.
#   Sum_sig_SNPs_per_trait  the old column, kept and renamed so it is obvious
#                           that it counts trait-marker pairs, not markers.
#   N_lead_SNPs             distinct lead markers across the gene's traits.
#   N_lead_intervals        distinct 250 kb intervals those lead markers fall
#                           in, which is what says whether a gene's support is
#                           one LD block or several.
#
# Input   <DIR>/{CTL,LIN}_gene_annotation.tsv
# Output  <DIR>/{CTL,LIN}_gene_counts.tsv
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })
DIR <- Sys.getenv("ANN_DIR", "data/gene_annotation_final")
INTERVAL <- 250000

for (cond in c("CTL", "LIN")) {
  f <- file.path(DIR, sprintf("%s_gene_annotation.tsv", cond))
  if (!file.exists(f)) { message("missing, skipped: ", f); next }
  a <- fread(f)
  cnt <- a[, .(
    N_Traits               = uniqueN(Trait),
    N_Traits_indiv         = uniqueN(Trait[grepl("individual", Layer)]),
    N_Traits_sumrat        = uniqueN(Trait[grepl("sum_ratio",  Layer)]),
    Max_sig_SNPs_per_trait = max(N_sig_SNPs),
    Sum_sig_SNPs_per_trait = sum(N_sig_SNPs),
    N_lead_SNPs            = uniqueN(Lead_SNP),
    N_lead_intervals       = uniqueN(paste(Lead_Chr, Lead_BP %/% INTERVAL)),
    Best_P                 = min(Best_P),
    Max_r2                 = max(Max_r2),
    Traits                 = paste(sort(unique(Trait)), collapse = " | ")
  ), by = .(Condition, GeneID, Gene_Chr, Gene_Start, Gene_End)]
  setorder(cnt, -N_Traits, Best_P)
  out <- file.path(DIR, sprintf("%s_gene_counts.tsv", cond))
  fwrite(cnt, out, sep = "\t")
  message(sprintf("%s: %d genes -> %s", cond, nrow(cnt), basename(out)))
}

cat("\n-- the two genes the manuscript quotes --\n")
k <- rbind(fread(file.path(DIR, "CTL_gene_counts.tsv")), fread(file.path(DIR, "LIN_gene_counts.tsv")))
print(k[GeneID %in% c("SORBI_3003G088800","SORBI_3003G088900","SORBI_3003G089000","SORBI_3006G244000"),
        .(Condition, GeneID, N_Traits, Max_sig_SNPs_per_trait, Sum_sig_SNPs_per_trait,
          N_lead_SNPs, N_lead_intervals, Best_P)], row.names = FALSE)
cat("\n-- how badly the old column overcounted --\n")
print(k[, .(genes = .N,
            median_inflation = round(median(Sum_sig_SNPs_per_trait / Max_sig_SNPs_per_trait), 2),
            worst_inflation  = round(max(Sum_sig_SNPs_per_trait / Max_sig_SNPs_per_trait), 1)),
        by = Condition])
