# =============================================================================
# Build the four gene-annotation tables from the work files written by
# 11_LD_annotate_genes.sh.
#
#   <OUT_DIR>/CTL_gene_annotation.tsv   one row per trait x gene
#   <OUT_DIR>/LIN_gene_annotation.tsv
#   <OUT_DIR>/CTL_gene_counts.tsv       one row per gene, recurrence across traits
#   <OUT_DIR>/LIN_gene_counts.tsv
#
# Both layers, individual lipids and class sums/ratios, go into the same file
# and are told apart by the Layer column, so a gene's recurrence can be read
# either overall or within a layer.
#
# Trait names are joined with " | " and never with ";", because some lipid
# species carry an internal semicolon, for example SPB 18:0;2OH, and a
# semicolon separator would silently split one trait into two.
#
# Usage   Rscript 12_make_annotation_tables.R
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })

BLUP    <- Sys.getenv("BLUP",    "/rsstu/users/r/rrellan/DOE_CAREER/SAP/results/spats_corrected/BLUP")
OUT_DIR <- Sys.getenv("OUT_DIR", file.path(BLUP, "gene_annotation_final"))
WORK    <- Sys.getenv("WORK",    file.path(OUT_DIR, "work"))

need <- file.path(WORK, c("trait_snp.tsv", "sig_snp_map.tsv", "snp_gene_LD.tsv"))
miss <- need[!file.exists(need)]
if (length(miss)) stop("run 11_LD_annotate_genes.sh first, missing: ", paste(basename(miss), collapse = ", "))

ts    <- fread(file.path(WORK, "trait_snp.tsv"), showProgress = FALSE)
map   <- fread(file.path(WORK, "sig_snp_map.tsv"), col.names = c("SNP","Chr","BP"), showProgress = FALSE)
sg    <- fread(file.path(WORK, "snp_gene_LD.tsv"), col.names = c("SNP","GeneID","R2"), showProgress = FALSE)
genes <- fread(file.path(OUT_DIR, "genes.range"),
               col.names = c("Gene_Chr","Gene_Start","Gene_End","GeneID"), showProgress = FALSE)

message("significant trait-SNP rows : ", nrow(ts))
message("resolved SNP positions     : ", nrow(map))
message("SNP-gene pairs             : ", nrow(sg))
message("genes in range file        : ", nrow(genes))

# ---- join -------------------------------------------------------------------
setkey(map, Chr, BP)
x <- merge(ts, map, by = c("Chr","BP"), all.x = TRUE)
unres <- sum(is.na(x$SNP))
if (unres) message("NOTE: ", unres, " significant rows have no bfile marker at that position and are dropped")
x <- x[!is.na(SNP)]

x <- merge(x, sg, by = "SNP", allow.cartesian = TRUE)
if (!nrow(x)) stop("no significant SNP mapped to a gene, check genes.range and the LD step")
x <- merge(x, genes, by = "GeneID", all.x = TRUE)

# ---- annotation, one row per condition x layer x trait x gene ---------------
ann <- x[, .(
  N_sig_SNPs = uniqueN(SNP),          # distinct markers for this trait x gene
  N_intervals = uniqueN(paste(Chr, BP %/% 250000)),
  Best_P     = min(P),
  Lead_SNP   = SNP[which.min(P)],
  Lead_Chr   = Chr[which.min(P)],
  Lead_BP    = BP[which.min(P)],
  Max_r2     = max(R2)
), by = .(Condition, Layer, Trait, GeneID, Gene_Chr, Gene_Start, Gene_End)]
setorder(ann, Condition, Layer, Best_P, Trait, GeneID)

# ---- counts, one row per condition x gene ------------------------------------
# Built from x, the row-level join, NOT from ann. Rolling ann up would have to
# sum its per-trait N_sig_SNPs, and the same SNP is significant for many traits,
# so that sum counts one marker once per trait it hits. For the chromosome 3
# block that turns 2 markers into 313. uniqueN over the SNP column is the
# number the text means when it says "across N significant SNPs".
cnt <- x[, .(
  N_Traits         = uniqueN(Trait),
  N_Traits_indiv   = uniqueN(Trait[grepl("individual", Layer)]),
  N_Traits_sumrat  = uniqueN(Trait[grepl("sum_ratio",  Layer)]),
  N_sig_SNPs       = uniqueN(SNP),
  N_intervals      = uniqueN(paste(Chr, BP %/% 250000)),
  Best_P           = min(P),
  Max_r2           = max(R2),
  Traits           = paste(sort(unique(Trait)), collapse = " | ")
), by = .(Condition, GeneID, Gene_Chr, Gene_Start, Gene_End)]
setorder(cnt, Condition, -N_Traits, Best_P)

# ---- write -------------------------------------------------------------------
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
for (cond in c("CTL","LIN")) {
  a <- ann[Condition == cond]; k <- cnt[Condition == cond]
  if (!nrow(a)) { message("no results for ", cond); next }
  fwrite(a, file.path(OUT_DIR, sprintf("%s_gene_annotation.tsv", cond)), sep = "\t")
  fwrite(k, file.path(OUT_DIR, sprintf("%s_gene_counts.tsv",     cond)), sep = "\t")
  message(sprintf("%s  annotation %d rows   counts %d genes", cond, nrow(a), nrow(k)))
}

# ---- summary -----------------------------------------------------------------
cat("\n-- candidate genes by condition and layer --\n")
print(ann[, .(traits = uniqueN(Trait), genes = uniqueN(GeneID),
              trait_gene_pairs = .N), by = .(Condition, Layer)])

cat("\n-- genes per condition, any layer --\n")
print(cnt[, .(genes = .N,
              genes_in_1_trait  = sum(N_Traits == 1),
              genes_in_5plus    = sum(N_Traits >= 5),
              max_traits_1_gene = max(N_Traits)), by = Condition])

cat("\n-- distinct significant markers per gene, sanity --\n")
print(head(cnt[order(-N_sig_SNPs), .(Condition, GeneID, N_Traits, N_sig_SNPs, N_intervals, Best_P)], 8), row.names = FALSE)

cat("\n-- top recurrent genes --\n")
for (cond in c("CTL","LIN")) {
  k <- cnt[Condition == cond]
  if (!nrow(k)) next
  print(head(k[, .(Condition, GeneID, Gene_Chr, Gene_Start, N_Traits, Best_P, Max_r2)], 10),
        row.names = FALSE)
}

cat("\nSaved to", OUT_DIR, "\n")
