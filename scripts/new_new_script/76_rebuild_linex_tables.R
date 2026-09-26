# ==============================================================================
# Supplementary Tables S12, S13, S14 and S15 -- rebuilt from current data.
#
#   Rscript scripts/new_new_script/76_rebuild_linex_tables.R
#
# WHY. Workbook S9 (LINEX reaction mapping) and workbook S10 (frozen species
# set) were the last supplementary tables still carrying pre-2026 numbers.
#
#   S12 reaction balance, written 2026-05-04. Its deltas disagreed with the
#       manuscript after the lyso rename moved PC(18:2/0:0) and PC(22:0/0:0)
#       into LPC, which changes the LPC and PC pools and therefore the LCAT*
#       and LRO1 scores.
#           branch    S12 (old)   recomputed
#           PNPLA3      +0.327      +0.314
#           LRO1        +0.232      +0.261
#           LCAT*       +0.238      +0.139
#           PNPLA1      -0.268      -0.256
#
#   S13/S14 GWAS gene support, same date. They named SORBI_3001G103800 and
#       SORBI_3001G448800 for the LCAT* branch, neither of which is a candidate
#       at Max_r2 >= 0.4 any more, and did not contain SORBI_3004G341900 (PSAT),
#       which the manuscript now names. Six rows still used _log10safe trait
#       names that the phenotype rebuild replaced with _log10ratio.
#
#   S15 frozen species set, written 2025-07-25. It carried five chlorophyll and
#       pheophytin features the annotation has since dropped, lacked LPC(18:2)
#       and LPC(22:0), and had no Category_Code column. It is now a straight
#       copy of data/metadata/final_lipid_classes.csv, which is what every other
#       script reads.
#
# S11 (the four reactions and their RHEA identifiers) is a hand-curated list
# with no data dependency and is left alone.
#
# Inputs
#   data/SPATS_fitted/non_normalized_intensities/Final_subset_{control,lowinput}_*.csv
#   data/gene_annotation_final/{CTL,LIN}_gene_annotation.tsv, GO_terms_all.txt
#   data/metadata/final_lipid_classes.csv
# Outputs
#   table/supp/SuppTable_S12_ReactionBalance.csv
#   table/supp/SuppTable_S13_LINEX_GWAS_GeneSupport.csv
#   table/supp/SuppTable_S14_LINEX_GWAS_BranchSummary.csv
#   table/supp/SuppTable_S15_final_lipid_classes.csv
# ==============================================================================
source("scripts/new_new_script/_common.R")
suppressPackageStartupMessages({ library(data.table); library(stringr) })

ANN    <- file.path(DATA_ROOT, "gene_annotation_final")
R2_CUT <- 0.4
N_BOOT <- 2000

# ---- the four branches -------------------------------------------------------
# substrates and products as in S11; `classes` is every lipid class the branch
# touches, used to decide whether a GWAS trait is on-branch.
BR <- list(
  list(id="LCAT*",  label="LCAT* (PC+DG <-> LPC+TG)",
       score="mean(log10(LPC),log10(TG)) - mean(log10(PC),log10(DG))",
       num=c("LPC","TG"), den=c("PC","DG"), classes=c("PC","DG","LPC","TG"),
       fam="lecithin-cholesterol|phospholipid--sterol", famlab="LCAT/PSAT superfamily"),
  list(id="PNPLA1", label="PNPLA1 branch (TG <-> DG)",
       score="log10(DG) - log10(TG)",
       num="DG", den="TG", classes=c("TG","DG"),
       fam="pnpla|patatin", famlab="PNPLA/patatin lipase"),
  list(id="LRO1",   label="LRO1 branch (PE+DG <-> TG+LPE)",
       score="mean(log10(TG),log10(LPE)) - mean(log10(PE),log10(DG))",
       num=c("TG","LPE"), den=c("PE","DG"), classes=c("PE","DG","TG","LPE"),
       fam="phospholipid:diacylglycerol|phospholipid diacylglycerol", famlab="PDAT"),
  list(id="PNPLA3", label="PNPLA3 branch (DG+MG <-> TG)",
       score="log10(TG) - mean(log10(DG),log10(MG))",
       num="TG", den=c("DG","MG"), classes=c("DG","MG","TG"),
       fam="^o-acyltransferase|diacylglycerol acyltransferase", famlab="DGAT")
)

# ---- per-genotype class sums -------------------------------------------------
class_sums <- function(path) {
  x <- read_trial(path)
  m <- as.matrix(x[, -1, drop = FALSE]); storage.mode(m) <- "numeric"
  m[!is.finite(m)] <- 0
  colnames(m) <- normalize_lipid_name(colnames(m))
  cl <- sub("\\(.*$", "", colnames(m))
  s <- sapply(c("PC","PE","DG","TG","MG","LPC","LPE"), function(k) {
    j <- which(cl == k); if (length(j)) rowSums(m[, j, drop = FALSE]) else rep(NA_real_, nrow(m))
  })
  rownames(s) <- as.character(x[[1]]); s
}
lg <- function(v) log10(pmax(v, 1e-12))
score_of <- function(d, b) rowMeans(sapply(b$num, function(k) lg(d[,k])), na.rm=TRUE) -
                           rowMeans(sapply(b$den, function(k) lg(d[,k])), na.rm=TRUE)

A <- class_sums(CTL_CSV); B <- class_sums(LIN_CSV)
g <- intersect(rownames(A), rownames(B))
message("genotype pairs: ", length(g))
A <- A[g,]; B <- B[g,]

# ---- S12, reaction balance ---------------------------------------------------
set.seed(1207)
s12 <- rbindlist(lapply(BR, function(b) {
  a <- score_of(A,b); z <- score_of(B,b)
  d <- z - a
  w <- wilcox.test(z, a, paired = TRUE)
  bo <- replicate(N_BOOT, { i <- sample(length(d), replace=TRUE); median(z[i]) - median(a[i]) })
  data.table(Reaction_branch = b$id, Score_definition = b$score, N_pairs = length(g),
             median_CTL = median(a), median_LIN = median(z),
             delta_median = median(z) - median(a),
             test = "paired Wilcoxon signed-rank",
             p_value = w$p.value,
             direction = if (median(z) > median(a)) "products_enriched_in_LIN" else "substrates_enriched_in_LIN",
             ci_low_2.5pct = unname(quantile(bo, .025)), ci_high_97.5pct = unname(quantile(bo, .975)))
}))
s12[, p_adj_BH := p.adjust(p_value, "BH")]
setcolorder(s12, c("Reaction_branch","Score_definition","N_pairs","median_CTL","median_LIN",
                   "delta_median","test","p_value","p_adj_BH","direction",
                   "ci_low_2.5pct","ci_high_97.5pct"))
fwrite(s12, file.path(TAB_SUPP, "SuppTable_S12_ReactionBalance.csv"))
cat("\n-- S12 --\n"); print(s12[, .(Reaction_branch, median_CTL=round(median_CTL,4),
  median_LIN=round(median_LIN,4), delta_median=round(delta_median,4),
  p_adj_BH=signif(p_adj_BH,3), direction)])

# ---- S13/S14, GWAS gene support ----------------------------------------------
gw <- rbindlist(lapply(c("CTL_gene_annotation.tsv","LIN_gene_annotation.tsv"),
                       function(f) fread(file.path(ANN,f), sep="\t", colClasses="character")))
gw[, `:=`(Best_P=as.numeric(Best_P), Max_r2=as.numeric(Max_r2))]
gw <- gw[Max_r2 >= R2_CUT]
gw[, Layer := ifelse(grepl("individual", Layer), "Individual", "SumRatio")]
fn <- unique(fread(file.path(ANN,"GO_terms_all.txt"), sep="\t", colClasses="character"),
             by="GeneID")[, .(GeneID, GeneName)]

# A trait is on-branch if it is a species of a branch class, or a sum/ratio
# naming one. Ratio names carry both sides, so a ratio spanning the branch and
# an off-branch class still counts; that is stated in the manuscript.
on_branch <- function(trait, classes) {
  sp <- str_match(normalize_lipid_name(trait), "^([A-Za-z0-9]+)\\(")[,2]
  if (!is.na(sp)) return(sp %in% classes)
  toks <- unlist(str_split(gsub("^Sum_|_log10(safe|ratio)$", "", trait), "_over_|_"))
  any(toks %in% classes)
}

s13 <- rbindlist(lapply(BR, function(b) {
  ids <- fn[grepl(b$fam, GeneName, ignore.case = TRUE), GeneID]
  d <- gw[GeneID %in% ids]
  if (!nrow(d)) return(NULL)
  d[, linked := vapply(Trait, on_branch, logical(1), classes = b$classes)]
  d[, .(Reaction_branch = b$label, Gene_family = b$famlab,
        annotation = fn[GeneID == .BY$GeneID, GeneName][1],
        Best_SNP = Lead_SNP[which.min(Best_P)], Best_P_Value = min(Best_P),
        Max_r2 = max(Max_r2), N_Phenotypes = .N, N_linked = sum(linked),
        linked_fraction = round(sum(linked)/.N, 3),
        Phenotypes = paste(sort(Trait), collapse = "; "),
        linked_examples = paste(sort(Trait[linked]), collapse = "; ")),
    by = .(GeneID, Condition, GWAS_type = Layer)]
}), fill = TRUE)
setorder(s13, Reaction_branch, Best_P_Value)
setcolorder(s13, c("Reaction_branch","GeneID","Gene_family","annotation","Condition",
                   "GWAS_type","Best_SNP","Best_P_Value","Max_r2","N_Phenotypes",
                   "N_linked","linked_fraction","Phenotypes","linked_examples"))
fwrite(s13, file.path(TAB_SUPP, "SuppTable_S13_LINEX_GWAS_GeneSupport.csv"))

s14 <- s13[, .(N_genes_with_hits = uniqueN(GeneID), Total_N_Phenotypes = sum(N_Phenotypes),
               Total_N_linked = sum(N_linked),
               median_linked_fraction = round(median(linked_fraction), 3),
               min_best_p = min(Best_P_Value),
               Genes = paste(sort(unique(GeneID)), collapse = "; ")),
           by = .(Reaction_branch, Condition)]
# branches with no candidate are reported as such rather than omitted
miss <- setdiff(sapply(BR, `[[`, "label"), unique(s14$Reaction_branch))
if (length(miss)) s14 <- rbind(s14, data.table(Reaction_branch = miss, Condition = "both",
  N_genes_with_hits = 0L, Total_N_Phenotypes = 0L, Total_N_linked = 0L,
  median_linked_fraction = NA_real_, min_best_p = NA_real_, Genes = ""), fill = TRUE)
setorder(s14, Reaction_branch, Condition)
fwrite(s14, file.path(TAB_SUPP, "SuppTable_S14_LINEX_GWAS_BranchSummary.csv"))
cat("\n-- S14 --\n"); print(s14[, .(Reaction_branch, Condition, N_genes_with_hits,
  Total_N_Phenotypes, Total_N_linked, min_best_p = signif(min_best_p,3))])

# ---- S15, frozen species set -------------------------------------------------
cls <- fread(file.path(DATA_ROOT, "metadata/final_lipid_classes.csv"))
fwrite(cls, file.path(TAB_SUPP, "SuppTable_S15_final_lipid_classes.csv"))
message("\nS15: ", nrow(cls), " species, ", ncol(cls), " columns")
