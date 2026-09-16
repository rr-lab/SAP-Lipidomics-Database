# ==============================================================================
# Supplementary Tables S30 and S31 -- the CTL/LIN candidate-gene overlap broken
# down by lipid class, and the shared genes ranked across both trait layers.
#
#   Rscript scripts/new_new_script/14_SuppTableS30S31_overlap_by_class_and_shared_ranked.R
#
# These two were computed in Python during the first pass
# (scripts/chapter2_addons/gwas_overlap.py) and never carried a supplementary
# number. This script re-derives both in R from the LD candidate master table,
# and checks the result against the Python output still on disk under
# table/overlap/, so the numbers in the paper and the numbers in the repository
# cannot drift apart silently.
#
# THE CLASS RULE, stated once.
#
# A candidate gene belongs to a lipid class if any phenotype it was called for
# resolves to that class. Resolution differs by trait layer:
#
#   individual  the species name is parsed for its leading abbreviation, so
#               TG(16:0_18:1_18:2) is TG and AEG(o-16:2/16:0) is AEG. A species
#               whose name carries no abbreviation -- the carotenoids, sterols,
#               tocopherols and free sphingoid bases -- is looked up in
#               final_lipid_classes.csv and mapped through SUBCLASS_FALLBACK
#               first, then CLASS_MAP, then falls back to its own Class label.
#
#   sumratio    Sum_PC_log10safe is PC. A ratio, Sum_PA_over_PC_log10safe,
#               belongs to BOTH numerator and denominator, so a gene called for
#               it is counted under PA and under PC. That is deliberate -- the
#               ratio carries no information about which side moved.
#
# A gene therefore appears in as many class rows as it has classes. The rows are
# not a partition of the candidate set and the columns do not sum to the totals
# in Supplementary Table S20.
#
# READ THE SUM/RATIO ROWS WITH THAT IN MIND. There are only 18 shared genes in
# that layer and 10 of them are called for ratios spanning most classes, so the
# same 10 genes carry nearly every sum/ratio row. The 25 rows are not 25
# independent tests and the fold enrichments in that block are not comparable to
# the individual-lipid block, where a gene lands in one class.
#
# Classes with fewer than three candidates in either trial are dropped, because
# the hypergeometric tail is meaningless there. q_BH is computed within a layer,
# not across both.
#
# The gene universe N is every gene with coordinates in genes.range (34,027),
# not the annotated subset -- the same universe the overall overlap in S20 uses.
#
# WHAT THE CHECK FOUND, recorded so it is not rediscovered.
#
# S30 reproduces the Python output to 1e-10 -- print formatting, nothing more.
#
# S31 does not. For 15 of the 354 shared genes the Python file carries a
# tot_phen exactly one higher than this script computes. Every other column, and
# every other gene, matches to the last digit. The cause is the species
# deduplication of 2026-09-03: all 15 of those genes are called for a species the
# dedup KEPT, and none of the 339 matching genes are, so each lost exactly one
# phenotype when its duplicate partner was merged away. The Python output
# predates the dedup; this script reads the post-dedup master, so the numbers
# here are the current ones and the file under table/overlap/ is stale.
#
# The same staleness reaches Supplementary Table S22, which was copied from that
# Python run: 3 of its N_Phenotypes_CTL and 11 of its N_Phenotypes_LIN values are
# one too high. This script rewrites S22 and S23 from the master for that reason,
# and says so when it does. Gene membership and every p-value were unaffected.
#
# Inputs
#   data/LD_mapped/candidate_tables/ALL_LD_candidate_genes_master.tsv
#   data/metadata/final_lipid_classes.csv
#   data/LD_mapped/genes_ranges/genes.range
#   table/overlap/gwas_overlap_by_class.csv            (reference, for the check)
#   table/overlap/shared_genes_all_layers_ranked.csv   (reference, for the check)
#
# Outputs
#   table/supp/SuppTable_S30_Overlap_by_lipid_class.csv
#   table/supp/SuppTable_S31_Shared_genes_ranked.csv
#   table/supp/SuppTable_S22_Shared_genes_individual.csv   (rewritten if stale)
#   table/supp/SuppTable_S23_Shared_genes_sumratio.csv     (rewritten if stale)
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr); library(vroom) })

REPO     <- Sys.getenv("SOLD_REPO", ".")
MASTER   <- file.path(REPO, "data/LD_mapped/candidate_tables/ALL_LD_candidate_genes_master.tsv")
CLASSES  <- file.path(REPO, "data/metadata/final_lipid_classes.csv")
RANGES   <- file.path(REPO, "data/LD_mapped/genes_ranges/genes.range")
OVL_DIR  <- file.path(REPO, "table/overlap")
TAB_SUPP <- file.path(REPO, "table/supp")
dir.create(TAB_SUPP, recursive = TRUE, showWarnings = FALSE)

MIN_PER_TRIAL <- 3

# ---- class vocabulary --------------------------------------------------------
ABBR <- c("MGDG","DGDG","SQDG","GalCer","DGTS","AEG","LPC","LPE","Cer","SM",
          "DG","TG","MG","PC","PE","PG","PA","PS","FA")

SUBCLASS_FALLBACK <- c(
  "Fatty acid and conjugates" = "FA", "Fatty acid amide" = "FA",
  "Fatty amide" = "FA", "Fatty acid ester" = "FA", "Oxidised fatty acid" = "FA",
  "Triterpenoid" = "Terpenoid", "Tetraterpenoid" = "Terpenoid",
  "Diterpenoid" = "Terpenoid", "Monoterpenoid" = "Terpenoid",
  "Sesquiterpenoid" = "Terpenoid", "Prenol" = "Terpenoid",
  "Cyclic monoterpenoid" = "Terpenoid", "Sterol" = "Sterol",
  "Vitamin" = "Other", "Vitmain" = "Other", "Coenzyme" = "Other",
  "Other sphingolipid" = "SPB", "Ceramide" = "Cer",
  "Galactosylceramide" = "GalCer", "Cardiolipin" = "CL",
  "Glycerophosphocholine" = "PC")

# The annotation table moved to the LIPID MAPS categories on 2026-09-16, so the
# fallback keys are the category names. A species whose own name carries no class
# abbreviation is labelled by its category here.
# The keys are the LIPID MAPS category names the annotation table carries since
# 2026-09-16. The values reproduce the labels this table used before that change,
# so the only difference from the previous version is the removed phenotypes.
# Note FA and "Fatty acid" remain separate groups here: FA comes from a species
# name that starts with the FA abbreviation, "Fatty acid" from a species with no
# abbreviation whose category is Fatty Acyls. That split predates this script.
CLASS_MAP <- c("Fatty Acyls" = "Fatty acid", "Prenol Lipids" = "Terpenoid",
               "Sterol Lipids" = "Sterol", "Sphingolipids" = "SPB")

meta <- vroom(CLASSES, show_col_types = FALSE) %>%
  transmute(Lipids = trimws(Lipids), Class = trimws(Class), SubClass = trimws(SubClass))
meta_cls <- setNames(meta$Class,    meta$Lipids)
meta_sub <- setNames(meta$SubClass, meta$Lipids)

# ---- gene universe -----------------------------------------------------------
univ <- vroom(RANGES, delim = "\t", col_names = FALSE, show_col_types = FALSE)
N_UNIVERSE <- length(unique(univ[[4]]))

# ---- phenotype -> class ------------------------------------------------------
pheno_classes <- function(p, layer) {
  p <- trimws(p)
  if (!nzchar(p)) return(character(0))
  if (layer == "sumratio") {
    if (p %in% ABBR) return(p)
    m <- regmatches(p, regexec("^Sum_(.+?)_over_(.+?)_log10safe$", p))[[1]]
    if (length(m) == 3) return(c(m[2], m[3]))
    m <- regmatches(p, regexec("^Sum_(.+?)(_log10safe)?$", p))[[1]]
    if (length(m) >= 2) return(m[2])
    return(character(0))
  }
  m <- regmatches(p, regexec("^([A-Za-z0-9]+)\\(", p))[[1]]
  if (length(m) == 2 && m[2] %in% ABBR) return(m[2])
  if (p %in% names(meta_cls)) {
    sub <- meta_sub[[p]]; cls <- meta_cls[[p]]
    if (!is.na(sub) && sub %in% names(SUBCLASS_FALLBACK)) return(unname(SUBCLASS_FALLBACK[[sub]]))
    if (!is.na(cls) && cls %in% names(CLASS_MAP))         return(unname(CLASS_MAP[[cls]]))
    return(cls)
  }
  "Other"
}

# ---- the candidate master ----------------------------------------------------
mst <- vroom(MASTER, delim = "\t", show_col_types = FALSE,
             col_types = cols(.default = col_character(),
                              Best_P_Value = col_double(),
                              N_Phenotypes = col_integer()))

# SPLIT ON "; ", NEVER ON ";". One species is named SPB 18:0;2OH -- the
# semicolon is LIPID MAPS shorthand for the 2-hydroxy sphingoid base and is part
# of the name. Splitting on a bare ";" tears it in two and manufactures a "2OH"
# phenotype that resolves to no class. The Python original
# (scripts/chapter2_addons/gwas_overlap.py) does exactly that, which is where the
# phantom 2OH in the old unresolved-phenotype list came from. It happened to be
# harmless there because the fragment fell through to "Other", which is not a
# reported class row, but the next species with a semicolon will not be.
gene_classes <- Map(function(ph, ly) {
    unique(unlist(lapply(strsplit(ph, "; ", fixed = TRUE)[[1]], pheno_classes, layer = ly)))
  }, mst$Phenotypes, mst$layer)

# Nothing should reach the class resolver with a bare semicolon still in it.
stopifnot(!any(grepl(";", unlist(strsplit(mst$Phenotypes, "; ", fixed = TRUE)),
                     fixed = TRUE) &
               !grepl("^SPB ", unlist(strsplit(mst$Phenotypes, "; ", fixed = TRUE)))))

long <- tibble(condition = mst$condition, layer = mst$layer, gene = mst$GeneID,
               cls = gene_classes) %>%
  tidyr::unnest_longer(cls) %>%
  distinct(condition, layer, gene, cls)

# ---- S30: per-class overlap --------------------------------------------------
# hyper_sf(k, N, K, n) = P(X >= k); phyper's q is P(X > q), hence k - 1.
overlap_row <- function(A, B, N) {
  a <- length(A); b <- length(B); ov <- length(intersect(A, B))
  exp <- a * b / N
  tibble(n_CTL = a, n_LIN = b, n_shared = ov,
         CTL_only = length(setdiff(A, B)), LIN_only = length(setdiff(B, A)),
         expected_shared = exp,
         fold_enrichment = if (exp > 0) ov / exp else NA_real_,
         jaccard = if (a + b > 0) ov / length(union(A, B)) else NA_real_,
         p_hypergeom = if (a > 0 && b > 0) phyper(ov - 1, b, N - b, a, lower.tail = FALSE) else NA_real_)
}

per <- lapply(c("individual", "sumratio"), function(ly) {
  d <- long %>% filter(layer == ly)
  lapply(sort(unique(d$cls)), function(cc) {
    A <- d %>% filter(condition == "CTL", cls == cc) %>% pull(gene) %>% unique()
    B <- d %>% filter(condition == "LIN", cls == cc) %>% pull(gene) %>% unique()
    if (length(A) < MIN_PER_TRIAL || length(B) < MIN_PER_TRIAL) return(NULL)
    bind_cols(tibble(Layer = ly, Class = cc), overlap_row(A, B, N_UNIVERSE))
  }) %>% bind_rows()
}) %>% bind_rows() %>%
  group_by(Layer) %>% mutate(q_BH = p.adjust(p_hypergeom, "BH")) %>% ungroup() %>%
  arrange(Layer, desc(n_CTL)) %>%
  select(Layer, Class, n_CTL, n_LIN, n_shared, CTL_only, LIN_only,
         expected_shared, fold_enrichment, jaccard, p_hypergeom, q_BH)

# ---- S31: shared genes ranked ------------------------------------------------
shared <- intersect(mst$GeneID[mst$condition == "CTL"], mst$GeneID[mst$condition == "LIN"])
ranked <- mst %>% filter(GeneID %in% shared) %>%
  group_by(GeneID) %>%
  summarise(GeneName  = dplyr::first(GeneName),
            n_records = dplyr::n(),
            tot_phen  = sum(N_Phenotypes),
            best_p    = min(Best_P_Value), .groups = "drop") %>%
  arrange(best_p, desc(tot_phen)) %>%
  select(GeneID, GeneName, n_records, tot_phen, best_p)

# ---- check against the Python output still on disk ---------------------------
compare <- function(new, ref_path, keys, label, tol = 1e-6) {
  if (!file.exists(ref_path)) { message("  ", label, ": no reference, skipped"); return(invisible(NULL)) }
  ref <- vroom(ref_path, show_col_types = FALSE)
  j <- inner_join(new, ref, by = keys, suffix = c(".new", ".ref"))
  if (nrow(j) != nrow(new) || nrow(j) != nrow(ref)) {
    # The reference under table/overlap/ predates the removal of the two phantom
    # phenotypes SM(d18:1_17:0) and DG(18:0_18:2_0:0) on 2026-09-16, so a
    # membership difference here is expected for the class-level table. Report it
    # rather than stopping, and check the rows named below.
    #
    # SPB 18:0;2OH is NOT a phantom. The semicolon is part of the LIPID MAPS
    # shorthand for the 2-hydroxy sphingoid base, not a phenotype separator, and
    # any code that splits the Phenotypes column on ';' rather than '; ' will
    # manufacture a bare 2OH entry out of it. Split on '; '.
    message(sprintf("  %-26s membership differs: new %d, reference %d, matched %d",
                    label, nrow(new), nrow(ref), nrow(j)))
    only_new <- dplyr::anti_join(new, ref, by = keys)
    only_ref <- dplyr::anti_join(ref, new, by = keys)
    if (nrow(only_new)) message("       only in new: ",
        paste(apply(only_new[keys], 1, paste, collapse = "/"), collapse = ", "))
    if (nrow(only_ref)) message("       only in reference: ",
        paste(apply(only_ref[keys], 1, paste, collapse = "/"), collapse = ", "))
  }
  num <- intersect(names(new), names(ref)) |> setdiff(keys)
  num <- num[vapply(new[num], is.numeric, logical(1))]
  d <- vapply(num, function(v) {
         a <- j[[paste0(v, ".new")]]; b <- j[[paste0(v, ".ref")]]
         max(abs(a - b) / pmax(abs(b), 1e-300), na.rm = TRUE) }, numeric(1))
  bad <- names(d)[d > tol]
  if (!length(bad)) {
    message(sprintf("  %-26s %4d rows  identical", label, nrow(new)))
  } else {
    message(sprintf("  %-26s %4d rows  DIFFERS on %s", label, nrow(new), paste(bad, collapse = ", ")))
    for (v in bad) message(sprintf("       %s: %d of %d rows differ",
      v, sum(abs(j[[paste0(v, ".new")]] - j[[paste0(v, ".ref")]]) > tol * pmax(abs(j[[paste0(v, ".ref")]]), 1)), nrow(j)))
  }
  invisible(bad)
}

# ---- S22 / S23: the shared-gene tables already in the supplement --------------
# These are corrections, not new tables, so they are written back in the row order
# and the number formatting the submitted files already use -- the only cells that
# change are the stale phenotype counts.
shared_layer <- function(ly, ref_path) {
  x <- mst %>% filter(layer == ly) %>%
    group_by(GeneID, condition) %>%
    summarise(bp = min(Best_P_Value), np = sum(N_Phenotypes), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = condition, values_from = c(bp, np)) %>%
    filter(!is.na(bp_CTL), !is.na(bp_LIN)) %>%
    left_join(distinct(mst, GeneID, GeneName), by = "GeneID") %>%
    transmute(GeneID, GeneName,
              Best_P_Value_CTL = bp_CTL, Best_P_Value_LIN = bp_LIN,
              N_Phenotypes_CTL = np_CTL, N_Phenotypes_LIN = np_LIN)
  if (file.exists(ref_path)) {
    ord <- vroom(ref_path, show_col_types = FALSE)$GeneID
    stopifnot(setequal(ord, x$GeneID))
    x <- x[match(ord, x$GeneID), ]
  }
  # The submitted files were written by Python, which prints the shortest
  # mantissa that round-trips; %.6e then stripping trailing zeros reproduces it,
  # so the corrected file differs from the submitted one only in the stale counts.
  pyfmt <- function(v) sub("\\.e", "e", sub("0+e", "e", sprintf("%.6e", v)))
  x %>% mutate(across(starts_with("Best_P_Value"), pyfmt))
}
s22 <- shared_layer("individual", file.path(OVL_DIR, "shared_genes_individual.csv"))
s23 <- shared_layer("sumratio",   file.path(OVL_DIR, "shared_genes_sumratio.csv"))

message("Reproducibility check against table/overlap/")
compare(per,    file.path(OVL_DIR, "gwas_overlap_by_class.csv"),          c("Layer", "Class"), "S30 overlap by class")
compare(ranked, file.path(OVL_DIR, "shared_genes_all_layers_ranked.csv"), "GeneID",            "S31 shared genes ranked")
compare(mutate(s22, across(starts_with("Best_P_Value"), as.numeric)),
        file.path(OVL_DIR, "shared_genes_individual.csv"), "GeneID", "S22 shared individual")
compare(mutate(s23, across(starts_with("Best_P_Value"), as.numeric)),
        file.path(OVL_DIR, "shared_genes_sumratio.csv"),   "GeneID", "S23 shared sum/ratio")

readr::write_csv(per,    file.path(TAB_SUPP, "SuppTable_S30_Overlap_by_lipid_class.csv"))
readr::write_csv(ranked, file.path(TAB_SUPP, "SuppTable_S31_Shared_genes_ranked.csv"))
readr::write_csv(s22,    file.path(TAB_SUPP, "SuppTable_S22_Shared_genes_individual.csv"), na = "")
readr::write_csv(s23,    file.path(TAB_SUPP, "SuppTable_S23_Shared_genes_sumratio.csv"), na = "")

message(sprintf("\nGene universe N = %d", N_UNIVERSE))
message(sprintf("S30  %d class rows (%d individual, %d sum/ratio)", nrow(per),
                sum(per$Layer == "individual"), sum(per$Layer == "sumratio")))
message(sprintf("S31  %d shared genes", nrow(ranked)))
message(sprintf("S22  %d shared genes, individual layer", nrow(s22)))
message(sprintf("S23  %d shared genes, sum/ratio layer",  nrow(s23)))
