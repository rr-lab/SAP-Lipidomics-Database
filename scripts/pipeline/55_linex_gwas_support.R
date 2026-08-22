#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
})

dir.create("table/supp", recursive = TRUE, showWarnings = FALSE)

message("Building LINEX-GWAS support tables...")

gwas_dir <- "/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database/table/supp"
ann_control <- "/Users/nirwantandukar/Downloads/control.txt"
ann_lowinput <- "/Users/nirwantandukar/Downloads/lowinput.txt"

files <- list(
  S8_control = file.path(gwas_dir, "SuppTable_S8_GWAS_candidate_genes_for_individual_lipid_traits_control.tsv"),
  S8_lowinput = file.path(gwas_dir, "SuppTable_S8_GWAS_candidate_genes_for_individual_lipid_traits_lowinput.tsv"),
  S9_control = file.path(gwas_dir, "SuppTable_S9_GWAS_candidate_genes_for_lipid_sum_ratio_traits_control.tsv"),
  S9_lowinput = file.path(gwas_dir, "SuppTable_S9_GWAS_candidate_genes_for_lipid_sum_ratio_traits_lowinput.tsv")
)

missing_files <- names(files)[!file.exists(unlist(files))]
if (length(missing_files) > 0) {
  stop("Missing GWAS input files: ", paste(missing_files, collapse = ", "))
}

gene_map <- tibble::tribble(
  ~GeneID, ~Reaction_branch, ~Role_label,
  "SORBI_3001G103800", "LCAT* (PC+DG <-> LPC+TG)", "LCAT-like 1",
  "SORBI_3006G214500", "LCAT* (PC+DG <-> LPC+TG)", "LCAT-like 4",
  "SORBI_3001G448800", "LCAT* (PC+DG <-> LPC+TG)", "LCAT-like 1 paralog",
  "SORBI_3001G041900", "PNPLA1 branch (TG <-> DG)", "PNPLA-domain TAG lipase-like",
  "SORBI_3010G270700", "PNPLA1 branch (TG <-> DG)", "GDSL esterase/lipase-like",
  "SORBI_3010G170000", "TG-DG synthesis arm", "DGAT1",
  "SORBI_3009G034600", "TG-DG synthesis arm", "DGAT-like"
)

branch_classes <- list(
  "LCAT* (PC+DG <-> LPC+TG)" = list(substrates = c("PC", "DG"), products = c("LPC", "TG")),
  "PNPLA1 branch (TG <-> DG)" = list(substrates = c("TG"), products = c("DG")),
  "TG-DG synthesis arm" = list(substrates = c("DG"), products = c("TG"))
)

all_classes <- c("TG", "DG", "MG", "PC", "PE", "DGDG", "MGDG", "SQDG", "LPC", "LPE", "PG", "PA", "PS", "FA", "Cer", "SM", "GalCer", "AEG")

parse_ratio_classes <- function(x) {
  # e.g., Sum_DG_over_TG_log10safe -> c("DG", "TG")
  m <- str_match(x, "^Sum_([A-Za-z0-9]+)_over_([A-Za-z0-9]+)_")
  if (!is.na(m[1, 1])) return(c(m[1, 2], m[1, 3]))
  character(0)
}

parse_single_class <- function(x) {
  # e.g., TG(18:1_18:2_18:3) -> TG
  m <- str_match(x, "^([A-Za-z0-9]+)\\(")
  if (!is.na(m[1, 2])) return(m[1, 2])
  # class-only names (rare)
  if (x %in% all_classes) return(x)
  character(0)
}

is_linked <- function(pheno, substrates, products) {
  ratio_cls <- parse_ratio_classes(pheno)
  union_cls <- unique(c(substrates, products))
  if (length(ratio_cls) == 2) return(any(ratio_cls %in% union_cls))
  single_cls <- parse_single_class(pheno)
  if (length(single_cls) == 1) return(single_cls %in% union_cls)
  FALSE
}

is_direct_ratio <- function(pheno, substrates, products) {
  ratio_cls <- parse_ratio_classes(pheno)
  if (length(ratio_cls) != 2) return(FALSE)
  a <- ratio_cls[1]
  b <- ratio_cls[2]
  (a %in% substrates && b %in% products) || (a %in% products && b %in% substrates)
}

process_gwas_file <- function(path, source_name) {
  df <- read_tsv(path, show_col_types = FALSE)
  df %>%
    mutate(
      Source = source_name,
      Condition = if_else(str_detect(source_name, "lowinput"), "LIN", "CTL"),
      GWAS_type = if_else(str_detect(source_name, "^S8"), "Individual", "SumRatio")
    ) %>%
    select(GeneID, Best_SNP, Best_P_Value, N_Phenotypes, Phenotypes, Source, Condition, GWAS_type)
}

gwas_all <- bind_rows(
  process_gwas_file(files$S8_control, "S8_control"),
  process_gwas_file(files$S8_lowinput, "S8_lowinput"),
  process_gwas_file(files$S9_control, "S9_control"),
  process_gwas_file(files$S9_lowinput, "S9_lowinput")
)

gwas_sel <- gwas_all %>%
  inner_join(gene_map, by = "GeneID")

if (nrow(gwas_sel) == 0) {
  stop("No LINEX-linked genes found in GWAS tables.")
}

scored <- gwas_sel %>%
  rowwise() %>%
  mutate(
    phenotype_list = list(str_split(Phenotypes %||% "", ";", simplify = FALSE)[[1]] %>% str_trim() %>% discard(~ .x == "")),
    branch_sub = list(branch_classes[[Reaction_branch]]$substrates),
    branch_pro = list(branch_classes[[Reaction_branch]]$products),
    N_linked = sum(map_lgl(phenotype_list, ~ is_linked(.x, branch_sub, branch_pro))),
    N_direct_ratio = sum(map_lgl(phenotype_list, ~ is_direct_ratio(.x, branch_sub, branch_pro))),
    linked_examples = {
      v <- phenotype_list[map_lgl(phenotype_list, ~ is_linked(.x, branch_sub, branch_pro))]
      paste(head(v, 5), collapse = "; ")
    },
    direct_ratio_examples = {
      v <- phenotype_list[map_lgl(phenotype_list, ~ is_direct_ratio(.x, branch_sub, branch_pro))]
      paste(head(v, 5), collapse = "; ")
    }
  ) %>%
  ungroup() %>%
  mutate(
    linked_fraction = if_else(N_Phenotypes > 0, N_linked / N_Phenotypes, NA_real_),
    direct_ratio_fraction = if_else(N_Phenotypes > 0, N_direct_ratio / N_Phenotypes, NA_real_)
  ) %>%
  arrange(Reaction_branch, desc(Condition), Best_P_Value)

read_annotation <- function(path, source_label) {
  if (!file.exists(path)) {
    return(tibble(GeneID = character(), annotation = character(), go_bp = character(), go_mf = character(), source = character()))
  }
  raw <- read_tsv(path, col_names = FALSE, show_col_types = FALSE, progress = FALSE)
  # expected columns from user file: V2=GeneID, V3=annotation/free text, V8/V9 GO fields (if present)
  out <- raw %>%
    transmute(
      GeneID = as.character(X2),
      annotation = as.character(X3),
      go_mf = as.character(if ("X8" %in% names(raw)) X8 else NA_character_),
      go_bp = as.character(if ("X9" %in% names(raw)) X9 else NA_character_),
      source = source_label
    ) %>%
    filter(!is.na(GeneID), GeneID != "") %>%
    distinct(GeneID, .keep_all = TRUE)
  out
}

ann <- bind_rows(
  read_annotation(ann_control, "control.txt"),
  read_annotation(ann_lowinput, "lowinput.txt")
) %>%
  group_by(GeneID) %>%
  summarise(
    annotation = first(na.omit(annotation)),
    go_mf = first(na.omit(go_mf)),
    go_bp = first(na.omit(go_bp)),
    .groups = "drop"
  )

out_gene <- scored %>%
  left_join(ann, by = "GeneID") %>%
  transmute(
    Reaction_branch,
    GeneID,
    Role_label,
    annotation,
    Condition,
    GWAS_type,
    Source,
    Best_SNP,
    Best_P_Value,
    N_Phenotypes,
    N_linked,
    linked_fraction = round(linked_fraction, 3),
    N_direct_ratio,
    direct_ratio_fraction = round(direct_ratio_fraction, 3),
    linked_examples,
    direct_ratio_examples
  ) %>%
  arrange(Reaction_branch, desc(Condition), Best_P_Value)

out_branch <- out_gene %>%
  group_by(Reaction_branch, Condition) %>%
  summarise(
    N_genes_with_hits = n_distinct(GeneID),
    Total_N_Phenotypes = sum(N_Phenotypes, na.rm = TRUE),
    Total_N_linked = sum(N_linked, na.rm = TRUE),
    Total_N_direct_ratio = sum(N_direct_ratio, na.rm = TRUE),
    median_linked_fraction = median(linked_fraction, na.rm = TRUE),
    min_best_p = min(Best_P_Value, na.rm = TRUE),
    Genes = paste(unique(GeneID), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(Reaction_branch, desc(Condition))

out_gene_env <- out_gene %>%
  group_by(Reaction_branch, GeneID, Role_label, annotation, Condition) %>%
  summarise(
    N_rows = n(),
    N_Phenotypes_total = sum(N_Phenotypes, na.rm = TRUE),
    N_linked_total = sum(N_linked, na.rm = TRUE),
    N_direct_ratio_total = sum(N_direct_ratio, na.rm = TRUE),
    best_p = min(Best_P_Value, na.rm = TRUE),
    best_snp = Best_SNP[which.min(Best_P_Value)],
    .groups = "drop"
  ) %>%
  mutate(
    linked_fraction_total = if_else(N_Phenotypes_total > 0, N_linked_total / N_Phenotypes_total, NA_real_),
    direct_ratio_fraction_total = if_else(N_Phenotypes_total > 0, N_direct_ratio_total / N_Phenotypes_total, NA_real_)
  ) %>%
  pivot_wider(
    names_from = Condition,
    values_from = c(N_rows, N_Phenotypes_total, N_linked_total, N_direct_ratio_total, linked_fraction_total, direct_ratio_fraction_total, best_p, best_snp),
    names_glue = "{.value}_{Condition}"
  ) %>%
  mutate(
    across(starts_with("N_"), ~ replace_na(.x, 0)),
    across(starts_with("linked_fraction_total_"), ~ replace_na(.x, 0)),
    across(starts_with("direct_ratio_fraction_total_"), ~ replace_na(.x, 0))
  ) %>%
  arrange(Reaction_branch, best_p_LIN, best_p_CTL)

write_csv(out_gene, "table/supp/SuppTable_tmp_LINEX_GWAS_GeneSupport.csv", na = "")
write_csv(out_branch, "table/supp/SuppTable_tmp_LINEX_GWAS_BranchSummary.csv", na = "")
write_csv(out_gene_env, "table/supp/SuppTable_tmp_LINEX_GWAS_GeneSummary_CTL_vs_LIN.csv", na = "")

message("Saved:")
message(" - table/supp/SuppTable_tmp_LINEX_GWAS_GeneSupport.csv")
message(" - table/supp/SuppTable_tmp_LINEX_GWAS_BranchSummary.csv")
message(" - table/supp/SuppTable_tmp_LINEX_GWAS_GeneSummary_CTL_vs_LIN.csv")
