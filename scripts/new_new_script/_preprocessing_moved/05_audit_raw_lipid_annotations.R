#!/usr/bin/env Rscript

# Audit final lipid labels against their original scan-level library annotations.
# This creates review tables only; it never overwrites analysis input matrices.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

project_root <- "/Users/nirwantandukar/Documents/Github/SoLD_paper"
annotation_root <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SAP_lipids_GWAS/data"
raw_root <- file.path(project_root, "data", "raw_lipid_intensities")
output_root <- file.path(raw_root, "annotation_audit")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
corrected_root <- file.path(raw_root, "corrected_lipid_names")
dir.create(corrected_root, recursive = TRUE, showWarnings = FALSE)
name_map_root <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SAP_lipids_GWAS/data/lipid_match"

# Normalize names only for comparison and reporting. The original name is retained.
normalize_lipid_name <- function(x) {
  x %>%
    str_trim() %>%
    str_replace("^Spectral Match to ", "") %>%
    str_replace(" from NIST14$", "") %>%
    str_replace("^Massbank[^ ]*:[^ ]*\\s+", "") %>%
    str_replace("^Massbank[^ ]*\\s+", "") %>%
    str_replace("^TAG\\s*", "TG") %>%
    str_replace("^DAG\\s*", "DG") %>%
    str_replace(";.*$", "") %>%
    str_replace(" - .*?$", "") %>%
    str_replace_all('"', "") %>%
    str_replace_all("_", "/") %>%
    str_replace_all("\\s+", "")
}

display_lipid_name <- function(x) {
  x %>%
    str_trim() %>%
    str_replace("^Spectral Match to ", "") %>%
    str_replace(" from NIST14$", "") %>%
    str_replace("^Massbank[^ ]*:[^ ]*\\s+", "") %>%
    str_replace("^Massbank[^ ]*\\s+", "") %>%
    str_replace("^TAG\\s*", "TG") %>%
    str_replace("^DAG\\s*", "DG") %>%
    str_replace(";.*$", "") %>%
    str_replace(" - .*?$", "") %>%
    str_replace_all('"', "") %>%
    str_replace_all("_", "/") %>%
    str_squish()
}

# The best source-library identity is selected once per raw scan, before names
# are compared to final labels. This avoids the original many-to-many join issue.
scan_lookup <- function(annotation) {
  annotation %>%
    mutate(
      scan_id = as.character(X.Scan.),
      source_name_raw = Compound_Name,
      source_name_display = display_lipid_name(Compound_Name),
      source_name_normalized = normalize_lipid_name(Compound_Name)
    ) %>%
    group_by(scan_id) %>%
    arrange(desc(MQScore), desc(SharedPeaks), MZErrorPPM, .by_group = TRUE) %>%
    summarise(
      source_name_best = first(source_name_raw),
      source_name_best_display = first(source_name_display),
      source_name_best_normalized = first(source_name_normalized),
      source_precursor_mz = first(Precursor_MZ),
      source_rt_seconds = first(RT_Query),
      source_mq_score = first(MQScore),
      source_shared_peaks = first(SharedPeaks),
      source_mz_error_ppm = first(MZErrorPPM),
      source_candidate_names = paste(sort(unique(source_name_normalized)), collapse = " | "),
      .groups = "drop"
    )
}

class_lookup <- read_csv(
  file.path(project_root, "data", "lipid_class", "final_lipid_classes.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    source_name_best_normalized = normalize_lipid_name(Lipids),
    corrected_class = Class,
    corrected_subclass = SubClass,
    corrected_sub_subclass = Sub_subclass
  ) %>%
  distinct(source_name_best_normalized, .keep_all = TRUE)

# Curated aliases standardize nomenclature such as full systematic PC/PG/PS/TG
# names to the shorthand used in the analysis. Entries with conflicting aliases
# are excluded and reported rather than resolved arbitrarily.
read_curated_aliases <- function(file_name) {
  read_csv(
    file.path(name_map_root, file_name),
    show_col_types = FALSE,
    name_repair = "minimal"
  ) %>%
    select(LipidsMatch, OtherName) %>%
    filter(!is.na(OtherName), OtherName != "", str_to_lower(OtherName) != "empty")
}

curated_aliases_raw <- bind_rows(
  read_curated_aliases("name_change_new.csv"),
  read_curated_aliases("remaining_lipid_name.csv")
)
ambiguous_alias_sources <- curated_aliases_raw %>%
  distinct(LipidsMatch, OtherName) %>%
  count(LipidsMatch, name = "n_aliases") %>%
  filter(n_aliases > 1) %>%
  pull(LipidsMatch)
curated_aliases <- curated_aliases_raw %>%
  filter(!LipidsMatch %in% ambiguous_alias_sources) %>%
  distinct(LipidsMatch, .keep_all = TRUE) %>%
  rename(manual_alias = OtherName)

audit_one_condition <- function(final_file, annotation_file, condition) {
  final <- read_csv(final_file, show_col_types = FALSE)
  annotation <- read.csv(annotation_file, stringsAsFactors = FALSE)

  final %>%
    mutate(
      final_row = row_number(),
      scan_id = as.character(X.Scan.),
      existing_name_normalized = normalize_lipid_name(Compound_Name)
    ) %>%
    left_join(scan_lookup(annotation), by = "scan_id") %>%
    left_join(
      curated_aliases %>% select(source_name_best = LipidsMatch, manual_alias_from_source = manual_alias),
      by = "source_name_best"
    ) %>%
    left_join(
      curated_aliases %>% select(Compound_Name = LipidsMatch, manual_alias_from_existing = manual_alias),
      by = "Compound_Name"
    ) %>%
    mutate(
      source_contains_existing_name = str_detect(
        paste0(" | ", source_candidate_names, " | "),
        fixed(paste0(" | ", existing_name_normalized, " | "))
      ),
      annotation_status = case_when(
        is.na(scan_id) ~ "untraceable_missing_scan",
        is.na(source_name_best) ~ "untraceable_scan_absent_from_annotation",
        source_contains_existing_name ~ "verified_label_matches_source_scan",
        TRUE ~ "corrected_label_from_source_scan"
      ),
      # A current name is not allowed to override a contradictory raw scan.
      manual_alias_from_existing = if_else(
        annotation_status == "corrected_label_from_source_scan",
        NA_character_,
        manual_alias_from_existing
      ),
      corrected_compound_name = case_when(
        annotation_status == "corrected_label_from_source_scan" ~ source_name_best_display,
        TRUE ~ Compound_Name
      ),
      corrected_compound_name = coalesce(
        manual_alias_from_source,
        manual_alias_from_existing,
        corrected_compound_name
      ),
      corrected_name_source = case_when(
        !is.na(manual_alias_from_source) ~ "curated_alias_from_source_name",
        !is.na(manual_alias_from_existing) ~ "curated_alias_from_existing_name",
        annotation_status == "corrected_label_from_source_scan" ~ "source_scan_best_annotation",
        TRUE ~ "existing_name_retained"
      ),
      corrected_name_normalized = normalize_lipid_name(corrected_compound_name),
      condition = condition
    ) %>%
    left_join(
      class_lookup %>% rename(corrected_name_normalized = source_name_best_normalized),
      by = "corrected_name_normalized"
    ) %>%
    mutate(
      corrected_class = case_when(
        !is.na(corrected_class) ~ corrected_class,
        TRUE ~ Class
      ),
      corrected_subclass = case_when(
        !is.na(corrected_subclass) ~ corrected_subclass,
        TRUE ~ SubClass
      ),
      corrected_sub_subclass = case_when(
        !is.na(corrected_sub_subclass) ~ corrected_sub_subclass,
        TRUE ~ Sub_subclass
      )
    )
}

ctl <- audit_one_condition(
  file.path(raw_root, "A_final_lipids.csv"),
  file.path(annotation_root, "lipid_class_A.csv"),
  "CTL"
)
lin <- audit_one_condition(
  file.path(raw_root, "B_final_lipids.csv"),
  file.path(annotation_root, "lipid_class_B.csv"),
  "LIN"
)

write_csv(ctl, file.path(output_root, "A_final_lipids_scan_annotation_audit.csv"))
write_csv(lin, file.path(output_root, "B_final_lipids_scan_annotation_audit.csv"))

# Full copies retain every original intensity column. Only rows with a direct
# source-scan conflict receive a replacement label/class; untraceable rows are
# deliberately kept and flagged rather than silently changed.
write_csv(
  ctl %>%
    mutate(
      Compound_Name = corrected_compound_name,
      Class = corrected_class,
      SubClass = corrected_subclass,
      Sub_subclass = corrected_sub_subclass
    ),
  file.path(corrected_root, "1_CTL_lipids_name_corrected.csv")
)
write_csv(
  lin %>%
    mutate(
      Compound_Name = corrected_compound_name,
      Class = corrected_class,
      SubClass = corrected_subclass,
      Sub_subclass = corrected_sub_subclass
    ),
  file.path(corrected_root, "2_LIN_lipids_name_corrected.csv")
)

# These rows need a separate raw-feature decision before any reaggregation.
# They are not summed or averaged by this audit script.
write_csv(
  lin %>%
    mutate(Compound_Name = corrected_compound_name) %>%
    add_count(Compound_Name, name = "n_rows_with_corrected_name") %>%
    filter(n_rows_with_corrected_name > 1) %>%
    select(
      Compound_Name, n_rows_with_corrected_name, final_row, X.Scan.,
      annotation_status, source_name_best, source_precursor_mz,
      source_rt_seconds, source_mq_score
    ) %>%
    arrange(Compound_Name, X.Scan.),
  file.path(output_root, "B_corrected_name_duplicate_groups_review.csv")
)
write_csv(
  bind_rows(ctl, lin) %>%
    filter(annotation_status == "corrected_label_from_source_scan") %>%
    select(
      condition, final_row, Compound_Name, corrected_compound_name,
      Class, corrected_class, SubClass, corrected_subclass,
      X.Scan., source_name_best, source_precursor_mz, source_rt_seconds,
      source_mq_score, source_shared_peaks, source_mz_error_ppm,
      source_candidate_names
    ),
  file.path(output_root, "final_lipid_label_corrections_required.csv")
)

summary <- bind_rows(ctl, lin) %>%
  count(condition, annotation_status, corrected_name_source, name = "n_rows")
write_csv(summary, file.path(output_root, "final_lipid_scan_annotation_audit_summary.csv"))
write_csv(
  curated_aliases %>% select(LipidsMatch, manual_alias),
  file.path(output_root, "curated_lipid_name_mappings_applied.csv")
)
write_csv(
  curated_aliases_raw %>% filter(LipidsMatch %in% ambiguous_alias_sources),
  file.path(output_root, "curated_lipid_name_mappings_excluded_ambiguous.csv")
)

print(summary)
cat("\nOutputs written to:", output_root, "\n")
