#!/usr/bin/env Rscript

# Map normalized SERRF feature IDs to their highest-scoring library annotation,
# retain only curated non-Other lipid classes, and reaggregate identical names.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

project_root <- "/Users/nirwantandukar/Documents/Github/SoLD_paper"
data_root <- "/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database/data/new_data"
annotation_root <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SAP_lipids_GWAS/data"
name_map_root <- file.path(annotation_root, "lipid_match")

settings <- tribble(
  ~condition, ~serrf_file, ~annotation_file,
  "CTL", file.path(data_root, "SERRF Result_CTL", "normalized by - SERRF.csv"),
  file.path(annotation_root, "lipid_class_A.csv"),
  "LIN", file.path(data_root, "SERRF Result_LIN", "normalized by - SERRF.csv"),
  file.path(annotation_root, "lipid_class_B.csv")
)

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
    str_squish()
}

read_aliases <- function(file_name) {
  read_csv(file.path(name_map_root, file_name), show_col_types = FALSE,
           name_repair = "minimal") %>%
    select(LipidsMatch, OtherName) %>%
    filter(!is.na(OtherName), OtherName != "", str_to_lower(OtherName) != "empty")
}

aliases_raw <- bind_rows(
  read_aliases("name_change_new.csv"),
  read_aliases("remaining_lipid_name.csv")
) %>% distinct(LipidsMatch, OtherName)
ambiguous_aliases <- aliases_raw %>% count(LipidsMatch, name = "n_aliases") %>%
  filter(n_aliases > 1L) %>% pull(LipidsMatch)
aliases <- aliases_raw %>% filter(!LipidsMatch %in% ambiguous_aliases) %>%
  distinct(LipidsMatch, .keep_all = TRUE) %>% rename(curated_name = OtherName)

# Combine the maintained class table with prior manual class assignments. The
# newest manually curated assignment is used when a normalized name recurs.
class_lookup <- bind_rows(
  read_csv(file.path(project_root, "data", "lipid_class", "final_lipid_classes.csv"),
           show_col_types = FALSE) %>%
    transmute(name_key = normalize_lipid_name(Lipids), Class, SubClass, Sub_subclass),
  read_csv(file.path(project_root, "data", "raw_lipid_intensities", "corrected_lipid_names",
                     "1_CTL_lipids_name_corrected.csv"), show_col_types = FALSE) %>%
    transmute(name_key = normalize_lipid_name(Compound_Name), Class, SubClass, Sub_subclass),
  read_csv(file.path(project_root, "data", "raw_lipid_intensities", "corrected_lipid_names",
                     "2_LIN_lipids_name_corrected.csv"), show_col_types = FALSE) %>%
    transmute(name_key = normalize_lipid_name(Compound_Name), Class, SubClass, Sub_subclass)
) %>%
  filter(!is.na(Class), Class != "", Class != "Other") %>%
  distinct(name_key, .keep_all = TRUE)

read_normalized <- function(path) {
  x <- read_csv(path, show_col_types = FALSE, name_repair = "minimal")
  if (!"label" %in% names(x)) stop("SERRF output has no label column: ", path)
  x %>% mutate(feature_id = as.character(label))
}

process_one <- function(condition, serrf_path, annotation_path) {
  serrf <- read_normalized(serrf_path)
  sample_columns <- grep("^S1_", names(serrf), value = TRUE)
  if (!length(sample_columns)) stop("No biological sample columns in: ", serrf_path)

  annotation <- read.csv(annotation_path, check.names = FALSE, stringsAsFactors = FALSE) %>%
    mutate(feature_id = as.character(`#Scan#`)) %>%
    filter(feature_id %in% serrf$feature_id) %>%
    group_by(feature_id) %>%
    arrange(desc(MQScore), desc(SharedPeaks), MZErrorPPM, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(
      feature_id,
      source_name = Compound_Name,
      MQScore, SharedPeaks, MZErrorPPM,
      library_precursor_mz = Precursor_MZ,
      library_rt_seconds = RT_Query
    ) %>%
    left_join(aliases, by = c("source_name" = "LipidsMatch")) %>%
    mutate(
      Compound_Name = coalesce(curated_name, normalize_lipid_name(source_name)),
      name_key = normalize_lipid_name(Compound_Name)
    ) %>%
    left_join(class_lookup, by = "name_key")

  audit <- serrf %>%
    select(feature_id, label) %>%
    left_join(annotation, by = "feature_id") %>%
    mutate(
      condition = condition,
      retain = !is.na(Class) & Class != "Other",
      exclusion_reason = case_when(
        is.na(source_name) ~ "no_library_annotation",
        is.na(Class) ~ "unclassified_annotation",
        Class == "Other" ~ "other_compound",
        TRUE ~ NA_character_
      )
    )

  named_features <- serrf %>%
    select(feature_id, all_of(sample_columns)) %>%
    inner_join(
      audit %>% filter(retain) %>%
        select(feature_id, Compound_Name, Class, SubClass, Sub_subclass,
               source_name, MQScore, SharedPeaks, MZErrorPPM,
               library_precursor_mz, library_rt_seconds),
      by = "feature_id"
    ) %>%
    select(Compound_Name, Class, SubClass, Sub_subclass, feature_id,
           source_name, MQScore, SharedPeaks, MZErrorPPM,
           library_precursor_mz, library_rt_seconds, all_of(sample_columns))

  # Different source features with the same curated lipid label represent
  # distinct signal contributions. Exact duplicated numeric values are counted once.
  reaggregated <- named_features %>%
    group_by(Compound_Name) %>%
    summarise(
      Class = first(Class),
      SubClass = first(SubClass),
      Sub_subclass = first(Sub_subclass),
      source_feature_count = n(),
      source_feature_ids = paste(feature_id, collapse = ";"),
      across(all_of(sample_columns), ~ sum(unique(.x), na.rm = TRUE)),
      .groups = "drop"
    )

  write_csv(audit, file.path(data_root, paste0("postSERRF_", condition, "_feature_annotation_audit.csv")))
  write_csv(named_features, file.path(data_root, paste0("postSERRF_", condition, "_named_features.csv")))
  write_csv(reaggregated, file.path(data_root, paste0("postSERRF_", condition, "_named_reaggregated.csv")))

  tibble(
    condition,
    serrf_normalized_features = nrow(serrf),
    library_annotated_features = sum(!is.na(audit$source_name)),
    retained_non_other_features = sum(audit$retain),
    retained_reaggregated_lipids = nrow(reaggregated),
    excluded_no_library_annotation = sum(audit$exclusion_reason == "no_library_annotation", na.rm = TRUE),
    excluded_unclassified_annotation = sum(audit$exclusion_reason == "unclassified_annotation", na.rm = TRUE),
    excluded_other_compound = sum(audit$exclusion_reason == "other_compound", na.rm = TRUE)
  )
}

summary <- bind_rows(lapply(seq_len(nrow(settings)), function(i) {
  process_one(settings$condition[i], settings$serrf_file[i], settings$annotation_file[i])
}))
write_csv(summary, file.path(data_root, "postSERRF_processing_summary.csv"))
print(summary)
