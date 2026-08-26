#!/usr/bin/env Rscript

# Build the current SERRF inputs from the raw feature matrices. Filtering occurs
# before zero replacement so pseudo-counts are never counted as detections.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

project_root <- "/Users/nirwantandukar/Documents/Github/SoLD_paper"
output_root <- "/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database/data/new_data"
old_root <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SAP_lipids_GWAS"

# A feature must be detected in more than two thirds of biological injections.
max_zero_fraction <- 1 / 3
blank_multiplier <- 10
min_sample_fraction_above_blank <- 0.70

settings <- tribble(
  ~condition, ~raw_file, ~template_file, ~preprocessed_file, ~serrf_file,
  "CTL", "SetA_lipid_FLO2019Control.csv", "A_SERRF_new.csv",
  "1_CTL_preSERRF_filtered_features.csv", "3_CTL_SERRF_input.csv",
  "LIN", "SetB_lipid_FLO2022_lowP.csv", "B_SERRF_new.csv",
  "2_LIN_preSERRF_filtered_features.csv", "4_LIN_SERRF_input.csv"
)

canonical_injection_name <- function(x) {
  type <- case_when(
    str_detect(x, "S1_Run") ~ "sample",
    str_detect(x, "QC_Run") ~ "qc",
    str_detect(x, "ISTD") ~ "validate",
    str_detect(x, "InjBL") ~ "blank",
    TRUE ~ NA_character_
  )
  run <- case_when(
    type == "sample" ~ str_match(x, "S1_Run(\\d+)")[, 2],
    type == "qc" ~ str_match(x, "QC_Run(\\d+)")[, 2],
    type == "validate" ~ str_match(x, "ISTD_Run(\\d+)")[, 2],
    type == "blank" ~ str_match(x, "InjBL[^_]*_Run(\\d+)")[, 2],
    TRUE ~ NA_character_
  )
  sample_id <- str_extract(x, "PI\\d+|CHECK\\d+_\\d+|Check_\\d+")
  case_when(
    type == "sample" ~ paste0("S1_", sample_id, "_Run", run),
    type == "qc" ~ paste0("QC_Run", run),
    type == "validate" ~ paste0("ISTD_Run", run),
    type == "blank" ~ paste0("InjBL_Run", run),
    TRUE ~ NA_character_
  )
}

read_template_metadata <- function(path) {
  x <- read.csv(path, header = FALSE, check.names = FALSE, stringsAsFactors = FALSE)
  metadata <- x[1:4, , drop = FALSE]
  colnames(metadata) <- as.character(metadata[4, ])
  metadata
}

build_one <- function(condition, raw_path, template_path, preprocessed_path, serrf_path) {
  raw <- read.csv(raw_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("row ID", "row retention time") %in% names(raw))) {
    stop("Required raw feature columns missing in ", raw_path)
  }
  # Raw vendor exports contain blank metadata headers, which dplyr rejects.
  raw <- raw[is.finite(raw[["row retention time"]]) & raw[["row retention time"]] >= 1, , drop = FALSE]

  peak_columns <- grep("Peak area$", names(raw), value = TRUE)
  injection_names <- canonical_injection_name(peak_columns)
  keep <- !is.na(injection_names)
  peak_columns <- peak_columns[keep]
  injection_names <- injection_names[keep]
  if (anyDuplicated(injection_names)) {
    stop("Duplicate canonical injection names in ", condition, ": ",
         paste(unique(injection_names[duplicated(injection_names)]), collapse = ", "))
  }

  values <- raw[, peak_columns, drop = FALSE]
  names(values) <- injection_names
  sample_columns <- grep("^S1_", names(values), value = TRUE)
  blank_columns <- grep("^InjBL_", names(values), value = TRUE)
  if (!length(sample_columns) || !length(blank_columns)) {
    stop("Could not identify biological samples and blanks for ", condition)
  }

  sample_matrix <- as.matrix(values[, sample_columns, drop = FALSE])
  storage.mode(sample_matrix) <- "numeric"
  blank_matrix <- as.matrix(values[, blank_columns, drop = FALSE])
  storage.mode(blank_matrix) <- "numeric"
  blank_mean <- rowMeans(blank_matrix, na.rm = TRUE)
  above_blank <- rowMeans(sample_matrix > blank_multiplier * blank_mean, na.rm = TRUE)
  zero_fraction <- rowMeans(!is.finite(sample_matrix) | sample_matrix <= 0, na.rm = TRUE)
  keep_feature <- is.finite(above_blank) & above_blank >= min_sample_fraction_above_blank &
    is.finite(zero_fraction) & zero_fraction < max_zero_fraction

  filtered <- bind_cols(
    raw[keep_feature, c("row ID", "row m/z", "row retention time"), drop = FALSE],
    as_tibble(values[keep_feature, setdiff(names(values), blank_columns), drop = FALSE])
  )
  names(filtered)[1:3] <- c("feature_id", "precursor_mz", "retention_time")
  filtered$feature_id <- as.character(filtered$feature_id)

  injection_columns <- setdiff(names(filtered), c("feature_id", "precursor_mz", "retention_time"))
  for (column in injection_columns) {
    z <- filtered[[column]]
    z[!is.finite(z) | z < 0] <- 0
    positive <- z[z > 0]
    if (!length(positive)) stop("All-zero retained feature/injection column: ", condition, " ", column)
    z[z == 0] <- (2 / 3) * min(positive)
    filtered[[column]] <- z
  }
  write_csv(filtered, preprocessed_path)

  metadata <- read_template_metadata(template_path)
  requested_columns <- colnames(metadata)[-(1:2)]
  if (!setequal(requested_columns, injection_columns)) {
    missing <- setdiff(requested_columns, injection_columns)
    extra <- setdiff(injection_columns, requested_columns)
    stop("Metadata and current raw injections differ for ", condition,
         "\nMissing: ", paste(missing, collapse = ", "),
         "\nExtra: ", paste(extra, collapse = ", "))
  }
  filtered <- filtered[, c("feature_id", requested_columns), drop = FALSE]
  output <- rbind(
    metadata,
    cbind(No = seq_len(nrow(filtered)), label = filtered$feature_id, filtered[, requested_columns, drop = FALSE])
  )
  # SERRF expects batch/sampleType/time/label as the first four rows, not a CSV header.
  write.table(output, serrf_path, sep = ",", row.names = FALSE, col.names = FALSE,
              quote = FALSE, na = "")

  tibble(
    condition = condition,
    raw_features_rt_ge_1 = nrow(raw),
    retained_features = nrow(filtered),
    removed_blank_filter = sum(above_blank < min_sample_fraction_above_blank, na.rm = TRUE),
    removed_zero_filter = sum(zero_fraction >= max_zero_fraction, na.rm = TRUE),
    n_samples = length(sample_columns),
    n_qc = sum(metadata[2, -(1:2)] == "qc"),
    n_validate = sum(metadata[2, -(1:2)] == "validate")
  )
}

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
summary <- bind_rows(lapply(seq_len(nrow(settings)), function(i) {
  build_one(
    settings$condition[i],
    file.path(project_root, "data", "raw_lipid_intensities", settings$raw_file[i]),
    file.path(old_root, "results", "SERRF", settings$template_file[i]),
    file.path(output_root, settings$preprocessed_file[i]),
    file.path(output_root, settings$serrf_file[i])
  )
}))
write_csv(summary, file.path(output_root, "preSERRF_filtering_summary.csv"))
print(summary)
