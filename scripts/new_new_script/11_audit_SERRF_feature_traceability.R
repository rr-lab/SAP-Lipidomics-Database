#!/usr/bin/env Rscript

# Verify that each curated lipid row can be traced to the exact feature sent to
# SERRF. This must pass before normalized values are assigned lipid identities.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

project_root <- "/Users/nirwantandukar/Documents/Github/SoLD_paper"
data_root <- "/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database/data/new_data"
curated_root <- file.path(project_root, "data", "raw_lipid_intensities", "corrected_lipid_names")

settings <- tribble(
  ~condition, ~curated_file, ~serrf_input,
  "CTL", "1_CTL_lipids_name_corrected.csv", "3_CTL_SERRF_input.csv",
  "LIN", "2_LIN_lipids_name_corrected.csv", "4_LIN_SERRF_input.csv"
)

read_serrf_input <- function(path) {
  x <- read.csv(path, header = FALSE, check.names = FALSE, stringsAsFactors = FALSE)
  labels <- as.character(x[-(1:4), 2])
  sample_names <- as.character(x[4, 3:ncol(x)])
  mat <- as.matrix(x[-(1:4), 3:ncol(x), drop = FALSE])
  storage.mode(mat) <- "numeric"
  colnames(mat) <- sample_names
  list(labels = labels, matrix = mat)
}

audit_one <- function(condition, curated_path, serrf_path) {
  curated <- read_csv(curated_path, show_col_types = FALSE)
  serrf <- read_serrf_input(serrf_path)
  sample_columns <- intersect(grep("^S1_", names(curated), value = TRUE), colnames(serrf$matrix))
  if (length(sample_columns) < 10L) stop("Insufficient common sample columns for ", condition)

  x <- as.matrix(curated[, sample_columns, drop = FALSE])
  storage.mode(x) <- "numeric"
  y <- serrf$matrix[, sample_columns, drop = FALSE]
  out <- vector("list", nrow(curated))

  for (i in seq_len(nrow(curated))) {
    v <- x[i, ]
    finite <- which(is.finite(v) & v > 0)
    anchors <- finite[seq_len(min(12L, length(finite)))]
    if (length(anchors) == 0L) {
      out[[i]] <- tibble(condition, final_row = curated$final_row[i], Compound_Name = curated$Compound_Name[i],
                         curated_scan_id = as.character(curated$scan_id[i]), n_common_samples = length(sample_columns),
                         match_status = "no_positive_values", matched_serrf_label = NA_character_, n_exact_matches = 0L)
      next
    }
    tol <- pmax(abs(v[anchors]) * 1e-10, 1e-6)
    candidate <- which(rowSums(abs(sweep(y[, anchors, drop = FALSE], 2, v[anchors], "-")) <= tol) == length(anchors))
    exact <- integer()
    if (length(candidate) > 0L) {
      tol_all <- pmax(abs(v[finite]) * 1e-10, 1e-6)
      exact <- candidate[rowSums(abs(sweep(y[candidate, finite, drop = FALSE], 2, v[finite], "-")) <= tol_all) == length(finite)]
    }
    status <- if (length(exact) == 1L) "exact_profile_match" else if (length(exact) > 1L) "ambiguous_exact_profile_match" else "no_exact_profile_match"
    out[[i]] <- tibble(condition, final_row = curated$final_row[i], Compound_Name = curated$Compound_Name[i],
                       curated_scan_id = as.character(curated$scan_id[i]), n_common_samples = length(sample_columns),
                       match_status = status,
                       matched_serrf_label = if (length(exact)) paste(serrf$labels[exact], collapse = ";") else NA_character_,
                       n_exact_matches = length(exact))
  }
  bind_rows(out)
}

audit <- bind_rows(lapply(seq_len(nrow(settings)), function(i) {
  audit_one(
    settings$condition[i],
    file.path(curated_root, settings$curated_file[i]),
    file.path(data_root, settings$serrf_input[i])
  )
}))

write_csv(audit, file.path(data_root, "SERRF_feature_traceability_audit.csv"))
summary <- audit %>% count(condition, match_status, name = "n_curated_rows")
write_csv(summary, file.path(data_root, "SERRF_feature_traceability_summary.csv"))
print(summary)
