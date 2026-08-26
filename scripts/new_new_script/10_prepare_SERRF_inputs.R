#!/usr/bin/env Rscript

# Copy the previously validated SERRF upload matrices without changing their
# batch, sampleType, time, QC, internal-standard, or zero-replacement rows.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

old_root <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SAP_lipids_GWAS"
output_dir <- "/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database/data/new_data"

inputs <- tribble(
  ~condition, ~source_file, ~output_file,
  "CTL", file.path(old_root, "results", "SERRF", "A_SERRF_new.csv"),
  file.path(output_dir, "3_CTL_SERRF_input.csv"),
  "LIN", file.path(old_root, "results", "SERRF", "B_SERRF_new.csv"),
  file.path(output_dir, "4_LIN_SERRF_input.csv")
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

inspect_serrf <- function(path, condition) {
  # SERRF's first four rows encode batch, sample type, injection order, and label.
  x <- read.csv(path, header = FALSE, check.names = FALSE, stringsAsFactors = FALSE)
  sample_types <- as.character(x[2, 3:ncol(x)])
  batches <- as.character(x[1, 3:ncol(x)])
  tibble(
    condition = condition,
    n_features = nrow(x) - 4L,
    n_injections = ncol(x) - 2L,
    n_qc = sum(sample_types == "qc"),
    n_samples = sum(sample_types == "sample"),
    n_validate = sum(sample_types == "validate"),
    n_batches = n_distinct(batches),
    batch_labels = paste(sort(unique(batches)), collapse = ";")
  )
}

for (i in seq_len(nrow(inputs))) {
  if (!file.exists(inputs$source_file[i])) {
    stop("Missing validated SERRF template: ", inputs$source_file[i])
  }
  ok <- file.copy(inputs$source_file[i], inputs$output_file[i], overwrite = TRUE)
  if (!ok) stop("Could not copy SERRF input for ", inputs$condition[i])
}

manifest <- bind_rows(lapply(seq_len(nrow(inputs)), function(i) {
  inspect_serrf(inputs$output_file[i], inputs$condition[i])
}))
write_csv(manifest, file.path(output_dir, "SERRF_input_manifest.csv"))
print(manifest)
