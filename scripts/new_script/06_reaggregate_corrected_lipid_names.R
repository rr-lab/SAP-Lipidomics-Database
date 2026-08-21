#!/usr/bin/env Rscript

# Reaggregate standardized lipid names without averaging intensities.
# Exact duplicate sample profiles are counted once; distinct profiles are summed.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

project_root <- "/Users/nirwantandukar/Documents/Github/SoLD_paper"
input_root <- file.path(project_root, "data", "raw_lipid_intensities", "corrected_lipid_names")

reaggregate_one <- function(input_file, output_file) {
  x <- read_csv(input_file, show_col_types = FALSE)
  sample_columns <- grep("^S1_", names(x), value = TRUE)

  if (length(sample_columns) == 0L) {
    stop("No sample-intensity columns beginning with 'S1_' found in: ", input_file)
  }

  out <- x %>%
    mutate(
      source_scan = as.character(X.Scan.),
      source_row = final_row
    ) %>%
    group_by(Compound_Name) %>%
    summarise(
      Class = first(na.omit(Class)),
      SubClass = first(na.omit(SubClass)),
      Sub_subclass = first(na.omit(Sub_subclass)),
      source_row_count = n(),
      source_rows = paste(source_row, collapse = ";"),
      source_scans = paste(na.omit(source_scan), collapse = ";"),
      across(all_of(sample_columns), ~ sum(unique(.x), na.rm = TRUE)),
      .groups = "drop"
    )

  write_csv(out, output_file)
  list(data = out, input_rows = nrow(x))
}

ctl <- reaggregate_one(
  file.path(input_root, "1_CTL_lipids_name_corrected.csv"),
  file.path(input_root, "1_CTL_lipids_name_corrected_reaggregated.csv")
)
lin <- reaggregate_one(
  file.path(input_root, "2_LIN_lipids_name_corrected.csv"),
  file.path(input_root, "2_LIN_lipids_name_corrected_reaggregated.csv")
)

summary <- tibble(
  condition = c("CTL", "LIN"),
  input_rows = c(ctl$input_rows, lin$input_rows),
  reaggregated_lipids = c(nrow(ctl$data), nrow(lin$data)),
  duplicate_rows_collapsed = c(
    ctl$input_rows - nrow(ctl$data),
    lin$input_rows - nrow(lin$data)
  )
)
write_csv(summary, file.path(input_root, "reaggregation_summary.csv"))
print(summary)
