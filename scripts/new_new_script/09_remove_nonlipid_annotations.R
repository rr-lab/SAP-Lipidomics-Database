#!/usr/bin/env Rscript

# Remove entries classified as Other after manual curation. The retained
# octadecylamine is an Amine/Fatty amines entry, not an Other compound.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

project_root <- "/Users/nirwantandukar/Documents/Github/SoLD_paper"
input_dir <- file.path(project_root, "data", "raw_lipid_intensities", "corrected_lipid_names")

filter_one <- function(file_name, condition) {
  path <- file.path(input_dir, file_name)
  x <- read_csv(path, show_col_types = FALSE)
  removed <- x %>%
    filter(Class == "Other") %>%
    transmute(condition, Compound_Name, Class, SubClass, Sub_subclass, final_row, X.Scan.)
  retained <- x %>% filter(Class != "Other")

  octadecylamine <- retained %>% filter(Compound_Name == "1-Octadecanamine")
  if (nrow(octadecylamine) != 1L || octadecylamine$Class != "Amine" ||
      octadecylamine$SubClass != "Fatty amines") {
    stop("1-Octadecanamine was not retained as Amine / Fatty amines in ", file_name)
  }

  write_csv(retained, path)
  removed
}

ctl_removed <- filter_one("1_CTL_lipids_name_corrected.csv", "CTL")
lin_removed <- filter_one("2_LIN_lipids_name_corrected.csv", "LIN")
removed <- bind_rows(ctl_removed, lin_removed)
write_csv(removed, file.path(input_dir, "excluded_other_compounds.csv"))
write_csv(
  removed %>% count(condition, name = "n_removed"),
  file.path(input_dir, "excluded_other_compound_summary.csv")
)
print(removed %>% count(condition, name = "n_removed"))
