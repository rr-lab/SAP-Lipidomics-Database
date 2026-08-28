#!/usr/bin/env Rscript

# Classify previously unassigned library annotations. Non-lipid chemicals and
# unresolved library IDs remain explicit as "Other" rather than forced into a lipid class.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

project_root <- "/Users/nirwantandukar/Documents/Github/SoLD_paper"
input_dir <- file.path(project_root, "data", "raw_lipid_intensities", "corrected_lipid_names")

curated_classes <- tribble(
  ~Compound_Name, ~Class, ~SubClass, ~Sub_subclass, ~classification_basis,
  "DG(12:0/0:0/12:0)", "Glycerolipid", "Diacylglycerol", "DG", "DG shorthand",
  "DG(8:0/8:0/0:0)", "Glycerolipid", "Diacylglycerol", "DG", "DG shorthand",
  "1-Octadecanamine", "Amine", "Fatty amines", "Long-chain primary amine", "C18 primary fatty amine",
  "Glyceryl monolinoleate", "Glycerolipid", "Monoacylglycerol", "MG", "monoacylglycerol",
  "Monoelaidin", "Glycerolipid", "Monoacylglycerol", "MG", "monoacylglycerol",
  "1-Hexadecanoyl-sn-glycero-3-phosphocholine", "Glycerophospholipid", "Lysophosphatidylcholine", "Lysophosphatidylcholine", "lysophosphatidylcholine nomenclature",
  "L-alpha-phosphatidyl choline", "Glycerophospholipid", "Phosphatidylcholine", "Phosphatidylcholine", "phosphatidylcholine nomenclature",
  "9,12-Octadecadiynoic acid", "Fatty acid", "Fatty acid", NA_character_, "fatty acid nomenclature",
  "cis-2-Decenoic acid", "Fatty acid", "Fatty acid", NA_character_, "fatty acid nomenclature",
  "9-Oxo-11-(3-pentyl-2-oxiranyl)-10E-undecenoic acid", "Fatty acid", "Fatty acid", "Oxidized fatty acid", "oxidized fatty acid",
  "13S-Hydroxy-9Z,11E,15Z-octadecatrienoic acid", "Fatty acid", "Fatty acid", "Oxylipin", "oxygenated fatty acid",
  "4.alpha./.beta.-Hydroxyprostanozol", "Steroid", "Androgen", NA_character_, "steroidal scaffold",
  "5(6)-Epoxy-8Z,11Z,14Z-eicosatrienoic acid, methyl ester", "Fatty acid derivative", "Epoxyeicosatrienoic acid", NA_character_, "eicosanoid ester",
  "Stearidonic acid methyl ester", "Fatty acid derivative", "Fatty acid ester", NA_character_, "fatty acid ester",
  "NCGC00386020-01_C18H28O3_8-{(1S,5R)-4-Oxo-5-[(2Z)-2-penten-1-yl]-2-cyclopenten-1-yl}octanoic acid", "Fatty acid derivative", "Oxylipin", NA_character_, "oxygenated fatty-acid derivative",
  "Erucamide", "Fatty acid amide", NA_character_, NA_character_, "fatty acid amide",
  "N-Cyclohexanecarbonylpentadecylamine", "Fatty acid amide derivative", NA_character_, NA_character_, "long-chain acyl amide",
  "(4R)-4-((5S,7R,9S,10S,13R,14S,17R)-7-hydroxy-10,13-dimethyl-3-oxohexadecahydro-1H-cyclopenta[a]phenanthren-17-yl)pentanoic acid", "Steroid", "Steroid derivative", NA_character_, "steroid scaffold",
  "5.alpha.-Pregnan-3.alpha.,17-diol-20-one 3-sulfate", "Steroid", "Progestogen", NA_character_, "pregnane steroid",
  "Nandrolone", "Steroid", "Androgen", NA_character_, "androgen steroid",
  "h_42_ethylestrenol_m", "Steroid", "Androgen", NA_character_, "androgen steroid",
  "Cholesta-4,6-dien-3-one", "Sterol", "Sterol", NA_character_, "sterol scaffold",
  "Cholesterol", "Sterol", "Sterol", NA_character_, "sterol nomenclature",
  "alpha-Carotene", "Terpenoid", "Tetraterpenoid", "Carotenoid", "carotenoid nomenclature",
  "Austinoneol", "Terpenoid", "Meroterpenoid", NA_character_, "published meroterpenoid annotation",
  "Dipterocarpol", "Terpenoid", "Triterpenoid", NA_character_, "triterpenoid nomenclature",
  "Glochidone", "Terpenoid", "Triterpenoid", NA_character_, "triterpenoid nomenclature",
  "Panaxatriol", "Terpenoid", "Triterpenoid", NA_character_, "triterpenoid nomenclature",
  "NCGC00384683-01_C30H50O_Lup-20(29)-en-3-ol,(3alpha)-", "Terpenoid", "Triterpenoid", NA_character_, "lupane triterpenoid",
  "NCGC00179780-02_C20H32O_(5xi,9xi)-Beyer-15-en-18-ol", "Terpenoid", "Diterpenoid", NA_character_, "beyerane diterpenoid"
)

classify_file <- function(file_name, condition) {
  path <- file.path(input_dir, file_name)
  x <- read_csv(path, show_col_types = FALSE)
  original_class <- if ("original_Class" %in% names(x)) x$original_Class else x$Class
  original_subclass <- if ("original_SubClass" %in% names(x)) x$original_SubClass else x$SubClass
  original_sub_subclass <- if ("original_Sub_subclass" %in% names(x)) x$original_Sub_subclass else x$Sub_subclass
  needs_class <- is.na(original_class) | original_class == ""

  # Start by explicitly retaining non-lipid or unresolved library assignments.
  x <- x %>%
    mutate(
      original_Class = original_class,
      original_SubClass = original_subclass,
      original_Sub_subclass = original_sub_subclass,
      classification_basis = if_else(needs_class, "Other: non-lipid, contaminant, or unresolved library annotation", NA_character_),
      Class = if_else(needs_class, "Other", Class),
      SubClass = if_else(needs_class, "Other", SubClass),
      Sub_subclass = if_else(needs_class, "Other", Sub_subclass)
    ) %>%
    left_join(curated_classes, by = "Compound_Name", suffix = c("", "_curated")) %>%
    mutate(
      Class = coalesce(Class_curated, Class),
      SubClass = coalesce(SubClass_curated, SubClass),
      Sub_subclass = coalesce(Sub_subclass_curated, Sub_subclass),
      classification_basis = coalesce(classification_basis_curated, classification_basis)
    ) %>%
    select(-ends_with("_curated"))

  assignments <- x %>%
    filter(needs_class) %>%
    transmute(
      condition = condition,
      Compound_Name,
      original_Class,
      assigned_Class = Class,
      assigned_SubClass = SubClass,
      assigned_Sub_subclass = Sub_subclass,
      classification_basis,
      final_row,
      X.Scan.
    )

  write_csv(
    x %>% select(-original_Class, -original_SubClass, -original_Sub_subclass, -classification_basis),
    path
  )
  assignments
}

ctl_assignments <- classify_file("1_CTL_lipids_name_corrected.csv", "CTL")
lin_assignments <- classify_file("2_LIN_lipids_name_corrected.csv", "LIN")
assignments <- bind_rows(ctl_assignments, lin_assignments)
write_csv(assignments, file.path(input_dir, "nonfocused_compound_class_assignments.csv"))
write_csv(
  assignments %>% count(condition, assigned_Class, assigned_SubClass, sort = TRUE),
  file.path(input_dir, "nonfocused_compound_class_assignment_summary.csv")
)
print(assignments %>% count(condition, assigned_Class, sort = TRUE))
