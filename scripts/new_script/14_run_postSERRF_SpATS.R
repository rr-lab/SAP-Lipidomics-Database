#!/usr/bin/env Rscript

# Fit the established spatial additive model separately to every post-SERRF
# reaggregated lipid. Outputs 5/6 are genotype BLUPs for GWAS; 7/8 retain
# fitted plot-level values and spatial coordinates for downstream diagnostics.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(SpATS)
  library(tibble)
})

project_root <- "/Users/nirwantandukar/Documents/Github/SoLD_paper"
data_root <- "/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database/data/new_data"

settings <- tribble(
  ~condition, ~input_file, ~fieldmap_file, ~blup_file, ~fitted_file,
  "CTL", "postSERRF_CTL_named_reaggregated.csv", "data/fieldmap/fieldmap_control.csv",
  "5_CTL_lipid_BLUPs_for_GWAS.csv", "7_CTL_lipid_SpATS_fitted.csv",
  "LIN", "postSERRF_LIN_named_reaggregated.csv", "data/fieldmap/fieldmap_lowinput.csv",
  "6_LIN_lipid_BLUPs_for_GWAS.csv", "8_LIN_lipid_SpATS_fitted.csv"
)

normalize_line <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^PI_", "PI", x)
  x
}

sample_to_line <- function(x) {
  normalize_line(sub("^S1_(.*)_Run[0-9]+$", "\\1", x))
}

read_fieldmap <- function(path) {
  field <- as.matrix(read.csv(path, header = FALSE, check.names = FALSE,
                              stringsAsFactors = FALSE))
  out <- expand.grid(row = seq_len(nrow(field)), col = seq_len(ncol(field))) %>%
    mutate(LineRaw = normalize_line(as.vector(t(field)))) %>%
    filter(!is.na(LineRaw), LineRaw != "", LineRaw != "X")
  if (anyDuplicated(out$LineRaw)) {
    # Repeated checks are legitimate; a biological PI accession must map once.
    duplicated_pi <- out %>% filter(grepl("^PI", LineRaw)) %>% count(LineRaw) %>% filter(n > 1L)
    if (nrow(duplicated_pi)) stop("Duplicate biological accession(s) in field map: ",
                                 paste(duplicated_pi$LineRaw, collapse = ", "))
  }
  out
}

fit_condition <- function(condition, input_path, fieldmap_path, blup_path, fitted_path) {
  wide <- read_csv(input_path, show_col_types = FALSE, name_repair = "minimal")
  sample_columns <- grep("^S1_", names(wide), value = TRUE)
  lipid_columns <- setdiff(names(wide), c("Compound_Name", "Class", "SubClass", "Sub_subclass",
                                           "source_feature_count", "source_feature_ids", sample_columns))
  if (length(lipid_columns)) stop("Unexpected non-sample columns: ", paste(lipid_columns, collapse = ", "))
  map <- read_fieldmap(fieldmap_path)

  # Assign a stable plot ID before pivoting, so every lipid observation from a
  # given injection joins back to the same field position in the fitted table.
  sample_layout <- tibble(Sample = sample_columns) %>%
    mutate(LineRaw = sample_to_line(Sample)) %>%
    inner_join(map, by = "LineRaw") %>%
    mutate(row = as.integer(row), col = as.integer(col), PlotID = row_number())

  phenotype <- wide %>%
    select(Compound_Name, all_of(sample_columns)) %>%
    pivot_longer(-Compound_Name, names_to = "Sample", values_to = "Value") %>%
    inner_join(sample_layout, by = "Sample")

  mapped_samples <- sample_layout %>% distinct(Sample) %>% nrow()
  if (mapped_samples != length(sample_columns)) {
    stop(condition, ": only ", mapped_samples, " of ", length(sample_columns),
         " biological samples map to the field layout")
  }

  traits <- unique(phenotype$Compound_Name)
  blups <- phenotype %>% distinct(LineRaw) %>% arrange(LineRaw)
  fitted <- phenotype %>% distinct(Sample, LineRaw, row, col, PlotID) %>% arrange(PlotID)
  diagnostics <- vector("list", length(traits))

  for (i in seq_along(traits)) {
    trait <- traits[i]
    dat <- phenotype %>% filter(Compound_Name == trait, is.finite(Value)) %>%
      transmute(LineRaw, row, col, PlotID, Trait = Value)
    fit <- tryCatch(
      SpATS(
        response = "Trait",
        spatial = ~ SAP(col, row, nseg = c(8, 2), degree = 3, pord = 2),
        genotype = "LineRaw",
        data = dat,
        control = list(tolerance = 1e-3, maxit = 500),
        genotype.as.random = TRUE
      ),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      diagnostics[[i]] <- tibble(condition, Compound_Name = trait, status = "failed",
                                 n_plots = nrow(dat), message = conditionMessage(fit))
      next
    }
    geno_names <- fit$terms$geno$geno_names
    trait_blup <- tibble(LineRaw = geno_names, value = as.numeric(fit$coeff[geno_names]))
    names(trait_blup)[2] <- trait
    blups <- left_join(blups, trait_blup, by = "LineRaw")

    trait_fitted <- tibble(PlotID = dat$PlotID, value = as.numeric(fit$fitted))
    names(trait_fitted)[2] <- trait
    fitted <- left_join(fitted, trait_fitted, by = "PlotID")
    diagnostics[[i]] <- tibble(condition, Compound_Name = trait, status = "ok",
                               n_plots = nrow(dat), message = NA_character_)
    message(condition, ": fitted ", i, "/", length(traits), " ", trait)
  }

  diagnostics <- bind_rows(diagnostics)
  write_csv(blups, blup_path)
  write_csv(fitted, fitted_path)
  write_csv(diagnostics, file.path(data_root, paste0("SpATS_", condition, "_fit_diagnostics.csv")))
  tibble(condition, n_lipids = length(traits), n_success = sum(diagnostics$status == "ok"),
         n_failed = sum(diagnostics$status == "failed"), n_genotypes = nrow(blups), n_plots = nrow(fitted))
}

summary <- bind_rows(lapply(seq_len(nrow(settings)), function(i) {
  fit_condition(
    settings$condition[i],
    file.path(data_root, settings$input_file[i]),
    file.path(project_root, settings$fieldmap_file[i]),
    file.path(data_root, settings$blup_file[i]),
    file.path(data_root, settings$fitted_file[i])
  )
}))
write_csv(summary, file.path(data_root, "SpATS_processing_summary.csv"))
print(summary)
