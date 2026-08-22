suppressPackageStartupMessages({library(readr); library(dplyr); library(stringr); library(purrr)})
input_dir <- '/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database/data/new_data'
out_dir <- file.path(input_dir, 'duplicate_profile_audit')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

canonicalize <- function(x) {
  # Collapse acyl-chain order only; retain class and all non-acyl descriptors.
  str_replace(x, '^([A-Za-z]+)\\(([^()]*)\\)$', function(m) {
    cls <- str_match(m, '^([A-Za-z]+)\\(([^()]*)\\)$')[, 2]
    body <- str_match(m, '^([A-Za-z]+)\\(([^()]*)\\)$')[, 3]
    if (is.na(body) || !str_detect(body, '^[0-9]+:[0-9]+([/_][0-9]+:[0-9]+)+$')) return(m)
    chains <- str_split(body, '[/_]', simplify = TRUE)
    paste0(cls, '(', paste(sort(chains[chains != '']), collapse = '_'), ')')
  })
}

audit_one <- function(path, condition) {
  x <- read_csv(path, show_col_types = FALSE)
  samples <- grep('^S1_', names(x), value = TRUE)
  stopifnot(length(samples) > 0)
  profile_key <- apply(as.matrix(x[, samples]), 1, function(z) paste(format(z, digits = 16, scientific = FALSE, trim = TRUE), collapse = '|'))
  x$profile_key <- profile_key
  x$canonical_unordered_species <- vapply(x$Compound_Name, canonicalize, character(1))

  exact <- x %>%
    group_by(profile_key) %>%
    filter(n() > 1) %>%
    arrange(Compound_Name, .by_group = TRUE) %>%
    mutate(
      condition = condition,
      duplicate_type = 'identical_profile',
      duplicate_group_size = n()
    ) %>%
    ungroup() %>%
    select(condition, duplicate_type, duplicate_group_size, Compound_Name, canonical_unordered_species,
           source_row_count, source_rows, source_scans, all_of(samples))

  ordered <- x %>%
    group_by(canonical_unordered_species) %>%
    filter(n() > 1, n_distinct(Compound_Name) > 1) %>%
    arrange(canonical_unordered_species, Compound_Name, .by_group = TRUE) %>%
    mutate(
      condition = condition,
      duplicate_type = 'reversed_acyl_order_nonidentical_profile',
      duplicate_group_size = n()
    ) %>%
    ungroup() %>%
    select(condition, duplicate_type, duplicate_group_size, Compound_Name, canonical_unordered_species,
           source_row_count, source_rows, source_scans, all_of(samples))

  list(exact = exact, ordered = ordered,
       summary = tibble(condition = condition, n_features = nrow(x),
                        n_exact_duplicate_rows = nrow(exact),
                        n_exact_duplicate_groups = n_distinct(profile_key[duplicated(profile_key) | duplicated(profile_key, fromLast = TRUE)]),
                        n_reversed_order_rows = nrow(ordered),
                        n_reversed_order_groups = n_distinct(ordered$canonical_unordered_species)))
}

ctl <- audit_one(file.path(input_dir, '1_CTL_lipids_reaggregated.csv'), 'CTL')
lin <- audit_one(file.path(input_dir, '2_LIN_lipids_reaggregated.csv'), 'LIN')
write_csv(bind_rows(ctl$summary, lin$summary), file.path(out_dir, 'duplicate_profile_summary.csv'))
write_csv(bind_rows(ctl$exact, lin$exact), file.path(out_dir, 'identical_intensity_profiles.csv'))
write_csv(bind_rows(ctl$ordered, lin$ordered), file.path(out_dir, 'reversed_acyl_order_profiles.csv'))
print(bind_rows(ctl$summary, lin$summary))
