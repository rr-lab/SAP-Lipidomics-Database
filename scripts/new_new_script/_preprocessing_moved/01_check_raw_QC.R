################################################################################
# Pre-SPATS QC Check
#
# Purpose:
# - Run QC diagnostics on raw/cleaned lipid intensity matrices BEFORE SpATS.
# - Works on Set A and Set B files under data/raw_lipid_intensities.
# - Treats CHECK injections as QC when explicit QC columns are absent.
#
# Outputs:
# - fig/supp/PreSPATS_QC_<Set>_TIC.png
# - fig/supp/PreSPATS_QC_<Set>_PCA.png
# - table/supp/PreSPATS_QC_<Set>_InjectionSummary.csv
# - table/supp/PreSPATS_QC_<Set>_QCFeatureStats.csv
# - table/supp/PreSPATS_QC_Summary.csv
################################################################################

suppressPackageStartupMessages({
  library(vroom)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(tibble)
})

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(file_arg) == 0) return(NULL)
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE))
}

find_project_root <- function() {
  starts <- unique(c(getwd(), get_script_dir()))
  starts <- starts[!is.na(starts) & nzchar(starts)]

  for (s in starts) {
    cur <- normalizePath(s, winslash = "/", mustWork = FALSE)
    repeat {
      if (dir.exists(file.path(cur, "data", "raw_lipid_intensities"))) {
        return(cur)
      }
      parent <- dirname(cur)
      if (identical(parent, cur)) break
      cur <- parent
    }
  }
  NULL
}

project_root <- find_project_root()
if (is.null(project_root)) {
  stop("Could not locate project root containing data/raw_lipid_intensities. ",
       "Current working directory: ", getwd())
}

fig_dir <- file.path(project_root, "fig", "supp")
table_dir <- file.path(project_root, "table", "supp")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

message("\n========================================")
message("Pre-SPATS QC Check")
message("========================================\n")
message("Project root: ", project_root)

pick_input_file <- function(set_id) {
  candidates <- c(
    file.path(project_root, "data", "raw_lipid_intensities", paste0(set_id, "_final_lipids.csv")),
    file.path(project_root, "data", "raw_lipid_intensities", paste0(set_id, "_cleaned_lipids_unsummed.csv")),
    file.path(project_root, "data", "raw_lipid_intensities", paste0(set_id, "_cleaned_lipids.csv"))
  )
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) return(NULL)
  hit
}

classify_injection <- function(x) {
  case_when(
    str_detect(x, regex("^QC_", ignore_case = TRUE)) ~ "QC",
    str_detect(x, regex("CHECK", ignore_case = TRUE)) ~ "QC",
    str_detect(x, regex("^InjBL", ignore_case = TRUE)) ~ "Blank",
    str_detect(x, regex("^ISTD", ignore_case = TRUE)) ~ "ISTD",
    str_detect(x, "^S\\d+_PI\\d+_Run\\d+") ~ "Sample",
    str_detect(x, "^S\\d+_Run\\d+") ~ "Sample",
    TRUE ~ "Other"
  )
}

extract_run <- function(x) {
  as.integer(str_extract(x, "(?<=Run)\\d+"))
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(0)
  stats::median(x)
}

load_set_matrix <- function(file_path) {
  dat <- suppressMessages(vroom(file_path, show_col_types = FALSE, delim = ","))

  inj_cols <- names(dat)[str_detect(names(dat), "Run\\d+")]
  if (length(inj_cols) == 0) {
    stop("No injection columns with pattern 'Run<digits>' found in: ", file_path)
  }

  mat <- dat %>%
    dplyr::select(all_of(inj_cols)) %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(.x)))) %>%
    as.matrix()

  mode(mat) <- "numeric"
  keep_rows <- rowSums(is.finite(mat) & mat != 0, na.rm = TRUE) > 0
  mat <- mat[keep_rows, , drop = FALSE]

  meta <- tibble(
    Injection = inj_cols,
    Type = classify_injection(inj_cols),
    Run = extract_run(inj_cols)
  ) %>%
    mutate(
      Run = ifelse(is.na(Run), row_number(), Run)
    ) %>%
    arrange(Run)

  mat <- mat[, meta$Injection, drop = FALSE]

  list(mat = mat, meta = meta)
}

make_tic_plot <- function(meta, mat, set_label) {
  tic_df <- meta %>%
    mutate(TIC = colSums(mat, na.rm = TRUE))

  p <- ggplot(tic_df, aes(x = Run, y = TIC, color = Type, group = Type)) +
    geom_line(alpha = 0.5, linewidth = 0.6) +
    geom_point(size = 2.3) +
    scale_y_continuous(labels = scales::scientific) +
    labs(
      title = paste0("Pre-SPATS TIC Across Injection Order (", set_label, ")"),
      x = "Injection Run Order",
      y = "Total Ion Current (TIC)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  list(plot = p, tic_df = tic_df)
}

make_pca_plot <- function(meta, mat, set_label) {
  x <- t(log10(mat + 1))
  x <- as.data.frame(x)

  for (j in seq_len(ncol(x))) {
    x[[j]][!is.finite(x[[j]])] <- NA_real_
    med <- safe_median(x[[j]])
    x[[j]][is.na(x[[j]])] <- med
  }

  sds <- apply(x, 2, stats::sd)
  keep <- is.finite(sds) & sds > 0
  x <- x[, keep, drop = FALSE]
  if (ncol(x) < 2 || nrow(x) < 3) {
    return(NULL)
  }

  pca <- stats::prcomp(x, center = TRUE, scale. = TRUE)
  ve <- (pca$sdev^2) / sum(pca$sdev^2)
  df <- bind_cols(
    meta,
    as_tibble(pca$x[, 1:2, drop = FALSE]) %>%
      setNames(c("PC1", "PC2"))
  )

  outliers <- df %>%
    group_by(Type) %>%
    slice_max(order_by = PC1^2 + PC2^2, n = 1, with_ties = FALSE) %>%
    ungroup()

  p <- ggplot(df, aes(PC1, PC2, color = Type)) +
    geom_point(size = 2.3, alpha = 0.85) +
    stat_ellipse(level = 0.8, linewidth = 0.7, linetype = "dashed", na.rm = TRUE) +
    ggrepel::geom_text_repel(
      data = outliers,
      aes(label = Injection),
      size = 2.8,
      max.overlaps = 8,
      show.legend = FALSE
    ) +
    labs(
      title = paste0("Pre-SPATS PCA by Injection Type (", set_label, ")"),
      x = paste0("PC1 (", round(ve[1] * 100, 1), "%)"),
      y = paste0("PC2 (", round(ve[2] * 100, 1), "%)")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  p
}

qc_feature_stats <- function(meta, mat) {
  qc_cols <- which(meta$Type == "QC")
  if (length(qc_cols) < 2) {
    return(tibble(
      Metric = c("N_QC_Injections", "Median_QC_CV_percent", "Features_abs_spearman_ge_0.6"),
      Value = c(length(qc_cols), NA_real_, NA_real_)
    ))
  }

  qc_mat <- mat[, qc_cols, drop = FALSE]
  qc_run <- meta$Run[qc_cols]

  mu <- rowMeans(qc_mat, na.rm = TRUE)
  sdv <- apply(qc_mat, 1, stats::sd, na.rm = TRUE)
  cv <- ifelse(mu > 0, 100 * sdv / mu, NA_real_)

  rho <- apply(qc_mat, 1, function(v) {
    suppressWarnings(cor(v, qc_run, method = "spearman", use = "pairwise.complete.obs"))
  })

  tibble(
    Metric = c("N_QC_Injections", "Median_QC_CV_percent", "Features_abs_spearman_ge_0.6"),
    Value = c(
      length(qc_cols),
      stats::median(cv, na.rm = TRUE),
      sum(abs(rho) >= 0.6, na.rm = TRUE)
    )
  )
}

run_set_qc <- function(set_id) {
  file_path <- pick_input_file(set_id)
  if (is.null(file_path)) {
    message("Skipping Set ", set_id, ": no input file found.")
    return(NULL)
  }

  message("Processing Set ", set_id, " using: ", file_path)
  payload <- load_set_matrix(file_path)
  mat <- payload$mat
  meta <- payload$meta

  message("  Features retained: ", nrow(mat))
  message("  Injections retained: ", ncol(mat))
  message("  Injection types: ", paste(names(table(meta$Type)), table(meta$Type), collapse = ", "))

  tic_out <- make_tic_plot(meta, mat, set_id)
  ggsave(
    filename = file.path(fig_dir, paste0("PreSPATS_QC_Set", set_id, "_TIC.png")),
    plot = tic_out$plot, width = 12, height = 6, dpi = 300, bg = "white"
  )

  pca_plot <- make_pca_plot(meta, mat, set_id)
  if (!is.null(pca_plot)) {
    ggsave(
      filename = file.path(fig_dir, paste0("PreSPATS_QC_Set", set_id, "_PCA.png")),
      plot = pca_plot, width = 8, height = 7, dpi = 300, bg = "white"
    )
  } else {
    message("  PCA skipped for Set ", set_id, " (insufficient data after filtering).")
  }

  inj_summary <- tic_out$tic_df %>%
    group_by(Type) %>%
    summarise(
      N = n(),
      Mean_TIC = mean(TIC, na.rm = TRUE),
      SD_TIC = sd(TIC, na.rm = TRUE),
      CV_TIC_percent = ifelse(Mean_TIC > 0, 100 * SD_TIC / Mean_TIC, NA_real_),
      .groups = "drop"
    ) %>%
    arrange(desc(N))

  feature_summary <- qc_feature_stats(meta, mat)

  write.csv(
    inj_summary,
    file.path(table_dir, paste0("PreSPATS_QC_Set", set_id, "_InjectionSummary.csv")),
    row.names = FALSE
  )
  write.csv(
    feature_summary,
    file.path(table_dir, paste0("PreSPATS_QC_Set", set_id, "_QCFeatureStats.csv")),
    row.names = FALSE
  )

  list(
    Set = set_id,
    File = file_path,
    Features = nrow(mat),
    Injections = ncol(mat),
    N_QC = sum(meta$Type == "QC"),
    N_Sample = sum(meta$Type == "Sample"),
    N_Blank = sum(meta$Type == "Blank"),
    N_ISTD = sum(meta$Type == "ISTD"),
    Median_QC_CV = suppressWarnings(feature_summary$Value[feature_summary$Metric == "Median_QC_CV_percent"])
  )
}

results <- list(
  run_set_qc("A"),
  run_set_qc("B")
)

results <- results[!vapply(results, is.null, logical(1))]

if (length(results) > 0) {
  summary_tbl <- bind_rows(results)
  out_sum <- file.path(table_dir, "PreSPATS_QC_Summary.csv")
  write.csv(summary_tbl, out_sum, row.names = FALSE)
  message("\nSaved: ", out_sum)
} else {
  message("\nNo Set A/B files found. No outputs generated.")
}

message("\nQC check complete.")
