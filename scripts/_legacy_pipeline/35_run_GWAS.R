# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5: GENOME-WIDE ASSOCIATION STUDIES OF LIPID TRAITS UNDER LOWINPUT
# ══════════════════════════════════════════════════════════════════════════════
#
# This script generates ONLY the figures and tables for Section 5:
# "Genome-wide Association Studies of Lipid Traits Under LowInput"
#
# OUTPUTS:
#   - Figure 3 (Main): GWAS Manhattan Plots by Functional Category (4x4 grid)
#   - Table 1 (Main): Representative GWAS candidate genes (already in LaTeX)
#   - Note: Supplementary Tables S8-S10 are generated separately via annotation scripts
#
# ══════════════════════════════════════════════════════════════════════════════

message("\n")
message("══════════════════════════════════════════════════════════════")
message("SECTION 5: GENOME-WIDE ASSOCIATION STUDIES")
message("══════════════════════════════════════════════════════════════")
message("\n")

# ──────────────────────────────────────────────────────────────────────────────
# SETUP
# ──────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(vroom)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ggstatsplot)
})

# Create output directories
dir.create("fig/main", recursive = TRUE, showWarnings = FALSE)
dir.create("fig/supp", recursive = TRUE, showWarnings = FALSE)
dir.create("table/supp", recursive = TRUE, showWarnings = FALSE)


# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
#
#            FIGURE 3 (MAIN): GWAS MANHATTAN PLOTS BY FUNCTIONAL CATEGORY
#
# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════

message("\n══════════════════════════════════════════════════════════════")
message("FIGURE 3 (MAIN): GWAS MANHATTAN PLOTS - FUNCTIONAL CATEGORIES (4x4 Grid)")
message("══════════════════════════════════════════════════════════════\n")

# ──────────────────────────────────────────────────────────────────────────────
# CONFIGURATION
# ──────────────────────────────────────────────────────────────────────────────

base_dir <- "results/gwas_raw"
out_dir <- "fig/main"

# Marker SNP positions for INDIVIDUAL phenotypes (rows 1-2)
markers_ind <- list(
  Lipid = list(chr = 10, pos = 50134510, gene = "SORBI_3010G170000"),      # DGAT1
  Nitrogen = list(chr = 1, pos = 4242712, gene = "SORBI_3001G056400"),     # Trp synthase
  Phosphorus = list(chr = 1, pos = 67183097, gene = "SORBI_3001G384300"),  # PHR1
  Cold = list(chr = 3, pos = 65113625, gene = "SORBI_3003G324400")         # AP37
)

# Marker SNP positions for SUM/RATIO phenotypes (rows 3-4)
markers_sr <- list(
  Lipid = list(chr = 1, pos = 55628277, gene = "SORBI_3001G283700"),       # GPAT
  Nitrogen = list(chr = 1, pos = 80515621, gene = "SORBI_3001G541900"),    # NRT1/PTR 6.3
  Phosphorus = list(chr = 2, pos = 50002465, gene = "SORBI_3002G161700"),  # PHL
  Cold = list(chr = 3, pos = 68071592, gene = "SORBI_3003G363400")         # ABI5 (bZIP)
)

# ──────────────────────────────────────────────────────────────────────────────
# FILE DEFINITIONS
# ──────────────────────────────────────────────────────────────────────────────

# Individual lipid files (2 per category)
ind_files <- list(
  Lipid = list.files(file.path(base_dir, "Lipid_ind_SORBI_3010G170000"),
                     pattern = "\\.txt$", full.names = TRUE)[1:2],
  Nitrogen = list.files(file.path(base_dir, "Nitrogen_ind_SORBI_3001G056400"),
                        pattern = "\\.txt$", full.names = TRUE)[1:2],
  Phosphorus = list.files(file.path(base_dir, "Phosphorus_ind_SORBI_3001G384300"),
                          pattern = "\\.txt$", full.names = TRUE)[1:2],
  Cold = list.files(file.path(base_dir, "Cold_ind_SORBI_3003G324400"),
                    pattern = "\\.txt$", full.names = TRUE)[1:2]
)

# Sum/Ratio files (2 per category)
sr_files <- list(
  Lipid = list.files(file.path(base_dir, "Lipid_SR_SORBI_3001G283700"),
                     pattern = "\\.txt$", full.names = TRUE)[1:2],
  Nitrogen = list.files(file.path(base_dir, "Nitrogen_SR_SORBI_3001G541900"),
                        pattern = "\\.txt$", full.names = TRUE)[1:2],
  Phosphorus = list.files(file.path(base_dir, "Phosphorus_SR_SORBI_3002G161700"),
                          pattern = "\\.txt$", full.names = TRUE)[1:2],
  Cold = list.files(file.path(base_dir, "Cold_SR_SORBI_3003G363400"),
                    pattern = "\\.txt$", full.names = TRUE)[1:2]
)

# Check that files exist
check_files <- function(file_list, name) {
  for (cat in names(file_list)) {
    files <- file_list[[cat]]
    if (length(files) < 2 || any(is.na(files))) {
      stop("Missing GWAS files for ", name, " - ", cat,
           "\nExpected directory: ", file.path(base_dir, paste0(cat, "_*")))
    }
  }
}

check_files(ind_files, "Individual")
check_files(sr_files, "Sum/Ratio")

message("All GWAS files found.")

# ──────────────────────────────────────────────────────────────────────────────
# HELPER FUNCTIONS
# ──────────────────────────────────────────────────────────────────────────────

read_assoc <- function(f) {
  message("  Reading: ", basename(f))
  df <- vroom::vroom(f, show_col_types = FALSE, col_select = c(chr, ps, p_wald))
  names(df) <- tolower(names(df))
  df %>%
    dplyr::transmute(
      chr = as.integer(chr),
      pos = as.numeric(ps),
      p = as.numeric(p_wald)
    ) %>%
    dplyr::filter(is.finite(chr) & is.finite(pos) & is.finite(p) & p > 0)
}

prep_manhattan <- function(df) {
  chr_info <- df %>%
    dplyr::group_by(chr) %>%
    dplyr::summarise(chr_len = max(pos, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(chr) %>%
    dplyr::mutate(chr_start = dplyr::lag(cumsum(chr_len), default = 0))

  df2 <- df %>%
    dplyr::left_join(chr_info, by = "chr") %>%
    dplyr::mutate(pos_cum = pos + chr_start, logp = -log10(p))

  axis_df <- chr_info %>%
    dplyr::mutate(center = chr_start + chr_len / 2)

  list(df = df2, axis = axis_df, chr_info = chr_info)
}

thin_for_plot <- function(df_plot, max_points = 400000L, cap_logp = 20) {
  df_plot <- df_plot %>% dplyr::filter(logp <= cap_logp)
  n <- nrow(df_plot)
  if (n <= max_points) return(df_plot)
  set.seed(42)
  keep <- sort(sample.int(n, max_points))
  df_plot[keep, , drop = FALSE]
}

label_from_file <- function(f) {
  lab <- basename(f)
  lab <- sub("_mod_sub_.*$", "", lab)
  lab <- sub("\\.txt$", "", lab)
  lab <- sub("_log10safe$", "", lab)
  lab <- gsub("_", "/", lab, fixed = TRUE)
  lab
}

prep_panel <- function(f, marker_ind, marker_sr, category, type) {
  df <- read_assoc(f)
  prep <- prep_manhattan(df)

  # Use correct marker based on type (Individual vs Sum/Ratio)
  marker <- if (type == "Individual") marker_ind else marker_sr

  # Calculate marker vline position
  vline_x <- marker$pos + prep$chr_info$chr_start[match(marker$chr, prep$chr_info$chr)]

  list(
    df = thin_for_plot(prep$df),
    axis = prep$axis,
    vline = vline_x,
    gene = marker$gene,
    title = label_from_file(f),
    category = category,
    type = type
  )
}

# ──────────────────────────────────────────────────────────────────────────────
# PREPARE ALL PANELS
# ──────────────────────────────────────────────────────────────────────────────

message("Preparing panels...")

categories <- c("Lipid", "Nitrogen", "Phosphorus", "Cold")
all_panels <- list()

for (cat in categories) {
  # Individual panels (rows 1-2) - use markers_ind
  for (i in 1:2) {
    panel <- prep_panel(ind_files[[cat]][i], markers_ind[[cat]], markers_sr[[cat]], cat, "Individual")
    all_panels[[length(all_panels) + 1]] <- panel
  }
  # Sum/Ratio panels (rows 3-4) - use markers_sr
  for (i in 1:2) {
    panel <- prep_panel(sr_files[[cat]][i], markers_ind[[cat]], markers_sr[[cat]], cat, "Sum/Ratio")
    all_panels[[length(all_panels) + 1]] <- panel
  }
}

# ──────────────────────────────────────────────────────────────────────────────
# DRAW FIGURE 3
# ──────────────────────────────────────────────────────────────────────────────

draw_grid <- function(outfile, device = c("png", "pdf")) {
  device <- match.arg(device)
  if (device == "png") {
    png(outfile, width = 4000, height = 3200, res = 300)
  } else {
    pdf(outfile, width = 16, height = 12)
  }

  # 4 rows x 4 columns
  par(mfrow = c(4, 4), mar = c(2.5, 3.5, 3.5, 0.8), oma = c(3, 3, 4, 1))

  for (row in 1:4) {
    for (col in 1:4) {
      idx <- (col - 1) * 4 + row  # Panel index
      p <- all_panels[[idx]]

      show_x <- row == 4
      show_y <- col == 1

      # Alternate chromosome colors
      pt_cols <- ifelse(p$df$chr %% 2 == 0, "grey65", "grey35")

      # Plot
      plot(
        p$df$pos_cum, p$df$logp,
        pch = 20, cex = 0.4, col = pt_cols,
        xlab = "", ylab = "",
        xaxt = "n", yaxt = "n",
        bty = "l",
        ylim = c(0, max(15, max(p$df$logp, na.rm = TRUE) + 1))
      )

      # Marker SNP vertical line
      abline(v = p$vline, lty = 2, lwd = 1.8, col = "#D95F02")

      # Significance thresholds
      abline(h = 7, lty = 2, lwd = 1.0, col = "red")
      abline(h = 5, lty = 3, lwd = 0.8, col = "#2C7FB8")

      # Gene label near SNP-line / threshold intersection
      x_rng <- range(p$df$pos_cum, na.rm = TRUE)
      x_span <- diff(x_rng)
      right_side <- is.finite(p$vline) && p$vline > (x_rng[1] + 0.70 * x_span)
      label_x <- if (right_side) p$vline - 0.01 * x_span else p$vline + 0.01 * x_span
      text(
        x = label_x, y = 7.25, labels = p$gene,
        cex = 1.35, font = 2, col = "#D95F02",
        adj = if (right_side) 1 else 0, xpd = NA
      )

      # Title (phenotype name)
      title(main = p$title, cex.main = 0.80, font.main = 1, line = 0.5)

      # X-axis (chromosome labels) - only bottom row
      if (show_x) {
        axis(1, at = p$axis$center, labels = p$axis$chr, cex.axis = 0.75, las = 1)
      }

      # Y-axis
      axis(2, cex.axis = 0.8, las = 1)
      if (show_y) {
        mtext(expression(-log[10](p)), side = 2, line = 2.5, cex = 0.7)
      }

      # Column headers (category names) - only top row
      if (row == 1) {
        mtext(categories[col], side = 3, line = 2.2, cex = 0.95, font = 2)
      }

      # Row type labels - only first column
      if (col == 1) {
        if (row <= 2) {
          mtext("Individual", side = 2, line = 4.5, cex = 0.7, font = 3)
        } else if (row == 3) {
          mtext("Sum/Ratio", side = 2, line = 4.5, cex = 0.7, font = 3)
        }
      }
    }
  }

  # Overall title
  mtext("GWAS Manhattan Plots by Functional Category",
        side = 3, outer = TRUE, line = 1.5, cex = 1.1, font = 2)

  # X-axis label
  mtext("Chromosome", side = 1, outer = TRUE, line = 1.5, cex = 0.9)

  dev.off()
  message("Saved: ", outfile)
}

# ──────────────────────────────────────────────────────────────────────────────
# GENERATE OUTPUT
# ──────────────────────────────────────────────────────────────────────────────

message("\nGenerating Figure 3...")

png_file <- file.path(out_dir, "Figure3_GWAS.png")

draw_grid(png_file, "png")

message("\n\u2713 Figure 3 (GWAS) complete!\n")

if (identical(Sys.getenv("ONLY_GWAS_MAIN", ""), "1")) {
  message("ONLY_GWAS_MAIN=1: completed main GWAS figure only.")
  quit(save = "no", status = 0)
}

# ══════════════════════════════════════════════════════════════════════════════
# SUPPLEMENTARY FIGURE S8: ALLELE EFFECTS FOR REPRESENTATIVE GWAS GENES
# ══════════════════════════════════════════════════════════════════════════════

message("\n══════════════════════════════════════════════════════════════")
message("OPTIONAL SUPPLEMENTARY FIGURE: ALLELE EFFECTS (BEST SNP x BEST PHENOTYPE)")
message("══════════════════════════════════════════════════════════════\n")

representative_genes <- tibble::tribble(
  ~Category,     ~GeneID,             ~GeneName,          ~PreferredSource,
  "Lipid",       "SORBI_3007G223801", "GPDH",             "Sum/Ratio",
  "Lipid",       "SORBI_3010G170000", "DGAT1",            "Individual",
  "Lipid",       "SORBI_3001G283700", "GPAT",             "Sum/Ratio",
  "Lipid",       "SORBI_3001G103800", "LCAT-like 1",      "Sum/Ratio",
  "Lipid",       "SORBI_3006G214500", "LCAT-like 4",      "Individual",
  "Nitrogen",    "SORBI_3001G541900", "NRT1/PTR 6.3",     "Sum/Ratio",
  "Nitrogen",    "SORBI_3001G056400", "Trp synthase",     "Individual",
  "Nitrogen",    "SORBI_3010G116400", "AA transporter",   "Individual",
  "Phosphorus",  "SORBI_3001G384300", "PHR1",             "Individual",
  "Phosphorus",  "SORBI_3001G186800", "PAP18",            "Individual",
  "Phosphorus",  "SORBI_3002G161700", "PHL",              "Sum/Ratio",
  "Cold",        "SORBI_3003G363400", "ABI5 (bZIP)",      "Sum/Ratio",
  "Cold",        "SORBI_3003G324400", "AP37",             "Individual"
) %>%
  dplyr::mutate(
    Category = factor(Category, levels = c("Lipid", "Nitrogen", "Phosphorus", "Cold")),
    GeneOrder = dplyr::row_number()
  )

ind_best <- vroom::vroom(
  "results/gene_count/individual/lowinput/gene_by_phenotype_bestP_p<=1e-07.tsv",
  show_col_types = FALSE
) %>%
  dplyr::mutate(Source = "Individual", Best_P_Value = as.numeric(Best_P_Value))

sr_best <- vroom::vroom(
  "results/gene_count/sum_ratio/lowinput/gene_by_phenotype_bestP_p<=1e-07.tsv",
  show_col_types = FALSE
) %>%
  dplyr::mutate(Source = "Sum/Ratio", Best_P_Value = as.numeric(Best_P_Value))

candidate_hits <- dplyr::bind_rows(ind_best, sr_best) %>%
  dplyr::inner_join(
    representative_genes %>% dplyr::select(GeneID, PreferredSource),
    by = "GeneID"
  ) %>%
  dplyr::filter(Source == PreferredSource) %>%
  dplyr::arrange(GeneID, Best_P_Value)

# Load phenotype matrices used in GWAS
pheno_ind <- vroom::vroom(
  "data/SPATS_fitted/BLUP_GWAS_phenotype/Final_lowinput_all_lipids_BLUPs.csv",
  show_col_types = FALSE
)
names(pheno_ind)[1] <- "Line"

pheno_sr <- vroom::vroom(
  "data/SPATS_fitted/BLUP_GWAS_phenotype/Final_lowinput_BLUPs_class_sums_and_ratios.csv",
  show_col_types = FALSE
)
names(pheno_sr)[1] <- "Line"

get_phenotype_values <- function(phenotype, source) {
  if (!is.character(phenotype) || length(phenotype) != 1 || is.na(phenotype)) {
    return(NULL)
  }
  mat <- if (identical(source, "Individual")) pheno_ind else pheno_sr

  pheno_col <- phenotype
  if (!(pheno_col %in% names(mat)) && identical(source, "Individual") && grepl("\\(", phenotype)) {
    inside <- sub("^.*\\((.*)\\)$", "\\1", phenotype)
    inside_slash <- gsub("_", "/", inside, fixed = TRUE)
    candidate <- sub("\\(.*\\)$", paste0("(", inside_slash, ")"), phenotype)
    if (candidate %in% names(mat)) {
      pheno_col <- candidate
    }
  }

  if (!(pheno_col %in% names(mat))) return(NULL)
  mat %>%
    dplyr::transmute(
      Line = as.character(Line),
      PhenotypeValue = as.numeric(.data[[pheno_col]])
    ) %>%
    dplyr::filter(is.finite(PhenotypeValue))
}

# SNP genotype extraction from HapMap files
genotype_dir <- "/Users/nirwantandukar/Documents/Research/data/SAP/genotype"

iupac_het <- list(
  R = c("A", "G"), Y = c("C", "T"), S = c("C", "G"),
  W = c("A", "T"), K = c("G", "T"), M = c("A", "C")
)

decode_call <- function(call, a1, a2) {
  call <- toupper(trimws(as.character(call)))
  if (!nzchar(call) || call %in% c("NA", "N", ".", "-", "0")) return(NA_character_)

  if (nchar(call) == 2) {
    cc <- strsplit(call, "", fixed = TRUE)[[1]]
    if (all(cc %in% c(a1, a2))) {
      if (all(cc == a1)) return("REF")
      if (all(cc == a2)) return("ALT")
      return("HET")
    }
    return(NA_character_)
  }

  if (nchar(call) == 1) {
    if (call == a1) return("REF")
    if (call == a2) return("ALT")
    if (!is.null(iupac_het[[call]]) &&
        setequal(sort(iupac_het[[call]]), sort(c(a1, a2)))) {
      return("HET")
    }
  }

  NA_character_
}

geno_cache <- new.env(parent = emptyenv())
pheno_cache <- new.env(parent = emptyenv())

extract_snp_genotypes <- function(snp_id, chr_num) {
  if (is.na(snp_id) || is.na(chr_num)) return(NULL)
  key <- paste0("chr", chr_num, "_", snp_id)
  if (exists(key, envir = geno_cache, inherits = FALSE)) {
    return(get(key, envir = geno_cache, inherits = FALSE))
  }

  hmp_file <- file.path(
    genotype_dir,
    sprintf("SAP_only_samples_bialleles_MAF_0.05_chr%s_filtered_80perc_hets_.hmp.txt", chr_num)
  )
  if (!file.exists(hmp_file)) return(NULL)

  header <- strsplit(readLines(hmp_file, n = 1), "\t", fixed = TRUE)[[1]]
  if (length(header) < 12) return(NULL)
  sample_ids <- header[12:length(header)]

  safe_snp <- gsub("[^A-Za-z0-9_]", "", snp_id)
  cmd <- sprintf("awk -F'\\t' '$1==\"%s\" {print; exit}' %s", safe_snp, shQuote(hmp_file))
  row_txt <- suppressWarnings(system(cmd, intern = TRUE))
  if (length(row_txt) == 0) return(NULL)

  fields <- strsplit(row_txt[[1]], "\t", fixed = TRUE)[[1]]
  if (length(fields) < 12) return(NULL)

  allele_field <- toupper(fields[2])
  allele_parts <- strsplit(allele_field, "/", fixed = TRUE)[[1]]
  if (length(allele_parts) < 2) return(NULL)
  a1 <- trimws(allele_parts[1])
  a2 <- trimws(allele_parts[2])
  if (!nzchar(a1) || !nzchar(a2)) return(NULL)

  calls <- fields[12:length(fields)]
  n <- min(length(sample_ids), length(calls))
  if (n == 0) return(NULL)

  out <- tibble::tibble(
    Line = sample_ids[seq_len(n)],
    Call = calls[seq_len(n)]
  ) %>%
    dplyr::mutate(
      GenotypeClass = vapply(Call, decode_call, character(1), a1 = a1, a2 = a2)
    ) %>%
    dplyr::filter(!is.na(GenotypeClass))

  assign(key, out, envir = geno_cache)
  out
}

build_panel_data <- function(phenotype, source, snp_id, chr_num) {
  pheno_key <- paste(source, phenotype, sep = "||")
  if (!exists(pheno_key, envir = pheno_cache, inherits = FALSE)) {
    assign(pheno_key, get_phenotype_values(phenotype, source), envir = pheno_cache)
  }
  pdat <- get(pheno_key, envir = pheno_cache, inherits = FALSE)
  gdat <- extract_snp_genotypes(snp_id, chr_num)
  if (is.null(pdat) || is.null(gdat)) return(NULL)

  merged <- gdat %>%
    dplyr::inner_join(pdat, by = "Line") %>%
    dplyr::filter(is.finite(PhenotypeValue))
  if (nrow(merged) == 0) return(NULL)

  merged %>%
    dplyr::mutate(GenotypeClass = factor(GenotypeClass, levels = c("REF", "HET", "ALT"))) %>%
    dplyr::filter(!is.na(GenotypeClass))
}

has_sig_pairwise <- function(dat, alpha = 0.05) {
  if (is.null(dat) || nrow(dat) < 3) return(FALSE)
  x <- dat %>% dplyr::filter(!is.na(GenotypeClass), is.finite(PhenotypeValue))
  if (nrow(x) < 3 || dplyr::n_distinct(x$GenotypeClass) < 2) return(FALSE)

  pw <- tryCatch(
    stats::pairwise.wilcox.test(
      x = x$PhenotypeValue,
      g = droplevels(x$GenotypeClass),
      p.adjust.method = "BH",
      exact = FALSE
    )$p.value,
    error = function(e) NULL
  )
  if (is.null(pw)) return(FALSE)
  any(is.finite(pw) & pw < alpha, na.rm = TRUE)
}

pick_gene_hit <- function(gene_rows) {
  candidates <- vector("list", nrow(gene_rows))
  for (i in seq_len(nrow(gene_rows))) {
    r <- gene_rows[i, ]
    dat <- build_panel_data(r$Phenotype, r$Source, r$Best_SNP, r$Chromosome)
    if (is.null(dat)) next

    n_ref <- sum(dat$GenotypeClass == "REF", na.rm = TRUE)
    n_het <- sum(dat$GenotypeClass == "HET", na.rm = TRUE)
    n_alt <- sum(dat$GenotypeClass == "ALT", na.rm = TRUE)
    n_groups <- dplyr::n_distinct(dat$GenotypeClass[!is.na(dat$GenotypeClass)])

    if (n_groups < 2) next

    alt_rank <- if (n_alt >= 5) 1L else if (n_alt >= 2) 2L else 3L
    sig_any <- has_sig_pairwise(dat, alpha = 0.05)
    candidates[[i]] <- list(
      hit = r,
      dat = dat,
      n_ref = n_ref,
      n_het = n_het,
      n_alt = n_alt,
      alt_rank = alt_rank,
      sig_any = sig_any
    )
  }
  candidates <- candidates[!vapply(candidates, is.null, logical(1))]
  if (length(candidates) == 0) return(NULL)

  sig_candidates <- candidates[vapply(candidates, function(x) isTRUE(x$sig_any), logical(1))]
  pool <- if (length(sig_candidates) > 0) sig_candidates else candidates

  best_idx <- order(
    vapply(pool, function(x) x$alt_rank, integer(1)),
    vapply(pool, function(x) x$hit$Best_P_Value, numeric(1))
  )[1]
  pool[[best_idx]]
}

selected_hits <- vector("list", nrow(representative_genes))
panel_data <- vector("list", nrow(representative_genes))

for (i in seq_len(nrow(representative_genes))) {
  g <- representative_genes[i, ]
  gene_rows <- candidate_hits %>%
    dplyr::filter(GeneID == g$GeneID) %>%
    dplyr::arrange(Best_P_Value) %>%
    dplyr::group_by(Best_SNP, Chromosome) %>%
    dplyr::slice_min(order_by = Best_P_Value, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::slice_head(n = 80)

  pick <- pick_gene_hit(gene_rows)
  if (is.null(pick)) {
    selected_hits[[i]] <- dplyr::bind_cols(
      g %>% dplyr::select(Category, GeneID, GeneName, PreferredSource, GeneOrder),
      tibble::tibble(
        Phenotype = NA_character_,
        Best_SNP = NA_character_,
        Best_P_Value = NA_real_,
        Source = NA_character_,
        Chromosome = NA_integer_,
        n_REF = NA_integer_,
        n_HET = NA_integer_,
        n_ALT = NA_integer_,
        ALT_selection = "No candidate"
      )
    )
    panel_data[[i]] <- NULL
  } else {
    selected_hits[[i]] <- dplyr::bind_cols(
      g %>% dplyr::select(Category, GeneID, GeneName, PreferredSource, GeneOrder),
      pick$hit %>% dplyr::select(Phenotype, Best_SNP, Best_P_Value, Source, Chromosome),
      tibble::tibble(
        n_REF = pick$n_ref,
        n_HET = pick$n_het,
        n_ALT = pick$n_alt,
        ALT_selection = dplyr::case_when(
          pick$n_alt >= 5 ~ "ALT>=5",
          pick$n_alt >= 2 ~ "ALT>=2 fallback",
          TRUE ~ "ALT<2 fallback"
        ),
        Has_sig_pairwise = pick$sig_any
      )
    )
    panel_data[[i]] <- pick$dat
  }
}

best_hits <- dplyr::bind_rows(selected_hits) %>%
  dplyr::arrange(Category, GeneOrder)

message("Selected allele-effect hits (significance-prioritized + ALT fallback):")
print(
  best_hits %>% dplyr::select(Category, GeneID, GeneName, Source, Phenotype, Best_SNP, Best_P_Value, n_ALT, ALT_selection, Has_sig_pairwise),
  n = nrow(best_hits)
)

blank_panel <- function(title_txt, subtitle_txt = NULL) {
  sub_txt <- if (is.null(subtitle_txt)) "" else subtitle_txt
  ggplot() +
    annotate("text", x = 0.5, y = 0.62, label = title_txt, fontface = "bold", size = 4) +
    annotate("text", x = 0.5, y = 0.45, label = sub_txt, size = 3.2, color = "grey30") +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "grey96", color = "grey80"),
      plot.margin = margin(6, 6, 6, 6)
    )
}

mk_panel <- function(hit, dat) {
  if (is.null(dat) || nrow(dat) < 5) {
    missing_reason <- if (!is.na(hit$ALT_selection) && grepl("No candidate", hit$ALT_selection, fixed = TRUE)) {
      hit$ALT_selection
    } else {
      "Insufficient genotype/phenotype overlap"
    }
    return(blank_panel(hit$GeneName, missing_reason))
  }

  # Keep only genotype levels present in this panel
  present <- c("REF", "HET", "ALT")[c("REF", "HET", "ALT") %in% as.character(unique(dat$GenotypeClass))]
  dat2 <- dat %>%
    dplyr::filter(as.character(GenotypeClass) %in% present) %>%
    dplyr::mutate(
      GenotypeClass = factor(as.character(GenotypeClass), levels = present)
    )

  if (nlevels(dat2$GenotypeClass) < 2) {
    return(blank_panel(hit$GeneName, "Only one genotype group"))
  }

  p <- tryCatch(
    ggstatsplot::ggbetweenstats(
      data = dat2,
      x = GenotypeClass,
      y = PhenotypeValue,
      type = "n",
      centrality.plotting = FALSE,
      pairwise.comparisons = TRUE,
      pairwise.display = "significant",
      p.adjust.method = "BH",
      messages = FALSE,
      bf.message = FALSE,
      xlab = "Alleles",
      ylab = "Phenotype (BLUP)",
      title = paste0(hit$GeneName, " | ", hit$Phenotype),
      subtitle = paste0(
        hit$Best_SNP, " | p=", format(hit$Best_P_Value, scientific = TRUE, digits = 2),
        " | ALT n=", hit$n_ALT
      ),
      results.subtitle = TRUE
    ),
    error = function(e) {
      ggplot(dat2, aes(x = GenotypeClass, y = PhenotypeValue, fill = GenotypeClass)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_jitter(width = 0.12, alpha = 0.6, size = 1.2) +
        labs(
          title = paste0(hit$GeneName, " | ", hit$Phenotype),
          subtitle = paste0(hit$Best_SNP, " | ALT n=", hit$n_ALT),
          x = "Alleles", y = "Phenotype (BLUP)"
        ) +
        theme_minimal(base_size = 9)
    }
  )

  p +
    theme(
      text = element_text(size = 6.8),
      plot.title = element_text(size = 8.3, face = "bold"),
      plot.subtitle = element_text(size = 6.4),
      axis.title = element_text(size = 7.0),
      axis.text = element_text(size = 6.6),
      plot.margin = margin(6, 6, 6, 6)
    )
}

# Layout with category columns
cats <- c("Lipid", "Nitrogen", "Phosphorus", "Cold")
n_rows <- max(table(representative_genes$Category))
n_cols <- length(cats)

plot_list <- list()
for (row in seq_len(n_rows)) {
  for (col in seq_len(n_cols)) {
    cat_name <- cats[col]
    idx_cat <- which(best_hits$Category == cat_name)
    if (row > length(idx_cat)) {
      plot_list[[length(plot_list) + 1L]] <- blank_panel(cat_name, "No panel")
    } else {
      idx <- idx_cat[row]
      plot_list[[length(plot_list) + 1L]] <- mk_panel(best_hits[idx, ], panel_data[[best_hits$GeneOrder[idx]]])
    }
  }
}

fig_s8 <- patchwork::wrap_plots(plot_list, ncol = n_cols, byrow = TRUE) +
  patchwork::plot_annotation(
    title = "Allele Effects for Representative GWAS Genes (Best SNP x Best Phenotype, LowInput)"
  )

supp_png <- "fig/supp/SuppFig_GWAS_Allele_Effects_optional.png"
ggsave(supp_png, fig_s8, width = 16, height = 18, dpi = 300, bg = "white")

message("Saved: ", supp_png)
message("\n\u2713 Optional GWAS allele-effects figure complete!\n")


# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════════════════

message("\n")
message("══════════════════════════════════════════════════════════════")
message("SECTION 5 COMPLETE!")
message("══════════════════════════════════════════════════════════════")
message("\n")
message("Generated outputs:")
message("  - fig/main/Figure3_GWAS.png")
message("  - fig/supp/SuppFig_GWAS_Allele_Effects_optional.png")
message("\n")
message("Note: Table 1 (GWAS candidate genes) is in the LaTeX document.")
message("Note: Supplementary Tables S8-S10 are generated by annotation scripts")
message("      (22.3_Annotation_table_GWAS.R and related scripts)")
message("\n")
