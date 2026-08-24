suppressPackageStartupMessages({
  library(vroom)
  library(dplyr)
})

message("\n==============================================================================")
message("STACKED GWAS MANHATTAN PLOTS - SEPARATE CONTROL AND LOWINPUT")
message("==============================================================================\n")

base_dir <- "/Users/nirwantandukar/Documents/Github/SoLD_paper/results/gwas_raw"
out_dir <- "/Users/nirwantandukar/Documents/Github/SoLD_paper/fig/main"
thesis_out_dir <- "/Users/nirwantandukar/Documents/Github/Thesis_NirwanTandukar_GeneticsGenomics/ncsuthesis-0.6/Figures/Chapter2/main"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(thesis_out_dir, recursive = TRUE, showWarnings = FALSE)

tested_snps <- c(
  CTL = 6106358,
  LIN = 6085245
)

thresholds <- lapply(tested_snps, function(n_snps) {
  strict_p <- 0.05 / n_snps
  list(
    strict_p = strict_p,
    strict_logp = -log10(strict_p),
    suggestive_logp = 7
  )
})

panel_groups <- list(
  "CTL individual 1" = list(
    folder = "CTL_individual_lipid_SORBI_3006G248300_SNP_58780498_chr6",
    gene_id = "SORBI_3006G248300",
    marker = list(chr = 6L, pos = 58780498, snp = "SNP_58780498")
  ),
  "CTL individual 2" = list(
    folder = "CTL_individual_lipid_SORBI_3001G276500_SNP_53871418_chr1",
    gene_id = "SORBI_3001G276500",
    marker = list(chr = 1L, pos = 53871418, snp = "SNP_53871418")
  ),
  "CTL sum/ratio 1" = list(
    folder = "CTL_SR_SORBI_3002G373700_SNP_73174964_chr2",
    gene_id = "SORBI_3002G373700",
    marker = list(chr = 2L, pos = 73174964, snp = "SNP_73174964")
  ),
  "CTL sum/ratio 2" = list(
    folder = "CTL_SR_SORBI_3003G141800_SNP_14021625_chr3",
    gene_id = "SORBI_3003G141800",
    marker = list(chr = 3L, pos = 14021625, snp = "SNP_14021625")
  ),
  "LIN individual 1" = list(
    folder = "LIN_individual_SORBI_3001G384300_SNP_67183097_chr1",
    gene_id = "SORBI_3001G384300",
    marker = list(chr = 1L, pos = 67183097, snp = "SNP_67183097")
  ),
  "LIN individual 2" = list(
    folder = "LIN_individual_SORBI_3004G107800_SNP_10366823_chr4",
    gene_id = "SORBI_3004G107800",
    marker = list(chr = 4L, pos = 10366823, snp = "SNP_10366823")
  ),
  "LIN sum/ratio 1" = list(
    folder = "LIN_SR_SORBI_3001G222700_SNP_21276112_chr1",
    gene_id = "SORBI_3001G222700",
    marker = list(chr = 1L, pos = 21276112, snp = "SNP_21276112")
  ),
  "LIN sum/ratio 2" = list(
    folder = "LIN_SR_SORBI_3001G256400_SNP_29811579_chr1",
    gene_id = "SORBI_3001G256400",
    marker = list(chr = 1L, pos = 29811579, snp = "SNP_29811579")
  )
)

group_files <- list(
  "CTL individual 1" = file.path(
    base_dir,
    panel_groups[["CTL individual 1"]]$folder,
    c(
      "TG(10:0_10:0_10:0)_mod_sub_Final_control_all_lipids_BLUPs.part13_SAP_bialleles_MAF_0.05_2.assoc.txt",
      "TG(12:0_18:1_18:2)_mod_sub_Final_control_all_lipids_BLUPs.part7_SAP_bialleles_MAF_0.05_2.assoc.txt",
      "TG(16:0_18:1_18:1)_mod_sub_Final_control_all_lipids_BLUPs.part24_SAP_bialleles_MAF_0.05_2.assoc.txt",
      "TG(16:0_18:1_18:2)_mod_sub_Final_control_all_lipids_BLUPs.part25_SAP_bialleles_MAF_0.05_2.assoc.txt"
    )
  ),
  "CTL individual 2" = file.path(
    base_dir,
    panel_groups[["CTL individual 2"]]$folder,
    c(
      "TG(12:0_16:0_18:2)_mod_sub_Final_control_all_lipids_BLUPs.part4_SAP_bialleles_MAF_0.05_2.assoc.txt",
      "TG(16:1_18:1_18:1)_mod_sub_Final_control_all_lipids_BLUPs.part10_SAP_bialleles_MAF_0.05_2.assoc.txt",
      "TG(18:0_18:2_18:2)_mod_sub_Final_control_all_lipids_BLUPs.part14_SAP_bialleles_MAF_0.05_2.assoc.txt",
      "TG(18:1_18:2_18:2)_mod_sub_Final_control_all_lipids_BLUPs.part18_SAP_bialleles_MAF_0.05_2.assoc.txt"
    )
  ),
  "CTL sum/ratio 1" = file.path(
    base_dir,
    panel_groups[["CTL sum/ratio 1"]]$folder,
    c(
      "Sum_DGDG_over_PC_log10safe_mod_sub_Final_control_BLUPs_class_sums_and_ratios_wide.part01.parta1_SAP_bialleles_MAF_0.05_11.assoc.txt",
      "Sum_MGDG_over_DGDG_log10safe_mod_sub_Final_control_BLUPs_class_sums_and_ratios_wide.part02.parta3_SAP_bialleles_MAF_0.05_12.assoc.txt",
      "Sum_MG_over_DG_log10safe_mod_sub_Final_control_BLUPs_class_sums_and_ratios_wide.part08.parta7_SAP_bialleles_MAF_0.05_4.assoc.txt",
      "Sum_SQDG_over_PC_log10safe_mod_sub_Final_control_BLUPs_class_sums_and_ratios_wide.part02.parta2_SAP_bialleles_MAF_0.05_12.assoc.txt"
    )
  ),
  "CTL sum/ratio 2" = file.path(
    base_dir,
    panel_groups[["CTL sum/ratio 2"]]$folder,
    c(
      "Sum_DGDG_over_PC_log10safe_mod_sub_Final_control_BLUPs_class_sums_and_ratios_wide.part01.parta1_SAP_bialleles_MAF_0.05_11.assoc.txt",
      "Sum_MG_over_LPC_log10safe_mod_sub_Final_control_BLUPs_class_sums_and_ratios_wide.part07.parta3_SAP_bialleles_MAF_0.05_3.assoc.txt",
      "Sum_PC_over_PE_log10safe_mod_sub_Final_control_BLUPs_class_sums_and_ratios_wide.part08.parta6_SAP_bialleles_MAF_0.05_4.assoc.txt",
      "Sum_SQDG_over_PS_log10safe_mod_sub_Final_control_BLUPs_class_sums_and_ratios_wide.part03.parta4_SAP_bialleles_MAF_0.05_13.assoc.txt"
    )
  ),
  "LIN individual 1" = file.path(
    base_dir,
    panel_groups[["LIN individual 1"]]$folder,
    c(
      "PC(16:0_20:3).txt",
      "PC(16:0_22:5)_mod_sub_Final_lowinput_all_lipids_BLUPs.part6_SAP_bialleles_MAF_0.05.assoc.txt",
      "PC(18:1_20:1)_mod_sub_Final_lowinput_all_lipids_BLUPs.part22_SAP_bialleles_MAF_0.05.assoc.txt",
      "PE(16:0_18:1).txt"
    )
  ),
  "LIN individual 2" = file.path(
    base_dir,
    panel_groups[["LIN individual 2"]]$folder,
    c(
      "DG(8:0_8:0)_mod_sub_Final_lowinput_all_lipids_BLUPs.part9_SAP_bialleles_MAF_0.05.assoc.txt",
      "FA(22:1)_mod_sub_Final_lowinput_all_lipids_BLUPs.part15_SAP_bialleles_MAF_0.05.assoc.txt",
      "PC(18:2_20:4)_mod_sub_Final_lowinput_all_lipids_BLUPs.part3_SAP_bialleles_MAF_0.05.assoc.txt",
      "SQDG(16:0_18:1)_mod_sub_Final_lowinput_all_lipids_BLUPs.part6_SAP_bialleles_MAF_0.05.assoc.txt"
    )
  ),
  "LIN sum/ratio 1" = file.path(
    base_dir,
    panel_groups[["LIN sum/ratio 1"]]$folder,
    c(
      "Sum_DG_over_PC_log10safe_mod_sub_Final_lowinput_BLUPs_class_sums_and_ratios_wide.part07.parta1_SAP_bialleles_MAF_0.05_7.assoc.txt",
      "Sum_LPE_over_TG_log10safe_mod_sub_Final_lowinput_BLUPs_class_sums_and_ratios_wide.part04.parta1_SAP_bialleles_MAF_0.05_4.assoc.txt",
      "Sum_PS_over_DGDG_log10safe_mod_sub_Final_lowinput_BLUPs_class_sums_and_ratios_wide.part10.parta1_SAP_bialleles_MAF_0.05_10.assoc.txt",
      "Sum_PS_over_MGDG_log10safe_mod_sub_Final_lowinput_BLUPs_class_sums_and_ratios_wide.part10.parta1_SAP_bialleles_MAF_0.05_10.assoc.txt"
    )
  ),
  "LIN sum/ratio 2" = file.path(
    base_dir,
    panel_groups[["LIN sum/ratio 2"]]$folder,
    c(
      "Sum_DG_over_DGDG_log10safe_mod_sub_Final_lowinput_BLUPs_class_sums_and_ratios_wide.part07.parta1_SAP_bialleles_MAF_0.05_7.assoc.txt",
      "Sum_LPC_over_PC_log10safe_mod_sub_Final_lowinput_BLUPs_class_sums_and_ratios_wide.part04.parta1_SAP_bialleles_MAF_0.05_4.assoc.txt",
      "Sum_LPE_over_PE_log10safe_mod_sub_Final_lowinput_BLUPs_class_sums_and_ratios_wide.part04.parta1_SAP_bialleles_MAF_0.05_4.assoc.txt",
      "Sum_PS_over_PA_log10safe_mod_sub_Final_lowinput_BLUPs_class_sums_and_ratios_wide.part10.parta1_SAP_bialleles_MAF_0.05_10.assoc.txt"
    )
  )
)

read_assoc <- function(f) {
  message("  Reading: ", basename(f))
  df <- vroom::vroom(f, show_col_types = FALSE, col_select = c(chr, ps, p_wald))
  names(df) <- tolower(names(df))
  df %>%
    transmute(
      chr = as.integer(chr),
      pos = as.numeric(ps),
      p = as.numeric(p_wald)
    ) %>%
    filter(is.finite(chr) & is.finite(pos) & is.finite(p) & p > 0)
}

prep_manhattan <- function(df) {
  chr_info <- df %>%
    group_by(chr) %>%
    summarise(chr_len = max(pos, na.rm = TRUE), .groups = "drop") %>%
    arrange(chr) %>%
    mutate(chr_start = lag(cumsum(chr_len), default = 0))

  df2 <- df %>%
    left_join(chr_info, by = "chr") %>%
    mutate(pos_cum = pos + chr_start, logp = -log10(p))

  axis_df <- chr_info %>%
    mutate(center = chr_start + chr_len / 2)

  list(df = df2, axis = axis_df, chr_info = chr_info)
}

label_from_file <- function(f) {
  lab <- basename(f)
  lab <- sub("_mod_sub_.*$", "", lab)
  lab <- sub("\\.assoc\\.txt$", "", lab)
  lab <- sub("\\.txt$", "", lab)
  lab <- sub("_log10safe$", "", lab)
  lab <- gsub("_over_", " / ", lab, fixed = TRUE)
  lab <- gsub("_", "/", lab, fixed = TRUE)
  lab
}

thin_for_plot <- function(df_plot, max_points = 400000L, cap_logp = 20) {
  df_plot <- df_plot %>% filter(logp <= cap_logp)
  n <- nrow(df_plot)
  if (n <= max_points) return(df_plot)
  set.seed(42)
  keep <- sort(sample.int(n, max_points))
  df_plot[keep, , drop = FALSE]
}

prep_panel <- function(f, marker, group_name) {
  df <- read_assoc(f)
  prep <- prep_manhattan(df)
  chr_idx <- match(marker$chr, prep$chr_info$chr)
  vline_x <- if (is.na(chr_idx)) NA_real_ else marker$pos + prep$chr_info$chr_start[chr_idx]

  list(
    df = thin_for_plot(prep$df),
    axis = prep$axis,
    vline = vline_x,
    marker_label = marker$snp,
    title = label_from_file(f),
    group = group_name
  )
}

message("Preparing panels...")

prepare_condition_panels <- function(group_names) {
  file_counts <- vapply(group_names, function(group_name) length(group_files[[group_name]]), integer(1))
  if (length(unique(file_counts)) != 1L) {
    stop("Groups within a condition must have the same number of .txt files. Counts: ",
         paste(group_names, file_counts, sep = "=", collapse = ", "))
  }

  n_rows <- as.integer(unique(file_counts)[1])
  panels <- vector("list", length = length(group_names) * n_rows)
  idx <- 1L
  for (group_name in group_names) {
    marker <- panel_groups[[group_name]]$marker
    for (f in group_files[[group_name]]) {
      panels[[idx]] <- prep_panel(f, marker, group_name)
      idx <- idx + 1L
    }
  }

  list(
    panels = panels,
    n_rows = n_rows
  )
}

draw_condition_grid <- function(outfile, group_names, condition_name, device = c("png", "pdf")) {
  device <- match.arg(device)
  prep <- prepare_condition_panels(group_names)
  all_panels <- prep$panels
  n_rows <- prep$n_rows
  thr <- thresholds[[condition_name]]

  # Global sizing tuned for manuscript readability in dense multi-panel layouts.
  point_cex <- 0.78
  snp_label_cex <- 1.55
  panel_title_cex <- 1.90
  axis_tick_cex <- 1.75
  ylab_cex <- 1.70
  header_cex <- 1.55
  header_gene_cex <- 1.35
  xlab_cex <- 1.80

  if (device == "png") {
    # Taller canvas makes each Manhattan panel less squat ("thinner" x-span look).
    png(outfile, width = 5600, height = 5200, res = 300)
  } else {
    pdf(outfile, width = 18, height = 16.5)
  }

  par(mfrow = c(as.integer(n_rows), as.integer(length(group_names))),
      mar = c(4.6, 6.0, 8.2, 1.5),
      oma = c(4.4, 4.2, 1.2, 1.0),
      cex.axis = axis_tick_cex,
      lwd = 1.5)

  group_headers <- vapply(
    group_names,
    function(group_name) {
      clean_group <- gsub(" [0-9]+$", "", group_name)
      clean_group
    },
    character(1)
  )
  gene_headers <- vapply(
    group_names,
    function(group_name) panel_groups[[group_name]]$gene_id,
    character(1)
  )

  panel_letters <- LETTERS[seq_along(group_names)]

  for (row in seq_len(n_rows)) {
    for (col in seq_along(group_names)) {
      panel_idx <- (col - 1L) * n_rows + row
      p <- all_panels[[panel_idx]]

      show_x <- row == n_rows
      show_y <- col == 1L
      pt_cols <- ifelse(p$df$chr %% 2 == 0, "grey65", "grey35")

      plot(
        p$df$pos_cum, p$df$logp,
        pch = 20, cex = point_cex, col = pt_cols,
        xlab = "", ylab = "",
        xaxt = "n", yaxt = "n",
        bty = "l",
        ylim = c(0, max(15, max(p$df$logp, na.rm = TRUE) + 1))
      )

      if (is.finite(p$vline)) {
        abline(v = p$vline, lty = 2, lwd = 2.6, col = "#D95F02")

        x_rng <- range(p$df$pos_cum, na.rm = TRUE)
        x_span <- diff(x_rng)
        right_side <- p$vline > (x_rng[1] + 0.70 * x_span)
        label_x <- if (right_side) p$vline - 0.01 * x_span else p$vline + 0.01 * x_span
        text(
          x = label_x, y = thr$strict_logp + 0.22, labels = p$marker_label,
          cex = snp_label_cex, font = 2, col = "#D95F02",
          adj = if (right_side) 1 else 0, xpd = NA
        )
      }

      abline(h = thr$strict_logp, lty = 2, lwd = 1.9, col = "red")
      abline(h = thr$suggestive_logp, lty = 3, lwd = 1.5, col = "#2C7FB8")

      title(main = p$title, cex.main = panel_title_cex, font.main = 1, line = 1.15)

      if (show_x) {
        axis(1, at = p$axis$center, labels = p$axis$chr, cex.axis = axis_tick_cex, las = 1, lwd = 1.4)
      }

      axis(2, cex.axis = axis_tick_cex, las = 1, lwd = 1.4)
      if (show_y) {
        mtext(expression(-log[10](p)), side = 2, line = 3.4, cex = ylab_cex)
      }

      if (row == 1L) {
        mtext(
          paste0(panel_letters[col], "  ", group_headers[col]),
          side = 3, line = 4.6, cex = header_cex, font = 2
        )
        mtext(
          gene_headers[col],
          side = 3, line = 3.0, cex = header_gene_cex, font = 2
        )
      }
    }
  }

  mtext("Chromosome", side = 1, outer = TRUE, line = 2.0, cex = xlab_cex)

  dev.off()
  message("Saved: ", outfile)
}

conditions <- list(
  CTL = c("CTL individual 1", "CTL individual 2", "CTL sum/ratio 1", "CTL sum/ratio 2"),
  LIN = c("LIN individual 1", "LIN individual 2", "LIN sum/ratio 1", "LIN sum/ratio 2")
)

outputs <- list(
  CTL = list(
    png = file.path(out_dir, "Figure2_GWAS_CTL.png"),
    thesis_png = file.path(thesis_out_dir, "Figure2_GWAS_CTL.png")
  ),
  LIN = list(
    png = file.path(out_dir, "Figure3_GWAS_LIN.png"),
    thesis_png = file.path(thesis_out_dir, "Figure3_GWAS_LIN.png")
  )
)

for (condition_name in names(conditions)) {
  message("\nGenerating ", condition_name, " stacked GWAS figure...")
  message(sprintf(
    "  Thresholds: strict p=%.3e (-log10=%.3f), suggestive -log10=%.1f",
    thresholds[[condition_name]]$strict_p,
    thresholds[[condition_name]]$strict_logp,
    thresholds[[condition_name]]$suggestive_logp
  ))
  draw_condition_grid(outputs[[condition_name]]$png, conditions[[condition_name]], condition_name, "png")
  file.copy(outputs[[condition_name]]$png, outputs[[condition_name]]$thesis_png, overwrite = TRUE)
  message("Copied to thesis figures: ", outputs[[condition_name]]$thesis_png)
}

message("\n✓ Condition-specific GWAS figures complete!")
for (condition_name in names(outputs)) {
  message("  - ", outputs[[condition_name]]$png)
}
