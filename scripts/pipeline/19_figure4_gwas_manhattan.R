# ==============================================================================
# Figure 4 -- Recurrent GWAS signals under CTL and LIN
#
# Four columns, one per condition x trait layer. Each column stacks four
# phenotypes that share the same candidate locus, so the recurring peak is
# visible down the column. A dashed vertical line marks the shared locus.
# Remaining GWAS results are explorable in the SoLD Shiny app.
#
#   A  CTL individual lipids   -- SORBI_3006G040500, phospholipid-transporting ATPase (chr6)
#   B  CTL class sums/ratios   -- SORBI_3010G151000, CMP-sialic acid transporter 1 (chr10)
#   C  LIN individual lipids   -- SORBI_3010G170000, DGAT1 / O-acyltransferase (chr10)
#   D  LIN class sums/ratios   -- SORBI_3006G195200, purine permease array (chr6)
#
# Input : GEMMA .assoc output (tab-delimited: chr rs ps ... p_wald)
# Output: fig/main/Figure4_GWAS_Manhattan.png
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

gwas_root <- "/Users/nirwantandukar/Documents/Research/data/SAP/GWAS_result"
out_file  <- "/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database/fig/main/Figure4_GWAS_Manhattan.png"

columns <- list(
  list(tag = "CTL", layer = "Individual lipids", dir = "CTL_ind",
       gene = "SORBI_3006G040500", gene_lab = "Phospholipid-transporting ATPase",
       chr = 6L, pos = 25522678,
       traits = c("TG(16:0_16:0_18:1)", "TG(16:0_18:0_18:2)",
                  "TG(18:0_18:2_18:2)", "TG(18:1_18:2_18:2)")),
  list(tag = "CTL", layer = "Class sums / ratios", dir = "CTL_sumratio",
       gene = "SORBI_3010G151000", gene_lab = "CMP-sialic acid transporter 1",
       chr = 10L, pos = 43625333,
       traits = c("Sum_AEG_over_DGDG", "Sum_PC_over_DG",
                  "Sum_PC_over_MGDG", "Sum_PC_over_SQDG")),
  list(tag = "LIN", layer = "Individual lipids", dir = "LIN_ind",
       gene = "SORBI_3010G170000", gene_lab = "DGAT1 (O-acyltransferase)",
       chr = 10L, pos = 50099836,
       traits = c("TG(18:1_18:3_22:0)", "TG(18:2_18:2_18:4)",
                  "TG(18:2_20:3_22:0)", "TG(18:3_18:3_18:3)")),
  list(tag = "LIN", layer = "Class sums / ratios", dir = "LIN_sumratio",
       gene = "SORBI_3006G195200", gene_lab = "Purine permease array",
       chr = 6L, pos = 54823646,
       traits = c("Sum_FA", "Sum_FA_over_DG", "Sum_FA_over_SM", "Sum_FA_over_SQDG"))
)

bonferroni <- 8.19e-9
suggestive <- 1e-5

# Rendering only: all SNPs below keep_all_below are drawn, the rest subsampled
# 1-in-thin_every. Set thin_every <- 1 to plot every SNP.
keep_all_below <- 0.01
thin_every     <- 25

condition_colors <- c(CTL = "#440154FF", LIN = "#FDE725FF")
marker_col       <- "#D55E00"   # Okabe-Ito vermillion, colourblind-safe

plot_theme <- theme_minimal(base_size = 13) +
  theme(axis.text = element_text(colour = "black", size = 7),
        axis.title = element_text(face = "bold", size = 9),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        plot.tag = element_text(face = "bold", size = 15),
        plot.margin = margin(3, 5, 3, 5))

# chromosome offsets are shared across panels so peaks line up down a column
read_gwas <- function(path) {
  d <- fread(path, select = c("chr", "ps", "p_wald"), showProgress = FALSE)
  setnames(d, c("chr", "ps", "p"))
  d <- d[is.finite(p) & p > 0]
  d[, chr := suppressWarnings(as.integer(chr))]
  d[!is.na(chr) & chr >= 1 & chr <= 10]
}

message("scanning chromosome lengths ...")
first <- read_gwas(file.path(gwas_root, columns[[1]]$dir, paste0(columns[[1]]$traits[1], ".txt")))
chrlen <- first[, .(len = max(ps)), by = chr][order(chr)]
chrlen[, offset := cumsum(as.numeric(len)) - as.numeric(len)]
axis_df <- chrlen[, .(chr, centre = offset + len / 2)]
rm(first); invisible(gc())

gpos_of <- function(ch, p) chrlen[chr == ch, offset] + p
genome_end <- chrlen[, max(offset + len)]
# keep the gene label inside the panel: right-align near the end of the genome,
# left-align near the start, centre otherwise
label_hjust <- function(g) { f <- g / genome_end
  if (f > 0.85) 1.02 else if (f < 0.15) -0.02 else 0.5 }

panel <- function(d, tag, trait, marker_g, show_x, gene_lab) {
  d <- merge(d, chrlen[, .(chr, offset)], by = "chr")
  d[, gpos := as.numeric(ps) + offset][, logp := -log10(p)]
  keep <- d$p < keep_all_below | (seq_len(nrow(d)) %% thin_every == 0L)
  d <- d[keep]
  base <- condition_colors[[tag]]
  alt  <- colorRampPalette(c(base, "grey30"))(3)[2]
  d[, band := factor(chr %% 2)]

  ggplot(d, aes(gpos, logp, colour = band)) +
    geom_vline(xintercept = marker_g, colour = marker_col, linewidth = .7, alpha = .9) +
    geom_point(size = .35, alpha = .75) +
    geom_hline(yintercept = -log10(bonferroni), linetype = "dashed",
               colour = "grey20", linewidth = .4) +
    geom_hline(yintercept = -log10(suggestive), linetype = "dotted",
               colour = "grey50", linewidth = .35) +
    scale_colour_manual(values = c("0" = base, "1" = alt)) +
    scale_x_continuous(breaks = axis_df$centre, labels = axis_df$chr,
                       expand = expansion(mult = .01)) +
    scale_y_continuous(expand = expansion(mult = c(0, .10))) +
    annotate("text", x = -Inf, y = Inf, label = trait, hjust = -0.06, vjust = 1.4,
             size = 2.5, fontface = "italic", colour = "grey20") +
    {if (!is.na(gene_lab))
       annotate("text", x = marker_g, y = Inf, label = gene_lab, vjust = -0.55,
                hjust = label_hjust(marker_g),
                size = 2.6, fontface = "bold", colour = marker_col) else NULL} +
    coord_cartesian(clip = "off") +
    labs(x = if (show_x) "Chromosome" else NULL,
         y = expression(bold(-log[10](italic(p))))) +
    plot_theme +
    theme(axis.text.x = if (show_x) element_text(size = 6.5) else element_blank(),
          plot.margin = margin(if (is.na(gene_lab)) 3 else 14, 5, 3, 5))
}

cols_built <- lapply(columns, function(cl) {
  mg <- gpos_of(cl$chr, cl$pos)
  ps <- lapply(seq_along(cl$traits), function(i) {
    tr <- cl$traits[i]
    f  <- file.path(gwas_root, cl$dir, paste0(tr, ".txt"))
    if (!file.exists(f)) stop("missing: ", f)
    message("  ", cl$dir, " / ", tr)
    panel(read_gwas(f), cl$tag, gsub("_", "/", tr, fixed = TRUE), mg,
          show_x = (i == length(cl$traits)),
          gene_lab = if (i == 1) cl$gene else NA_character_)
  })
  wrap_plots(ps, ncol = 1)
})

fig <- wrap_plots(cols_built, nrow = 1) + plot_annotation(tag_levels = "A")
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
ggsave(out_file, fig, width = 17, height = 11, dpi = 300, bg = "white", limitsize = FALSE)
message("Saved: ", out_file)
