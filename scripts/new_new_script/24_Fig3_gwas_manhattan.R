# ==============================================================================
# Figure 3 -- The same candidate locus recurs across multiple lipid phenotypes
#             within each trial.
#
# Four columns, one per trial x trait layer. Each column stacks four phenotypes
# that share the same candidate locus, so the recurring peak reads straight down
# the column.
#
# Panel letters run A to D, one per column, not one per plot: sixteen letters
# would be a lookup table rather than a label. Each column therefore carries its
# tag on its top plot only, and the four phenotypes inside a column are named by
# the trait label in the corner of each plot.
#
# The shared locus is marked by a short arrow sitting just above the lead SNP in
# every plot, and the candidate gene is named once, above the arrow in the top
# plot of each column. The dashed vertical line that used to run the full height
# of every plot has been removed: it drew the eye down the whole y-axis and
# implied the locus mattered at every -log10(p), which it does not.
#
#   A  CTL individual lipids   -- SORBI_3006G040500, phospholipid-transporting ATPase (chr6)
#   B  CTL class sums/ratios   -- SORBI_3010G151000, CMP-sialic acid transporter 1 (chr10)
#   C  LIN individual lipids   -- SORBI_3010G170000, DGAT1 / O-acyltransferase (chr10)
#   D  LIN class sums/ratios   -- SORBI_3006G195200, purine permease array (chr6)
#
# plot_theme supplies the house styling. Sixteen plots on one page need smaller
# type than a standalone figure, so axis text and titles are scaled down here
# and the vertical grid is dropped, since a grid line at each chromosome centre
# competes with the chromosome banding.
#
# Input : GEMMA .assoc output (tab-delimited: chr rs ps ... p_wald)
# Output: fig/main/Figure3_GWAS_Manhattan.png
# ==============================================================================
suppressPackageStartupMessages({ library(data.table) })
source("scripts/new_new_script/_common.R")

gwas_root <- Sys.getenv("GWAS_ROOT",
  "/Users/nirwantandukar/Documents/Research/data/SAP/GWAS_result")
out_file  <- Sys.getenv("FIG3_OUT", file.path(FIG_MAIN, "Figure3_GWAS_Manhattan.png"))

columns <- list(
  list(tag = "CTL", layer = "Individual lipids", dir = "CTL_ind",
       gene = "SORBI_3006G040500", gene_lab = "Phospholipid ATPase",
       chr = 6L, pos = 25522678,
       traits = c("TG(16:0_16:0_18:1)", "TG(16:0_18:0_18:2)",
                  "TG(18:0_18:2_18:2)", "TG(18:1_18:2_18:2)")),
  list(tag = "CTL", layer = "Class sums / ratios", dir = "CTL_sumratio",
       gene = "SORBI_3010G151000", gene_lab = "CMP-sialic acid transporter 1",
       chr = 10L, pos = 43625333,
       traits = c("Sum_AEG_over_DGDG", "Sum_PC_over_DG",
                  "Sum_PC_over_MGDG", "Sum_PC_over_SQDG")),
  list(tag = "LIN", layer = "Individual lipids", dir = "LIN_ind",
       gene = "SORBI_3010G170000", gene_lab = "DGAT1",
       chr = 10L, pos = 50099836,
       traits = c("TG(18:1_18:3_22:0)", "TG(18:2_18:2_18:4)",
                  "TG(18:2_20:3_22:0)", "TG(18:3_18:3_18:3)")),
  list(tag = "LIN", layer = "Class sums / ratios", dir = "LIN_sumratio",
       gene = "SORBI_3006G195200", gene_lab = "Purine permease array",
       chr = 6L, pos = 54823646,
       traits = c("Sum_FA", "Sum_FA_over_DG", "Sum_FA_over_SM", "Sum_FA_over_SQDG"))
)

bonferroni <- 8.19e-9

# Rendering only: all SNPs below keep_all_below are drawn, the rest subsampled
# 1-in-thin_every. Set thin_every <- 1 to plot every SNP.
keep_all_below <- 0.01
thin_every     <- 25

# The arrow is anchored on the strongest SNP within this distance of the
# candidate gene, which is the LD window the candidate mapping itself used.
peak_window <- 250e3

band_dark  <- "#3B3B3B"
band_light <- "#9E9E9E"

theme_gwas <- plot_theme +
  theme(axis.text          = element_text(size = 8, colour = "black"),
        axis.title         = element_text(size = 10, face = "bold"),
        panel.grid.major.x = element_blank(),
        legend.position    = "none",
        plot.tag           = element_text(face = "bold", size = 18),
        plot.tag.position  = "topleft",
        plot.margin        = margin(4, 6, 4, 6))

read_gwas <- function(path) {
  d <- fread(path, select = c("chr", "ps", "p_wald"), showProgress = FALSE)
  setnames(d, c("chr", "ps", "p"))
  d <- d[is.finite(p) & p > 0]
  d[, chr := suppressWarnings(as.integer(chr))]
  d[!is.na(chr) & chr >= 1 & chr <= 10]
}

message("scanning chromosome lengths ...")
first  <- read_gwas(file.path(gwas_root, columns[[1]]$dir,
                              paste0(columns[[1]]$traits[1], ".txt")))
chrlen <- first[, .(len = max(ps)), by = chr][order(chr)]
chrlen[, offset := cumsum(as.numeric(len)) - as.numeric(len)]
axis_df <- chrlen[, .(chr, centre = offset + len / 2)]
rm(first); invisible(gc())

gpos_of    <- function(ch, p) chrlen[chr == ch, offset] + p
genome_end <- chrlen[, max(offset + len)]

# Keep the gene label inside the plot: right-align near the end of the genome,
# left-align near the start, centre otherwise.
label_hjust <- function(g) { f <- g / genome_end
  if (f > 0.85) 1.02 else if (f < 0.15) -0.02 else 0.5 }

panel <- function(d, trait, marker_g, show_x, gene_lab, col_tag = NULL) {
  d <- merge(d, chrlen[, .(chr, offset)], by = "chr")
  d[, gpos := as.numeric(ps) + offset][, logp := -log10(p)]
  keep <- d$p < keep_all_below | (seq_len(nrow(d)) %% thin_every == 0L)
  d <- d[keep]
  d[, band := factor(chr %% 2)]

  # Anchor the arrow on the lead SNP at the shared locus, so it points at the
  # peak this plot actually has rather than at a fixed coordinate.
  near <- d[abs(gpos - marker_g) <= peak_window]
  if (nrow(near)) {
    px <- near$gpos[which.max(near$logp)]
    py <- max(near$logp)
  } else {
    px <- marker_g; py <- 0
  }
  ytop <- max(d$logp)
  a_lo <- py + 0.05 * ytop          # arrow tip, just clear of the point
  a_hi <- py + 0.20 * ytop          # arrow tail
  head_room <- if (is.na(gene_lab)) 0.26 else 0.34

  ggplot(d, aes(gpos, logp, colour = band)) +
    geom_point(size = .9, alpha = .85) +
    geom_hline(yintercept = -log10(bonferroni), linetype = "dashed",
               colour = bonf_col, linewidth = .6) +
    annotate("segment", x = px, xend = px, y = a_hi, yend = a_lo,
             colour = marker_col, linewidth = .7,
             arrow = arrow(type = "closed", length = unit(0.07, "inches"))) +
    scale_colour_manual(values = c("0" = band_dark, "1" = band_light)) +
    scale_x_continuous(breaks = axis_df$centre, labels = axis_df$chr,
                       expand = expansion(mult = .01)) +
    scale_y_continuous(expand = expansion(mult = c(0, head_room))) +
    annotate("text", x = -Inf, y = Inf, label = trait, hjust = -0.06, vjust = 1.5,
             size = 2.9, fontface = "bold", colour = "grey15") +
    {if (!is.na(gene_lab))
       annotate("text", x = px, y = a_hi, label = gene_lab, vjust = -0.45,
                hjust = label_hjust(px), size = 3, fontface = "bold",
                colour = marker_col) else NULL} +
    coord_cartesian(clip = "off") +
    labs(x = if (show_x) "Chromosome" else NULL,
         y = expression(bold(-log[10](italic(p)))),
         tag = col_tag) +
    theme_gwas +
    theme(axis.text.x = if (show_x) element_text(size = 7.5) else element_blank(),
          plot.margin = margin(if (is.null(col_tag)) 4 else 14, 6, 4, 6))
}

# One tag per column, carried by that column's top plot.
cols_built <- lapply(seq_along(columns), function(ci) {
  cl <- columns[[ci]]
  mg <- gpos_of(cl$chr, cl$pos)
  ps <- lapply(seq_along(cl$traits), function(i) {
    tr <- cl$traits[i]
    f  <- file.path(gwas_root, cl$dir, paste0(tr, ".txt"))
    if (!file.exists(f)) stop("missing: ", f)
    message("  ", cl$dir, " / ", tr)
    panel(read_gwas(f), gsub("_", "/", tr, fixed = TRUE), mg,
          show_x   = (i == length(cl$traits)),
          gene_lab = if (i == 1) sprintf("%s  %s", cl$gene_lab, cl$gene) else NA_character_,
          col_tag  = if (i == 1) LETTERS[ci] else NULL)
  })
  wrap_plots(ps, ncol = 1)
})

fig <- wrap_plots(cols_built, nrow = 1)
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
ggsave(out_file, fig, width = 17, height = 11, dpi = 300, bg = "white", limitsize = FALSE)
message("Saved: ", out_file)
