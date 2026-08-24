# ==============================================================================
# Figure 6, panel C -- GWAS support for the LINEX-supported DG-TG branches
#
# Replaces the previous panel C, which was labelled with SORBI_3001G103800 --
# a gene that is NOT a candidate (zero rows in the LD-mapped candidate table).
# The real LCAT-like candidate is SORBI_3006G214500 (LIN individual,
# p = 3.86e-10, r2 = 1.00), which associates with DG(18:2_18:3) and MG(18:2).
#
# All three genes are LIN-exclusive candidates, so every panel is LIN.
#
# Columns  DGAT1      SORBI_3010G170000  chr10  3 traits
#          LCAT-like  SORBI_3006G214500  chr6   2 traits
#          SDP1-like  SORBI_3001G041900  chr1   2 traits
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

gwas_root <- "/Users/nirwantandukar/Documents/Research/data/SAP/GWAS_result"
out_file  <- "/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database/fig/main/Figure6C_LINEX_GWAS.png"

columns <- list(
  list(gene = "SORBI_3010G170000", lab = "DGAT1", chr = 10L, pos = 50099836, dir = "LIN_ind",
       traits = c("TG(18:1_18:3_22:0)", "TG(18:2_18:2_18:4)")),
  list(gene = "SORBI_3006G214500", lab = "LCAT-like", chr = 6L, pos = 56286921, dir = "LINEX_gwas",
       traits = c("DG(18:2_18:3)", "MG(18:2)")),
  list(gene = "SORBI_3001G041900", lab = "SDP1-like", chr = 1L, pos = 3125810, dir = "LINEX_gwas",
       traits = c("TG(16:1_20:1_20:2)", "TG(18:1_18:2_22:0)"))
)

bonferroni <- 8.19e-9
suggestive <- 1e-5

keep_all_below <- 0.01   # rendering only; set thin_every <- 1 to plot every SNP
thin_every     <- 25

band_dark  <- "#3B3B3B"
band_light <- "#9E9E9E"
marker_col <- "#D55E00"   # Okabe-Ito vermillion, colourblind-safe
bonf_col   <- "#B2182B"
sugg_col   <- "#2166AC"

plot_theme <- theme_minimal(base_size = 13) +
  theme(axis.text = element_text(colour = "black", size = 8),
        axis.title = element_text(face = "bold", size = 10),
        axis.line = element_line(colour = "black", linewidth = .5),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        plot.margin = margin(4, 6, 4, 6))

read_gwas <- function(path) {
  stopifnot(file.exists(path))
  d <- fread(path, select = c("chr", "ps", "p_wald"), showProgress = FALSE)
  setnames(d, c("chr", "ps", "p"))
  d <- d[is.finite(p) & p > 0]
  d[, chr := suppressWarnings(as.integer(chr))]
  d[!is.na(chr) & chr >= 1 & chr <= 10]
}

message("scanning chromosome lengths ...")
f0 <- file.path(gwas_root, columns[[1]]$dir, paste0(columns[[1]]$traits[1], ".txt"))
chrlen <- read_gwas(f0)[, .(len = max(ps)), by = chr][order(chr)]
chrlen[, offset := cumsum(as.numeric(len)) - as.numeric(len)]
axis_df <- chrlen[, .(chr, centre = offset + len / 2)]
genome_end <- chrlen[, max(offset + len)]
gpos_of <- function(ch, p) chrlen[chr == ch, offset] + p
label_hjust <- function(g) { f <- g / genome_end
  if (f > 0.85) 1.02 else if (f < 0.15) -0.02 else 0.5 }

panel <- function(d, trait, marker_g, show_x, gene_lab) {
  d <- merge(d, chrlen[, .(chr, offset)], by = "chr")
  d[, gpos := as.numeric(ps) + offset][, logp := -log10(p)]
  d <- d[p < keep_all_below | (seq_len(nrow(d)) %% thin_every == 0L)]
  d[, band := factor(chr %% 2)]

  ggplot(d, aes(gpos, logp, colour = band)) +
    geom_vline(xintercept = marker_g, colour = marker_col,
               linetype = "dashed", linewidth = .8) +
    geom_point(size = .9, alpha = .85) +
    geom_hline(yintercept = -log10(bonferroni), linetype = "dashed",
               colour = bonf_col, linewidth = .6) +
    geom_hline(yintercept = -log10(suggestive), linetype = "dashed",
               colour = sugg_col, linewidth = .5) +
    scale_colour_manual(values = c("0" = band_dark, "1" = band_light)) +
    scale_x_continuous(breaks = axis_df$centre, labels = axis_df$chr,
                       expand = expansion(mult = .01)) +
    scale_y_continuous(expand = expansion(mult = c(0, .12))) +
    annotate("text", x = -Inf, y = Inf, label = trait, hjust = -0.06, vjust = 1.5,
             size = 2.9, fontface = "bold", colour = "grey15") +
    {if (!is.na(gene_lab))
       annotate("text", x = marker_g, y = Inf, label = gene_lab, vjust = -0.6,
                hjust = label_hjust(marker_g), size = 3.0,
                fontface = "bold", colour = marker_col) else NULL} +
    coord_cartesian(clip = "off") +
    labs(x = if (show_x) "Chromosome" else NULL,
         y = expression(bold(-log[10](italic(p))))) +
    plot_theme +
    theme(axis.text.x = if (show_x) element_text(size = 7.5) else element_blank(),
          plot.margin = margin(if (is.na(gene_lab)) 4 else 18, 6, 4, 6))
}

nmax <- max(vapply(columns, function(cl) length(cl$traits), integer(1)))

cols_built <- lapply(columns, function(cl) {
  mg <- gpos_of(cl$chr, cl$pos)
  ps <- lapply(seq_along(cl$traits), function(i) {
    tr <- cl$traits[i]
    f  <- file.path(gwas_root, cl$dir, paste0(tr, ".txt"))
    message("  ", cl$lab, " / ", tr)
    panel(read_gwas(f), gsub("_", "/", tr, fixed = TRUE), mg,
          show_x = (i == length(cl$traits)),
          gene_lab = if (i == 1) sprintf("%s  %s", cl$lab, cl$gene) else NA_character_)
  })
  while (length(ps) < nmax) ps <- c(ps, list(patchwork::plot_spacer()))
  wrap_plots(ps, ncol = 1)
})

fig <- wrap_plots(cols_built, nrow = 1)
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
ggsave(out_file, fig, width = 15, height = 6.6, dpi = 300, bg = "white", limitsize = FALSE)
message("Saved: ", out_file)
