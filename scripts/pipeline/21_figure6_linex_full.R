# ==============================================================================
# Figure 6 -- LINEX-supported lipid reaction remodeling between LIN and CTL
#
#   A  reaction subnetwork: substrates -> candidate reactions -> products
#   B  reaction-balance scores, LIN - CTL, paired genotype-matched Wilcoxon
#   C  GWAS support for the three branch-supporting candidate genes (LIN only)
#
# Panel C replaces the previous version, which was labelled SORBI_3001G103800 --
# a gene with zero rows in the LD-mapped candidate table. The real LCAT-like
# candidate is SORBI_3006G214500 (LIN individual, p = 3.86e-10, r2 = 1.00).
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})

sold_data <- Sys.getenv("SOLD_DATA", "/Users/nirwantandukar/Documents/Github/SoLD_paper/data")
gwas_root <- Sys.getenv("GWAS_ROOT", "/Users/nirwantandukar/Documents/Research/data/SAP/GWAS_result")
out_file  <- Sys.getenv("FIG6_OUT",
  "/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database/fig/main/Figure6_LINEX_Full.png")
spats <- file.path(sold_data, "SPATS_fitted/non_normalized_intensities")

marker_col <- "#D55E00"; bonf_col <- "#B2182B"
band_dark  <- "#3B3B3B"; band_light <- "#9E9E9E"
prod_col   <- "#1B9E77"; subs_col   <- "#D95F02"
# panel A: border colour = lipid class (palette from 32_Linex_network_clean.R)
lipid_class_colors <- c(PC="#00441B", LPC="#41AB5D", PE="#1B7837", LPE="#78C679",
                        DG="#54278F", MG="#8941ED", TG="#ED804A")
lip_class <- function(x) sub("\\(.*$", "", x)

plot_theme <- theme_minimal(base_size = 13) +
  theme(axis.text = element_text(colour = "black", size = 8),
        axis.title = element_text(face = "bold", size = 10),
        axis.line = element_line(colour = "black", linewidth = .5),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        legend.title = element_blank(), plot.tag = element_text(face = "bold", size = 16))

# ---------------- A: reaction subnetwork -------------------------------------
# Node/edge table with explicit direction, so PNPLA1 (TG -> DG) is drawn the way
# 32_Linex_network_clean.R defines it and the way panel B scores it.
rxn <- data.frame(
  lab = c("LCAT*", "PNPLA1", "LRO1", "PNPLA3"),
  eq  = c("PC+DG <-> LPC+TG", "TG <-> DG", "PE+DG <-> TG+LPE", "DG+MG <-> TG"),
  x   = 2, y = c(4, 3, 2, 1))

lipL <- data.frame(lab = c("PC(16:0_18:1)","PE(16:0_18:1)","DG(16:0_18:1)","DG(16:0_18:2)",
                           "DG(18:1_18:2)","MG(16:0)","MG(18:1)","MG(18:2)"),
                   x = 1, y = seq(4.7, 0.3, length.out = 8))
lipR <- data.frame(lab = c("TG(16:0_18:1_18:2)","LPC(16:0)","LPE(16:0)"),
                   x = 3, y = c(4.1, 2.6, 1.2))
lipids <- rbind(lipL, lipR)
nodes  <- rbind(lipids[, c("lab","x","y")], rxn[, c("lab","x","y")])

DGs <- lipL$lab[grepl("^DG", lipL$lab)]
MGs <- lipL$lab[grepl("^MG", lipL$lab)]
PCs <- lipL$lab[grepl("^PC", lipL$lab)]
PEs <- lipL$lab[grepl("^PE", lipL$lab)]
TGp <- lipR$lab[grepl("^TG", lipR$lab)]

# R1 LCAT*  : PC -> R1, DG -> R1, R1 -> LPC, R1 -> TG
# R2 PNPLA1 : TG -> R2, R2 -> DG
# R3 LRO1   : PE -> R3, DG -> R3, R3 -> TG, R3 -> LPE
# R4 PNPLA3 : DG -> R4, MG -> R4, R4 -> TG
edges <- rbind(
  data.frame(from = c(PCs, DGs),        to = "LCAT*"),
  data.frame(from = "LCAT*",            to = c("LPC(16:0)", TGp)),
  data.frame(from = TGp,                to = "PNPLA1"),
  data.frame(from = "PNPLA1",           to = DGs),
  data.frame(from = c(PEs, DGs),        to = "LRO1"),
  data.frame(from = "LRO1",             to = c(TGp, "LPE(16:0)")),
  data.frame(from = c(DGs, MGs),        to = "PNPLA3"),
  data.frame(from = "PNPLA3",           to = TGp))
edges <- edges |>
  left_join(nodes, by = c("from" = "lab")) |> rename(x1 = x, y1 = y) |>
  left_join(nodes, by = c("to"   = "lab")) |> rename(x2 = x, y2 = y)
# trim both ends so arrowheads sit clear of the boxes
trim <- 0.16
edges <- within(edges, {
  dx <- x2 - x1; dy <- y2 - y1; len <- sqrt(dx^2 + dy^2)
  xa <- x1 + dx/len * trim; ya <- y1 + dy/len * trim
  xb <- x2 - dx/len * trim; yb <- y2 - dy/len * trim
})
edges$rev <- edges$x2 < edges$x1   # back-edges (PNPLA1)

pA <- ggplot() +
  geom_segment(data = edges[!edges$rev, ], aes(x = xa, xend = xb, y = ya, yend = yb),
               colour = "grey68", linewidth = .35,
               arrow = arrow(length = unit(4.5, "pt"), type = "closed")) +
  geom_segment(data = edges[edges$rev, ], aes(x = xa, xend = xb, y = ya, yend = yb),
               colour = "grey35", linewidth = .45, linetype = "22",
               arrow = arrow(length = unit(4.5, "pt"), type = "closed")) +
  geom_label(data = lipids, aes(x, y, label = lab, colour = lip_class(lab)), size = 2.5,
             fill = "white", label.size = .6, label.r = unit(2, "pt")) +
  geom_label(data = rxn, aes(x, y, label = paste0(lab, "\n", eq)), size = 2.6,
             colour = "grey20", fontface = "bold", fill = "grey96",
             label.size = .5, label.r = unit(2, "pt")) +
  scale_colour_manual(values = lipid_class_colors, guide = "none") +
  scale_x_continuous(limits = c(0.55, 3.45)) +
  scale_y_continuous(limits = c(0.05, 5.05)) +
  theme_void(base_size = 13) + theme(plot.tag = element_text(face = "bold", size = 16))

# ---------------- B: reaction-balance scores ---------------------------------
read_cond <- function(f) as.data.table(read.csv(file.path(spats, f), check.names = FALSE))
C <- read_cond("Final_subset_control_all_lipids_fitted_phenotype_non_normalized.csv")
L <- read_cond("Final_subset_lowinput_all_lipids_fitted_phenotype_non_normalized.csv")

class_sum <- function(D, k) {
  cols <- grep(paste0("^", k, "\\("), names(D), value = TRUE)
  if (!length(cols)) return(rep(NA_real_, nrow(D)))
  rowSums(sapply(D[, ..cols], function(v) suppressWarnings(as.numeric(v))), na.rm = TRUE)
}
scores <- function(D) {
  eps <- 0.5 * min(unlist(lapply(c("PC","DG","TG","LPC","LPE","PE","MG"),
        function(k) { s <- class_sum(D, k); s[s > 0] })), na.rm = TRUE)
  l <- function(k) log10(class_sum(D, k) + eps)
  data.table(Line = D[[1]],
    `LCAT*`  = rowMeans(cbind(l("LPC"), l("TG"))) - rowMeans(cbind(l("PC"), l("DG"))),
    PNPLA1   = l("DG") - l("TG"),
    LRO1     = rowMeans(cbind(l("TG"), l("LPE"))) - rowMeans(cbind(l("PE"), l("DG"))),
    PNPLA3   = l("TG") - rowMeans(cbind(l("DG"), l("MG"))))
}
sc <- merge(melt(scores(C), id.vars = "Line", variable.name = "rxn", value.name = "CTL"),
            melt(scores(L), id.vars = "Line", variable.name = "rxn", value.name = "LIN"),
            by = c("Line", "rxn"))
sc <- sc[is.finite(CTL) & is.finite(LIN)]
# Report the MEDIAN paired difference with empirical 2.5/97.5 percentile bounds,
# matching SuppTable_S12_ReactionBalance and the values quoted in the Results.
# (Plotting mean(d) here previously put slightly different numbers on the bars
#  than the text beside them.)
res <- sc[, {
  w <- wilcox.test(LIN, CTL, paired = TRUE)
  d <- LIN - CTL
  q <- quantile(d, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
  .(delta = median(d), lo = q[1], hi = q[2],
    med_CTL = median(CTL), med_LIN = median(LIN),
    p = w$p.value, n = .N)
}, by = rxn]
res[, q := p.adjust(p, "BH")]
res[, stars := fifelse(q < .001, "***", fifelse(q < .01, "**", fifelse(q < .05, "*", "")))]
res[, dir := fifelse(delta > 0, "products enriched in LIN", "substrates enriched in LIN")]
res[, rxn := factor(rxn, levels = rev(c("LCAT*","PNPLA1","LRO1","PNPLA3")))]

pB <- ggplot(res, aes(rxn, delta, fill = dir)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45", linewidth = .45) +
  geom_col(width = .62, colour = "black", linewidth = .3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .18, linewidth = .5) +
  geom_text(aes(label = stars, y = fifelse(delta > 0, hi + .03, lo - .03)),
            size = 4, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = c("products enriched in LIN" = prod_col,
                               "substrates enriched in LIN" = subs_col)) +
  labs(x = NULL, y = "LIN - CTL reaction score (median, log10)") +
  plot_theme + theme(legend.position = "bottom", panel.grid.major.y = element_blank())

# ---------------- C: GWAS support --------------------------------------------
columns <- list(
  list(gene="SORBI_3010G170000", lab="DGAT1",     chr=10L, pos=50099836, dir="LIN_ind",
       traits=c("TG(18:1_18:3_22:0)","TG(18:2_18:2_18:4)")),
  list(gene="SORBI_3006G214500", lab="LCAT-like", chr=6L,  pos=56286921, dir="LINEX_gwas",
       traits=c("DG(18:2_18:3)","MG(18:2)")),
  list(gene="SORBI_3001G041900", lab="SDP1-like", chr=1L,  pos=3125810,  dir="LINEX_gwas",
       traits=c("TG(16:1_20:1_20:2)","TG(18:1_18:2_22:0)")))

bonferroni <- 8.19e-9
keep_all_below <- 0.01; thin_every <- 25

read_gwas <- function(p) {
  d <- fread(p, select=c("chr","ps","p_wald"), showProgress=FALSE)
  setnames(d, c("chr","ps","p")); d <- d[is.finite(p) & p>0]
  d[, chr := suppressWarnings(as.integer(chr))]; d[!is.na(chr) & chr>=1 & chr<=10]
}
message("scanning chromosome lengths ...")
chrlen <- read_gwas(file.path(gwas_root, columns[[1]]$dir,
          paste0(columns[[1]]$traits[1],".txt")))[, .(len=max(ps)), by=chr][order(chr)]
chrlen[, offset := cumsum(as.numeric(len)) - as.numeric(len)]
axis_df <- chrlen[, .(chr, centre = offset + len/2)]
genome_end <- chrlen[, max(offset+len)]
gpos_of <- function(ch,p) chrlen[chr==ch, offset] + p
lab_h <- function(g){f<-g/genome_end; if(f>.85) 1.02 else if(f<.15) -0.02 else .5}

mpanel <- function(d, trait, mg, show_x, gene_lab) {
  d <- merge(d, chrlen[, .(chr, offset)], by="chr")
  d[, gpos := as.numeric(ps)+offset][, logp := -log10(p)]
  d <- d[p < keep_all_below | (seq_len(nrow(d)) %% thin_every == 0L)]
  d[, band := factor(chr %% 2)]
  ggplot(d, aes(gpos, logp, colour=band)) +
    geom_vline(xintercept=mg, colour=marker_col, linetype="dashed", linewidth=.8) +
    geom_point(size=.9, alpha=.85) +
    geom_hline(yintercept=-log10(bonferroni), linetype="dashed", colour=bonf_col, linewidth=.6) +
    scale_colour_manual(values=c("0"=band_dark,"1"=band_light)) +
    scale_x_continuous(breaks=axis_df$centre, labels=axis_df$chr, expand=expansion(mult=.01)) +
    scale_y_continuous(expand=expansion(mult=c(0,.12))) +
    annotate("text", x=-Inf, y=Inf, label=trait, hjust=-0.06, vjust=1.5,
             size=2.9, fontface="bold", colour="grey15") +
    {if(!is.na(gene_lab)) annotate("text", x=mg, y=Inf, label=gene_lab, vjust=-0.6,
        hjust=lab_h(mg), size=3, fontface="bold", colour=marker_col) else NULL} +
    coord_cartesian(clip="off") +
    labs(x = if(show_x) "Chromosome" else NULL, y = expression(bold(-log[10](italic(p))))) +
    plot_theme + theme(legend.position="none",
      axis.text.x = if(show_x) element_text(size=7.5) else element_blank(),
      plot.margin = margin(if(is.na(gene_lab)) 4 else 18, 6, 4, 6))
}

cols_built <- lapply(columns, function(cl) {
  mg <- gpos_of(cl$chr, cl$pos)
  wrap_plots(lapply(seq_along(cl$traits), function(i) {
    tr <- cl$traits[i]
    message("  ", cl$lab, " / ", tr)
    mpanel(read_gwas(file.path(gwas_root, cl$dir, paste0(tr,".txt"))),
           gsub("_","/",tr,fixed=TRUE), mg, show_x=(i==length(cl$traits)),
           gene_lab = if(i==1) sprintf("%s  %s", cl$lab, cl$gene) else NA_character_)
  }), ncol = 1)
})
pC <- wrap_plots(cols_built, nrow = 1)

fig <- (pA | pB) / pC + plot_layout(heights = c(1, 1.35)) + plot_annotation(tag_levels = "A")
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
ggsave(out_file, fig, width = 15, height = 12, dpi = 300, bg = "white", limitsize = FALSE)
# NOTE: an S12 writer was drafted here and deliberately disabled. This script
# assigns lipids to classes by name pattern (^CLASS\(), whereas the archived
# SuppTable_S12 was evidently built using data/lipid_class/final_lipid_classes.csv.
# The two disagree for any species named in IUPAC form rather than CLASS(chains)
# notation, which shifts the LCAT* and LRO1 branches (both use LPC/LPE) while
# leaving PNPLA1 and PNPLA3 unchanged. Resolve the class-assignment question
# before regenerating S12 from this script.
message("Saved: ", out_file)
print(res[, .(rxn, delta = round(delta,4), q = signif(q,3), stars, n)])
