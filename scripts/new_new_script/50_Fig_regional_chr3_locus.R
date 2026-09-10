# ==============================================================================
# Figure 7. The SoLD Shiny application and the chromosome 3 candidate locus
#
#   Rscript scripts/new_new_script/50_Fig_regional_chr3_locus.R
#
# Panel A  the Shiny homepage, fig/main/individual_figs/Fig7A.png
# Panel B  the Gene Hits module at r2 >= 0.9 under LIN,
#          fig/main/individual_figs/Fig7B.png
# Panel C  the chr3 locus that the top of that Gene Hits table points to,
#          drawn here as a regional association plot in three tiers that share
#          one x axis and one r2 scale
#            top     every GWAS SNP in chr3 7.40-7.90 Mb for the LIN fatty-acid
#                    class-sum trait, as -log10(p) from the GEMMA output,
#                    coloured by its r2 with the lead SNP. Dashed line, the LIN
#                    Bonferroni threshold 0.05 / 6,085,245. The lead SNP is the
#                    diamond; a second SNP at 3:7,739,242 ties it at p = 4.46e-20.
#            middle  the 49 annotated genes in the window, coloured by the
#                    highest r2 any variant inside them reaches with the lead
#                    SNP, which is exactly how candidate genes are assigned
#                    throughout the paper. Grey genes carry no variant at
#                    r2 >= 0.4 and are not candidates, however close they sit.
#            bottom  pairwise r2 between every variant pair in the window less
#                    than 150 kb apart, as the usual hanging triangle. Each
#                    pixel takes the r2 of the closest variant pair, about
#                    360 bp per pixel.
#
# LD throughout comes from the same SAP LD panel used for candidate assignment,
# so all three tiers are on one r2 scale and the single colour bar reads across
# the whole of panel C.
#
# plot_theme overrides for this figure, all local, _common.R is untouched
#   legend.position moved from the inside top-right to the bottom of the figure
#   and the legend box dropped, since one colour bar is shared by three tiers
#   axis text on the two stacked tiers that carry no x axis is switched off
#   panel grid dropped on the gene track and the LD triangle
#
# Inputs   GWAS_result/LIN_sumratio/Sum_FA.txt              (GEMMA, not in repo)
#          fig/main/individual_figs/Fig7A.png
#          fig/main/individual_figs/Fig7B.png
#          data/LD_mapped/genes_ranges/genes.range
#          data/LD_mapped/candidate_tables/gene_trait_LD.csv
#          data/LD_mapped/LD_chr3/thinned/chr3_7.4_7.9.thin3.ld
#          data/LD_mapped/LD_chr3/thinned/chr3_7.4_7.9.thin3.snplist
#          data/LD_mapped/LD_chr3/thinned/lead_SNP_7732730_r2.tsv
# Output   fig/main/Figure7_Shiny_App.png
#
# The three thinned LD files are derived from the full PLINK run, which is
# 8,347 x 8,347 and 714 MB and is not carried in the repo. To rebuild them
#
#   plink --bfile SAP_ldpanel --chr 3 --from-bp 7400000 --to-bp 7900000 \
#         --r2 square --write-snplist --out chr3_7.4_7.9
#   cd data/LD_mapped/LD_chr3 && mkdir -p thinned
#   awk 'FNR%3==1' chr3_7.4_7.9.ld | cut -f "$(seq -s, 1 3 8347)" \
#       > thinned/chr3_7.4_7.9.thin3.ld
#   awk 'FNR%3==1' chr3_7.4_7.9.snplist > thinned/chr3_7.4_7.9.thin3.snplist
#   LN=$(grep -n '^SNP_7732730$' chr3_7.4_7.9.snplist | cut -d: -f1)
#   sed -n "${LN}p" chr3_7.4_7.9.ld | tr '\t' '\n' > /tmp/leadrow.txt
#   paste chr3_7.4_7.9.snplist /tmp/leadrow.txt > thinned/lead_SNP_7732730_r2.tsv
#
# Every third variant is kept for the triangle only, which still leaves 2,783
# variants over 500 kb, about one every 180 bp, finer than the pixel grid. The
# scatter and the gene track use the full, unthinned LD.
# ==============================================================================
source("scripts/new_new_script/_common.R")
suppressPackageStartupMessages({
  library(ggrepel); library(data.table); library(png); library(grid)
})

GWAS  <- Sys.getenv("GWAS_DIR", "~/mnt/GWAS_result/LIN_sumratio")
LD    <- file.path(REPO, "data/LD_mapped/LD_chr3/thinned")
IND   <- file.path(FIG_MAIN, "individual_figs")
BONF  <- 0.05 / 6085245                 # LIN genome-wide threshold
X0    <- 7.40; X1 <- 7.90               # window, Mb
MAXD  <- 0.15                           # deepest pair distance drawn, Mb
NX    <- 1400; NY <- 340                # triangle raster
BLUES <- c("#F7FBFF","#DEEBF7","#C6DBEF","#9ECAE1","#6BAED6",
           "#4292C6","#2171B5","#08519C","#08306B")
FOCUS <- c("SORBI_3003G088800","SORBI_3003G088900","SORBI_3003G089000")
LAB   <- c(SORBI_3003G088800 = "SORBI_3003G088800\nplastocyanin-like,\nphosphate starvation",
           SORBI_3003G088900 = "SORBI_3003G088900\nplastocyanin-like,\nphosphate starvation",
           SORBI_3003G089000 = "SORBI_3003G089000\nubiquitin-conjugating E2,\nclosest homolog PHO2")

## ---- association and LD with the lead SNP ----------------------------------
gw <- read.delim(file.path(GWAS, "Sum_FA.txt")) |>
  filter(chr == 3, ps >= X0*1e6, ps <= X1*1e6) |>
  transmute(pos_mb = ps/1e6, logp = -log10(p_wald)) |>
  filter(is.finite(logp))

ldlead <- read.delim(file.path(LD, "lead_SNP_7732730_r2.tsv"), header = FALSE,
                     col.names = c("snp","r2lead")) |>
  transmute(pos_mb = as.numeric(sub("^SNP_", "", snp))/1e6, r2lead)

gw   <- left_join(gw, ldlead, by = "pos_mb")
lead <- gw |> slice_max(logp, n = 1, with_ties = FALSE)

## ---- genes, coloured by the r2 that made them candidates -------------------
r2g <- read.csv(file.path(REPO, "data/LD_mapped/candidate_tables/gene_trait_LD.csv")) |>
  filter(condition == "LIN") |>
  group_by(GeneID) |> summarise(r2 = max(max_r2), .groups = "drop")

genes <- read.delim(file.path(REPO, "data/LD_mapped/genes_ranges/genes.range"),
                    header = FALSE, col.names = c("chr","start","end","GeneID")) |>
  filter(chr == 3, start >= X0*1e6, start <= X1*1e6) |>
  left_join(r2g, by = "GeneID") |>
  mutate(s = start/1e6, e = end/1e6, mid = (s+e)/2, candidate = !is.na(r2))

## ---- pairwise LD triangle --------------------------------------------------
pos <- as.numeric(sub("^SNP_", "",
        readLines(file.path(LD, "chr3_7.4_7.9.thin3.snplist"))))/1e6
M   <- as.matrix(fread(file.path(LD, "chr3_7.4_7.9.thin3.ld")))
stopifnot(nrow(M) == length(pos), ncol(M) == length(pos))

xg <- seq(X0, X1, length.out = NX)
yg <- seq(0, MAXD/2, length.out = NY)
g  <- expand.grid(xi = seq_len(NX), yi = seq_len(NY))
pa <- xg[g$xi] - yg[g$yi]                       # left member of the pair
pb <- xg[g$xi] + yg[g$yi]                       # right member of the pair
ok <- pa >= X0 & pb <= X1
ia <- round(approx(pos, seq_along(pos), xout = pa, rule = 2)$y)
ib <- round(approx(pos, seq_along(pos), xout = pb, rule = 2)$y)
val <- rep(NA_real_, nrow(g)); val[ok] <- M[cbind(ia[ok], ib[ok])]
tri <- data.frame(x = xg[g$xi], y = -yg[g$yi], r2 = val)
rm(M, g, pa, pb, ia, ib, val); invisible(gc())

## ---- shared scales ---------------------------------------------------------
XSC <- scale_x_continuous(limits = c(X0, X1), expand = expansion(mult = c(.01,.04)))
FSC <- scale_fill_gradientn(
  colours = BLUES, limits = c(0,1), na.value = "transparent",
  name = expression(italic(r)^2),
  guide = guide_colourbar(title.position = "left", title.vjust = .85,
                          barwidth = unit(14,"cm"), barheight = unit(.55,"cm"),
                          ticks.colour = "grey30", frame.colour = "grey30"))
LEADLINE <- geom_vline(xintercept = lead$pos_mb, colour = "grey60",
                       linetype = "dotted", linewidth = .45)
BARE <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              axis.line.x = element_blank())
# plot_theme sets panel.grid.major explicitly, so a blanket panel.grid = blank
# is ignored; both children have to be switched off by name.
NOGRID <- theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())

## ---- panel C, three tiers --------------------------------------------------
cTop <- ggplot(gw, aes(pos_mb, logp)) + LEADLINE +
  geom_hline(yintercept = -log10(BONF), linetype = "dashed",
             colour = bonf_col, linewidth = .5) +
  geom_point(data = filter(gw, is.na(r2lead)), colour = "grey75",
             size = 1.2, alpha = .5) +
  geom_point(data = arrange(filter(gw, !is.na(r2lead)), r2lead),
             aes(fill = r2lead), shape = 21, size = 1.9, stroke = .12,
             colour = "grey65", alpha = .9) +
  geom_point(data = lead, shape = 23, fill = "#08306B", colour = "black",
             stroke = .5, size = 4.6) +
  geom_text(data = lead,
            aes(label = sprintf("lead SNP  3:%s",
                                format(round(pos_mb*1e6), big.mark = ","))),
            hjust = -0.12, vjust = .3, size = 5.2, colour = "#08306B",
            fontface = "bold") +
  annotate("text", x = X0, y = -log10(BONF), vjust = -0.55, hjust = 0,
           label = "genome-wide threshold", size = 4.6, colour = bonf_col) +
  FSC + XSC + labs(x = NULL, y = expression(bold(-log[10]~(italic(p))))) +
  plot_theme + BARE

cMid <- ggplot(genes) + LEADLINE +
  geom_rect(data = filter(genes, !candidate),
            aes(xmin = s, xmax = e, ymin = -.09, ymax = .09), fill = "grey80") +
  geom_rect(data = filter(genes, candidate),
            aes(xmin = s, xmax = e, ymin = -.09, ymax = .09, fill = r2),
            colour = "grey30", linewidth = .2) +
  geom_text_repel(data = filter(genes, GeneID %in% FOCUS),
                  aes(x = mid, y = -.09, label = LAB[GeneID]), size = 4.6,
                  lineheight = .95, nudge_y = -.42, direction = "x",
                  segment.size = .35, segment.colour = "grey45",
                  box.padding = .35, min.segment.length = 0, max.overlaps = Inf) +
  FSC + XSC + scale_y_continuous(limits = c(-.62,.18), breaks = NULL) +
  labs(x = NULL, y = NULL) + plot_theme + BARE + NOGRID +
  theme(axis.line = element_blank())

cBot <- ggplot(tri, aes(x, y, fill = r2)) + geom_raster() + FSC + XSC +
  scale_y_continuous(breaks = NULL, expand = expansion(mult = c(.02,.02))) +
  labs(x = "Chromosome 3 position (Mb)", y = NULL) + plot_theme + NOGRID +
  theme(axis.line.y = element_blank())

panelC <- cTop / cMid / cBot + plot_layout(heights = c(2.3, 0.85, 1.45))

## ---- panels A and B, the two app screenshots -------------------------------
# Wrapped as ggplots rather than raw grobs so they take patchwork tags and keep
# their own aspect ratio whatever the panel they land in.
img_panel <- function(path) {
  im  <- readPNG(path)
  asp <- dim(im)[1] / dim(im)[2]
  ggplot() +
    annotation_custom(rasterGrob(im, interpolate = TRUE),
                      xmin = 0, xmax = 1, ymin = 0, ymax = asp) +
    scale_x_continuous(limits = c(0,1),   expand = c(0,0)) +
    scale_y_continuous(limits = c(0,asp), expand = c(0,0)) +
    coord_fixed() + theme_void() + theme(plot.margin = margin(2,2,2,2))
}

## ---- assemble --------------------------------------------------------------
fig <- (img_panel(file.path(IND,"Fig7A.png")) + labs(tag = "A") |
        img_panel(file.path(IND,"Fig7B.png")) + labs(tag = "B")) /
       (panelC & labs(tag = NULL)) +
  plot_layout(heights = c(1, 1.12), guides = "collect")

# the tag for panel C rides on its first tier, which is where it belongs
fig[[2]][[1]] <- fig[[2]][[1]] + labs(tag = "C")

fig <- fig & TAG_THEME &
  theme(legend.position   = "bottom",
        legend.direction  = "horizontal",
        legend.background = element_blank(),
        legend.title      = element_text(size = 18, face = "bold"),
        legend.text       = element_text(size = 15))

save_fig(fig, "Figure7_Shiny_App.png", width = 16, height = 16.8, subdir = "main")

cat("GWAS SNPs matched to LD panel:", sum(!is.na(gw$r2lead)), "of", nrow(gw),
    "| genes:", nrow(genes), "| candidates:", sum(genes$candidate),
    "| lead:", lead$pos_mb, "logp", round(lead$logp,2),
    "| tied leads:", sum(gw$logp == max(gw$logp)), "\n")
