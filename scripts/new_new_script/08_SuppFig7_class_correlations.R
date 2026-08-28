# ==============================================================================
# Supp Figure S7 -- Correlation structure among lipid classes, within each trial.
#
# Two panels, CTL and LIN. The former panel C, a delta-r map of
# r_LIN - r_CTL with sign-reversal boxes, has been removed: it is a CTL-LIN
# difference in compositional space, which the manuscript no longer reports.
# Each matrix here is computed inside one trial and the two are never
# differenced.
#
# The transform is the original one. Class sums are closed to a constant sum
# per sample, zeros are replaced by half the smallest positive value in that
# sample, and the centered log-ratio is taken. Unlike the superseded version,
# classes are assigned with lipid_class() from _common.R, so a species written
# PC(18:2/0:0) counts as LPC here and everywhere else in the paper rather than
# as PC in this one figure.
#
# plot_theme is used unmodified except for the legend, which sits outside on the
# right: a colour bar inside a coord_fixed correlation matrix covers cells.
#
# Inputs
#   data/SPATS_fitted/non_normalized_intensities/Final_subset_{control,lowinput}_*.csv
#
# Output
#   fig/supp/SuppFig_S6_CLR_Correlations.png
#     (filename kept so the manuscript's \includegraphics path does not move;
#      the printed number is S7)
# ==============================================================================
source("scripts/new_new_script/_common.R")
suppressPackageStartupMessages({ library(tidyr); library(purrr); library(tibble) })

CLASSES <- c("DG", "DGDG", "LPC", "LPE", "MG", "MGDG",
             "PA", "PC", "PE", "PG", "PS", "SQDG", "TG")

# ---- class-level CLR matrix --------------------------------------------------
class_clr <- function(path) {
  d <- read_trial(path)
  feats <- names(d)[-1]
  feats <- feats[grepl("\\(", feats)]
  m <- as.matrix(d[, feats, drop = FALSE]); storage.mode(m) <- "numeric"
  m[!is.finite(m)] <- 0
  labs <- lipid_class(feats)
  cm <- vapply(CLASSES, function(cl) {
    j <- which(labs == cl)
    if (length(j)) rowSums(m[, j, drop = FALSE]) else rep(0, nrow(m))
  }, numeric(nrow(m)))
  # Close each sample to a constant sum, replacing zeros with half the smallest
  # positive class value in that sample so the log is defined.
  for (i in seq_len(nrow(cm))) {
    z <- cm[i, ]; pos <- z[z > 0]
    if (length(pos)) { z[z == 0] <- min(pos) * .5; cm[i, ] <- z / sum(z) }
  }
  log(cm) - rowMeans(log(cm))
}

cor_stats <- function(x) {
  r <- cor(x, method = "pearson", use = "pairwise.complete.obs")
  p <- matrix(NA_real_, ncol(x), ncol(x), dimnames = dimnames(r))
  diag(p) <- 0
  for (i in seq_len(ncol(x) - 1L)) for (j in (i + 1L):ncol(x)) {
    p[i, j] <- p[j, i] <- cor.test(x[, i], x[, j], method = "pearson")$p.value
  }
  list(r = r, p = p)
}

ctl <- cor_stats(class_clr(CTL_CSV))
lin <- cor_stats(class_clr(LIN_CSV))

# ---- one panel ---------------------------------------------------------------
# Lower triangle carries the coefficient, upper triangle the significance stars,
# so one square holds both without a second figure.
plot_cor <- function(cs) {
  d <- expand_grid(X = CLASSES, Y = CLASSES) %>%
    mutate(ix = match(X, CLASSES), iy = match(Y, CLASSES),
           R  = map2_dbl(X, Y, ~ cs$r[.x, .y]),
           P  = map2_dbl(X, Y, ~ cs$p[.x, .y]),
           diag = ix == iy, upper = ix < iy,
           sig = case_when(P < .001 ~ "***", P < .01 ~ "**", P < .05 ~ "*", TRUE ~ ""),
           label = case_when(diag ~ sprintf("%.2f", R),
                             upper ~ sig,
                             TRUE ~ sprintf("%.2f", R)),
           X = factor(X, levels = CLASSES),
           Y = factor(Y, levels = rev(CLASSES)))
  ggplot(d, aes(X, Y, fill = R)) +
    geom_tile(colour = "white", linewidth = .45) +
    geom_tile(data = filter(d, diag), fill = "black", colour = "white", linewidth = .45) +
    geom_text(aes(label = label, colour = ifelse(diag, "white", "black")), size = 4.1) +
    scale_colour_identity() +
    scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                         midpoint = 0, limits = c(-1, 1), name = "CLR r") +
    coord_fixed() + labs(x = NULL, y = NULL) +
    plot_theme +
    theme(panel.grid = element_blank(),   # tiles fill the panel; a grid only shows
                                          # in the whitespace coord_fixed leaves
          axis.text.x  = element_text(angle = 45, hjust = 1, colour = "black", size = 14),
          axis.text.y  = element_text(colour = "black", size = 14),
          legend.title = element_text(size = 14, face = "bold"),
          legend.background = element_blank(),
          legend.position = "right")
}

fig <- (plot_cor(ctl) | plot_cor(lin)) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A", theme = TAG_THEME)

save_fig(fig, "SuppFig_S6_CLR_Correlations.png", width = 17, height = 8.5, subdir = "supp")

# ---- console summary ---------------------------------------------------------
up <- upper.tri(ctl$r)
cat("\n-- class pairs, within-trial CLR correlation --\n")
cat(sprintf("  CTL  n pairs = %d  significant (p<0.05) = %d  |r| range %.2f to %.2f\n",
            sum(up), sum(ctl$p[up] < .05), min(abs(ctl$r[up])), max(abs(ctl$r[up]))))
cat(sprintf("  LIN  n pairs = %d  significant (p<0.05) = %d  |r| range %.2f to %.2f\n",
            sum(up), sum(lin$p[up] < .05), min(abs(lin$r[up])), max(abs(lin$r[up]))))
