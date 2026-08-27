# ==============================================================================
# Shared setup for every figure and table script in scripts/figures/.
#
#   source("scripts/figures/_common.R")
#
# Everything each script needs in common lives here -- paths, palettes, the
# house plot theme and the save helper -- so a change to any of them happens
# once rather than in twenty files.
#
# Paths can be overridden with environment variables, which keeps the scripts
# runnable on a machine where the two repositories sit somewhere else.
# ==============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(patchwork); library(vroom)
})

# ---- paths -------------------------------------------------------------------
REPO      <- Sys.getenv("SOLD_REPO", ".")
DATA_ROOT <- Sys.getenv("SOLD_DATA", file.path(REPO, "data"))
SPATS     <- file.path(DATA_ROOT, "SPATS_fitted/non_normalized_intensities")
FIG_MAIN  <- file.path(REPO, "fig/main")
FIG_SUPP  <- file.path(REPO, "fig/supp")
TAB_SUPP  <- file.path(REPO, "table/supp")
for (d in c(FIG_MAIN, FIG_SUPP, TAB_SUPP)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

CTL_CSV <- file.path(SPATS, "Final_subset_control_all_lipids_fitted_phenotype_non_normalized.csv")
LIN_CSV <- file.path(SPATS, "Final_subset_lowinput_all_lipids_fitted_phenotype_non_normalized.csv")

# ---- palettes ----------------------------------------------------------------
condition_colors <- c(CTL = "#440154FF", LIN = "#FDE725FF")

class_colors <- c(PC = "#00441B", PA = "#1B7837", PE = "#41AB5D", LPC = "#66C2A4",
                  LPE = "#2CA25F", PG = "#78C679", PS = "#C2E699", DG = "#54278F",
                  DGDG = "#F768A1", MG = "#8941ED", MGDG = "#FBB4D9",
                  SQDG = "#9D4D6C", TG = "#ED804A")

# Okabe-Ito, colourblind-safe, for markers and annotations
marker_col <- "#D55E00"
bonf_col   <- "#B2182B"

# ---- theme ---------------------------------------------------------------
# plot_theme, as Nirwan defines it. Do not edit this block without being asked.
# Anything a figure needs to change is overridden in that figure's own script
# and listed in its header, so this stays the single source of truth.
#
# One change from the original, requested 2026-08-26: panel.grid was
# element_blank(); it is now a very light grey major grid on both axes with the
# minor grid still off. Change GRID_COLOUR / GRID_WIDTH below to tune it.
GRID_COLOUR <- "grey92"   # very light grey; "grey88" is a touch stronger
GRID_WIDTH  <- 0.35

plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title     = element_text(
      size   = 14,
      face   = "bold",
      hjust  = 0.5,
      margin = margin(b = 10)
    ),
    axis.title.x   = element_text(
      size = 16,      # X-axis title size
      face = "bold"
    ),
    axis.title.y   = element_text(
      size = 16,      # Y-axis title size
      face = "bold"
    ),
    axis.text.x    = element_text(
      size = 16,      # X-axis tick label size
      color = "black"
    ),
    axis.text.y    = element_text(
      size = 16,      # Y-axis tick label size
      color = "black"
    ),
    axis.line      = element_line(color = "black"),
    panel.grid.major = element_line(color = GRID_COLOUR, linewidth = GRID_WIDTH),
    panel.grid.minor = element_blank(),

    legend.position      = c(0.95, 0.95),
    legend.justification = c("right","top"),
    legend.background    = element_rect(fill="white", color="grey70", size=0.4),
    legend.direction     = "vertical",
    legend.spacing.y     = unit(0.2,"cm"),
    legend.title         = element_blank(),
    legend.text          = element_text(size=16),

    plot.margin    = margin(15, 15, 15, 15)
  )

# Panel letters are patchwork tags, since the paper puts no titles inside
# panels. Styled per figure via plot_annotation(), never by editing plot_theme.
TAG_THEME <- theme(plot.tag = element_text(face = "bold", size = 20))

# ---- helpers -----------------------------------------------------------------
read_trial <- function(path) {
  vroom(path, show_col_types = FALSE) %>% dplyr::select(-c(2, 3, 4))
}

# The two trial matrices name lyso species differently: CTL writes them
# CLASS(x:y/0:0) (a headgroup with an empty second chain), LIN writes them
# LCLASS(x:y). Left alone, the same chemical species is counted as PC in one
# trial and LPC in the other. Normalise before doing anything else.
# Note the 3-position guard: DG(18:0/18:2/0:0) is a diacylglycerol, not a lyso.
normalize_lipid_name <- function(x) {
  x <- sub("^(PC|PE|PG|PS|PA)\\(([^/()]+)/0:0\\)$", "L\\1(\\2)", x)
  x <- sub("^(LPC|LPE)\\(([^/()]+)/0:0\\)$", "\\1(\\2)", x)
  x
}

lipid_class <- function(x) {
  cls <- sub("\\(.*$", "", normalize_lipid_name(x))
  ifelse(cls %in% names(class_colors), cls, "Other")
}

save_fig <- function(plot, filename, width, height, subdir = "main", dpi = 300) {
  dir <- if (subdir == "main") FIG_MAIN else FIG_SUPP
  path <- file.path(dir, filename)
  ggsave(path, plot, width = width, height = height, dpi = dpi,
         bg = "white", limitsize = FALSE)
  message("Saved: ", path)
  invisible(path)
}

save_table <- function(df, filename) {
  path <- file.path(TAB_SUPP, filename)
  write.csv(df, path, row.names = FALSE)
  message("Saved: ", path, "  (", nrow(df), " rows)")
  invisible(path)
}

# ---- shared data builders ----------------------------------------------------
# These are here, not in a figure script, because more than one figure needs
# them and they must agree exactly between figures. Each figure still owns its
# own plotting code in its own file.

CLASS_ORDER <- names(class_colors)
PURE        <- c("Bicolor", "Caudatum", "Durra", "Guinea", "Kafir")
RACE_ORDER  <- c(PURE, "Mixed")

race_colors <- c(Bicolor = "#E69F00", Caudatum = "#0072B2", Durra = "#009E73",
                 Guinea = "#CC79A7", Kafir = "#D55E00", Mixed = "#999999")
k_colors    <- c("1" = "#E69F00", "2" = "#56B4E9", "3" = "#009E73",
                 "4" = "#F0E442", "5" = "#0072B2", "6" = "#CC79A7")

# The five classical races, with every intermediate designation collapsed into
# Mixed. S. verticilliflorum is a wild relative and is dropped.
race_group <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x %in% c("", "NA")] <- NA_character_
  x[grepl("verticilliflorum", x, ignore.case = TRUE)] <- NA_character_
  ifelse(is.na(x), NA_character_, ifelse(x %in% PURE, x, "Mixed"))
}

# Per-feature share of total annotated signal, taken per sample then averaged,
# so a few high-signal samples cannot dominate the mean composition.
pct_tic <- function(path, label) {
  d <- read_trial(path)
  feats <- names(d)[-1]
  m <- as.matrix(d[, feats]); storage.mode(m) <- "numeric"
  m[!is.finite(m)] <- 0
  share <- sweep(m, 1, pmax(rowSums(m), 1e-12), "/") * 100
  tibble::tibble(Feature = feats, pct = colMeans(share, na.rm = TRUE), Condition = label)
}

# Per-accession class composition. TotalLipid is the summed signal over every
# annotated feature, so the 13 focal classes do not add up to 100%; the
# remainder is sterols, carotenoids and the other non-focal categories.
class_percent_tic <- function(path, condition) {
  x <- read_trial(path)
  geno <- x[[1]]
  m <- as.matrix(x[, -1, drop = FALSE]); storage.mode(m) <- "numeric"
  m[!is.finite(m)] <- 0
  cls <- lipid_class(colnames(m))
  total <- rowSums(m)
  out <- tibble::tibble(LineRaw = geno, Condition = condition, TotalLipid = total)
  for (cl in CLASS_ORDER) {
    j <- which(cls == cl)
    out[[cl]] <- if (length(j)) 100 * rowSums(m[, j, drop = FALSE]) / total else 0
  }
  out
}

# Class composition for both trials, joined to botanical race and marker-based
# genetic cluster. Used by the population-structure figure and its PCA.
population_table <- function() {
  geoloc_csv <- Sys.getenv("SAP_GEOLOC", file.path(DATA_ROOT, "SAP_geoloc.csv"))
  stopifnot(file.exists(CTL_CSV), file.exists(LIN_CSV), file.exists(geoloc_csv))
  geo <- vroom(geoloc_csv, show_col_types = FALSE)
  names(geo)[1] <- "Taxa"                      # the file carries a UTF-8 BOM
  geo <- geo %>% dplyr::select(Taxa, K.Cluster, Original_Race) %>%
    distinct(Taxa, .keep_all = TRUE)
  bind_rows(class_percent_tic(CTL_CSV, "CTL"), class_percent_tic(LIN_CSV, "LIN")) %>%
    left_join(geo, by = c("LineRaw" = "Taxa")) %>%
    mutate(Condition        = factor(Condition, c("CTL", "LIN")),
           RaceGroup        = factor(race_group(Original_Race), RACE_ORDER),
           KCluster         = factor(K.Cluster, levels = 1:6),
           TotalLipid_log10 = log10(TotalLipid))
}

# Kruskal-Wallis with epsilon^2 = (H - k + 1) / (n - k), the proportion of rank
# variance explained, truncated at zero. Groups smaller than min_per_group are
# dropped before testing.
kw_eps2 <- function(values, groups, min_per_group = 5) {
  ok <- !is.na(values) & !is.na(groups)
  values <- values[ok]; groups <- droplevels(factor(groups[ok]))
  big <- names(which(table(groups) >= min_per_group))
  ok <- groups %in% big
  values <- values[ok]; groups <- droplevels(factor(groups[ok]))
  k <- nlevels(groups); n <- length(values)
  if (k < 3L) return(NULL)
  kt <- kruskal.test(values, groups)
  H <- unname(kt$statistic)
  tibble::tibble(n = n, k_groups = k, H = H, p = kt$p.value,
                 epsilon2 = max(0, (H - k + 1) / (n - k)))
}

# PCA of the 13 class percentages within one trial, on log10 values scaled to
# unit variance so no single dominant class sets the axes.
trial_class_pca <- function(dat, cond) {
  sub <- dat %>% filter(Condition == cond)
  m <- as.matrix(sub[, CLASS_ORDER, drop = FALSE])
  m <- log10(pmax(m, 0) + 1e-4)
  keep <- apply(m, 2, sd, na.rm = TRUE) > 0
  p <- prcomp(m[, keep, drop = FALSE], center = TRUE, scale. = TRUE)
  list(scores = sub %>% mutate(PC1 = p$x[, 1], PC2 = p$x[, 2]),
       ve = p$sdev^2 / sum(p$sdev^2))
}
