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
