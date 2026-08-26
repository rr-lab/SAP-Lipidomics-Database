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

# ---- theme -------------------------------------------------------------------
# No titles or subtitles anywhere. Panel identity comes from the tag letter and
# the caption, which is how journals expect it.
plot_theme <- theme_minimal(base_size = 13) +
  theme(axis.text        = element_text(colour = "black", size = 9),
        axis.title       = element_text(face = "bold", size = 11),
        axis.line        = element_line(colour = "black", linewidth = .5),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.title     = element_blank(),
        plot.tag         = element_text(face = "bold", size = 16))

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
