# ==============================================================================
# Supplementary Figure 1 -- Field-to-database workflow.
#
# A single hand-made schematic (FigureLabs.ai), not a plotted panel, so this
# script only places it and writes it out at the right size. It replaces the
# former analytical-comparability figure, whose QC panels now live in
# Supplementary Figure S1.
#
# Input : fig/main/individual_figs/Fig1A_working_pipeline.png
# Output: fig/supp/SuppFig_S1_Workflow.png
# ==============================================================================
source("scripts/new_new_script/_common.R")
suppressPackageStartupMessages(library(png))

schematic <- Sys.getenv("SUPPFIG1_SCHEMATIC",
  file.path(FIG_MAIN, "individual_figs/Fig1A_working_pipeline.png"))

if (!file.exists(schematic))
  stop("Workflow schematic not found at: ", schematic,
       "\n  Export it from FigureLabs.ai and save it there, or set SUPPFIG1_SCHEMATIC.")

raster <- png::readPNG(schematic)
dims   <- dim(raster)                      # rows (height) x cols (width) x channels
aspect <- dims[2] / dims[1]

fig <- ggplot() +
  annotation_custom(grid::rasterGrob(raster, interpolate = TRUE),
                    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  coord_cartesian(expand = FALSE) +
  theme_void()

# Height follows the image's own aspect ratio, so the schematic fills the width
# instead of being shrunk inside an over-wide panel.
w <- 13.5
save_fig(fig, "SuppFig_S1_Workflow.png", width = w, height = w / aspect, subdir = "supp")
message(sprintf("schematic %d x %d px (aspect %.2f)", dims[2], dims[1], aspect))
