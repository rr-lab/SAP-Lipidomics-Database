# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6: LINEX-BASED LIPID REACTION NETWORK ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════
#
# This script generates ONLY the figures and tables for Section 6:
# "LINEX-based lipid reaction network analysis"
#
# OUTPUTS:
#   - Figure 4 (Main): LINEX Enriched Subnetwork
#   - Table 2 (Main): LINEX GWAS support (in LaTeX document)
#   - Supplementary Table S11: LINEX Reactions with RHEA annotations
#
# ══════════════════════════════════════════════════════════════════════════════

message("\n")
message("══════════════════════════════════════════════════════════════")
message("SECTION 6: LINEX-BASED LIPID REACTION NETWORK ANALYSIS")
message("══════════════════════════════════════════════════════════════")
message("\n")

# ──────────────────────────────────────────────────────────────────────────────
# SETUP
# ──────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(igraph)
  library(stringr)
  library(dplyr)
  library(tibble)
})

# Create output directories
dir.create("fig/main", recursive = TRUE, showWarnings = FALSE)
dir.create("table/supp", recursive = TRUE, showWarnings = FALSE)


# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
#
#            FIGURE 4 (MAIN): LINEX ENRICHED SUBNETWORK
#
# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════

message("\n══════════════════════════════════════════════════════════════")
message("FIGURE 4 (MAIN): LINEX ENRICHED SUBNETWORK")
message("══════════════════════════════════════════════════════════════\n")

# ──────────────────────────────────────────────────────────────────────────────
# 1) Define lipid species (individual molecules)
# ──────────────────────────────────────────────────────────────────────────────

lipid_species <- c(
  # Monoglycerides
  "MG(16:0)", "MG(18:1)", "MG(18:2)",
  # Diglycerides
  "DG(16:0_18:1)", "DG(16:0_18:2)", "DG(18:1_18:2)",
  # Triglyceride
  "TG(16:0_18:1_18:2)",
  # Phospholipids involved in reactions
  "PC(16:0_18:1)", "PE(16:0_18:1)",
  # Lysophospholipids (products)
  "LPC(16:0)", "LPE(16:0)"
)

# Helper function to extract lipid class
lip_class <- function(x) str_extract(x, "^[A-Za-z]+")

# ──────────────────────────────────────────────────────────────────────────────
# 2) Reaction nodes with enzyme names AND reversible reactions
# ──────────────────────────────────────────────────────────────────────────────

R1_short <- "LCAT*\nPC+DG <-> LPC+TG"
R2_short <- "PNPLA1\nTG <-> DG"
R3_short <- "LRO1\nPE+DG <-> TG+LPE"
R4_short <- "PNPLA3\nDG+MG <-> TG"

reactions_short <- c(R1_short, R2_short, R3_short, R4_short)
enzyme_names <- c("LCAT*", "PNPLA1", "LRO1", "PNPLA3")

# ──────────────────────────────────────────────────────────────────────────────
# 3) Create Supplementary Table S11 with full RHEA info
# ──────────────────────────────────────────────────────────────────────────────

message("Creating Supplementary Table S11 (LINEX Reactions)...")

reaction_table <- tibble(
  Reaction_ID = c("R1", "R2", "R3", "R4"),
  Gene_Symbol = c("LCAT", "PNPLA1", "LRO1", "PNPLA3"),
  Reaction_Type = rep("L_FAdelete", 4),
  Substrates = c("PC, DG", "TG", "PE, DG", "DG, MG"),
  Products = c("LPC, TG", "DG", "TG, LPE", "TG"),
  RHEA_IDs = c(
    "RHEA:32843, RHEA:44236",
    "RHEA:61528, RHEA:44676, RHEA:40615, RHEA:40619",
    "RHEA:44252, RHEA:32847, RHEA:44232",
    "RHEA:44436, RHEA:44440"
  ),
  Notes = c(
    "Hypothesized based on reaction specificity (not from LINEX)",
    "From LINEX database",
    "From LINEX database",
    "From LINEX database"
  )
)

# Save table as S11
write.csv(reaction_table, "table/supp/SuppTable_S11_LINEX_Reactions.csv", row.names = FALSE)
message("Saved: table/supp/SuppTable_S11_LINEX_Reactions.csv")

# ──────────────────────────────────────────────────────────────────────────────
# 4) Node styling (neutral lipid circles, labeled species)
# ──────────────────────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────────────────────
# 5) Build DIRECTED edges showing conversions
# ──────────────────────────────────────────────────────────────────────────────

edges <- data.frame(from = character(), to = character(), stringsAsFactors = FALSE)

# Get species by class
DG_species <- lipid_species[lip_class(lipid_species) == "DG"]
MG_species <- lipid_species[lip_class(lipid_species) == "MG"]
TG_species <- lipid_species[lip_class(lipid_species) == "TG"]
PC_species <- lipid_species[lip_class(lipid_species) == "PC"]
PE_species <- lipid_species[lip_class(lipid_species) == "PE"]
LPC_species <- lipid_species[lip_class(lipid_species) == "LPC"]
LPE_species <- lipid_species[lip_class(lipid_species) == "LPE"]

# R1 (LCAT*): PC + DG → LPC + TG
for (pc in PC_species) edges <- rbind(edges, data.frame(from = pc, to = R1_short))
for (dg in DG_species) edges <- rbind(edges, data.frame(from = dg, to = R1_short))
for (lpc in LPC_species) edges <- rbind(edges, data.frame(from = R1_short, to = lpc))
for (tg in TG_species) edges <- rbind(edges, data.frame(from = R1_short, to = tg))

# R2 (PNPLA1): TG → DG (+ FA released)
for (tg in TG_species) edges <- rbind(edges, data.frame(from = tg, to = R2_short))
for (dg in DG_species) edges <- rbind(edges, data.frame(from = R2_short, to = dg))

# R3 (LRO1): PE + DG → TG + LPE
for (pe in PE_species) edges <- rbind(edges, data.frame(from = pe, to = R3_short))
for (dg in DG_species) edges <- rbind(edges, data.frame(from = dg, to = R3_short))
for (tg in TG_species) edges <- rbind(edges, data.frame(from = R3_short, to = tg))
for (lpe in LPE_species) edges <- rbind(edges, data.frame(from = R3_short, to = lpe))

# R4 (PNPLA3): DG + MG → TG (transacylation)
for (dg in DG_species) edges <- rbind(edges, data.frame(from = dg, to = R4_short))
for (mg in MG_species) edges <- rbind(edges, data.frame(from = mg, to = R4_short))
for (tg in TG_species) edges <- rbind(edges, data.frame(from = R4_short, to = tg))

# ──────────────────────────────────────────────────────────────────────────────
# 6) Build igraph (DIRECTED)
# ──────────────────────────────────────────────────────────────────────────────

all_nodes <- unique(c(lipid_species, reactions_short))
type <- ifelse(all_nodes %in% reactions_short, "reaction", "lipid")

g <- graph_from_data_frame(edges, directed = TRUE, vertices = data.frame(name = all_nodes))

V(g)$type <- type[match(V(g)$name, all_nodes)]
V(g)$shape <- ifelse(V(g)$type == "reaction", "square", "circle")

# Use neutral fills (no class coloring) so species labels drive interpretation
V(g)$lipid_class <- lip_class(V(g)$name)
V(g)$color <- ifelse(
  V(g)$type == "reaction",
  "#F7F7F7",
  "white"
)
V(g)$color[is.na(V(g)$color)] <- "grey85"

V(g)$frame.color <- "grey35"
V(g)$size <- ifelse(V(g)$type == "reaction", 42, 28)

# Labels
V(g)$label <- V(g)$name
V(g)$label.color <- "black"
V(g)$label.cex <- ifelse(V(g)$type == "reaction", 1.0, 1.0)

# Font: 1=plain, 2=bold, 3=italic, 4=bold italic
# LCAT* gets italic, other reactions get bold, lipids get plain
V(g)$label.font <- ifelse(grepl("^LCAT", V(g)$name), 3,
                          ifelse(V(g)$type == "reaction", 2, 1))

# Edge styling
E(g)$color <- "grey55"
E(g)$width <- 1.0
E(g)$arrow.size <- 0.3
E(g)$arrow.width <- 1.0

# ──────────────────────────────────────────────────────────────────────────────
# 7) Layout - Force-directed with enzyme row in middle
# ──────────────────────────────────────────────────────────────────────────────

set.seed(42)
lay <- layout_with_fr(g, niter = 2000)

# Identify node indices by type
lip_idx <- which(V(g)$type == "lipid")
rx_idx <- which(V(g)$type == "reaction")

# Group lipids by class for positioning
substrate_classes <- c("PC", "PE", "DG", "MG")
product_classes <- c("LPC", "LPE", "TG")

# Position enzymes in a column in the middle
enzyme_order <- reactions_short
for (i in seq_along(rx_idx)) {
  enz_name <- V(g)$name[rx_idx[i]]
  enz_pos <- match(enz_name, enzyme_order)
  if (!is.na(enz_pos)) {
    lay[rx_idx[i], 1] <- 0
    lay[rx_idx[i], 2] <- 2.5 - (enz_pos - 1) * 1.8
  }
}

# Position substrate lipids on the left (grouped by class)
substrate_idx <- lip_idx[V(g)$lipid_class[lip_idx] %in% substrate_classes]
if (length(substrate_idx) > 0) {
  sub_classes <- V(g)$lipid_class[substrate_idx]
  sub_order <- order(factor(sub_classes, levels = c("PC", "PE", "DG", "MG")), V(g)$name[substrate_idx])
  substrate_idx <- substrate_idx[sub_order]

  y_positions <- seq(from = 1.5, to = -3.5, length.out = length(substrate_idx))
  for (i in seq_along(substrate_idx)) {
    lay[substrate_idx[i], 1] <- -3.5
    lay[substrate_idx[i], 2] <- y_positions[i]
  }
}

# Position product lipids on the right
product_idx <- lip_idx[V(g)$lipid_class[lip_idx] %in% product_classes]
if (length(product_idx) > 0) {
  prod_classes <- V(g)$lipid_class[product_idx]
  prod_order <- order(factor(prod_classes, levels = c("TG", "LPC", "LPE")))
  product_idx <- product_idx[prod_order]

  y_positions <- seq(from = 1.5, to = -2.5, length.out = length(product_idx))
  for (i in seq_along(product_idx)) {
    lay[product_idx[i], 1] <- 3.5
    lay[product_idx[i], 2] <- y_positions[i]
  }
}

# ──────────────────────────────────────────────────────────────────────────────
# 8) Plot Figure 4
# ──────────────────────────────────────────────────────────────────────────────

out_png <- "fig/main/Figure4_LINEX_Enriched_Subnetwork.png"
out_pdf <- "fig/main/Figure4_LINEX_Enriched_Subnetwork.pdf"

plot_network <- function() {
  plot(
    g,
    layout = lay,
    main = "",
    vertex.label = V(g)$label,
    vertex.label.color = V(g)$label.color,
    vertex.label.cex = V(g)$label.cex,
    vertex.label.font = V(g)$label.font,
    vertex.label.family = "sans",
    edge.curved = 0.15,
    asp = 0.7
  )

  # Notes
  legend(
    "bottomright",
    legend = c("* Hypothesized enzyme", "Arrows show reaction direction", "Lipid species shown as labeled circles"),
    bty = "n",
    cex = 1.0,
    text.font = c(3, 1, 1)
  )
}

# PNG
png(out_png, width = 3600, height = 2600, res = 300)
par(mar = c(2.5, 2.5, 2.5, 2.5), cex.main = 1.5)
plot_network()
dev.off()
message("Saved: ", out_png)

# PDF
pdf(out_pdf, width = 11, height = 8)
par(mar = c(2.5, 2.5, 2.5, 2.5), cex.main = 1.5)
plot_network()
dev.off()
message("Saved: ", out_pdf)

message("\n\u2713 Figure 4 (LINEX Network) complete!\n")


# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════════════════

message("\n")
message("══════════════════════════════════════════════════════════════")
message("SECTION 6 COMPLETE!")
message("══════════════════════════════════════════════════════════════")
message("\n")
message("Generated outputs:")
message("  - fig/main/Figure4_LINEX_Enriched_Subnetwork.png")
message("  - fig/main/Figure4_LINEX_Enriched_Subnetwork.pdf")
message("  - table/supp/SuppTable_S11_LINEX_Reactions.csv")
message("\n")
message("Note: Table 2 (LINEX GWAS support) is in the LaTeX document.")
message("\n")

# Print reaction table for reference
message("=== REACTION TABLE (S11) ===")
print(reaction_table)
