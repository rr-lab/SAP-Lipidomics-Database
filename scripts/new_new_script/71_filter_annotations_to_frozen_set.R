# ==============================================================================
# Reconcile the GWAS gene annotations with the frozen lipid species set.
#
#   Rscript scripts/new_new_script/71_filter_annotations_to_frozen_set.R
#
# Why this exists. The GWAS was run in March on the phenotype tables as they stood
# then. The lipid species set was curated afterwards, on 2026-09-16, and that
# curation did two things the annotation tables never saw
#
#   1. It merged annotation pairs that point at the same MS feature, keeping one
#      name per feature. 13-Keto-9Z,11E-octadecadienoic acid and 13S-HOTrE are the
#      same scan at the same m/z with the same formula; the curated set keeps
#      13S-HOTrE. Both names still carry GWAS hits, over an identical 420-gene set
#      in LIN, so every one of those genes is credited with two phenotypes where it
#      should have one.
#   2. It settled which species count as detected in which trial, so a handful of
#      names carry hits in a trial whose frozen set does not contain them.
#
# The rule applied here is the one the manuscript states: a candidate gene is
# reported for a lipid species only if that species is in that trial's frozen set
# (214 under CTL, 216 under LIN, from data/final_species_set/species_inventory.csv).
# Class sums and ratios are not species and pass through untouched.
#
# Two naming conventions have to be reconciled before the comparison is meaningful
#   - GWAS trait names replace "/" with "_" so they are safe as file names
#   - lyso species were written with an explicit 0:0 tail before curation;
#     PC(18:2/0:0) is LPC(18:2), and DG(18:0/18:2/0:0) is DG(18:0/18:2)
# Traits are matched on the reconciled name and written out under the curated name,
# so the candidate tables and the species tables use one vocabulary.
#
# Inputs and outputs both live in data/gene_annotation_final/. The untouched files
# are moved once to data/gene_annotation_final/raw/ and read from there on every
# run, so the script is idempotent and the originals are never lost.
#
# Everything downstream of the annotation tables -- Supplementary Tables S7-S10 via
# 72_rebuild_S7toS10.R and Table 4 via 73_top_loci_table.R -- must be rebuilt after
# this runs.
# ==============================================================================
suppressPackageStartupMessages(library(data.table))

REPO   <- Sys.getenv("SOLD_REPO", ".")
ANNDIR <- file.path(REPO, "data/gene_annotation_final")
RAWDIR <- file.path(ANNDIR, "raw")
INV    <- file.path(REPO, "data/final_species_set/species_inventory.csv")
FILES  <- c(CTL = "CTL_gene_annotation.tsv", LIN = "LIN_gene_annotation.tsv")

dir.create(RAWDIR, recursive = TRUE, showWarnings = FALSE)
for (f in FILES) {
  if (!file.exists(file.path(RAWDIR, f))) {
    stopifnot(file.exists(file.path(ANNDIR, f)))
    file.copy(file.path(ANNDIR, f), file.path(RAWDIR, f))
    message("Archived unfiltered annotations: raw/", f)
  }
}

inv <- fread(INV)

# GWAS trait name -> curated species name -------------------------------------
reconcile <- function(x) {
  s <- gsub("_", "/", x)
  # PC(18:2/0:0) and PE(16:0/0:0) are the lyso species
  s <- sub("^([A-Z]+)\\(([^/()]+)/0:0\\)$", "L\\1(\\2)", s)
  # DG(18:0/18:2/0:0) is a two-chain DG written with an explicit empty sn-3
  s <- sub("^DG\\(([^/()]+)/([^/()]+)/0:0\\)$", "DG(\\1/\\2)", s)
  s
}

report <- rbindlist(lapply(names(FILES), function(cond) {
  g <- fread(file.path(RAWDIR, FILES[[cond]]), sep = "\t", colClasses = "character")
  frozen <- inv[[if (cond == "CTL") "In_CTL" else "In_LIN"]]
  keep_sp <- inv$Species[frozen]

  g[, Species := reconcile(Trait)]
  is_ind <- g$Layer == "individual_lipids_final"
  ok     <- !is_ind | g$Species %in% keep_sp

  dropped <- g[is_ind & !ok, .(Genes = .N), by = .(Trait, Species)]
  dropped[, `:=`(Condition = cond,
                 Reason = fifelse(Species %in% inv$Species,
                                  "species not detected in this trial after curation",
                                  "species not in the frozen set"))]

  # write the curated species name into Trait for the individual layer
  g[is_ind, Trait := Species]
  g[, Species := NULL]
  out <- g[ok]
  fwrite(out, file.path(ANNDIR, FILES[[cond]]), sep = "\t", quote = FALSE)

  message(sprintf("%s  %d -> %d rows, %d individual traits -> %d",
                  FILES[[cond]], nrow(g), nrow(out),
                  uniqueN(g[is_ind, Trait]), uniqueN(out[Layer == "individual_lipids_final", Trait])))
  dropped[, .(Condition, Trait, Species, Genes, Reason)]
}))

cat("\nTraits removed\n")
print(as.data.frame(report))
fwrite(report, file.path(ANNDIR, "traits_removed_vs_frozen_set.tsv"), sep = "\t", quote = FALSE)
message("\nSaved: data/gene_annotation_final/traits_removed_vs_frozen_set.tsv")
