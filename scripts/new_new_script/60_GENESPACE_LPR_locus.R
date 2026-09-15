# ==============================================================================
# Comparative genomics of the chromosome 3 LPR locus across five grasses
#
#   Rscript scripts/new_new_script/60_GENESPACE_LPR_locus.R
#
# Asks what the LOW PHOSPHATE ROOT locus looks like outside sorghum. Rice
# carries five LPR genes in a 65 kb tandem array on chromosome 1 (Cao et al.
# 2016); sorghum contributes two adjacent candidates, SORBI_3003G088800 and
# SORBI_3003G088900, 1.8 kb apart on chromosome 3. This run asks whether those
# two sit in the block syntenic with the rice array, and how many copies maize,
# Setaria and Brachypodium carry.
#
# Adapted from scripts/new_script/15.GENESPACE.R, which ran maize against
# sorghum for a different locus. Three changes from that script. Five genomes
# rather than two, the per-species GFF3 parsing folded into one function with
# a hard check that the BED and FASTA actually intersect, and the region of
# interest set to the LPR pair with query_pangenes rather than query_hits, so
# the output is orthologues per genome rather than raw blast hits.
#
# THIS DOES NOT RUN IN THE CLOUD. It needs OrthoFinder and a compiled MCScanX
# and roughly 20-40 GB of scratch, so run it locally.
#
#   conda create -n orthofinder -c bioconda orthofinder
#   conda activate orthofinder
#   git clone https://github.com/wyp1125/MCScanX.git && cd MCScanX && make
#   open -na rstudio
#
# Inputs, all from Phytozome, one directory per genome under GENOME_REPO, each
# holding the gene GFF3 and the primary-transcript peptide FASTA
#
#   sorghum      Sbicolor_454_v3.1.1                  (v3, matches SORBI_3003G ids)
#   maize        Zmays_833_Zm-B73-REFERENCE-NAM-5.0.55
#   rice         Osativa_323_v7.0
#   setaria      Sitalica_312_v2.2
#   brachy       Bdistachyon_314_v3.1
#
#   <genome>/*.gene.gff3.gz
#   <genome>/*.protein_primaryTranscriptOnly.fa.gz
#
# Sorghum v3.1.1 is not optional. The candidate ids are v3; any other assembly
# and none of them resolve.
#
# Output   results/genespace/  GENESPACE working directory
#          fig/supp/SuppFig_LPR_riparian.png
#          table/supp/SuppTable_LPR_pangenes.tsv
# ==============================================================================
suppressPackageStartupMessages({
  library(GENESPACE); library(Biostrings); library(seqinr)
  library(dplyr); library(tidyr); library(ggplot2)
})

GENOME_REPO <- path.expand("~/Downloads/GENESPACE_LPR/genomeRepo")
WD          <- path.expand("~/Downloads/GENESPACE_LPR/workingDirectory")
MCSCANX     <- path.expand("~/GitHub/MCScanX")
OUT_FIG     <- "fig/supp"
OUT_TAB     <- "table/supp"

GENOMES <- c("sorghum", "maize", "rice", "setaria", "brachy")

# The locus. Sorghum chromosome 3, the two adjacent LPR candidates plus a
# little flanking sequence so the syntenic block is not clipped.
ROI_CHR   <- "chr3"
ROI_START <- 7700000
ROI_END   <- 7780000
FOCAL     <- c("SORBI_3003G088800", "SORBI_3003G088900")

## ---- GFF3 and peptide FASTA to the BED and FASTA GENESPACE wants -----------
# Phytozome names the peptide sequences by transcript, and each assembly spells
# the transcript slightly differently from its own GFF3. The fix is per genome
# and it is the only species-specific part of this script. If a fix is wrong
# the intersect comes back empty, which the stopifnot below catches rather than
# letting the run continue on nothing.
ID_FIX <- list(
  sorghum = function(x) sub("v3\\.2", "p", x),
  maize   = function(x) sub("T([0-9]+)$", "P\\1", x),
  rice    = function(x) x,
  setaria = function(x) x,
  brachy  = function(x) x
)

CHR_KEEP <- list(
  sorghum = paste0("chr", 1:10), maize   = paste0("chr", 1:10),
  rice    = paste0("chr", 1:12), setaria = paste0("chr", 1:9),
  brachy  = paste0("chr", 1:5)
)

prepare_genome <- function(g) {
  dir  <- file.path(GENOME_REPO, g)
  gff  <- list.files(dir, pattern = "gene\\.gff3(\\.gz)?$", full.names = TRUE)[1]
  faa  <- list.files(dir, pattern = "protein.*\\.fa(\\.gz)?$", full.names = TRUE)[1]
  stopifnot(!is.na(gff), !is.na(faa))

  x <- read.table(gff, sep = "\t", quote = "", comment.char = "#",
                  stringsAsFactors = FALSE)
  x <- x[x$V3 == "mRNA", c(1, 4, 5, 9)]

  # the attribute column is key=value pairs; the first value is the mRNA id
  x$id <- vapply(strsplit(x$V9, ";"), function(a)
    sub("^[^=]+=", "", a[1]), character(1))
  x$id <- ID_FIX[[g]](x$id)

  x$V1 <- sub("^Chr0?", "chr", x$V1)   # Chr01 and Chr1 both become chr1
  x$V1 <- sub("^chr0", "chr", x$V1)
  x <- x[x$V1 %in% CHR_KEEP[[g]], c("V1", "V4", "V5", "id")]
  x <- x[order(factor(x$V1, levels = CHR_KEEP[[g]]), x$V4), ]

  pep    <- read.fasta(faa, seqtype = "AA")
  names(pep) <- sub("\\s.*$", "", names(pep))
  shared <- intersect(x$id, names(pep))

  message(sprintf("%-8s gff mRNA %6d | peptides %6d | shared %6d",
                  g, nrow(x), length(pep), length(shared)))
  stopifnot(length(shared) > 0.5 * nrow(x))   # a bad ID_FIX fails here

  dir.create(file.path(WD, "bed"),     recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(WD, "peptide"), recursive = TRUE, showWarnings = FALSE)
  write.table(x[x$id %in% shared, ], file.path(WD, "bed", paste0(g, ".bed")),
              row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  write.fasta(pep[shared], names = shared,
              file.out = file.path(WD, "peptide", paste0(g, ".fa")))
}

invisible(lapply(GENOMES, prepare_genome))

## ---- run ------------------------------------------------------------------
gpar <- init_genespace(wd = WD, path2mcscanx = MCSCANX)
out  <- run_genespace(gpar, overwrite = TRUE)

## ---- the locus ------------------------------------------------------------
roi <- data.frame(genome = "sorghum", chr = ROI_CHR,
                  start = ROI_START, end = ROI_END)

pg <- query_pangenes(gsParam = out, bed = roi)
write.table(pg, file.path(OUT_TAB, "SuppTable_LPR_pangenes.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# how many copies each genome carries in the syntenic block
cat("\n--- orthologues per genome in the sorghum chr3 LPR block ---\n")
print(pg)

cat("\n--- the two focal genes ---\n")
print(pg[apply(pg, 1, function(r) any(grepl(paste(FOCAL, collapse = "|"), r))), ])

## ---- figure ---------------------------------------------------------------
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
png(file.path(OUT_FIG, "SuppFig_LPR_riparian.png"),
    width = 12, height = 7, units = "in", res = 300, bg = "white")
plot_riparian(gsParam = out, refGenome = "sorghum", useRegions = FALSE,
              highlightBed = roi[, c("genome", "chr", "start", "end")])
dev.off()
