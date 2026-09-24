#!/usr/bin/env Rscript
# =============================================================================
# Rebuild the Shiny app's GWAS objects from the post-filter, post-deduplication
# association files. Run this on the server, where the .assoc.txt files live.
#
#   Rscript 60_rebuild_shiny_gwas_rds.R
#
# Environment overrides
#   ROOT_IND   parent of <control|lowinput>/<GWAS_SUBDIR> for individual lipids
#   ROOT_SR    the same for class sums and ratios
#   GWAS_SUBDIR   folder holding the .assoc.txt files (default raw_gwas_hetfiltered,
#                 falls back to raw_gwas when that folder is absent)
#   GFF        Sorghum GFF3
#   OUTDIR     repository root, files land under OUTDIR/data/GWAS_RDS/
#   REMOVED    high-het removal list, used only to verify the inputs are filtered
#   FLANK_BP   physical window each SNP is given, default 25000
#
# What it writes, per condition and layer, for each of three p thresholds
#   data/GWAS_RDS/<individual|sum_ratio>/all_annotations_<cond>_<layer>_plog10N.rds
#   data/GWAS_RDS/<individual|sum_ratio>/<layer>_gwas_manifest_<cond>_plog10N.csv
#   data/GWAS_RDS/<individual|sum_ratio>/unique_genes_<cond>_<layer>_plog10N.txt
#
# The RDS is a named list, one element per trait, each a data.frame of
# GeneID, Chromosome, `log(p)`. That is the shape app.R expects. Trait names are
# the bare names, matching gene_trait_LD.csv, so the Gene Hits module joins.
#
# Each file is read once and annotated once at p <= 0.05; the three thresholds
# are then taken as subsets of the same gene table, so the expensive overlap
# step runs a third as often as running the script three times.
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(IRanges); library(rtracklayer)
})

ROOT_IND <- Sys.getenv('ROOT_IND',
  '/rsstu/users/r/rrellan/DOE_CAREER/SAP/results/spats_corrected/BLUP/individual_lipids_final')
ROOT_SR  <- Sys.getenv('ROOT_SR',
  '/rsstu/users/r/rrellan/DOE_CAREER/SAP/results/spats_corrected/BLUP/sum_ratio_lipids_final')
GWAS_SUBDIR <- Sys.getenv('GWAS_SUBDIR', 'raw_gwas_hetfiltered')
GFF      <- Sys.getenv('GFF',
  '/rsstu/users/r/rrellan/DOE_CAREER/SAP/ref/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.54.gff3')
OUTDIR   <- Sys.getenv('OUTDIR', '.')
REMOVED  <- Sys.getenv('REMOVED', '')
FLANK    <- as.integer(Sys.getenv('FLANK_BP', '25000'))
ANNOT_P  <- 0.05
THRESH   <- c(plog105 = 1e-5, plog106 = 1e-6, plog107 = 1e-7)

LAYERS <- list(
  individual = list(root = ROOT_IND, dir = 'individual', stem = 'individual'),
  sum_ratio  = list(root = ROOT_SR,  dir = 'sum_ratio',  stem = 'sum_ratio'))
CONDS <- c('control', 'lowinput')

gwas_dir <- function(root, cond) {
  d <- file.path(root, cond, GWAS_SUBDIR)
  if (dir.exists(d)) return(d)
  alt <- file.path(root, cond, 'raw_gwas')
  if (dir.exists(alt)) { message('  ', GWAS_SUBDIR, ' absent, using raw_gwas'); return(alt) }
  stop('no association folder under ', file.path(root, cond))
}

# Trait name from the file stem. Everything the GWAS runner appended is dropped,
# so the names match gene_trait_LD.csv and the app's ld_normalise_trait().
trait_of <- function(path) {
  n <- basename(path)
  n <- sub('\\.assoc\\.txt$', '', n); n <- sub('\\.txt$', '', n)
  n <- sub('_mod_sub_.*$', '', n)
  n <- sub('_SAP_bialleles.*$', '', n)
  n
}
norm_chr <- function(x) { x <- sub('^[Cc]hr_?', '', as.character(x)); suppressWarnings(as.integer(x)) }

message('reading gene models from ', GFF)
gr <- rtracklayer::import(GFF)
genes <- gr[mcols(gr)$type == 'gene']
mcols(genes)$GeneID <- sub('^gene:', '', mcols(genes)$ID)
seqlevels(genes) <- sub('^[Cc]hr_?', '', seqlevels(genes))
message('  ', length(genes), ' genes')

bad_keys <- NULL
if (nzchar(REMOVED) && file.exists(REMOVED)) {
  b <- fread(REMOVED, showProgress = FALSE)
  bad_keys <- paste(as.integer(b$chr), b$pos, sep = ':')
  message('removal list loaded, ', length(bad_keys), ' markers')
}

annotate_one <- function(f) {
  d <- fread(f, select = c('chr','rs','ps','p_wald'), showProgress = FALSE)
  setnames(d, c('chr','rs','ps','p_wald'))
  d[, chr := norm_chr(chr)]
  d <- d[is.finite(chr) & is.finite(ps) & is.finite(p_wald)]
  n_snps <- nrow(d)
  n_bad <- if (is.null(bad_keys)) NA_integer_ else
             sum(paste(d$chr, d$ps, sep = ':') %chin% bad_keys)
  s <- d[p_wald <= ANNOT_P]
  if (!nrow(s)) return(list(tab = data.table(GeneID = character(), Chromosome = character(),
                                             minp = numeric()), n_snps = n_snps, n_bad = n_bad))
  snp <- GRanges(seqnames = as.character(s$chr),
                 ranges = IRanges(start = s$ps, end = s$ps), pvalue = s$p_wald)
  win <- resize(snp, width = 2L * FLANK + 1L, fix = 'center')
  ov  <- findOverlaps(genes, win)
  if (!length(ov)) return(list(tab = data.table(GeneID = character(), Chromosome = character(),
                                                minp = numeric()), n_snps = n_snps, n_bad = n_bad))
  tab <- data.table(GeneID = mcols(genes)$GeneID[queryHits(ov)],
                    Chromosome = as.character(seqnames(genes)[queryHits(ov)]),
                    p = mcols(win)$pvalue[subjectHits(ov)])
  tab <- tab[, .(minp = min(p)), by = .(GeneID, Chromosome)]
  list(tab = tab, n_snps = n_snps, n_bad = n_bad)
}

for (ln in names(LAYERS)) {
  L <- LAYERS[[ln]]
  outdir <- file.path(OUTDIR, 'data', 'GWAS_RDS', L$dir)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  for (cond in CONDS) {
    dir_in <- gwas_dir(L$root, cond)
    fs <- list.files(dir_in, pattern = '\\.assoc\\.txt$|\\.txt$', full.names = TRUE)
    message('\n== ', ln, ' ', cond, '  ', length(fs), ' files in ', dir_in)
    if (!length(fs)) { warning('no files, skipping'); next }

    res <- vector('list', length(fs)); nm <- character(length(fs))
    meta <- data.table(trait = character(), n_snps = numeric(), n_high_het = numeric())
    for (k in seq_along(fs)) {
      r <- annotate_one(fs[k])
      nm[k] <- trait_of(fs[k]); res[[k]] <- r$tab
      meta <- rbind(meta, data.table(trait = nm[k], n_snps = r$n_snps, n_high_het = r$n_bad))
      if (k %% 25 == 0) message('   ', k, ' / ', length(fs))
    }
    names(res) <- nm
    if (!is.null(bad_keys)) {
      tot <- sum(meta$n_high_het, na.rm = TRUE)
      if (tot > 0) warning('!! ', tot, ' high-het markers still present in ', ln, ' ', cond,
                           ' -- point GWAS_SUBDIR at the filtered files')
      else message('   verified, no high-het markers in these inputs')
    }

    for (tn in names(THRESH)) {
      thr <- THRESH[[tn]]
      lst <- lapply(res, function(t) {
        if (!nrow(t)) return(data.frame(GeneID = character(), Chromosome = character(),
                                        `log(p)` = numeric(), check.names = FALSE))
        u <- t[minp <= thr]
        data.frame(GeneID = as.character(u$GeneID),
                   Chromosome = as.character(u$Chromosome),
                   `log(p)` = as.numeric(-log10(u$minp)),
                   check.names = FALSE, stringsAsFactors = FALSE)
      })
      names(lst) <- names(res)
      saveRDS(lst, file.path(outdir,
        sprintf('all_annotations_%s_%s_%s.rds', cond, L$stem, tn)))

      ug <- sort(unique(unlist(lapply(lst, `[[`, 'GeneID'))))
      write.table(ug, file.path(outdir,
        sprintf('unique_genes_%s_%s_%s.txt', cond, L$stem, tn)),
        row.names = FALSE, col.names = FALSE, quote = FALSE)

      mf <- copy(meta)
      mf[, n_genes_in_rds := vapply(lst, nrow, integer(1))[trait]]
      fwrite(mf[, .(trait, n_snps, n_genes_in_rds)], file.path(outdir,
        sprintf('%s_gwas_manifest_%s_%s.csv', L$stem, cond, tn)))

      message(sprintf('   %s  traits %d  genes %d', tn, length(lst), length(ug)))
    }
  }
}
message('\ndone. Copy data/GWAS_RDS/ and data/LD_mapped/candidate_tables/gene_trait_LD.csv into the app bundle.')
