#!/usr/bin/env Rscript
# =============================================================================
# Build the Shiny GWAS-module annotation objects by LINKAGE DISEQUILIBRIUM,
# at r2 = 0.4 0.5 0.6 0.7 0.8 0.9. Physical distance is not used anywhere.
#
#   Rscript scripts/new_new_script/62_build_LD_rds.R
#
# Source of the genes   data/gene_annotation_final/{CTL,LIN}_gene_annotation.tsv
#   one row per trait x gene, carrying the largest r2 any variant inside the
#   gene reaches with that trait's lead SNP. A gene enters a level when that r2
#   is at or above the level.
#
# Source of the trait list, so that traits with no significant SNP still appear
# in the app's picker
#   individual   the trait names of the superseded distance-based RDS
#   sum / ratio  the column names of the class-sum phenotype file
#   plus any trait present in the annotation table but missing from those
#
# Writes, per condition x layer x level
#   data/GWAS_RDS/<layer>/all_annotations_<cond>_<layer>_LD<nn>.rds
#   data/GWAS_RDS/<layer>/<layer>_gwas_manifest_<cond>_LD<nn>.csv
#   data/GWAS_RDS/<layer>/unique_genes_<cond>_<layer>_LD<nn>.txt
#
# Each RDS is a named list, one element per trait, each a data.frame of
# GeneID, Chromosome, `log(p)`, r2.
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })
REPO <- Sys.getenv('SOLD_REPO', '.')
ANN  <- Sys.getenv('SOLD_ANN', file.path(REPO,'data/gene_annotation_final'))
PHE  <- Sys.getenv('SOLD_PHE', file.path(REPO,'data/SPATS_fitted/BLUP_GWAS_phenotype'))
RDS  <- Sys.getenv('SOLD_RDS', file.path(REPO,'data/GWAS_RDS'))
OUT  <- Sys.getenv('SOLD_OUT', file.path(REPO,'data/GWAS_RDS'))
LEVELS <- c(0.4,0.5,0.6,0.7,0.8,0.9)

A <- rbind(fread(file.path(ANN,'CTL_gene_annotation.tsv')),
           fread(file.path(ANN,'LIN_gene_annotation.tsv')))
tagc <- c(CTL='control', LIN='lowinput')
lays <- c(individual_lipids_final='individual', sum_ratio_lipids_final='sum_ratio')
empty <- function() data.frame(GeneID=character(), Chromosome=character(),
                               `log(p)`=numeric(), r2=numeric(), check.names=FALSE)

trait_list <- function(cond, st, seen) {
  tag <- tagc[[cond]]
  base <- character(0)
  if (st == 'individual') {
    f <- file.path(RDS,'individual', sprintf('all_annotations_%s_individual_plog105.rds', tag))
    if (file.exists(f)) base <- names(readRDS(f))
  } else {
    f <- file.path(PHE, sprintf('Final_%s_BLUPs_class_sums_and_logratios.csv', tag))
    if (file.exists(f)) base <- setdiff(names(fread(f, nrows=0)), c('Line','LineRaw'))
  }
  sort(union(base, seen))
}

for (cond in c('CTL','LIN')) for (ly in names(lays)) {
  st <- unname(lays[[ly]])
  d  <- A[Condition==cond & Layer==ly]
  tl <- trait_list(cond, st, unique(d$Trait))
  od <- file.path(OUT, st); dir.create(od, recursive=TRUE, showWarnings=FALSE)
  for (thr in LEVELS) {
    sfx <- sprintf('LD%02d', round(thr*10))
    g <- d[Max_r2 >= thr, .(lp = max(-log10(Best_P)), r2 = max(Max_r2)),
           by = .(Trait, GeneID, Chromosome = as.character(Gene_Chr))]
    setkey(g, Trait)
    lst <- setNames(lapply(tl, function(t) {
      u <- g[.(t), nomatch = 0L]
      if (!nrow(u)) return(empty())
      data.frame(GeneID = u$GeneID, Chromosome = u$Chromosome,
                 `log(p)` = u$lp, r2 = u$r2, check.names = FALSE, stringsAsFactors = FALSE)
    }), tl)
    saveRDS(lst, file.path(od, sprintf('all_annotations_%s_%s_%s.rds', tagc[[cond]], st, sfx)))
    ug <- sort(unique(unlist(lapply(lst, `[[`, 'GeneID'))))
    write.table(ug, file.path(od, sprintf('unique_genes_%s_%s_%s.txt', tagc[[cond]], st, sfx)),
                row.names=FALSE, col.names=FALSE, quote=FALSE)
    fwrite(data.table(trait = tl, n_genes_in_rds = vapply(lst, nrow, integer(1))),
           file.path(od, sprintf('%s_gwas_manifest_%s_%s.csv', st, tagc[[cond]], sfx)))
    cat(sprintf('%s %-10s r2>=%.1f   traits %3d   with genes %3d   genes %5d\n',
                cond, st, thr, length(tl), sum(vapply(lst, nrow, integer(1)) > 0), length(ug)))
  }
}
