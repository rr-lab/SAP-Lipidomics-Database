#!/usr/bin/env Rscript
# Rebuild data/LD_mapped/candidate_tables/gene_trait_LD.csv, the table behind the
# Shiny Gene Hits module, from the post-filter LD annotation tables.
#
# Schema the app expects, read at app.R line 676
#   condition, layer, GeneID, trait, max_r2, best_p
#   condition is CTL / LIN, layer is individual / sumratio
suppressPackageStartupMessages({library(data.table)})
IN  <- Sys.getenv('SOLD_ANN', '/mnt/user-data/uploads/SAP-Lipidomics-Database/data/gene_annotation_final')
OUT <- Sys.getenv('SOLD_OUT', '/home/claude/appgen')

A <- rbind(fread(file.path(IN,'CTL_gene_annotation.tsv')),
           fread(file.path(IN,'LIN_gene_annotation.tsv')))
lay <- c(individual_lipids_final = 'individual', sum_ratio_lipids_final = 'sumratio')
out <- A[, .(condition = Condition,
             layer     = unname(lay[Layer]),
             GeneID,
             trait     = Trait,
             max_r2    = Max_r2,
             best_p    = Best_P)]
stopifnot(!any(is.na(out$layer)))
# one row per condition x layer x gene x trait, keeping the strongest evidence
out <- out[, .(max_r2 = max(max_r2), best_p = min(best_p)),
           by = .(condition, layer, GeneID, trait)]
setorder(out, condition, layer, GeneID, trait)
fwrite(out, file.path(OUT, 'gene_trait_LD.csv'))
cat(sprintf('rows %d   genes %d   traits %d\n', nrow(out),
            uniqueN(out$GeneID), uniqueN(out$trait)))
print(out[, .(rows = .N, genes = uniqueN(GeneID), traits = uniqueN(trait)),
          by = .(condition, layer)])

## reproduce the two numbers the manuscript quotes
cat('\n-- Gene Hits, LIN, all trait sources, r2 >= 0.9 --\n')
d <- out[condition == 'LIN' & max_r2 >= 0.9]
g <- d[, .(Phenotypes = uniqueN(trait), Highest_logp = round(-log10(min(best_p)), 2),
           Max_r2 = round(max(max_r2), 3)), by = GeneID][order(-Phenotypes, -Highest_logp)]
cat(sprintf('candidate genes output  %d\n', nrow(g)))
print(head(g, 8))
