# ==============================================================================
# 40b_relabel_candidate_phenotypes.R
#
# Relabel duplicate phenotype names in the LD-mapping master candidate table.
#
#   Rscript scripts/new_script/40b_relabel_candidate_phenotypes.R
#
# 40_deduplicate_species.R removed species that were duplicate annotations of a
# retained species. An association reported against a removed name is the same
# association against the name that was kept, because the two columns were the
# same measurement. So the removed names are renamed to their retained partner
# and de-duplicated within each gene, rather than deleted.
#
# Nothing is lost. No gene leaves the table, no p-value changes, and no locus
# changes position. Only Phenotypes and N_Phenotypes are rewritten, which lowers
# a few recurrence counts because two names that were counted separately are now
# correctly counted once.
#
# This runs on the master table so that everything derived from it agrees:
# Supplementary Tables S7-S10 (38_export_supp_tables_S7toS10.R) and the main
# recurrent-locus table (new_new_script/06e_MainTable_top_recurrent_loci.R).
# Run it after 40 and before either of those.
#
# Input   data/LD_mapped/candidate_tables/ALL_LD_candidate_genes_master.tsv
#         table/supp/SuppTable_SpeciesDeduplication_Decisions.csv
# Output  the master table, relabelled, original archived under
#         data/LD_mapped/candidate_tables/_pre_dedup/
#         table/supp/SuppTable_PhenotypeRelabel_Summary.csv
#
# Run from the repository root. Idempotent: running it twice changes nothing.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table) })

REPO   <- Sys.getenv("SOLD_DB_REPO", ".")
MASTER <- file.path(REPO, "data/LD_mapped/candidate_tables/ALL_LD_candidate_genes_master.tsv")
DEC    <- file.path(REPO, "table/supp/SuppTable_SpeciesDeduplication_Decisions.csv")
ARCH   <- file.path(REPO, "data/LD_mapped/candidate_tables/_pre_dedup")
OUT    <- file.path(REPO, "table/supp")
for (d in c(ARCH, OUT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(MASTER), file.exists(DEC))

# The Phenotypes field is delimited by "; " (semicolon SPACE), not by ";" alone.
# That matters because one species is named "SPB 18:0;2OH" and carries a
# semicolon inside its own name. Splitting on ";" tears it into "SPB 18:0" and
# "2OH", which silently invents a phenotype and destroys a real one. Splitting on
# "; " is unambiguous, because no phenotype name contains a semicolon followed by
# a space. The guard below enforces that assumption rather than trusting it.
SEP <- "; "

m <- fread(MASTER, sep = "\t", quote = "", colClasses = "character")
if (!file.exists(file.path(ARCH, basename(MASTER))))
  invisible(file.copy(MASTER, file.path(ARCH, basename(MASTER))))

# candidate tables write acyl separators as "_", the species tables as "/"
u  <- function(x) gsub("/", "_", x, fixed = TRUE)
dd <- unique(fread(DEC)[, .(trial, kept, dropped)])
maps <- lapply(split(dd, dd$trial), function(z) setNames(u(z$kept), u(z$dropped)))

bad <- unlist(lapply(maps, function(z) c(names(z), unname(z))))
bad <- unique(bad[grepl(SEP, bad, fixed = TRUE)])
if (length(bad))
  stop("these species names contain the list separator \"; \" and would corrupt ",
       "the Phenotypes field:\n  ", paste(bad, collapse = "\n  "), call. = FALSE)

# row-wise: GeneID is not unique across the individual and sumratio layers,
# so the before/after comparison has to stay positional.
before_ph <- m$Phenotypes
before_n  <- as.integer(m$N_Phenotypes)

# ---- relabel, as plain vectors so nothing is recycled ------------------------
new_ph <- vapply(seq_len(nrow(m)), function(j) {
  if (is.na(before_ph[j])) return(NA_character_)
  v  <- trimws(strsplit(before_ph[j], SEP, fixed = TRUE)[[1]]); v <- v[nzchar(v)]
  mp <- maps[[ m$condition[j] ]]
  if (!is.null(mp)) { h <- v %in% names(mp); if (any(h)) v[h] <- unname(mp[v[h]]) }
  paste(unique(v), collapse = SEP)
}, character(1))

new_n   <- vapply(seq_along(new_ph), function(j)
  if (is.na(new_ph[j]) || !nzchar(new_ph[j])) 0L
  else length(strsplit(new_ph[j], SEP, fixed = TRUE)[[1]]), integer(1))
# NA-safe, so an empty phenotype list can never be counted as a change
changed <- !is.na(new_ph) & !is.na(before_ph) & new_ph != before_ph

chg <- data.table(condition = m$condition, layer = m$layer, GeneID = m$GeneID,
                  phenotypes_before = before_n, phenotypes_after = new_n)[changed]

m[, Phenotypes   := new_ph]
m[, N_Phenotypes := as.character(new_n)]
stopifnot(!anyNA(m$N_Phenotypes))

fwrite(m, MASTER, sep = "\t", quote = FALSE)
fwrite(chg[order(condition, layer, -phenotypes_before)],
       file.path(OUT, "SuppTable_PhenotypeRelabel_Summary.csv"))

cat("\nRelabelled", nrow(chg), "rows in", basename(MASTER), "\n\n")
rep <- m[, .(genes = .N,
             max_phenotypes = max(as.integer(N_Phenotypes)),
             genes_ge7 = sum(as.integer(N_Phenotypes) >= 7)),
         by = .(condition, layer)]
rep <- merge(rep, chg[, .(relabelled = .N), by = .(condition, layer)],
             by = c("condition", "layer"), all.x = TRUE)
rep[is.na(relabelled), relabelled := 0L]
print(rep[order(condition, layer)])

cat("\nOriginal archived in", ARCH, "\n")
