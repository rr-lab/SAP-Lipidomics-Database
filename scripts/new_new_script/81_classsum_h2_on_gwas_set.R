# ==============================================================================
# Class-sum heritability refitted on the GWAS genotype sets.
#
#   Rscript scripts/new_new_script/81_classsum_h2_on_gwas_set.R
#
# WHY. Table 3 reports class-sum heritabilities fitted on every accession of a
# trial that is present in the relationship matrix: 389 under CTL and 357 under
# LIN. The matrix holds 401 accessions, more than the 357 genotypes that survive
# into GWAS, so the CTL estimate rests on 38 accessions the GWAS never sees. This
# script asks what Table 3 would say if CTL were restricted to the 351 accessions
# common to both trials and to the matrix, the set the rest of the manuscript
# quotes. LIN is already 357 either way and is refitted only as a control.
#
# NOTHING IS OVERWRITTEN. 21_heritability_compare_grm.R and the shipped
# SuppTable_S26 are left alone; this writes its own file and prints a comparison.
#
# The estimator, the class list and the scaling are copied verbatim from
# 21_heritability_compare_grm.R so the only difference between the two runs is
# which accessions are in the fit. The full-set run is a reproduction check: if
# it does not return the shipped values the restricted run means nothing.
#
# One detail is not in 21_heritability_compare_grm.R and had to be recovered.
# No script in the repository rebuilds SuppTable_S26, so the phenotype it was
# fitted on was found by sweeping the construction against the shipped numbers:
# genotype means per species, class sum, log10(x + 1), then scale(). That
# reproduces all 34 shipped values to the grid step and gives the same n and
# N_species. On the raw scale only 17 of 34 agree, so the log is not optional.
#
# Input   data/SPATS_fitted/non_normalized_intensities/Final_subset_*.csv
#         data/kinship/sap_grm.rel and .rel.id
#         table/supp/SuppTable_S26_Heritability_Class_Sums.csv   (for the check)
# Output  table/new_table/classsum_h2_full_vs_gwasset.csv
# ==============================================================================
suppressPackageStartupMessages({ library(data.table) })

ROOT  <- Sys.getenv("SOLD_REPO", ".")
SPATS <- file.path(ROOT, "data/SPATS_fitted/non_normalized_intensities")
GRM   <- file.path(ROOT, "data/kinship/sap_grm.rel")
SHIP  <- file.path(ROOT, "table/supp/SuppTable_S26_Heritability_Class_Sums.csv")
OUT   <- file.path(ROOT, "table/new_table/classsum_h2_full_vs_gwasset.csv")
GRID  <- seq(0, 0.999, by = 0.001)
CLASSES <- c("DGDG","MGDG","SQDG","GalCer","LPC","LPE","Cer","AEG","SM","FA",
             "TG","DG","MG","PC","PE","PG","PA","PS")
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)

K <- as.matrix(fread(GRM, header = FALSE, data.table = FALSE))
idf <- fread(paste0(GRM, ".id"), header = FALSE, data.table = FALSE)
idf <- idf[!startsWith(trimws(as.character(idf[[1]])), "#"), , drop = FALSE]
ids <- as.character(idf[[ncol(idf)]]); stopifnot(nrow(K) == length(ids))
K <- (K + t(K)) / 2; dimnames(K) <- list(ids, ids)
message("relationship matrix: ", length(ids), " accessions")

prof_h2 <- function(y, K) {
  n <- length(y); p <- 1
  e <- eigen(K, symmetric = TRUE); d <- pmax(e$values, 0)
  yt <- as.numeric(crossprod(e$vectors, y))
  xt <- as.numeric(crossprod(e$vectors, rep(1, n)))
  ll <- vapply(GRID, function(h) {
    v <- h * d + (1 - h); if (any(v <= 1e-12)) return(-Inf)
    w <- 1 / v; XtWX <- sum(w * xt * xt); if (!is.finite(XtWX) || XtWX <= 0) return(-Inf)
    b <- sum(w * xt * yt) / XtWX; r <- yt - xt * b
    s2 <- sum(w * r * r) / (n - p); if (!is.finite(s2) || s2 <= 0) return(-Inf)
    -0.5 * (sum(log(v)) + (n - p) * log(s2) + log(XtWX))
  }, numeric(1))
  i <- which.max(ll); keep <- which(ll >= ll[i] - qchisq(0.95, 1) / 2)
  c(h2 = GRID[i], lo = GRID[min(keep)], hi = GRID[max(keep)])
}

norm <- function(x) { x <- sub("^(PC|PE|PG|PS|PA)\\(([^/()]+)/0:0\\)$", "L\\1(\\2)", x)
                      sub("^(LPC|LPE)\\(([^/()]+)/0:0\\)$", "\\1(\\2)", x) }

read_sums <- function(cond) {
  f <- file.path(SPATS, sprintf("Final_subset_%s_all_lipids_fitted_phenotype_non_normalized.csv",
                                if (cond == "CTL") "control" else "lowinput"))
  x <- fread(f, data.table = FALSE, check.names = FALSE)
  lip <- setdiff(names(x), c("LineRaw","PlotID","row","col")); lip <- lip[grepl("\\(", lip)]
  dt <- as.data.table(x)[, c("LineRaw", lip), with = FALSE]
  ag <- dt[, lapply(.SD, function(v) mean(v[is.finite(v)])), by = LineRaw, .SDcols = lip]
  m  <- as.matrix(ag[, ..lip]); m[!is.finite(m)] <- NA_real_
  rownames(m) <- as.character(ag$LineRaw)
  cls <- sub("\\(.*$", "", norm(colnames(m)))
  n_sp <- sapply(CLASSES, function(cl) sum(cls == cl))
  sums <- sapply(CLASSES, function(cl) {
    j <- which(cls == cl); if (!length(j)) return(rep(NA_real_, nrow(m)))
    rowSums(m[, j, drop = FALSE], na.rm = TRUE) })
  colnames(sums) <- CLASSES
  list(sums = sums, n_species = n_sp)
}

ctl <- read_sums("CTL"); lin <- read_sums("LIN")
paired <- intersect(intersect(rownames(ctl$sums), rownames(lin$sums)), ids)
message("CTL in matrix: ", length(intersect(rownames(ctl$sums), ids)),
        " | LIN in matrix: ", length(intersect(rownames(lin$sums), ids)),
        " | common to both: ", length(paired))

fit <- function(d, cond, keep, label) {
  M <- d$sums
  rbindlist(lapply(CLASSES, function(cl) {
    if (d$n_species[[cl]] == 0) return(NULL)
    ln <- intersect(rownames(M), keep)
    y <- M[ln, cl]; ok <- is.finite(y); ln <- ln[ok]; y <- y[ok]
    if (length(ln) < 50 || sd(y) == 0) return(NULL)
    y <- as.numeric(scale(log10(y + 1)))
    h <- prof_h2(y, K[ln, ln])
    data.table(Set = label, Condition = cond, Class = cl, N_species = d$n_species[[cl]],
               n = length(ln), h2 = round(h[["h2"]], 3),
               lo = round(h[["lo"]], 3), hi = round(h[["hi"]], 3),
               NonZero = h[["lo"]] > 0)
  }), use.names = TRUE)
}

res <- rbindlist(list(
  fit(ctl, "CTL", ids,    "full"),       # 389, reproduces Table 3
  fit(lin, "LIN", ids,    "full"),       # 357
  fit(ctl, "CTL", paired, "gwas_set"),   # 351
  fit(lin, "LIN", paired, "gwas_set")    # 351, for symmetry
), use.names = TRUE)
fwrite(res, OUT)

cat("\n== reproduction check against the shipped SuppTable_S26 ==\n")
if (file.exists(SHIP)) {
  s <- fread(SHIP)[, .(Condition, Class, N_ship = N_species, n_ship = n, h2_ship = h2,
                       lo_ship = CI_lo, hi_ship = CI_hi)]
  m <- merge(res[Set == "full"], s, by = c("Condition","Class"))
  m[, same := abs(h2 - h2_ship) < 0.0015 & n == n_ship & N_species == N_ship]
  cat(sprintf("  %d of %d class sums reproduce (h2 within 0.0015 and same n)\n",
              sum(m$same), nrow(m)))
  if (any(!m$same)) print(m[!(same), .(Condition, Class, N_species, N_ship, n, n_ship, h2, h2_ship)])
} else cat("  shipped table not found, skipping\n")

cat("\n== CTL: 389 accessions (Table 3) vs 351 (GWAS set) ==\n")
w <- dcast(res[Condition == "CTL"], Class + N_species ~ Set, value.var = c("n","h2","lo","hi","NonZero"))
w[, shift := h2_gwas_set - h2_full]
print(w[order(-abs(shift)),
        .(Class, N = N_species, n_389 = n_full, h2_389 = h2_full, lo_389 = lo_full, hi_389 = hi_full,
          n_351 = n_gwas_set, h2_351 = h2_gwas_set, lo_351 = lo_gwas_set, hi_351 = hi_gwas_set,
          shift, nz_389 = NonZero_full, nz_351 = NonZero_gwas_set)], row.names = FALSE)
cat("\n  classes whose non-zero call flips:",
    paste(w[NonZero_full != NonZero_gwas_set, Class], collapse = ", "), "\n")
cat("\nSaved:", OUT, "\n")
