# =============================================================================
# Does excluding the high-heterozygosity markers change heritability?
#
# Runs the same profile-likelihood h2 used in 22f_heritability_species_paired.R
# twice, once under the GRM built on all markers and once under the GRM built
# after the heterozygosity filter, on the same species and the same genotypes.
# The estimator is copied from that script unchanged so the two runs differ only
# in the relationship matrix.
#
# Phenotypes are the SpATS fitted intensities, which is what the heritability
# scripts use. They are not the BLUPs: fitting h2 on BLUPs would be circular,
# since a BLUP is already shrunk by a genetic model.
#
# Both individual species and lipid-class sums are run. Class sums are included
# because the shipped SuppTable_S26 reports them and nothing in the repository
# rebuilds it. Class ratios are NOT run: h2 of a ratio is not a property of a
# lipid class and none of the manuscript's heritability claims rest on one.
#
# Input   data/SPATS_fitted/non_normalized_intensities/Final_subset_*.csv
#         <GRM_DIR>/sap_grm_{allmarkers,hetfiltered}.rel and .rel.id
# Output  <OUT>/heritability_grm_comparison.csv
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })

ROOT    <- Sys.getenv("SOLD_ROOT", ".")
GRM_DIR <- Sys.getenv("GRM_DIR", "/rsstu/users/r/rrellan/DOE_CAREER/SAP/results/spats_corrected/BLUP/kinship_hetfiltered")
OUT     <- Sys.getenv("OUT_DIR", GRM_DIR)
SPATS   <- file.path(ROOT, "data/SPATS_fitted/non_normalized_intensities")
GRID    <- seq(0, 0.999, by = 0.001)
CLASSES <- c("DGDG","MGDG","SQDG","GalCer","LPC","LPE","Cer","AEG","SM","FA",
             "TG","DG","MG","PC","PE","PG","PA","PS")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

load_k <- function(tag) {
  f <- file.path(GRM_DIR, sprintf("sap_grm_%s.rel", tag))
  if (!file.exists(f)) stop("missing GRM, run 20_build_filtered_grm.sh first: ", f)
  K   <- as.matrix(fread(f, header = FALSE, data.table = FALSE))
  idf <- fread(paste0(f, ".id"), header = FALSE, data.table = FALSE)
  idf <- idf[!startsWith(trimws(as.character(idf[[1]])), "#"), , drop = FALSE]
  ids <- as.character(idf[[ncol(idf)]]); stopifnot(nrow(K) == length(ids))
  K <- (K + t(K)) / 2; dimnames(K) <- list(ids, ids)
  list(K = K, ids = ids)
}

# profile likelihood, copied from 22f so the two runs are comparable
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

read_spats <- function(cond) {
  f <- file.path(SPATS, sprintf("Final_subset_%s_all_lipids_fitted_phenotype_non_normalized.csv",
                                if (cond == "CTL") "control" else "lowinput"))
  x <- fread(f, data.table = FALSE, check.names = FALSE)
  lip <- setdiff(names(x), c("LineRaw","PlotID","row","col"))
  lip <- lip[grepl("\\(", lip)]
  # genotype means, so a genotype with more plots is not weighted more heavily
  dt <- as.data.table(x)[, c("LineRaw", lip), with = FALSE]
  ag <- dt[, lapply(.SD, function(v) mean(v[is.finite(v)])), by = LineRaw, .SDcols = lip]
  m  <- as.matrix(ag[, ..lip]); m[!is.finite(m)] <- NA_real_
  rownames(m) <- as.character(ag$LineRaw)
  cls <- sub("\\(.*$", "", norm(colnames(m)))
  sums <- sapply(CLASSES, function(cl) {
    j <- which(cls == cl); if (!length(j)) return(rep(NA_real_, nrow(m)))
    rowSums(m[, j, drop = FALSE], na.rm = TRUE) })
  colnames(sums) <- paste0("Sum_", CLASSES)
  list(species = m, sums = sums)
}

Kall <- load_k("allmarkers"); Kflt <- load_k("hetfiltered")
message("GRM individuals: all ", length(Kall$ids), " | filtered ", length(Kflt$ids))
message("max |K_all - K_filtered| over the shared block: ",
        signif(max(abs(Kall$K[Kflt$ids, Kflt$ids] - Kflt$K[Kflt$ids, Kflt$ids])), 3))
message("correlation of off-diagonal entries: ",
        signif(cor(Kall$K[lower.tri(Kall$K)], Kflt$K[Kflt$ids, Kflt$ids][lower.tri(Kflt$K)]), 6))

res <- rbindlist(lapply(c("CTL","LIN"), function(cond) {
  d <- read_spats(cond)
  rbindlist(lapply(c("species","class_sum"), function(kind) {
    M <- if (kind == "species") d$species else d$sums
    rbindlist(lapply(colnames(M), function(tr) {
      rbindlist(lapply(c(allmarkers = "allmarkers", hetfiltered = "hetfiltered"), function(tag) {
        KK <- if (tag == "allmarkers") Kall else Kflt
        ln <- intersect(rownames(M), KK$ids)
        y  <- M[ln, tr]; ok <- is.finite(y); ln <- ln[ok]; y <- y[ok]
        if (length(ln) < 50 || sd(y) == 0) return(NULL)
        y <- as.numeric(scale(y))
        h <- prof_h2(y, KK$K[ln, ln])
        data.table(Condition = cond, Kind = kind, Trait = tr, GRM = tag,
                   n = length(ln), h2 = h[["h2"]], lo = h[["lo"]], hi = h[["hi"]])
      }), use.names = TRUE)
    }), use.names = TRUE)
  }), use.names = TRUE)
}), use.names = TRUE)

w <- dcast(res, Condition + Kind + Trait + n ~ GRM, value.var = c("h2","lo","hi"))
w[, delta := h2_hetfiltered - h2_allmarkers]
fwrite(w, file.path(OUT, "heritability_grm_comparison.csv"))

cat("\n-- does the filter change h2? --\n")
print(w[, .(traits = .N,
            median_h2_all      = round(median(h2_allmarkers), 4),
            median_h2_filtered = round(median(h2_hetfiltered), 4),
            median_abs_delta   = round(median(abs(delta)), 5),
            max_abs_delta      = round(max(abs(delta)), 4),
            n_flipped_at_0.05  = sum((h2_allmarkers > 0.05) != (h2_hetfiltered > 0.05))),
          by = .(Condition, Kind)])
cat("\n-- largest movers --\n")
print(head(w[order(-abs(delta)), .(Condition, Kind, Trait, n, h2_allmarkers, h2_hetfiltered, delta)], 12),
      row.names = FALSE)
cat("\nSaved:", file.path(OUT, "heritability_grm_comparison.csv"), "\n")
