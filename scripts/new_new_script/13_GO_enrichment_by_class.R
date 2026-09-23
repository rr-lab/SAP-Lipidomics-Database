# =============================================================================
# GO over-representation per lipid class, with an LD-aware permutation null.
#
# Strata are Condition x Layer x lipid Class, matching the published table.
#
# BACKGROUND. Set GO_UNIVERSE to a genome-wide gene-to-GO file and the test is
# run against the genome, which is what the published table did and what these
# numbers need to mean.
#
# Without it the background falls back to the candidate pool, and the script
# marks the output diagnostic-only. That fallback cannot be used for results.
# One lipid class supplies most of the candidate pool in every stratum -- TG is
# 89% of the CTL individual-species intervals and FA is 83% of the LIN class
# trait intervals -- so for exactly the classes the paper discusses, the
# foreground and the background are nearly the same genes and no term can be
# distinguishable from the pool. The classes that would still test cleanly are
# the small ones nobody is interested in.
#
# TWO TESTS, as before.
#   q_gene  Fisher exact on genes, BH corrected. Treats genes as independent,
#           which they are not, since one locus contributes its whole LD block.
#   q_LD    the test that matters, and the one the Methods describes. The
#           observed loci are replaced by the same number of 250 kb intervals
#           drawn from the genome and matched on gene density, and from each
#           the same number of genes the real interval contributed. So the null
#           has the same number of loci, the same gene count, and genes that
#           are just as physically clustered. A term carried by one big LD
#           block cannot pass, because a single random interval reproduces it
#           just as easily.
#
# Input   <IN_DIR>/{CTL,LIN}_gene_annotation.tsv
#         <IN_DIR>/GO_terms_all.txt
# Output  <OUT_DIR>/Table_GO_enrichment_all.tsv          same columns as before
#         <OUT_DIR>/GO_gene_term_map.tsv                 parsed annotation
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })

IN_DIR   <- Sys.getenv("GO_IN",  "data/gene_annotation_final")
OUT_DIR  <- Sys.getenv("GO_OUT", "table/go_enrichment")
UNIVERSE <- Sys.getenv("GO_UNIVERSE", "")
# The permutation p cannot go below 1/(N_PERM+1). After BH over several hundred
# terms a q below 0.05 needs a raw p near 1e-4, so a few thousand permutations
# make passing arithmetically impossible whatever the data says. Two stages are
# used: a cheap screen over every term, then a deep run only for the terms that
# survived it, which is where the resolution is actually needed. The BH family
# is still every term tested, so the correction is not weakened.
N_PERM   <- as.integer(Sys.getenv("N_PERM", "2000"))        # stage 1, all terms
N_PERM2  <- as.integer(Sys.getenv("N_PERM2", "100000"))     # stage 2, survivors
SCREEN_HITS <- 20L      # a term goes to stage 2 if stage 1 gave it fewer hits
INTERVAL <- 250000
MIN_GENES <- 2                       # a term needs at least this many candidates in the stratum
set.seed(20260923)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CLASSES <- c("DGDG","MGDG","SQDG","GalCer","LPC","LPE","Cer","AEG","SM","FA",
             "TG","DG","MG","PC","PE","PG","PA","PS")

# ---- 1. annotation ----------------------------------------------------------
a <- rbind(fread(file.path(IN_DIR, "CTL_gene_annotation.tsv")),
           fread(file.path(IN_DIR, "LIN_gene_annotation.tsv")))
a[, Layer := ifelse(grepl("individual", Layer), "Individual lipid species", "Class sums and ratios")]

go <- fread(file.path(IN_DIR, "GO_terms_all.txt"), quote = "")
parse_go <- function(dt, col, onto) {
  x <- dt[nzchar(get(col)), .(GeneID, cell = get(col))]
  if (!nrow(x)) return(NULL)
  y <- x[, .(item = trimws(unlist(strsplit(cell, ";", fixed = TRUE)))), by = GeneID]
  y[, GO   := sub("^.*\\((GO:[0-9]+)\\)$", "\\1", item)]
  y[, Term := sub("\\s*\\(GO:[0-9]+\\)$", "", item)]
  y <- y[grepl("^GO:[0-9]+$", GO)]
  unique(y[, .(GeneID, GO, Term, Ontology = onto)])
}
gmap <- rbindlist(list(parse_go(go, "GO_BP", "BP"),
                       parse_go(go, "GO_MF", "MF"),
                       parse_go(go, "GO_CC", "CC")), use.names = TRUE)

# genome-wide annotation for the background, when supplied
if (nzchar(UNIVERSE)) {
  u <- fread(UNIVERSE, quote = "")
  umap <- rbindlist(list(parse_go(u, "GO_BP", "BP"),
                         parse_go(u, "GO_MF", "MF"),
                         parse_go(u, "GO_CC", "CC")), use.names = TRUE)
  message("genome universe: ", uniqueN(umap$GeneID), " genes, ", uniqueN(umap$GO), " terms")
  if (uniqueN(umap$GeneID) < 10000)
    warning("the universe has fewer than 10,000 genes, is it really genome-wide?")
  DIAGNOSTIC <- FALSE
} else {
  umap <- gmap
  DIAGNOSTIC <- TRUE
  message("\n!! No GO_UNIVERSE supplied. Background falls back to the candidate pool.")
  message("!! Output is DIAGNOSTIC ONLY and must not be used as a result.\n")
}
fwrite(gmap, file.path(OUT_DIR, "GO_gene_term_map.tsv"), sep = "\t")
message("GO annotation: ", uniqueN(gmap$GeneID), " genes, ", uniqueN(gmap$GO), " terms")

# ---- 2. trait -> lipid class -------------------------------------------------
norm <- function(x) { x <- sub("^(PC|PE|PG|PS|PA)\\(([^/()]+)/0:0\\)$", "L\\1(\\2)", x)
                      sub("^(LPC|LPE)\\(([^/()]+)/0:0\\)$", "\\1(\\2)", x) }
trait_class <- function(tr) {
  # a ratio names two classes and is counted under both, since the trait is the
  # contrast between them and neither class is the subject on its own
  m <- regmatches(tr, regexec("^Sum_(.+?)_over_(.+?)_log10ratio$", tr))[[1]]
  if (length(m)) return(m[2:3])
  if (tr %in% CLASSES) return(tr)                       # class sum
  p <- sub("\\(.*$", "", norm(tr))                      # individual species
  if (p %in% CLASSES) return(p)
  character(0)                                          # sterols, carotenoids, etc.
}
tc <- unique(a[, .(Trait)])[, .(Class = trait_class(Trait)), by = Trait]
message("traits mapped to a class: ", uniqueN(tc$Trait), " of ", uniqueN(a$Trait),
        "   (", uniqueN(a$Trait) - uniqueN(tc$Trait), " non-lipid-class compounds dropped)")

x <- merge(a, tc, by = "Trait", allow.cartesian = TRUE)
x[, Interval := paste(Lead_Chr, Lead_BP %/% INTERVAL, sep = ":")]

# ---- genome intervals, for the density-matched null --------------------------
gr <- fread(file.path(IN_DIR, "genes.range"),
            col.names = c("chr", "start", "end", "GeneID"))
gr[, Interval := paste(chr, start %/% INTERVAL, sep = ":")]
# genes are carried as integer indices throughout the null, because matching
# gene names on every draw is what makes a permutation loop slow
gene_lvl <- sort(unique(gr$GeneID))
gr[, gidx := match(GeneID, gene_lvl)]
giv <- gr[, .(n_genes = .N, gidx = list(gidx)), by = Interval]
# ten density strata, so a gene-dense observed locus is replaced by a gene-dense
# random one and enrichment cannot come from simply landing in a crowded region
qs <- unique(quantile(giv$n_genes, probs = seq(0, 1, 0.1)))
giv[, dens := if (length(qs) > 2) cut(n_genes, qs, include.lowest = TRUE, labels = FALSE) else 1L]
giv_by_dens <- split(seq_len(nrow(giv)), giv$dens)
message("genome intervals: ", nrow(giv), " of ", INTERVAL/1000, " kb, ",
        length(giv_by_dens), " density strata")

# ---- 3. enrichment -----------------------------------------------------------
res <- list()
for (cond in unique(x$Condition)) {
  bg_genes <- if (DIAGNOSTIC) intersect(unique(x[Condition == cond]$GeneID), gmap$GeneID)
              else unique(umap$GeneID)
  gm   <- gmap[GeneID %in% intersect(unique(x[Condition == cond]$GeneID), gmap$GeneID)]
  bgm  <- umap[GeneID %in% bg_genes]
  bg_n <- length(bg_genes)
  message("\n", cond, ": background = ", bg_n,
          if (DIAGNOSTIC) " annotated candidate genes (DIAGNOSTIC)" else " genes in the genome universe")
  # one gene -> term index for the whole condition, so the null costs a single
  # tabulate per draw instead of rebuilding the map for every class
  glob_key <- unique(paste(bgm$Ontology, bgm$GO))
  n_glob   <- length(glob_key)
  bg2      <- bgm[GeneID %in% gene_lvl]
  gt_glob  <- vector("list", length(gene_lvl))
  tmp      <- split(match(paste(bg2$Ontology, bg2$GO), glob_key), match(bg2$GeneID, gene_lvl))
  gt_glob[as.integer(names(tmp))] <- tmp

  for (lay in unique(x[Condition == cond]$Layer)) {
    d <- x[Condition == cond & Layer == lay & GeneID %in% bg_genes]
    if (!nrow(d)) next
    # interval -> genes, and the observed (interval, class) assignments
    iv_gene <- unique(d[, .(Interval, GeneID)])
    setkey(iv_gene, Interval)

    for (cls in unique(d$Class)) {
      fg <- unique(d[Class == cls]$GeneID)
      if (length(fg) < MIN_GENES) next
      # Intervals is per TERM: how many distinct 250 kb loci carry it. A term
      # sitting in one interval is one LD block, whatever its gene count says.
      giv_of <- unique(d[Class == cls, .(GeneID, Interval)])
      obs <- merge(gm[GeneID %in% fg], giv_of, by = "GeneID", allow.cartesian = TRUE)[
               , .(Genes = uniqueN(GeneID), Intervals = uniqueN(Interval),
                   Overlap_genes = paste(sort(unique(GeneID)), collapse = ";")),
               by = .(Ontology, GO, Term)][Genes >= MIN_GENES]
      if (!nrow(obs)) next
      bgc <- bgm[, .(Background = uniqueN(GeneID)), by = .(Ontology, GO, Term)]
      obs <- merge(obs, bgc, by = c("Ontology","GO","Term"))
      obs[, `:=`(Condition = cond, Layer = lay, Class = cls,
                 Fold = round((Genes / length(fg)) / (Background / bg_n), 2))]
      obs[, p_gene := phyper(Genes - 1, Background, bg_n - Background, length(fg), lower.tail = FALSE)]

      # LD-aware null, density-matched genome intervals
      key <- paste(obs$Ontology, obs$GO)

      # what the observed loci actually contributed: per interval, how many
      # candidate genes of this class, and how gene-dense that interval is
      oc <- d[Class == cls, .(n_cand = uniqueN(GeneID)), by = Interval]
      oc <- merge(oc, giv[, .(Interval, n_genes, dens)], by = "Interval", all.x = TRUE)
      oc[is.na(dens), dens := 1L]

      gi      <- match(key, glob_key)        # this class's terms in the global space
      n_terms <- length(key); obs_n <- obs$Genes

      # the density pools are flattened into one vector with offsets, so all
      # nr interval draws are one vectorised lookup instead of nr list lookups
      pools     <- lapply(oc$dens, function(dd) giv_by_dens[[as.character(dd)]])
      pool_flat <- unlist(pools, use.names = FALSE)
      pool_len  <- lengths(pools)
      pool_off  <- c(0L, cumsum(pool_len)[-length(pool_len)])
      ncand     <- oc$n_cand
      nr        <- nrow(oc)
      draw <- function() {
        js <- pool_flat[pool_off + ceiling(runif(nr) * pool_len)]
        picked <- vector("list", nr)
        for (r in seq_len(nr)) {
          g <- giv$gidx[[ js[r] ]]
          k <- if (ncand[r] < length(g)) ncand[r] else length(g)
          picked[[r]] <- if (k == length(g)) g else g[sample.int(length(g), k)]
        }
        # genes with no GO annotation contribute nothing, exactly as in the
        # observed set, so they stay in the draw rather than being resampled
        idx <- unlist(gt_glob[unique(unlist(picked, use.names = FALSE))], use.names = FALSE)
        if (is.null(idx)) return(integer(n_terms))
        tabulate(idx, nbins = n_glob)[gi]
      }
      hits <- integer(n_terms)
      for (i in seq_len(N_PERM)) hits <- hits + (draw() >= obs_n)
      deep <- which(hits < SCREEN_HITS)
      obs[, p_LD := (hits + 1) / (N_PERM + 1)]
      obs[, N_perm := N_PERM]
      if (length(deep) && N_PERM2 > 0) {
        hits2 <- integer(n_terms)
        for (i in seq_len(N_PERM2)) hits2 <- hits2 + (draw() >= obs_n)
        obs$p_LD[deep]   <- (hits2[deep] + 1) / (N_PERM2 + 1)
        obs$N_perm[deep] <- N_PERM2
      }
      res[[length(res) + 1L]] <- obs
    }
  }
}
out <- rbindlist(res, use.names = TRUE)
out[, q_gene := p.adjust(p_gene, "BH"), by = .(Condition, Layer, Ontology)]
out[, q_LD   := p.adjust(p_LD,   "BH"), by = .(Condition, Layer, Ontology)]
setcolorder(out, c("Ontology","Layer","Condition","Class","Term","GO","Genes","Background",
                   "Intervals","Fold","q_gene","q_LD","Overlap_genes"))
setorder(out, Ontology, Layer, Condition, Class, q_LD, -Fold)
fn <- if (DIAGNOSTIC) "Table_GO_enrichment_DIAGNOSTIC_candidate_background.tsv" else "Table_GO_enrichment_all.tsv"
fwrite(out[, .(Ontology,Layer,Condition,Class,Term,GO,Genes,Background,Intervals,Fold,
               q_gene = signif(q_gene,3), q_LD = signif(q_LD,3), Overlap_genes)],
       file.path(OUT_DIR, fn), sep = "\t")
if (DIAGNOSTIC) cat("\n!! candidate-pool background, diagnostic only, not a result\n")

cat("\n-- tested terms --\n"); print(out[, .(terms = .N), by = .(Condition, Layer, Ontology)])
cat("\n-- passing q_LD < 0.05 --\n"); print(out[q_LD < 0.05, .(terms = .N), by = .(Condition, Layer, Ontology)])
cat("\n-- top by fold among those passing --\n")
print(head(out[q_LD < 0.05][order(-Fold), .(Condition, Layer, Class, Term, Genes, Background, Intervals, Fold, q_LD)], 20),
      row.names = FALSE)
cat("\nSaved:", file.path(OUT_DIR, fn), " (", nrow(out), "rows )\n")
