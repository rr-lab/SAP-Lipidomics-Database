# ==============================================================================
# 45_two_environment_heritability.R
#
#   Rscript scripts/new_script/45_two_environment_heritability.R
#
# Heritability using both trials at once, rather than one trial at a time.
#
# With the same genotypes measured in two trials, each genotype has a mean
# across trials and a difference between them. Standardising within trial first,
# so the two are on one scale, these two derived traits separate what the
# single-trial fits cannot:
#
#   mean   = (z_CTL + z_LIN) / 2   genotype main effect, the part of a genotype's
#                                  lipid phenotype that is the same in both trials
#   diff   =  z_LIN - z_CTL        the response, how much a genotype changes
#                                  between trials; its heritability is genotype x
#                                  environment variance in single-kernel form
#
# From the four variance components the genetic correlation between trials
# follows, since Vg(mean) = (Vg1 + Vg2 + 2C)/4 and Vg(diff) = Vg1 + Vg2 - 2C,
# so C = Vg(mean) - Vg(diff)/4 and rg = C / sqrt(Vg1 * Vg2). An rg near one means
# the same genotypes rank the same way in both trials; near zero means the
# trials are measuring different genetic variation.
#
# WHAT THIS CANNOT SAY. The two trials differ in year, field, planting date and
# analytical batch as well as in nutrient supply, and neither trial is
# replicated. "Environment" here is the whole trial contrast, not the nutrient
# treatment, so a heritable response means genotypes differ in how they change
# between these two trials, not that they differ in response to low input.
#
# Input   data/SPATS_fitted/non_normalized_intensities/Final_subset_{control,lowinput}_*.csv
#         data/kinship/sap_grm.rel, sap_grm.rel.id
# Output  table/new_table/SuppTable_TwoEnvironment_Heritability.csv
#
# Run from the repository root.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(dplyr); library(readr); library(stringr) })

root <- Sys.getenv('SOLD_DB_REPO','.'); kdir <- file.path(root,'data','kinship')
tabdir <- file.path(root,'table','new_table'); dir.create(tabdir, recursive=TRUE, showWarnings=FALSE)
GRID <- seq(0,1,by=0.002)

load_k <- function(){
  K <- as.matrix(fread(file.path(kdir,'sap_grm.rel'), header=FALSE, data.table=FALSE))
  idf <- fread(file.path(kdir,'sap_grm.rel.id'), header=FALSE, data.table=FALSE, sep='\t')
  idf <- idf[!grepl('^#', idf[[1]]), , drop=FALSE]
  ids <- as.character(idf[[2]]); dimnames(K) <- list(ids,ids); list(K=K, ids=ids) }
read_spats <- function(cond){
  f <- file.path(root,'data','SPATS_fitted','non_normalized_intensities',
        sprintf('Final_subset_%s_all_lipids_fitted_phenotype_non_normalized.csv',
                if(cond=='CTL') 'control' else 'lowinput'))
  x <- fread(f, data.table=FALSE, check.names=FALSE)
  lip <- setdiff(names(x), c('LineRaw','PlotID','row','col')); lip <- lip[grepl('\\(', lip)]
  m <- as.matrix(x[,lip,drop=FALSE]); m[!is.finite(m)] <- NA_real_
  rownames(m) <- as.character(x$LineRaw); list(mat=m, lipids=lip) }

prof <- function(yt, xt, d){
  n <- length(yt)
  ll <- vapply(GRID, function(h){ v <- h*d+(1-h); if(any(v<=1e-12)) return(-Inf)
    w <- 1/v; XtWX <- sum(w*xt*xt); if(!is.finite(XtWX)||XtWX<=0) return(-Inf)
    b <- sum(w*xt*yt)/XtWX; r <- yt-xt*b; s2 <- sum(w*r*r)/(n-1)
    if(!is.finite(s2)||s2<=0) return(-Inf)
    -0.5*(sum(log(v))+(n-1)*log(s2)+log(XtWX)) }, numeric(1))
  i <- which.max(ll); keep <- which(ll >= ll[i]-qchisq(0.95,1)/2)
  c(h2=GRID[i], lo=GRID[min(keep)], hi=GRID[max(keep)]) }

C <- read_spats('CTL'); L <- read_spats('LIN'); KIN <- load_k()
sp    <- intersect(C$lipids, L$lipids)
lines <- Reduce(intersect, list(rownames(C$mat), rownames(L$mat), KIN$ids))
lines <- lines[!duplicated(lines)]
message(sprintf('shared species %d, shared genotypes %d', length(sp), length(lines)))

Ksub <- KIN$K[lines, lines]; e <- eigen(Ksub, symmetric=TRUE)
d <- pmax(e$values,0); V <- e$vectors
xt <- as.numeric(crossprod(V, rep(1,length(lines))))
z  <- function(v) (v-mean(v))/sd(v)

res <- rbindlist(lapply(sp, function(s){
  a <- log10(C$mat[lines,s]+1); b <- log10(L$mat[lines,s]+1)
  if(!all(is.finite(c(a,b))) || sd(a)==0 || sd(b)==0) return(NULL)
  za <- z(a); zb <- z(b)
  mn <- prof(as.numeric(crossprod(V,(za+zb)/2)), xt, d)
  df <- prof(as.numeric(crossprod(V, zb-za)), xt, d)
  data.table(Species=s, Class=str_match(s,'^([A-Za-z]+)\\(')[,2],
             h2_mean=mn[['h2']], lo_mean=mn[['lo']], hi_mean=mn[['hi']],
             h2_diff=df[['h2']], lo_diff=df[['lo']], hi_diff=df[['hi']]) }))
res[, `:=`(mean_nonzero = lo_mean>0, diff_nonzero = lo_diff>0)]
write_csv(res, file.path(tabdir,'SuppTable_TwoEnvironment_Heritability.csv'))

cat('\n== across-trial decomposition,', nrow(res), 'species ==\n')
cat(sprintf('  genotype main effect (trial mean) : median h2 %.3f, non-zero in %d\n',
            median(res$h2_mean), sum(res$mean_nonzero)))
cat(sprintf('  response       (trial difference) : median h2 %.3f, non-zero in %d\n',
            median(res$h2_diff), sum(res$diff_nonzero)))
cat('\n== by class ==\n')
print(as.data.frame(res[, .(n=.N, h2_mean=round(median(h2_mean),3), mean_nz=sum(mean_nonzero),
       h2_diff=round(median(h2_diff),3), diff_nz=sum(diff_nonzero)), by=Class][order(-n)]), row.names=FALSE)
