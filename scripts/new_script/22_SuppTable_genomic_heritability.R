suppressPackageStartupMessages({library(data.table);library(dplyr);library(tidyr);library(stringr);library(readr);library(rrBLUP);library(ggplot2);library(patchwork);library(jsonlite)})
root <- '/Users/nirwantandukar/Documents/Github/SoLD_paper'
out <- '/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database'
tabdir <- file.path(out,'table/new_table'); figdir <- file.path(out,'fig/new_figures')
dir.create(tabdir, recursive=TRUE, showWarnings=FALSE); dir.create(figdir, recursive=TRUE, showWarnings=FALSE)
plot_theme <- theme_minimal(base_size=13)+theme(axis.text=element_text(colour='black',size=9),axis.title=element_text(face='bold',size=11),axis.line=element_line(colour='black'),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),legend.title=element_blank())
core <- c('AEG','Cer','SM','GalCer','FA','MG','DG','TG','MGDG','DGDG','SQDG','PA','PC','PE','PG','PS','LPC','LPE','DGTS')
cols <- c(AEG='#D55E00',Cer='#CC79A7',SM='#7B61A8',GalCer='#8C6BB1',FA='#E6AB02',MG='#8A2BE2',DG='#B8860B',TG='#D53E9D',MGDG='#129DA5',DGDG='#B59A00',SQDG='#A15DD7',PA='#1B7837',PC='#0099C7',PE='#2CA25F',PG='#78C679',PS='#C2E699',LPC='#66C2A4',LPE='#2CA25F',DGTS='#7FCDBB','Fatty acid'='#E9C46A',Sphingolipid='#6D597A',Terpenoid='#E76F51',Sterol='#2A9D8F','Ether lipid'='#457B9D',Prenol='#F4A261','Betaine lipid'='#9B2226')
line_id <- function(x) sub('^S1_','',sub('_Run[0-9]+$','',x))
class_id <- function(lipid,broad) {p <- str_match(lipid,'^([A-Za-z]+)\\(')[,2]; z <- ifelse(!is.na(p)&p%in%core,p,broad); z[z=='Glycerolipid'] <- 'Other glycerolipid'; z[z=='Glycerophospholipid'] <- 'Other glycerophospholipid'; z[is.na(z)|z==''] <- 'Other annotated'; z}
read_final <- function(cond) {
 f <- file.path(root,'data/raw_lipid_intensities',if(cond=='CTL')'A_final_lipids.csv' else 'B_final_lipids.csv'); x <- fread(f,data.table=FALSE,check.names=FALSE)
 smp <- names(x)[grep('^S1_PI',names(x))]; list(meta=x[,c('Compound_Name','Class','ambiguous')],mat=as.matrix(x[,smp,drop=FALSE]),samples=smp)
}
read_serrf <- function(cond) {
 f <- file.path(root,'data/SEERF_normalized_intensities',if(cond=='CTL')'A_normalized_SERRF.csv' else 'B_normalized_SERRF.csv'); x <- fread(f,data.table=FALSE,check.names=FALSE); names(x)[1] <- 'Scan'; x$Scan <- as.character(x$Scan)
 smp <- names(x)[grep('^S1_PI',names(x))]; list(x=x[,c('Scan',smp),drop=FALSE],samples=smp)
}
map_serrf <- function(cond) {
 f <- read_final(cond); s <- read_serrf(cond); common <- intersect(f$samples,s$samples)
 a <- log1p(pmax(as.matrix(f$mat[,match(common,f$samples),drop=FALSE]),0)); b <- log1p(pmax(as.matrix(s$x[,common,drop=FALSE]),0))
 a <- t(scale(t(a))); b <- t(scale(t(b))); a[!is.finite(a)] <- 0; b[!is.finite(b)] <- 0
 cc <- tcrossprod(a,b)/max(1,ncol(a)-1); ix <- max.col(cc); r <- cc[cbind(seq_len(nrow(cc)),ix)]
 audit <- tibble(Condition=cond,Lipid=f$meta$Compound_Name,BroadClass=f$meta$Class,Ambiguous=f$meta$ambiguous,Scan=as.character(s$x$Scan[ix]),ProfileCorrelation=as.numeric(r),MappingStatus=ifelse(r>=.95,'retained_high_confidence','excluded_low_confidence'))
 keep <- audit |> filter(MappingStatus=='retained_high_confidence') |> arrange(desc(ProfileCorrelation),Ambiguous) |> group_by(Scan) |> slice(1) |> ungroup()
 list(s=s,keep=keep,audit=audit)
}
# Genomic relationship matrix built from 60,000 LD-thinned HapMap SNPs (one SNP per 10 kb bin).
# It is intentionally unweighted and is used only for genomic h2 estimation.
load_grm <- function() {
  gd <- file.path(out, 'data/genomic_grm_serrf')
  ids <- read_tsv(file.path(gd, 'samples.tsv'), show_col_types=FALSE)$sample_id
  K <- as.matrix(fread(file.path(gd, 'K_weighted.tsv'), header=FALSE, data.table=FALSE))
  list(K=(K+t(K))/2, ids=ids, nmarkers=as.integer(jsonlite::fromJSON(file.path(gd, 'grm_summary.json'))$counts$snps_kept))
}
fit_h2 <- function(y,ids,grm) {
  i <- match(ids,grm$ids); ok <- !is.na(i)&is.finite(y); y <- y[ok]; i <- i[ok]
  if(length(y)<30||sd(y)==0)return(tibble(n_genotypes=length(y),n_markers=grm$nmarkers,Vu=NA_real_,Ve=NA_real_,h2=NA_real_,model_status='insufficient_variation'))
  q <- tryCatch(mixed.solve(y,K=grm$K[i,i,drop=FALSE]),error=function(e)e)
  if(inherits(q,'error'))return(tibble(n_genotypes=length(y),n_markers=grm$nmarkers,Vu=NA_real_,Ve=NA_real_,h2=NA_real_,model_status=conditionMessage(q)))
  tibble(n_genotypes=length(y),n_markers=grm$nmarkers,Vu=q$Vu,Ve=q$Ve,h2=q$Vu/(q$Vu+q$Ve),model_status='ok')
}
S1 <- read_serrf('CTL'); S2 <- read_serrf('LIN'); grm <- load_grm(); message('GRM: ',length(grm$ids),' genotypes, ',grm$nmarkers,' markers')
estimate <- function(cond) {
 o <- map_serrf(cond); smp <- o$s$samples; x <- o$s$x |> filter(Scan%in%o$keep$Scan) |> pivot_longer(all_of(smp),names_to='Sample',values_to='Intensity') |> inner_join(o$keep |> select(Scan,Lipid,BroadClass),by='Scan') |> mutate(Line=line_id(Sample),Intensity=as.numeric(Intensity),Class=class_id(Lipid,BroadClass)) |> filter(Line%in%grm$ids) |> group_by(Lipid,BroadClass,Class,Line) |> summarise(Intensity=sum(Intensity,na.rm=TRUE),.groups='drop')
 ind <- x |> group_by(Lipid,BroadClass,Class) |> group_modify(~fit_h2(log10(.x$Intensity+1),.x$Line,grm)) |> ungroup() |> mutate(Condition=cond,FeatureType='Individual lipid',.before=1)
 sum <- x |> group_by(Line,Class) |> summarise(Intensity=sum(Intensity,na.rm=TRUE),.groups='drop') |> group_by(Class) |> group_modify(~fit_h2(log10(.x$Intensity+1),.x$Line,grm)) |> ungroup() |> mutate(Condition=cond,FeatureType='Lipid class sum',Lipid=NA_character_,BroadClass=NA_character_,.before=1)
 list(audit=o$audit,ind=ind,sum=sum)
}
ctl <- estimate('CTL'); lin <- estimate('LIN'); audit <- bind_rows(ctl$audit,lin$audit); ind <- bind_rows(ctl$ind,lin$ind); sums <- bind_rows(ctl$sum,lin$sum); all <- bind_rows(ind,sums)
write_csv(audit,file.path(tabdir,'SuppTable_SERRF_Scan_to_Lipid_Mapping_Audit.csv'));write_csv(ind,file.path(tabdir,'SuppTable_Genomic_Heritability_Individual_Lipids.csv'));write_csv(sums,file.path(tabdir,'SuppTable_Genomic_Heritability_Lipid_Class_Sums.csv'));write_csv(all,file.path(tabdir,'SuppTable_Genomic_Heritability_All.csv'))
pind <- function(x,t) {d <- x |> filter(is.finite(h2),model_status=='ok') |> slice_max(h2,n=20,with_ties=FALSE) |> arrange(h2) |> mutate(Lab=factor(str_wrap(Lipid,35),levels=str_wrap(Lipid,35)));ggplot(d,aes(h2,Lab,colour=Class))+geom_segment(aes(x=0,xend=h2,y=Lab,yend=Lab),colour='grey70')+geom_point(size=2.8)+scale_colour_manual(values=cols,na.value='grey50')+labs(x=expression('Genomic heritability ('*h^2*')'),y=NULL,tag=t)+plot_theme+theme(axis.text.y=element_text(size=7),legend.position='right')}
psum <- function(x,t) {d <- x |> filter(is.finite(h2),model_status=='ok') |> arrange(h2) |> mutate(Class=factor(Class,levels=Class));ggplot(d,aes(h2,Class,fill=Class))+geom_col(colour='black',linewidth=.2)+geom_text(aes(label=sprintf('%.2f',h2)),hjust=-.2,size=3)+scale_fill_manual(values=cols,na.value='grey50')+scale_x_continuous(expand=expansion(mult=c(0,.15)))+labs(x=expression('Genomic heritability ('*h^2*')'),y=NULL,tag=t)+plot_theme+theme(legend.position='none')}
fig <- (pind(filter(ind,Condition=='CTL'),'A')|pind(filter(ind,Condition=='LIN'),'B'))/(psum(filter(sums,Condition=='CTL'),'C')|psum(filter(sums,Condition=='LIN'),'D'))+plot_layout(heights=c(1.15,.85))+plot_annotation(theme=theme(plot.tag=element_text(face='bold',size=16)))
ggsave(file.path(figdir,'SuppFig_Lipid_Heritability.png'),fig,width=15,height=13,dpi=350,bg='white')
