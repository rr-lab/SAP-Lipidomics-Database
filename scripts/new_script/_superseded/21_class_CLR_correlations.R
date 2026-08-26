# Original CLR-space Pearson class-correlation workflow, rendered as a standalone supplementary figure.
suppressPackageStartupMessages({library(vroom);library(dplyr);library(tidyr);library(stringr);library(ggplot2);library(patchwork);library(purrr)})
root <- '/Users/nirwantandukar/Documents/Github/SoLD_paper'
out <- '/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database/fig/new_figures/SuppFig_Class_CLR_Correlations.png'
classes <- c('DG','DGDG','LPC','LPE','MG','MGDG','PA','PC','PE','PG','PS','SQDG','TG')
read_x <- function(path, condition) {x <- vroom(path,show_col_types=FALSE) |> select(-c(2,3,4)); names(x)[1] <- 'Sample'; x |> mutate(Condition=condition)}
class_clr <- function(d) {
  lipid_cols <- names(d)[grepl('\\(', names(d))]
  labels <- str_extract(lipid_cols, '^[A-Za-z0-9]+(?=\\()')
  mat <- as.matrix(d[, lipid_cols, drop=FALSE]); storage.mode(mat) <- 'numeric'
  cm <- sapply(classes, function(cl) {idx <- which(labels == cl); if(length(idx)) rowSums(mat[,idx,drop=FALSE],na.rm=TRUE) else rep(0,nrow(mat))})
  for(i in seq_len(nrow(cm))) {z <- cm[i,]; pos <- z[z>0]; if(length(pos)) {z[z==0] <- min(pos)*.5; cm[i,] <- z/sum(z)}}
  log(cm) - rowMeans(log(cm))
}
cor_stats <- function(x) {
  r <- cor(x, method='pearson', use='pairwise.complete.obs')
  p <- matrix(0,ncol(x),ncol(x),dimnames=dimnames(r))
  for(i in seq_len(ncol(x)-1)) for(j in (i+1):ncol(x)) {p[i,j] <- p[j,i] <- cor.test(x[,i],x[,j],method='pearson')$p.value}
  list(r=r,p=p)
}
ctl <- read_x(file.path(root,'data/SPATS_fitted/non_normalized_intensities/Final_subset_control_all_lipids_fitted_phenotype_non_normalized.csv'),'Control')
lin <- read_x(file.path(root,'data/SPATS_fitted/non_normalized_intensities/Final_subset_lowinput_all_lipids_fitted_phenotype_non_normalized.csv'),'LowInput')
c_ctl <- cor_stats(class_clr(ctl)); c_lin <- cor_stats(class_clr(lin)); delta <- c_lin$r-c_ctl$r
keys <- outer(classes,classes,function(a,b) ifelse(a<b,paste(a,b,sep='|'),paste(b,a,sep='|')))
upper <- upper.tri(keys)
flips <- unique(keys[upper][c_ctl$p[upper] < .05 & c_lin$p[upper] < .05 & c_ctl$r[upper] * c_lin$r[upper] < 0])
plot_cor <- function(r,p=NULL,show_diag=TRUE) {
  d <- expand_grid(X=classes,Y=classes) |> mutate(ix=match(X,classes),iy=match(Y,classes),R=map2_dbl(X,Y,~r[.x,.y]),P=if(is.null(p)) NA_real_ else map2_dbl(X,Y,~p[.x,.y]),diag=ix==iy,upper=ix<iy,key=map2_chr(X,Y,~ifelse(.x<.y,paste(.x,.y,sep='|'),paste(.y,.x,sep='|'))),sig=case_when(P<.001~'***',P<.01~'**',P<.05~'*',TRUE~''),label=case_when(diag & show_diag~sprintf('%.2f',R),upper & !is.null(p)~sig,!upper & !diag~sprintf('%.2f',R),TRUE~''),flip=!is.null(p) & key %in% flips, X=factor(X,levels=classes),Y=factor(Y,levels=rev(classes)))
  q <- ggplot(d,aes(X,Y,fill=R))+geom_tile(color='white',linewidth=.45)
  if(show_diag) q <- q+geom_tile(data=filter(d,diag),fill='black',color='white',linewidth=.45)
  if(!is.null(p)) q <- q+geom_point(data=filter(d,flip),shape=0,size=7.8,stroke=1.3,color='#C51B7D')
  q+geom_text(aes(label=label,color=ifelse(diag & show_diag,'white','black')),size=2.7)+scale_color_identity()+scale_fill_gradient2(low='#0072B2',mid='white',high='#D55E00',midpoint=0,limits=c(-1,1),name=if(is.null(p))expression(Delta*r) else 'CLR r')+coord_fixed()+labs(x=NULL,y=NULL)+theme_minimal(base_size=12)+theme(axis.text.x=element_text(angle=45,hjust=1,color='black',size=9),axis.text.y=element_text(color='black',size=9),panel.grid=element_blank(),axis.line=element_line(color='black',linewidth=.5),legend.position='right')
}
fig <- plot_cor(c_ctl$r,c_ctl$p) + plot_cor(c_lin$r,c_lin$p) + plot_cor(delta,NULL,FALSE) + plot_layout(ncol=3,guides='collect') + plot_annotation(tag_levels='A',theme=theme(plot.tag=element_text(face='bold',size=16)))
dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
ggsave(out,fig,width=18,height=7,units='in',dpi=300,bg='white')
