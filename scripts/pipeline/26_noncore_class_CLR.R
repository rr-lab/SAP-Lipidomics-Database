# CLR contrasts for lipid classes outside the core class-ratio/LINEX framework.
suppressPackageStartupMessages({library(vroom);library(dplyr);library(tidyr);library(ggplot2);library(patchwork);library(scales);library(grid)})
root <- '/Users/nirwantandukar/Documents/Github/SoLD_paper'
out <- '/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database/fig/new_figures/SuppFig_Noncore_Lipid_Class_Composition_CLR.png'
core <- c('Glycerolipid','Glycerophospholipid')
# Colors distinguish the non-core lipid superclasses without reusing the focused-class palette.
class_cols <- c('Terpenoid'='#E76F51','Fatty acid'='#E9C46A','Sphingolipid'='#6D597A','Ether lipid'='#457B9D','Sterol'='#2A9D8F','Prenol'='#F4A261','Steroid'='#A8DADC','Tetrapyrrole'='#8AB17D','Flavonoid'='#CDB4DB','Betaine lipid'='#9B2226')
read_matrix <- function(path, condition) { x <- vroom(path,show_col_types=FALSE); x <- x |> select(-c(2,3,4)); names(x)[1] <- 'Sample'; x |> mutate(Condition=condition) }
ctl <- read_matrix(file.path(root,'data/SPATS_fitted/non_normalized_intensities/Final_subset_control_all_lipids_fitted_phenotype_non_normalized.csv'),'Control')
lin <- read_matrix(file.path(root,'data/SPATS_fitted/non_normalized_intensities/Final_subset_lowinput_all_lipids_fitted_phenotype_non_normalized.csv'),'LowInput')
classes <- vroom(file.path(root,'data/lipid_class/final_lipid_classes.csv'),show_col_types=FALSE) |> transmute(Lipid=Lipids,Class)
make_class <- function(x) {
  lipid_cols <- intersect(names(x), classes$Lipid)
  x |> select(Sample,Condition,all_of(lipid_cols)) |> pivot_longer(-c(Sample,Condition),names_to='Lipid',values_to='Intensity') |> left_join(classes,by='Lipid') |> filter(!is.na(Class),Class %in% names(class_cols)) |> group_by(Sample,Condition,Class) |> summarise(class_sum=sum(as.numeric(Intensity),na.rm=TRUE),.groups='drop')
}
# CLR is calculated across all annotated lipid classes, then non-core contrasts are displayed.
all_sums <- bind_rows(ctl,lin) |> select(Sample,Condition,any_of(classes$Lipid)) |> pivot_longer(-c(Sample,Condition),names_to='Lipid',values_to='Intensity') |> left_join(classes,by='Lipid') |> filter(!is.na(Class)) |> group_by(Sample,Condition,Class) |> summarise(value=sum(as.numeric(Intensity),na.rm=TRUE),.groups='drop') |> pivot_wider(names_from=Class,values_from=value,values_fill=0)
mat <- as.matrix(all_sums[,setdiff(names(all_sums),c('Sample','Condition'))]); for(i in seq_len(nrow(mat))){z<-mat[i,]; z[z<0|!is.finite(z)]<-0; if(any(z==0)){p<-z[z>0];z[z==0]<-min(p)*.5};mat[i,]<-z/sum(z)}; clr<-log(mat)-rowMeans(log(mat)); keep<-intersect(names(class_cols),colnames(clr)); idx_c<-which(all_sums$Condition=='Control');idx_l<-which(all_sums$Condition=='LowInput'); eff<-colMeans(clr[idx_l,keep,drop=FALSE])-colMeans(clr[idx_c,keep,drop=FALSE]);set.seed(2026);boot<-replicate(500,{c0<-sample(idx_c,length(idx_c),TRUE);l0<-sample(idx_l,length(idx_l),TRUE);colMeans(clr[l0,keep,drop=FALSE])-colMeans(clr[c0,keep,drop=FALSE])}); ct<-tibble(Class=keep,Effect=eff,CI_Low=apply(boot,1,quantile,.025),CI_High=apply(boot,1,quantile,.975)) |> mutate(Class=factor(Class,levels=rev(Class[order(abs(Effect),decreasing=TRUE)])))
p_a <- ggplot(ct,aes(Class,Effect,color=Class))+geom_hline(yintercept=0,linetype='dashed',color='grey45')+geom_errorbar(aes(ymin=CI_Low,ymax=CI_High),width=.2)+geom_point(size=2.5)+coord_flip()+scale_color_manual(values=class_cols)+labs(x=NULL,y='CLR contrast (LowInput - Control)')+theme_minimal(base_size=12)+theme(legend.position='none',axis.text=element_text(color='black',size=11),axis.title=element_text(face='bold'),panel.grid.major.y=element_blank(),axis.line=element_line(color='black',linewidth=.5))

# Evaluate molecular species in the same all-lipid CLR space and retain the largest shifts.
species_long <- bind_rows(ctl, lin) |> select(Sample, Condition, any_of(classes$Lipid)) |> pivot_longer(-c(Sample, Condition), names_to='Lipid', values_to='Intensity') |> left_join(classes, by='Lipid') |> filter(Class %in% names(class_cols))
species_wide <- species_long |> select(Sample, Condition, Lipid, Intensity) |> pivot_wider(names_from=Lipid, values_from=Intensity, values_fill=0)
sp_mat <- as.matrix(species_wide[, setdiff(names(species_wide), c('Sample','Condition'))])
for (i in seq_len(nrow(sp_mat))) { z <- sp_mat[i,]; z[!is.finite(z) | z < 0] <- 0; pos <- z[z > 0]; if (length(pos) == 0) next; z[z == 0] <- min(pos) * 0.5; sp_mat[i,] <- z / sum(z) }
sp_clr <- log(sp_mat) - rowMeans(log(sp_mat))
sp_keep <- colSums(is.finite(sp_clr)) == nrow(sp_clr)
sp_clr <- sp_clr[, sp_keep, drop=FALSE]
sp_eff <- colMeans(sp_clr[idx_l,,drop=FALSE]) - colMeans(sp_clr[idx_c,,drop=FALSE])
set.seed(2026)
sp_boot <- replicate(500, { c0 <- sample(idx_c,length(idx_c),TRUE); l0 <- sample(idx_l,length(idx_l),TRUE); colMeans(sp_clr[l0,,drop=FALSE])-colMeans(sp_clr[c0,,drop=FALSE]) })
top_species <- tibble(Lipid=names(sp_eff), Effect=unname(sp_eff), CI_Low=apply(sp_boot,1,quantile,.025), CI_High=apply(sp_boot,1,quantile,.975)) |> left_join(classes, by='Lipid') |> arrange(desc(abs(Effect))) |> slice_head(n=12) |> mutate(Lipid=factor(Lipid,levels=rev(Lipid)))
wrap_lipid <- function(x) stringr::str_wrap(x, width=33)
p_b <- ggplot(top_species, aes(Lipid, Effect, color=Class)) + geom_hline(yintercept=0,linetype='dashed',color='grey45') + geom_errorbar(aes(ymin=CI_Low,ymax=CI_High),width=.2) + geom_point(size=2.5) + coord_flip() + scale_color_manual(values=class_cols) + scale_y_continuous(labels=scales::label_number(accuracy=.1)) + scale_x_discrete(labels=wrap_lipid) + labs(x=NULL,y='CLR contrast (LowInput - Control)') + theme_minimal(base_size=12) + theme(legend.position='right',axis.text=element_text(color='black',size=10),axis.title=element_text(face='bold'),panel.grid.major.y=element_blank(),axis.line=element_line(color='black',linewidth=.5))
fig <- (p_a + p_b + plot_layout(widths=c(1,1.45))) + plot_annotation(tag_levels='A',theme=theme(plot.tag=element_text(face='bold',size=16)))
dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
ggsave(out,fig,width=15,height=8,units='in',dpi=300,bg='white')
