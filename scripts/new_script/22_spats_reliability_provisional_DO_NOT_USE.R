# SpATS genotype-reliability estimates from the validated legacy pre-SpATS inputs.
# A_final_summed_lipids.csv contains the 394 CTL field entries; B contains 363 LIN entries.
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(ggplot2); library(patchwork); library(SpATS); library(scales)
})

root <- '/Users/nirwantandukar/Documents/Github/SoLD_paper'
script_root <- '/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database'
table_dir <- file.path(script_root, 'table/new_table')
figure_dir <- file.path(script_root, 'fig/new_figures')
dir.create(table_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(figure_dir, recursive=TRUE, showWarnings=FALSE)

class_map <- read_csv(file.path(root, 'data/lipid_class/final_lipid_classes.csv'), show_col_types=FALSE) |>
  transmute(Lipid=Lipids, BroadClass=Class)
core_classes <- c('AEG','Cer','SM','GalCer','FA','MG','DG','TG','MGDG','DGDG','SQDG','PA','PC','PE','PG','PS','LPC','LPE','DGTS')
class_colors <- c(AEG='#D55E00', Cer='#CC79A7', SM='#7B61A8', GalCer='#8C6BB1', FA='#E6AB02', MG='#8A2BE2', DG='#54278F', TG='#E76F51', MGDG='#F2A6CF', DGDG='#F768A1', SQDG='#A64D79', PA='#1B7837', PC='#00441B', PE='#41AB5D', PG='#78C679', PS='#C2E699', LPC='#66C2A4', LPE='#2CA25F', DGTS='#7FCDBB', Terpenoid='#E76F51', 'Fatty acid'='#E9C46A', Sphingolipid='#6D597A', 'Ether lipid'='#457B9D', Sterol='#2A9D8F', Prenol='#F4A261', 'Betaine lipid'='#9B2226')

plot_theme <- theme_minimal(base_size=13) + theme(
  axis.text=element_text(colour='black', size=10), axis.title=element_text(face='bold', size=12),
  axis.line=element_line(colour='black', linewidth=.5), panel.grid.major.y=element_blank(),
  panel.grid.minor=element_blank(), legend.title=element_blank(), legend.text=element_text(size=9),
  plot.margin=margin(10,12,10,12)
)

normalize_line <- function(x) sub('^PI_', 'PI', trimws(as.character(x)))

read_legacy_condition <- function(file, condition, fieldmap) {
  raw <- read_csv(file, show_col_types=FALSE, name_repair='minimal')
  sample_cols <- names(raw)[14:ncol(raw)]
  keep_samples <- sample_cols[str_detect(sample_cols, '^PI')]
  # Keep only features represented in the final legacy phenotype set.
  final_file <- file.path(root, 'data/SPATS_fitted/non_normalized_intensities',
    if (condition == 'CTL') 'Final_subset_control_all_lipids_fitted_phenotype_non_normalized.csv' else 'Final_subset_lowinput_all_lipids_fitted_phenotype_non_normalized.csv')
  final_features <- names(read_csv(final_file, show_col_types=FALSE))
  feature_info <- raw |>
    transmute(Lipid=Compound_Name, BroadClass=Class) |>
    filter(Lipid %in% final_features) |>
    distinct(Lipid, .keep_all=TRUE) |>
    left_join(class_map, by='Lipid', suffix=c('_raw','_final')) |>
    mutate(BroadClass=coalesce(BroadClass_final, BroadClass_raw)) |>
    select(Lipid, BroadClass)
  vals <- raw |>
    select(Compound_Name, all_of(keep_samples)) |>
    filter(Compound_Name %in% feature_info$Lipid) |>
    distinct(Compound_Name, .keep_all=TRUE) |>
    pivot_longer(-Compound_Name, names_to='Sample', values_to='Value') |>
    mutate(LineRaw=normalize_line(Sample), Value=as.numeric(Value))
  map <- as.matrix(read.csv(fieldmap, header=FALSE, check.names=FALSE, stringsAsFactors=FALSE))
  layout <- expand.grid(row=seq_len(nrow(map)), col=seq_len(ncol(map))) |>
    mutate(LineRaw=normalize_line(as.vector(t(map)))) |>
    filter(str_detect(LineRaw, '^PI')) |>
    distinct(LineRaw, .keep_all=TRUE) |>
    mutate(row=as.integer(row), col=as.integer(col))
  vals <- vals |> inner_join(layout, by='LineRaw') |> inner_join(feature_info, by=c('Compound_Name'='Lipid'))
  list(condition=condition, values=vals, features=feature_info)
}

# Individual molecular classes retain their explicit lipid class; non-traditional
# annotations are summarized at the broad class level.
display_class <- function(lipid, broad) {
  prefix <- str_match(lipid, '^([A-Za-z]+)\\(')[,2]
  out <- ifelse(!is.na(prefix) & prefix %in% core_classes, prefix, broad)
  ifelse(is.na(out) | out == '', 'Other annotated', out)
}

fit_reliability <- function(dat) {
  dat <- dat |> filter(is.finite(Trait), !is.na(LineRaw), !is.na(row), !is.na(col))
  if (nrow(dat) < 30L || dplyr::n_distinct(dat$LineRaw) < 30L || sd(dat$Trait) == 0) {
    return(tibble(n_field_entries=nrow(dat), n_genotypes=dplyr::n_distinct(dat$LineRaw), SpATS_genotype_reliability=NA_real_, model_status='insufficient_variation'))
  }
  mod <- tryCatch(suppressMessages(SpATS(response='Trait', spatial=~SAP(col,row,nseg=c(8,2),degree=3,pord=2), genotype='LineRaw', data=dat, control=list(tolerance=1e-3,maxit=500), genotype.as.random=TRUE)), error=function(e)e)
  if(inherits(mod,'error')) return(tibble(n_field_entries=nrow(dat), n_genotypes=dplyr::n_distinct(dat$LineRaw), SpATS_genotype_reliability=NA_real_, model_status=conditionMessage(mod)))
  term <- mod$model$geno$genotype
  h2 <- tryCatch(as.numeric(mod$eff.dim[term] / mod$dim.nom[term]), error=function(e) NA_real_)
  tibble(n_field_entries=nrow(dat), n_genotypes=dplyr::n_distinct(dat$LineRaw), SpATS_genotype_reliability=h2, model_status='ok')
}

estimate_individual <- function(obj) {
  feats <- obj$features$Lipid
  out <- vector('list', length(feats))
  for(i in seq_along(feats)) {
    f <- feats[i]
    d <- obj$values |> filter(Compound_Name == f) |> transmute(LineRaw,row,col,Trait=Value)
    out[[i]] <- fit_reliability(d) |> mutate(Condition=obj$condition, FeatureType='Individual lipid', Lipid=f, Class=display_class(f, obj$features$BroadClass[match(f,obj$features$Lipid)]), .before=1)
    if(i %% 25L == 0L || i == length(feats)) message(obj$condition, ': individual ',i,'/',length(feats))
  }
  bind_rows(out)
}

estimate_sums <- function(obj) {
  sums <- obj$values |> mutate(Class=display_class(Compound_Name, BroadClass)) |>
    group_by(LineRaw,row,col,Class) |> summarise(Trait=sum(Value,na.rm=TRUE), .groups='drop')
  cls <- sort(unique(sums$Class))
  out <- vector('list', length(cls))
  for(i in seq_along(cls)) {
    cl <- cls[i]
    out[[i]] <- fit_reliability(sums |> filter(Class==cl) |> select(LineRaw,row,col,Trait)) |>
      mutate(Condition=obj$condition, FeatureType='Lipid class sum', Lipid=NA_character_, Class=cl, .before=1)
  }
  bind_rows(out)
}

ctl <- read_legacy_condition(file.path(root,'data/summed_lipid_intensities/A_final_summed_lipids.csv'), 'CTL', file.path(root,'data/fieldmap/fieldmap_control.csv'))
lin <- read_legacy_condition(file.path(root,'data/summed_lipid_intensities/B_final_summed_lipids.csv'), 'LIN', file.path(root,'data/fieldmap/fieldmap_lowinput.csv'))
individual <- bind_rows(estimate_individual(ctl), estimate_individual(lin))
class_sums <- bind_rows(estimate_sums(ctl), estimate_sums(lin))
all_h2 <- bind_rows(individual, class_sums)
write_csv(individual, file.path(table_dir,'SuppTable_SpATS_Genotype_Reliability_Individual_Lipids.csv'))
write_csv(class_sums, file.path(table_dir,'SuppTable_SpATS_Genotype_Reliability_Lipid_Class_Sums.csv'))
write_csv(all_h2, file.path(table_dir,'SuppTable_SpATS_Genotype_Reliability_All.csv'))

top20 <- individual |> filter(model_status=='ok') |> group_by(Condition) |> slice_max(SpATS_genotype_reliability,n=20,with_ties=FALSE) |> ungroup()
write_csv(top20, file.path(table_dir,'SuppTable_Top20_SpATS_Genotype_Reliability_Individual_Lipids.csv'))

plot_lollipop <- function(d) {
  d <- d |> arrange(SpATS_genotype_reliability) |> mutate(Label=factor(str_wrap(Lipid, 34), levels=str_wrap(Lipid,34)))
  ggplot(d,aes(SpATS_genotype_reliability,Label,colour=Class)) +
    geom_segment(aes(x=0,xend=SpATS_genotype_reliability,yend=Label),colour='grey70',linewidth=.55) + geom_point(size=2.8) +
    scale_colour_manual(values=class_colors, na.value='grey35') + scale_x_continuous(limits=c(0,1),breaks=seq(0,1,.2)) +
    labs(x='SpATS genotype reliability',y=NULL) + plot_theme + theme(legend.position='right',axis.text.y=element_text(size=8))
}
plot_sums <- function(d) {
  d <- d |> arrange(SpATS_genotype_reliability) |> mutate(Class=factor(Class,levels=Class))
  ggplot(d,aes(SpATS_genotype_reliability,Class,fill=Class)) + geom_col(width=.7,colour='black',linewidth=.2) +
    geom_text(aes(label=ifelse(is.na(SpATS_genotype_reliability),'NA',sprintf('%.2f',SpATS_genotype_reliability))),hjust=-.15,size=3) +
    scale_fill_manual(values=class_colors,na.value='grey35') + scale_x_continuous(limits=c(0,1.08),breaks=seq(0,1,.2)) +
    labs(x='SpATS genotype reliability',y=NULL) + plot_theme + theme(legend.position='none',axis.text.y=element_text(size=9))
}
p_a <- plot_lollipop(top20 |> filter(Condition=='CTL'))
p_b <- plot_lollipop(top20 |> filter(Condition=='LIN'))
p_c <- plot_sums(class_sums |> filter(Condition=='CTL',model_status=='ok'))
p_d <- plot_sums(class_sums |> filter(Condition=='LIN',model_status=='ok'))
fig <- (p_a + p_b) / (p_c + p_d) + plot_annotation(tag_levels='A',theme=theme(plot.tag=element_text(face='bold',size=16)))
ggsave(file.path(figure_dir,'SuppFig_Lipid_Heritability.png'),fig,width=16,height=16,units='in',dpi=300,bg='white')
message('Saved figure and tables to: ', script_root)
