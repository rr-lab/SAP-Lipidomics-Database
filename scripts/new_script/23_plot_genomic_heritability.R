suppressPackageStartupMessages({library(readr);library(dplyr);library(stringr);library(ggplot2);library(patchwork)})
out <- '/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database'
tab <- file.path(out,'table/new_table')
plot_theme <- theme_minimal(base_size=13)+theme(axis.text=element_text(colour='black',size=9),axis.title=element_text(face='bold',size=11),axis.line=element_line(colour='black'),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),legend.title=element_blank())
cols <- c(AEG='#D55E00',Cer='#CC79A7',SM='#7B61A8',GalCer='#8C6BB1',FA='#E6AB02',MG='#8A2BE2',DG='#B8860B',TG='#D53E9D',MGDG='#129DA5',DGDG='#B59A00',SQDG='#A15DD7',PA='#1B7837',PC='#0099C7',PE='#2CA25F',PG='#78C679',PS='#C2E699',LPC='#66C2A4',LPE='#2CA25F',DGTS='#7FCDBB','Fatty acid'='#E9C46A','Fatty acid amide'='#C2A15A',Sphingolipid='#6D597A',Terpenoid='#E76F51',Sterol='#2A9D8F','Ether lipid'='#457B9D',Prenol='#F4A261','Betaine lipid'='#9B2226')
ind <- read_csv(file.path(tab,'SuppTable_Genomic_Heritability_Individual_Lipids.csv'),show_col_types=FALSE) |> filter(model_status=='ok',is.finite(h2),Class%in%names(cols))
short <- function(x) ifelse(nchar(x)>42,paste0(substr(x,1,39),'...'),x)
top_plot <- function(x,t) {d <- x |> slice_max(h2,n=20,with_ties=FALSE) |> arrange(h2) |> mutate(Label=factor(short(Lipid),levels=short(Lipid)));ggplot(d,aes(h2,Label,colour=Class))+geom_segment(aes(x=0,xend=h2,y=Label,yend=Label),colour='grey70',linewidth=.45)+geom_point(size=2.8)+scale_colour_manual(values=cols)+scale_x_continuous(expand=expansion(mult=c(0,.08)))+labs(x=expression('Genomic heritability ('*h^2*')'),y=NULL,tag=t)+plot_theme+theme(axis.text.y=element_text(size=8),legend.position='right')}
fig <- top_plot(filter(ind,Condition=='CTL'),'A') | top_plot(filter(ind,Condition=='LIN'),'B')
ggsave(file.path(out,'fig/new_figures/SuppFig_Lipid_Heritability.png'),fig,width=15,height=9,dpi=350,bg='white')
