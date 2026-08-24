suppressPackageStartupMessages({library(readr);library(dplyr);library(stringr);library(ggplot2);library(patchwork)})
repo <- '/Users/nirwantandukar/Documents/Github/SAP-Lipidomics-Database'
infile <- file.path(repo,'data/LD_mapped/go_enrichment/GO_BP_enrichment_by_lipid_class_sumratio_LD.tsv')
outfile <- file.path(repo,'fig/new_figures/SuppFig_GO_Enrichment_Class_Sums_Ratios_By_Lipid.png')
plot_theme <- theme_minimal(base_size=13)+theme(axis.text=element_text(colour='black',size=9),axis.title=element_text(face='bold',size=11),axis.line=element_line(colour='black'),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),legend.title=element_blank(),plot.margin=margin(12,90,12,12),plot.tag=element_text(face='bold',size=16),plot.tag.position=c(0,1))
d <- read_tsv(infile,show_col_types=FALSE) |> filter(p_value<.05) |> mutate(Condition=factor(condition,levels=c('CTL','LIN')),Term=str_wrap(go_term,42),BH=p_adj_bh<.05,score=-log10(p_adj_bh))
make_panel <- function(x,tag) {
  x <- x |> arrange(lipid_class,desc(BH),p_adj_bh,go_term) |> mutate(y_id=row_number())
  blocks <- x |> group_by(lipid_class) |> summarise(ymin=min(y_id)-.45,ymax=max(y_id)+.45,ymid=mean(c(min(y_id),max(y_id))),.groups='drop')
  xmax <- max(x$score,na.rm=TRUE)*1.45
  ggplot(x,aes(score,y_id))+
    geom_rect(data=blocks,aes(xmin=xmax*1.04,xmax=xmax*1.20,ymin=ymin,ymax=ymax),inherit.aes=FALSE,fill='grey96',colour='grey72',linewidth=.3)+
    geom_segment(aes(x=0,xend=score,y=y_id,yend=y_id),colour='grey78',linewidth=.45)+
    geom_point(aes(size=term_test_count,fill=BH),shape=21,colour='black',stroke=.25)+
    geom_text(data=blocks,aes(x=xmax*1.12,y=ymid,label=lipid_class),inherit.aes=FALSE,fontface='bold',size=3.3)+
    geom_vline(xintercept=-log10(.05),linetype='dashed',colour='#0072B2',linewidth=.7)+
    scale_fill_manual(values=c('FALSE'='#E69F00','TRUE'='#009E73'),labels=c('Nominal only','BH-significant'))+
    scale_size_continuous(name='Candidate genes',range=c(2.8,7))+
    scale_y_continuous(breaks=x$y_id,labels=x$Term,expand=expansion(mult=c(.02,.02)))+
    scale_x_continuous(expand=expansion(mult=c(.02,.45)))+
    labs(x=expression(-log[10]*'(BH-adjusted p)'),y=NULL,fill=NULL,tag=tag)+
    plot_theme+coord_cartesian(clip='off')
}
fig <- make_panel(filter(d,Condition=='CTL'),'A') | make_panel(filter(d,Condition=='LIN'),'B')
ggsave(outfile,fig,width=17,height=max(11,.40*max(table(d$Condition))+4),dpi=350,bg='white')
