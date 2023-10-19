# Load R libraries

library(data.table)
library(stringr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)
library(UpSetR)
library(RColorBrewer)
library(cowplot)
library(scales)
library(ape)
library(vegan)
library(ComplexUpset)

# Custom functions

plot_table <- function(my_table) {
  table <- tableGrob(my_table)
  grid.newpage()
  gt <- gTree(children=gList(table))
  grid.draw(gt)
}
                                             
set_size = function(w, h, factor=1.5) {
    s = 1 * factor
    options(
        repr.plot.width=w * s,
        repr.plot.height=h * s,
        repr.plot.res=100 / factor,
        jupyter.plot_mimetypes='image/png',
        jupyter.plot_scale=1
    )
}                                             


# Define colors

colors <- c('0'='#8DD3C7','HOMO'='#BEBADA','STR'='#FB8072','TR'='#80B1D3','LTR'='#FDB462','LINE'='#B3DE69','SINE'='#FCCDE5','Retroposon'='#D9D9D9','DNA'='#FFFFB3')
colors_mobile <- c('0'='#8DD3C7','STR'='#FB8072','TR'='#80B1D3','ME - fragment'='#E0BBE4','ME - complete'='#957DAD')
colors_mobile_2 <- c('0'='#8DD3C7','STR'='#FB8072','TR'='#80B1D3','ME'='#957DAD')
colors_sharedness <- c('singleton'='#E2CFFF','polymorphic'='#7ADEA8','major'='#E995BA','shared'='#3DC1CB')
colors_range <- c('Aus'='#9292D1','Aus-absent'='#D1F0A4','Global'='#C6E8EE')
colors_novelty <- c('Found - both'='#E0D0F5','Found - Gnomad only'='#FCE4BA','Found - DECODE only'='#E1B894','Not found'='#F3AFCC','Not lifted'='#B7CAB4')
colors_community <- c('NCIGP1'='#F07DD1','NCIGP2'='#986BE3','NCIGP3'='#8DBAF3','NCIGP4'='#6DB3AC','non-NCIG'='#FCB48D')
colors_aus <- c('Community-specific'='#F4C573','Widespread'='#D95050','Private'='#7ECBE6','Ancestral'='#82C493')                                             
colors_discovery <- c('nonredundant'='#F7ECB2','shared'='#3DC1CB')
colors_community_extended <- c('NCIGP1'='#F07DD1','NCIGP2'='#986BE3','NCIGP3'='#8DBAF3','NCIGP4'='#6DB3AC','non-NCIG'='#FCB48D','NCIG'='#64369C')
colors_repeats_summary <- c('Non-repetitive'='#8DD3C7','Tandem repeats'='#E88794','Mobile elements'='#AA8FC4')                                       
colors_gene_overlap <- c('CDS'='#9BAB65','UTR/Intron/2kb flanking'='#F0AA4F')    
colors_aus_extended <- c('Community-specific'='#F4C573','Widespread'='#D95050','Private'='#7ECBE6','Ancestral'='#82C493','Aus-absent'='#D1F0A4','Global'='#C6E8EE')   

# Load master table

vars <- fread('/path/to/Figures/all.filtered.joint_call.v3.master_table.tsv')

# Get intersection info   
                                             
vars_anno <- unique(vars[FILTER=="PASS" & SVTYPE %in% c("INS","DEL"),.(ID,LEN_LABEL,SVLEN,SVTYPE,REPEAT_STATUS,REPEAT_LABEL,INTERSPERSED_COMPLETENESS,EXON_INTERSECT_GENE_ID,INTRON_INTERSECT_GENE_ID,SHAREDNESS_FULL_LABEL)])
vars_anno[EXON_INTERSECT_GENE_ID!="0" & INTRON_INTERSECT_GENE_ID!="0",INTRON_INTERSECT_GENE_ID:="0"]
vars_anno_long <- data.table(melt(vars_anno,id.vars=c("ID","LEN_LABEL","SVLEN","SVTYPE","REPEAT_STATUS","REPEAT_LABEL","INTERSPERSED_COMPLETENESS","SHAREDNESS_FULL_LABEL"),variable.name="FEATURE",value.name="ANNOTATION"))
vars_anno_long[ANNOTATION=="0",OVERLAP:="No"]
vars_anno_long[ANNOTATION!="0",OVERLAP:="Yes"]

vars_anno_long[FEATURE=="EXON_INTERSECT_GENE_ID",FEATURE:="CDS"]                                             
vars_anno_long[FEATURE=="INTRON_INTERSECT_GENE_ID",FEATURE:="UTR/Intron/2kb flanking"]
                                             
vars_anno_long <- rbind(vars_anno_long[OVERLAP=="Yes"],vars_anno[ID %in% vars_anno_long[OVERLAP=="No",.N,by=list(LEN_LABEL,ID)][N==2,ID],.(ID,LEN_LABEL,SVLEN,SVTYPE,REPEAT_STATUS,REPEAT_LABEL,INTERSPERSED_COMPLETENESS,FEATURE="Intergenic",ANNOTATION="0",OVERLAP="Yes",SHAREDNESS_FULL_LABEL)])                                             
                                             
vars_anno_count <- vars_anno_long[,.N,by=list(LEN_LABEL,FEATURE,REPEAT_STATUS,REPEAT_LABEL,INTERSPERSED_COMPLETENESS,OVERLAP)]
                                             
 
vars_anno_count[REPEAT_LABEL=="Non-repetitive",REPEAT_LABEL:="0"]
vars_anno_count[REPEAT_LABEL=="Tandem repeats" & REPEAT_STATUS=="STR",REPEAT_LABEL:="STR"]
vars_anno_count[REPEAT_LABEL=="Tandem repeats" & REPEAT_STATUS=="HOMO",REPEAT_LABEL:="STR"]
vars_anno_count[REPEAT_LABEL=="Tandem repeats" & REPEAT_STATUS=="TR",REPEAT_LABEL:="TR"]
vars_anno_count[REPEAT_LABEL=="Mobile elements" & INTERSPERSED_COMPLETENESS=="INCOMPLETE",REPEAT_LABEL:="ME - fragment"]
vars_anno_count[REPEAT_LABEL=="Mobile elements" & INTERSPERSED_COMPLETENESS=="COMPLETE",REPEAT_LABEL:="ME - complete"]

vars_anno_count_2 <- vars_anno_count[,.(N=sum(N)),by=list(LEN_LABEL,FEATURE,REPEAT_LABEL,INTERSPERSED_COMPLETENESS,OVERLAP)]

vars_anno_count_2 <- merge(vars_anno_count_2,vars_anno_count_2[,.(TOTAL=sum(N)),by=list(LEN_LABEL,FEATURE)],by=c("LEN_LABEL","FEATURE")) 

pdf('/path/to/Figures/main_figures_v2/figure5/overlap_genome_features_count.pdf',height=8,width=12)
                                             
ggplot(vars_anno_count_2,aes(FEATURE,N/TOTAL*100,fill=REPEAT_LABEL))+
    geom_bar(stat="identity")+
    geom_text(aes(x=FEATURE,y=102,label=TOTAL),check_overlap=TRUE,vjust = 0, size=5)+
    scale_fill_manual(values=colors_mobile)+
    scale_y_continuous(breaks=c(0,50,100))+
    facet_wrap(~LEN_LABEL)+
    theme_bw(base_size=20)+
    theme(aspect.ratio=1)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.text.x=element_text(angle=20,hjust=1))+
    theme(legend.position="bottom")+                                         
    labs(x=NULL,y="Percentage (%)",fill="Repeat status")

dev.off()

# Statistics

# Number of variants in protein-coding loci - indels and SV
unique(vars_anno_count_2[,.(LEN_LABEL,FEATURE,TOTAL)])[FEATURE %in% c("CDS","UTR/Intron/2kb flanking"),.(TOTAL=sum(TOTAL)),by=LEN_LABEL]

# Total number of variants in protein-coding loci - indels and SV
sum(unique(vars_anno_count_2[,.(LEN_LABEL,FEATURE,TOTAL)])[FEATURE %in% c("CDS","UTR/Intron/2kb flanking"),.(TOTAL=sum(TOTAL)),by=LEN_LABEL][,TOTAL])

# Variants in CDS
sum(unique(vars_anno_count_2[,.(LEN_LABEL,FEATURE,TOTAL)])[FEATURE %in% c("CDS"),.(TOTAL=sum(TOTAL)),by=LEN_LABEL][,TOTAL])

# CDS/Intron per individual
vars_feature_per_individual <- data.table(melt(vars[FILTER=="PASS" & SVTYPE %in% c("INS","DEL"),.(ID,Index,EXON_INTERSECT_GENE_ID,INTRON_INTERSECT_GENE_ID)],id.vars=c("ID","Index"),variable.name="FEATURE",value.name="ANNOTATION"))
vars_feature_per_individual[FEATURE=="EXON_INTERSECT_GENE_ID",FEATURE:="CDS"]                                           
vars_feature_per_individual[FEATURE=="INTRON_INTERSECT_GENE_ID",FEATURE:="UTR/Intron/2kb flanking"]
vars_feature_per_individual[ANNOTATION=="0",OVERLAP:="No"]
vars_feature_per_individual[ANNOTATION!="0",OVERLAP:="Yes"]
vars_feature_per_individual[OVERLAP=="Yes",.N,by=list(Index,FEATURE)][,.(Mean=mean(N),SD=sd(N)),by=FEATURE]

# Private and polymorphic variants CDS
vars_anno_count <- vars_anno_long[,.N,by=list(LEN_LABEL,FEATURE,REPEAT_STATUS,REPEAT_LABEL,INTERSPERSED_COMPLETENESS,OVERLAP,SHAREDNESS_FULL_LABEL)]
vars_anno_count_stats <- vars_anno_count[,.(N=sum(N)),by=list(FEATURE,SHAREDNESS_FULL_LABEL)]
vars_anno_count_stats <- merge(vars_anno_count_stats,vars_anno_count_stats[,.(TOTAL=sum(N)),by=FEATURE])
vars_anno_count_stats[,PERCENT:=N/TOTAL*100]

# Repeat status in CDS
vars_anno_count_repeat_stat <- vars_anno_count[,.(N=sum(N)),by=list(FEATURE,REPEAT_STATUS)]
vars_anno_count_repeat_stat <- merge(vars_anno_count_repeat_stat,vars_anno_count_repeat_stat[,.(TOTAL=sum(N)),by=FEATURE])
vars_anno_count_repeat_stat[,PERCENT:=N/TOTAL*100]

#############################################
                                             
# LOEUF 

vars_anno <- unique(vars[FILTER=="PASS" & SVTYPE %in% c("INS","DEL"),.(ID,LEN_LABEL,SVLEN,SVTYPE,REPEAT_STATUS,EXON_INTERSECT_GENE_ID,INTRON_INTERSECT_GENE_ID)])
vars_anno[EXON_INTERSECT_GENE_ID!="0" & INTRON_INTERSECT_GENE_ID!="0",INTRON_INTERSECT_GENE_ID:="0"]
vars_anno_long <- data.table(melt(vars_anno,id.vars=c("ID","LEN_LABEL","SVLEN","SVTYPE","REPEAT_STATUS"),variable.name="FEATURE",value.name="ANNOTATION"))
vars_anno_long[ANNOTATION=="0",OVERLAP:="No"]
vars_anno_long[ANNOTATION!="0",OVERLAP:="Yes"]

vars_anno_long[FEATURE=="EXON_INTERSECT_GENE_ID",FEATURE:="CDS"]                                             
vars_anno_long[FEATURE=="INTRON_INTERSECT_GENE_ID",FEATURE:="UTR/Intron/2kb flanking"]
                                             
vars_anno_long <- rbind(vars_anno_long[OVERLAP=="Yes"],vars_anno[ID %in% vars_anno_long[OVERLAP=="No",.N,by=list(LEN_LABEL,ID)][N==2,ID],.(ID,LEN_LABEL,SVLEN,SVTYPE,REPEAT_STATUS,FEATURE="Intergenic",ANNOTATION="0",OVERLAP="Yes")])                                                                                             
vars_genes <- vars_anno_long[,.(gene_id=tstrsplit(ANNOTATION,",",fixed=TRUE)),by=list(ID,LEN_LABEL,SVLEN,SVTYPE,REPEAT_STATUS,FEATURE)]  
                                             
vars_genes$gene_id <- as.character(vars_genes$gene_id)
                                             
vars_genes_uniq <- unique(vars_genes)
                                             
vars_genes_uniq[,gene_id:=tstrsplit(gene_id,".",fixed=TRUE)[1]]                                             
                                             
gnomad_scores <- fread('/path/to/Figures/main_figures_v2/figure5/gnomad.v2.1.1.lof_metrics.by_gene.simple.tsv')
gnomad_sizes <- fread('/path/to/references/chm13_reference_dir/chm13.draft_v1.1.gene_annotation.v4.cds_and_non-cds_sizes.tsv')
gnomad_sizes[,gene_id:=tstrsplit(gene_id,".",fixed=TRUE)[1]]
gnomad_scores <- merge(gnomad_scores,gnomad_sizes,by="gene_id")

vars_genes_uniq <- merge(vars_genes_uniq,gnomad_scores,by="gene_id")         
                                             
vars_genes_uniq <- unique(vars_genes_uniq[,.(gene,gene_id,FEATURE,gene,pLI,oe_lof,oe_lof_upper_rank,oe_lof_upper_bin)])      
                                             
vars_genes_uniq_count <- vars_genes_uniq[!is.na(oe_lof_upper_bin),.N,by=list(FEATURE,oe_lof_upper_bin)]

#vars_genes_uniq_count <- merge(vars_genes_uniq_count,vars_genes_uniq_count[,.(TOTAL=sum(N)),by=FEATURE],by="FEATURE") 

vars_genes_uniq_count <- merge(vars_genes_uniq_count,gnomad_scores[,.(TOTAL=.N),by=oe_lof_upper_bin],by="oe_lof_upper_bin") 
                                             
vars_genes_uniq_count[,PERCENT:=N/TOTAL*100]
vars_genes_uniq_count[FEATURE=="UTR/Intron/2kb flanking",FEATURE:="UTR_Intron_2kbflanking"]
vars_genes_uniq_wide <- data.table(dcast(vars_genes_uniq_count,oe_lof_upper_bin~FEATURE,value.var="PERCENT"))
sf <- max(vars_genes_uniq_wide$CDS)/max(vars_genes_uniq_wide$UTR_Intron_2kbflanking)
vars_genes_uniq_wide[,UTR_Intron_2kbflanking:=UTR_Intron_2kbflanking*sf]
vars_genes_uniq_rescaled <- data.table(melt(vars_genes_uniq_wide,id.vars = "oe_lof_upper_bin"))
vars_genes_uniq_rescaled[variable=="UTR_Intron_2kbflanking",variable:="UTR/Intron/2kb flanking"]

pdf('/path/to/Figures/main_figures_v2/figure5/gnomad_loeuf_scores_distribution.pdf',height=8,width=12)

ggplot(vars_genes_uniq_rescaled,aes(oe_lof_upper_bin,value,color=variable))+
    geom_point(alpha=0.6,size=3)+
    scale_y_continuous(name = "CDS (%)",sec.axis = sec_axis(~./sf, name="UTR/Intron/2kb flanking (%)"))+
    scale_color_manual(values=colors_gene_overlap)+
    scale_x_continuous(breaks=0:9,labels=1:10)+                                         
    theme_bw(base_size=20)+
    theme(aspect.ratio=1/2)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(legend.position="bottom")+
    guides(color="none")+
    labs(x="LOEUF decile",y="Percentage (%)",fill=NULL)

dev.off()

# SV density per LOEUF

vars_anno <- unique(vars[FILTER=="PASS" & SVTYPE %in% c("INS","DEL"),.(ID,LEN_LABEL,SVLEN,SVTYPE,REPEAT_STATUS,EXON_INTERSECT_GENE_ID,INTRON_INTERSECT_GENE_ID)])
vars_anno[EXON_INTERSECT_GENE_ID!="0" & INTRON_INTERSECT_GENE_ID!="0",INTRON_INTERSECT_GENE_ID:="0"]
vars_anno_long <- data.table(melt(vars_anno,id.vars=c("ID","LEN_LABEL","SVLEN","SVTYPE","REPEAT_STATUS"),variable.name="FEATURE",value.name="ANNOTATION"))
vars_anno_long[ANNOTATION=="0",OVERLAP:="No"]
vars_anno_long[ANNOTATION!="0",OVERLAP:="Yes"]

vars_anno_long[FEATURE=="EXON_INTERSECT_GENE_ID",FEATURE:="CDS"]                                             
vars_anno_long[FEATURE=="INTRON_INTERSECT_GENE_ID",FEATURE:="UTR/Intron/2kb flanking"]
                                             
vars_anno_long <- rbind(vars_anno_long[OVERLAP=="Yes"],vars_anno[ID %in% vars_anno_long[OVERLAP=="No",.N,by=list(LEN_LABEL,ID)][N==2,ID],.(ID,LEN_LABEL,SVLEN,SVTYPE,REPEAT_STATUS,FEATURE="Intergenic",ANNOTATION="0",OVERLAP="Yes")])                                                                                             

vars_genes <- vars_anno_long[,.(gene_id=tstrsplit(ANNOTATION,",",fixed=TRUE)),by=list(ID,LEN_LABEL,SVLEN,SVTYPE,REPEAT_STATUS,FEATURE)]  

vars_genes$gene_id <- as.character(vars_genes$gene_id)
                                             
vars_genes_uniq <- unique(vars_genes)
                                             
vars_genes_uniq[,gene_id:=tstrsplit(gene_id,".",fixed=TRUE)[1]]  

gnomad_scores <- fread('/path/to/Figures/main_figures_v2/figure5/gnomad.v2.1.1.lof_metrics.by_gene.simple.tsv')
gnomad_sizes <- fread('/path/to/references/chm13_reference_dir/chm13.draft_v1.1.gene_annotation.v4.cds_and_non-cds_sizes.tsv')
gnomad_sizes[,gene_id:=tstrsplit(gene_id,".",fixed=TRUE)[1]]
gnomad_scores <- merge(gnomad_scores,gnomad_sizes,by="gene_id")

gnomad_scores_sum <- gnomad_scores[!is.na(oe_lof_upper_bin),.(nCDS_bases=sum(nCDS_bases),nNon_CDS_bases=sum(nNon_CDS_bases)),by=oe_lof_upper_bin][order(oe_lof_upper_bin)]
gnomad_scores_sum_long <- data.table(melt(gnomad_scores_sum,id.vars = "oe_lof_upper_bin",variable.name="FEATURE",value.name="SIZE"))
gnomad_scores_sum_long[FEATURE=="nCDS_bases",FEATURE:="CDS"]
gnomad_scores_sum_long[FEATURE=="nNon_CDS_bases",FEATURE:="UTR/Intron/2kb flanking"]

vars_genes_uniq <- merge(vars_genes_uniq,gnomad_scores,by="gene_id") 

vars_loef_density <- vars_genes_uniq[!is.na(oe_lof_upper_bin),.N,by=list(FEATURE,oe_lof_upper_bin)]

vars_loef_density <- merge(vars_loef_density,gnomad_scores_sum_long,by=c("oe_lof_upper_bin","FEATURE"))
vars_loef_density[,DENSITY:=N/SIZE]
#vars_loef_density[FEATURE=="UTR/Intron/2kb flanking",FEATURE:="UTR_Intron_2kbflanking"]

#sf <- max(vars_loef_density_wide$CDS)/max(vars_loef_density_wide$UTR_Intron_2kbflanking)
#vars_loef_density_wide[,UTR_Intron_2kbflanking:=UTR_Intron_2kbflanking*sf]
#vars_loef_density_rescaled <- data.table(melt(vars_loef_density_wide,id.vars = "oe_lof_upper_bin"))
#vars_loef_density_rescaled[variable=="UTR_Intron_2kbflanking",variable:="UTR/Intron/2kb flanking"]

pdf('/path/to/Figures/main_figures_v2/figure5/gnomad_loeuf_scores_var_density.pdf',height=8,width=12)

ggplot(vars_loef_density,aes(oe_lof_upper_bin,DENSITY*10000,color=FEATURE))+
    geom_point(alpha=0.6,size=5)+
    scale_color_manual(values=colors_gene_overlap)+
    #scale_y_continuous(name = "Var density (CDS)",sec.axis = sec_axis(~./sf, name="Var density (UTR/Intron/2kb flanking)"))+
    scale_x_continuous(breaks=0:9,labels=1:10)+                                         
    theme_bw(base_size=20)+
    theme(aspect.ratio=1/2)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(legend.position="bottom")+
    guides(color="none")+
    labs(x="LOEUF decile",y="Var density",fill=NULL)

dev.off()

# Statistics

# Number of variants in LOEUF 1/2
dim(vars_genes_uniq[oe_lof_upper_bin<=1 & FEATURE=="CDS" & LEN_LABEL==">= 50bp"])[1]

# LOEUF 1/2 novel variants
vars_range <- vars[FILTER=="PASS" & LEN_LABEL==">= 50bp" & SVTYPE %in% c("INS","DEL"),.(RunID,Index,ID,LEN_LABEL,REPEAT_STATUS,SHAREDNESS_FULL,SHAREDNESS_FULL_LABEL,SHAREDNESS_NCIG,SHAREDNESS_NCIG_LABEL,GLOBAL_RANGE,UCSC_CHAIN_LIFTOVER,UCSC_GNOMAD,UCSC_DECODE,SVLEN)]

vars_range[UCSC_CHAIN_LIFTOVER=="Lifted" & UCSC_GNOMAD=="Known",UCSC_GNOMAD:="Found"]
vars_range[UCSC_CHAIN_LIFTOVER=="Lifted" & UCSC_GNOMAD=="Novel",UCSC_GNOMAD:="Not found"]
vars_range[UCSC_CHAIN_LIFTOVER!="Lifted" & is.na(UCSC_GNOMAD),UCSC_GNOMAD:="Not lifted"]

vars_range[UCSC_CHAIN_LIFTOVER=="Lifted" & UCSC_DECODE=="Overlap",UCSC_DECODE:="Found"]
vars_range[UCSC_CHAIN_LIFTOVER=="Lifted" & is.na(UCSC_DECODE),UCSC_DECODE:="Not found"]
vars_range[UCSC_CHAIN_LIFTOVER!="Lifted" & is.na(UCSC_DECODE),UCSC_DECODE:="Not lifted"]

vars_range[UCSC_GNOMAD=="Not lifted" & UCSC_DECODE=="Not lifted",UCSC_BOTH:="Not lifted"]
vars_range[UCSC_GNOMAD=="Not found" & UCSC_DECODE=="Not found",UCSC_BOTH:="Not found"]
vars_range[UCSC_GNOMAD=="Found" & UCSC_DECODE=="Found",UCSC_BOTH:="Found - both"]
vars_range[UCSC_GNOMAD=="Found" & UCSC_DECODE=="Not found",UCSC_BOTH:="Found - Gnomad only"]
vars_range[UCSC_GNOMAD=="Not found" & UCSC_DECODE=="Found",UCSC_BOTH:="Found - DECODE only"]

vars_genes_uniq_range <- merge(vars_genes_uniq[oe_lof_upper_bin<=1 & FEATURE=="CDS" & LEN_LABEL==">= 50bp",.(ID,FEATURE,oe_lof_upper_bin)],unique(vars_range[,.(ID,UCSC_BOTH)]),by="ID")

dim(vars_genes_uniq_range)[1]

dim(vars_genes_uniq_range[UCSC_BOTH=="Not found"])[1]

#######################

# Statistics

# CDS var density decile 1
vars_loef_density[FEATURE=="CDS" & oe_lof_upper_bin==0,DENSITY]

# CDS var density decile 10
vars_loef_density[FEATURE=="CDS" & oe_lof_upper_bin==9,DENSITY]

# CDS Ratio decile 10/decile 1
vars_loef_density[FEATURE=="CDS" & oe_lof_upper_bin==9,DENSITY]/vars_loef_density[FEATURE=="CDS" & oe_lof_upper_bin==0,DENSITY]

# non-CDS Ratio decile 10/decile 1
vars_loef_density[FEATURE=="UTR/Intron/2kb flanking" & oe_lof_upper_bin==9,DENSITY]/vars_loef_density[FEATURE=="UTR/Intron/2kb flanking" & oe_lof_upper_bin==0,DENSITY]

# SV density per LOEUF - broken down by repeat

vars_anno <- unique(vars[FILTER=="PASS" & SVTYPE %in% c("INS","DEL"),.(ID,LEN_LABEL,SVLEN,SVTYPE,REPEAT_STATUS,EXON_INTERSECT_GENE_ID,INTRON_INTERSECT_GENE_ID,INTERSPERSED_COMPLETENESS)])
vars_anno[EXON_INTERSECT_GENE_ID!="0" & INTRON_INTERSECT_GENE_ID!="0",INTRON_INTERSECT_GENE_ID:="0"]
vars_anno_long <- data.table(melt(vars_anno,id.vars=c("ID","LEN_LABEL","SVLEN","SVTYPE","REPEAT_STATUS","INTERSPERSED_COMPLETENESS"),variable.name="FEATURE",value.name="ANNOTATION"))
vars_anno_long[ANNOTATION=="0",OVERLAP:="No"]
vars_anno_long[ANNOTATION!="0",OVERLAP:="Yes"]

vars_anno_long[FEATURE=="EXON_INTERSECT_GENE_ID",FEATURE:="CDS"]                                             
vars_anno_long[FEATURE=="INTRON_INTERSECT_GENE_ID",FEATURE:="UTR/Intron/2kb flanking"]
                                             
vars_anno_long <- rbind(vars_anno_long[OVERLAP=="Yes"],vars_anno[ID %in% vars_anno_long[OVERLAP=="No",.N,by=list(LEN_LABEL,ID)][N==2,ID],.(ID,LEN_LABEL,SVLEN,SVTYPE,REPEAT_STATUS,INTERSPERSED_COMPLETENESS,FEATURE="Intergenic",ANNOTATION="0",OVERLAP="Yes")])                                                                                             

vars_genes <- vars_anno_long[,.(gene_id=tstrsplit(ANNOTATION,",",fixed=TRUE)),by=list(ID,LEN_LABEL,SVLEN,SVTYPE,REPEAT_STATUS,INTERSPERSED_COMPLETENESS,FEATURE)]  

vars_genes$gene_id <- as.character(vars_genes$gene_id)
                                             
vars_genes_uniq <- unique(vars_genes)
                                             
vars_genes_uniq[,gene_id:=tstrsplit(gene_id,".",fixed=TRUE)[1]]  

gnomad_scores <- fread('/path/to/Figures/main_figures_v2/figure5/gnomad.v2.1.1.lof_metrics.by_gene.simple.tsv')
gnomad_sizes <- fread('/path/to/references/chm13_reference_dir/chm13.draft_v1.1.gene_annotation.v4.cds_and_non-cds_sizes.tsv')
gnomad_sizes[,gene_id:=tstrsplit(gene_id,".",fixed=TRUE)[1]]
gnomad_scores <- merge(gnomad_scores,gnomad_sizes,by="gene_id")

gnomad_scores_sum <- gnomad_scores[!is.na(oe_lof_upper_bin),.(nCDS_bases=sum(nCDS_bases),nNon_CDS_bases=sum(nNon_CDS_bases)),by=oe_lof_upper_bin][order(oe_lof_upper_bin)]
gnomad_scores_sum_long <- data.table(melt(gnomad_scores_sum,id.vars = "oe_lof_upper_bin",variable.name="FEATURE",value.name="SIZE"))
gnomad_scores_sum_long[FEATURE=="nCDS_bases",FEATURE:="CDS"]
gnomad_scores_sum_long[FEATURE=="nNon_CDS_bases",FEATURE:="UTR/Intron/2kb flanking"]

vars_genes_uniq <- merge(vars_genes_uniq,gnomad_scores,by="gene_id")

vars_genes_uniq[REPEAT_STATUS=="0",REPEAT_LABEL:="0"]
vars_genes_uniq[REPEAT_STATUS=="STR",REPEAT_LABEL:="STR"]
vars_genes_uniq[REPEAT_STATUS=="HOMO",REPEAT_LABEL:="STR"]
vars_genes_uniq[REPEAT_STATUS=="TR",REPEAT_LABEL:="TR"]
vars_genes_uniq[REPEAT_STATUS %in% c('STR','LINE','LTR','SINE','Retroposon','DNA') & INTERSPERSED_COMPLETENESS=="INCOMPLETE",REPEAT_LABEL:="ME"]
vars_genes_uniq[REPEAT_STATUS %in% c('STR','LINE','LTR','SINE','Retroposon','DNA') & INTERSPERSED_COMPLETENESS=="COMPLETE",REPEAT_LABEL:="ME"]

vars_loef_density <- vars_genes_uniq[!is.na(oe_lof_upper_bin),.N,by=list(FEATURE,oe_lof_upper_bin,REPEAT_LABEL)]

vars_loef_density <- merge(vars_loef_density,gnomad_scores_sum_long,by=c("oe_lof_upper_bin","FEATURE"))
vars_loef_density[,DENSITY:=N/SIZE]

vars_loef_density <- merge(vars_loef_density,vars_loef_density[,.(MAX=max(DENSITY)),by=list(FEATURE,REPEAT_LABEL)],by=c("FEATURE","REPEAT_LABEL"))

pdf('/path/to/Figures/main_figures_v2/figure5/gnomad_loeuf_scores_var_density.intro_repeat.relative.pdf',height=8,width=12)

ggplot(vars_loef_density[FEATURE!="CDS"],aes(oe_lof_upper_bin,(DENSITY/MAX),color=REPEAT_LABEL))+
    geom_point(alpha=0.6,size=5)+
    scale_color_manual(values=colors_mobile_2)+
    scale_x_continuous(breaks=0:9,labels=1:10)+    
    scale_y_continuous(breaks=c(0.5,0.75,1))+
    theme_bw(base_size=20)+
    theme(aspect.ratio=1/2)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(legend.position="bottom")+
    labs(x="LOEUF decile",y="Relative var density",color=NULL)

dev.off()


# Statistics

# 
1/vars_loef_density[FEATURE!="CDS" & REPEAT_LABEL == "ME" & oe_lof_upper_bin==0,DENSITY]
vars_loef_density[FEATURE!="CDS" & REPEAT_LABEL == "STR" & oe_lof_upper_bin==0,MAX/DENSITY]

# Breakdown of Aus-specific variants based on LOEUF                                             
vars_anno_individual <- vars[FILTER=="PASS" & SVTYPE %in% c("INS","DEL"),.(ID,RunID,Index,LEN_LABEL,SVLEN,SVTYPE,GLOBAL_RANGE,Community_Country_Level,COMMUNITY_RANGE_LABEL,SHAREDNESS_FULL_LABEL,EXON_INTERSECT_GENE_ID,INTRON_INTERSECT_GENE_ID)]     
vars_anno_individual[EXON_INTERSECT_GENE_ID!="0" & INTRON_INTERSECT_GENE_ID!="0",INTRON_INTERSECT_GENE_ID:="0"]                                        
vars_anno_individual_long <- data.table(melt(vars_anno_individual,id.vars=c("ID","RunID","Index","LEN_LABEL","SVLEN","SVTYPE","GLOBAL_RANGE","Community_Country_Level","COMMUNITY_RANGE_LABEL","SHAREDNESS_FULL_LABEL"),variable.name="FEATURE",value.name="ANNOTATION"))
vars_anno_individual_long[ANNOTATION=="0",OVERLAP:="No"]
vars_anno_individual_long[ANNOTATION!="0",OVERLAP:="Yes"]

vars_anno_individual_long[FEATURE=="EXON_INTERSECT_GENE_ID",FEATURE:="CDS"]                                             
vars_anno_individual_long[FEATURE=="INTRON_INTERSECT_GENE_ID",FEATURE:="UTR/Intron/2kb flanking"]                                             
                                             
vars_anno_individual_long <- vars_anno_individual_long[OVERLAP=="Yes"]                                                                                                                                          
vars_genes_individual <- vars_anno_individual_long[,.(gene_id=tstrsplit(ANNOTATION,",",fixed=TRUE)),by=list(ID,RunID,Index,LEN_LABEL,SVLEN,SVTYPE,GLOBAL_RANGE,Community_Country_Level,COMMUNITY_RANGE_LABEL,SHAREDNESS_FULL_LABEL,FEATURE)]
                                             
vars_genes_individual$gene_id <- as.character(vars_genes_individual$gene_id)
                                             
vars_genes_individual <- unique(vars_genes_individual)
                                             
vars_genes_individual[,gene_id:=tstrsplit(gene_id,".",fixed=TRUE)[1]]
                                             
gnomad_scores <- fread('/path/to/Figures/main_figures_v2/figure5/gnomad.v2.1.1.lof_metrics.by_gene.simple.tsv')
                                             
vars_genes_individual <- merge(vars_genes_individual,gnomad_scores,by="gene_id")
                                             
vars_genes_individual <- unique(vars_genes_individual[,.(gene_id,ID,GLOBAL_RANGE,COMMUNITY_RANGE_LABEL,SHAREDNESS_FULL_LABEL,FEATURE,oe_lof_upper_bin)])
                                             
vars_genes_individual[is.na(COMMUNITY_RANGE_LABEL) & GLOBAL_RANGE=="Global",COMMUNITY_RANGE_LABEL:="Global"]           
vars_genes_individual[is.na(COMMUNITY_RANGE_LABEL) & GLOBAL_RANGE=="Aus-absent",COMMUNITY_RANGE_LABEL:="Aus-absent"]                                     
vars_genes_individual_range_count <- vars_genes_individual[,.N,by=list(FEATURE,oe_lof_upper_bin,COMMUNITY_RANGE_LABEL)] 
vars_genes_individual_range_count <- merge(vars_genes_individual_range_count,vars_genes_individual_range_count[,.(TOTAL=sum(N)),by=list(FEATURE,oe_lof_upper_bin)])                                              

pdf('/path/to/Figures/main_figures_v2/figure5/gnomad_loeuf_scores_distribution.var_range.pdf',height=8,width=12)                                                                                    
                                             
ggplot(vars_genes_individual_range_count,aes(oe_lof_upper_bin,N/TOTAL*100,fill=COMMUNITY_RANGE_LABEL))+
    geom_bar(stat="identity")+
    scale_fill_manual(values=colors_aus_extended)+
    scale_y_continuous(limits=c(0,110),breaks=c(0,50,100))+
    geom_text(aes(y=108,label=TOTAL),check_overlap=TRUE,size=5)+                                         
    facet_wrap(~FEATURE)+
    scale_x_continuous(breaks=0:9,labels=1:10)+
    theme_bw(base_size=20)+
    #theme(aspect.ratio=1)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(legend.position="bottom")+    
    coord_flip()+
    labs(x="LOEUF decile",y="Percentage (%)",fill=NULL)               
                                             
dev.off()  
