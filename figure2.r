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

# Custom functions

plot_table <- function(my_table) {
  table <- tableGrob(my_table)
  grid.newpage()
  gt <- gTree(children=gList(table))
  grid.draw(gt)
}

# Define colors

colors <- c('0'='#8DD3C7','HOMO'='#BEBADA','STR'='#FB8072','TR'='#80B1D3','LTR'='#FDB462','LINE'='#B3DE69','SINE'='#FCCDE5','Retroposon'='#D9D9D9','DNA'='#FFFFB3')
colors_mobile <- c('0'='#8DD3C7','STR'='#FB8072','TR'='#80B1D3','ME - fragment'='#E0BBE4','ME - complete'='#957DAD')
colors_sharedness <- c('singleton'='#E2CFFF','polymorphic'='#7ADEA8','major'='#E995BA','shared'='#3DC1CB')
colors_range <- c('Aus'='#9292D1','Aus-absent'='#D1F0A4','Global'='#C6E8EE')
colors_novelty <- c('High overlap'='#E0D0F5','Moderate overlap'='#FCE4BA','Low overlap'='#E1B894','No overlap'='#F3AFCC','Not lifted'='#B7CAB4')
colors_novelty_extended <- c('High overlap'='#F3AFCC','Moderate overlap'='#E1B894','Low overlap'='#FCE4BA','No overlap'='#E0D0F5','Deleted'='#7ECCC5','Partially deleted'='#91F2FF','Split'='#83A1FF')
colors_community <- c('NCIGP1'='#F07DD1','NCIGP2'='#986BE3','NCIGP3'='#8DBAF3','NCIGP4'='#6DB3AC','non-NCIG'='#FCB48D')
colors_aus <- c('Community-specific'='#F4C573','Widespread'='#D95050','Private'='#7ECBE6','Ancestral'='#82C493')                                             
colors_discovery <- c('nonredundant'='#F7ECB2','shared'='#3DC1CB')
colors_community_extended <- c('NCIGP1'='#F07DD1','NCIGP2'='#986BE3','NCIGP3'='#8DBAF3','NCIGP4'='#6DB3AC','non-NCIG'='#FCB48D','NCIG'='#64369C')
colors_repeats_summary <- c('Non-repetitive'='#8DD3C7','Tandem repeats'='#E88794','Mobile elements'='#AA8FC4')                                       
colors_gene_overlap <- c('CDS'='#9BAB65','UTR/Intron/2kb flanking'='#F0AA4F')    
colors_aus_extended <- c('Community-specific'='#F4C573','Widespread'='#D95050','Private'='#7ECBE6','Ancestral'='#82C493','Aus-absent'='#D1F0A4','Global'='#C6E8EE') 
colors_genome <- c('HG38'='#789CBB','CHM13'='#D9ADB6')

# Load master table

vars <- fread('/path/to/Figures/all.filtered.joint_call.v5.master_table.tsv')
vars_nr <- unique(vars[,.(UNIQUE_ID,CHROM,POS,ID,QUAL,FILTER,SVTYPE,SVLEN,LEN_LABEL,TANDEM_REPEAT_STATUS,INTERSPERSED_REPEAT_STATUS,INTERSPERSED_COMPLETENESS,REPEAT_STATUS,REPEAT_LABEL,GLOBAL_RANGE,COMMUNITY_RANGE,COMMUNITY_RANGE_LABEL,SHAREDNESS_FULL,SHAREDNESS_FULL_LABEL,SHAREDNESS_NCIG,SHAREDNESS_FULL_LABEL,SHAREDNESS_NON_NCIG,SHAREDNESS_NON_NCIG_LABEL,CHROM_SIZE,DIST_TELOMERE,DIST_TELOMERE_BIN,UCSC_CHAIN_LIFTOVER,UCSC_GNOMAD,UCSC_HGSVC2,UCSC_DECODE,UCSC_ALL_DATABASES,CHM13_SPECIFIC_REGION)])

vars_nr$SVTYPE <- factor(vars_nr$SVTYPE,levels=c("INS","DEL","INV","DUP","TRA")) 
vars_nr$REPEAT_STATUS <- factor(vars_nr$REPEAT_STATUS,levels=c("0","HOMO","STR","TR","LTR","LINE","SINE","Retroposon","DNA","NA"))
vars_nr$SHAREDNESS_FULL_LABEL <- factor(vars_nr$SHAREDNESS_FULL_LABEL,levels=c("singleton","polymorphic","major","shared"))

# Min QUAL we filtered
vars_nr[,.(MIN=min(QUAL),MAX=max(QUAL)),by=FILTER]

#################################################
# Statistics
# avg number of SVs per individual
mean(vars[LEN_LABEL==">= 50bp" & SVTYPE %in% c("INS","DEL"),.N,by=RunID][,N]) 

# avg number of SVs per individual after filtering
mean(vars[LEN_LABEL==">= 50bp" & FILTER=="PASS",.N,by=RunID][,N]) 

# avg insertion and deletion per individual
vars[LEN_LABEL==">= 50bp" & FILTER=="PASS",.N,by=list(RunID,SVTYPE)][SVTYPE %in% c("INS","DEL"),.(N=mean(N)),by=SVTYPE] 

# number of non-redundant SVs
dim(vars_nr[LEN_LABEL==">= 50bp" & FILTER=="PASS" & SVTYPE %in% c("INS","DEL")])[1] 

# number of larger indels
dim(vars_nr[LEN_LABEL=="< 50bp" & FILTER=="PASS" & SVTYPE %in% c("INS","DEL")])[1] 

# number of INS and DEL SVs
vars_nr[LEN_LABEL==">= 50bp" & FILTER=="PASS",.N,by=SVTYPE][SVTYPE %in% c("INS","DEL"),.(N=mean(N)),by=SVTYPE] 

# number of larger INS and DEL
vars_nr[LEN_LABEL=="< 50bp" & FILTER=="PASS",.N,by=SVTYPE][SVTYPE %in% c("INS","DEL"),.(N=mean(N)),by=SVTYPE] 

# Repeat label counts 
vars_nr_stats <- vars_nr[FILTER=="PASS" & LEN_LABEL==">= 50bp" & SVTYPE %in% c("INS","DEL"),.N,by=REPEAT_STATUS]
vars_nr_stats[,TOTAL:=sum(vars_nr_stats$N)]

# Proportion of STR/TR/Mobile element (Abstract)
(vars_nr_stats[REPEAT_STATUS=="STR",N]+vars_nr_stats[REPEAT_STATUS=="TR",N]+sum(vars_nr_stats[REPEAT_STATUS %in% c("LTR","SINE","DNA","LINE","Retroposon"),N]))/sum(vars_nr_stats$N)*100

# Repeat label counts (Results)
vars_nr_stats <- vars_nr[FILTER=="PASS" & SVTYPE %in% c("INS","DEL"),.N,by=REPEAT_STATUS]
vars_nr_stats[,TOTAL:=sum(vars_nr_stats$N)]

(vars_nr_stats[REPEAT_STATUS=="STR",N]+vars_nr_stats[REPEAT_STATUS=="TR",N]+sum(vars_nr_stats[REPEAT_STATUS %in% c("LTR","SINE","DNA","LINE","Retroposon"),N]))/sum(vars_nr_stats$N)*100

# Number of mobile elements
sum(vars_nr_stats[REPEAT_STATUS %in% c("LTR","SINE","DNA","LINE","Retroposon"),N])

# Proportion of mobile elements
sum(vars_nr_stats[REPEAT_STATUS %in% c("LTR","SINE","DNA","LINE","Retroposon"),N])/sum(vars_nr_stats$N)*100

# Number of homopolymers
sum(vars_nr_stats[REPEAT_STATUS %in% c("HOMO"),N])

# Proportion of homopolymers
sum(vars_nr_stats[REPEAT_STATUS %in% c("HOMO"),N])/sum(vars_nr_stats$N)*100

# Number of non-repetitive
sum(vars_nr_stats[REPEAT_STATUS %in% c("0"),N])

# Proportion of non-repetitive
sum(vars_nr_stats[REPEAT_STATUS %in% c("0"),N])/sum(vars_nr_stats$N)*100

# Number of complete transposition events
vars_nr[LEN_LABEL==">= 50bp" & FILTER=="PASS" & SVTYPE %in% c("INS","DEL")][INTERSPERSED_COMPLETENESS=="COMPLETE",.(UNIQUE_ID,REPEAT_STATUS)][,.N,by=REPEAT_STATUS]

# Number of very large SVs
vars_nr[LEN_LABEL==">= 50bp" & FILTER=="PASS" & SVTYPE %in% c("INS","DEL")][SVLEN>50000,.(UNIQUE_ID,SVLEN,SVTYPE,REPEAT_STATUS)][,.N,by=list(SVTYPE,REPEAT_STATUS)]

# Max deletion size
max(vars_nr[LEN_LABEL==">= 50bp" & FILTER=="PASS" & SVTYPE %in% c("INS","DEL")][SVLEN>50000,.(UNIQUE_ID,SVLEN,SVTYPE,REPEAT_STATUS)][,SVLEN])

# proportion of (STR+TR+ME) per individual after filtering (Discussion)
vars_repeat_proportion <- vars[LEN_LABEL==">= 50bp" & FILTER=="PASS" & Community_Country_Level=="NCIG" & SVTYPE %in% c("INS","DEL"),.N,by=list(RunID,REPEAT_STATUS)]
vars_repeat_proportion <- merge(vars_repeat_proportion,vars_repeat_proportion[,.(TOTAL=sum(N)),by=RunID])
vars_repeat_proportion[,PERCENT:=N/TOTAL*100]

#################################################


# Number of SVs per repeat status

vars_nr_repeat_count <- vars_nr[FILTER=="PASS" & SVTYPE %in% c("INS","DEL"),.N,by=list(SVTYPE,LEN_LABEL,REPEAT_STATUS)]
vars_nr_repeat_count[SVTYPE=="INS",ALPHA:=TRUE]
vars_nr_repeat_count[SVTYPE=="DEL",ALPHA:=FALSE]
vars_nr_repeat_count$SVTYPE <- factor(vars_nr_repeat_count$SVTYPE,levels=c("DEL","INS"))

pdf('/path/to/Figures/Review/figure2/repeat_status_count.pdf',height=8,width=8)

ggplot(vars_nr_repeat_count,aes(REPEAT_STATUS,N/1000,fill=REPEAT_STATUS,by=SVTYPE,alpha=ALPHA))+
    geom_bar(stat="identity",position=position_dodge())+
    scale_fill_manual(values=colors)+
    scale_alpha_discrete(range = c(0.6, 1))+
    facet_grid(~LEN_LABEL,drop=TRUE)+
    theme_bw(base_size=20)+
    theme(aspect.ratio=1)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.text.x = element_text(angle=45,hjust=1))+
    theme(legend.position="bottom")+
    guides(alpha="none")+
    labs(x="Repeat status",y="Count (x1,000)",fill=NULL)

dev.off()

# Location in the chromosome

vars_nr_location <- vars_nr[FILTER=="PASS" & SVTYPE %in% c("INS","DEL"),.(UNIQUE_ID,SVTYPE,REPEAT_STATUS,LEN_LABEL,CHROM,CHROM_SIZE,DIST_TELOMERE,DIST_TELOMERE_BIN)]

dist_to_telomere_count <- vars_nr_location[,.N,by=list(CHROM,DIST_TELOMERE_BIN)]

dist_to_telomere_count[CHROM %in% c("chr13","chr14","chr15","chr21","chr22"),CHROM_CLASSIFICATION:="Acrocentric"]
dist_to_telomere_count[is.na(CHROM_CLASSIFICATION) & !(CHROM %in% c("chrX","chrY")),CHROM_CLASSIFICATION:="Other autosomes"]

pdf('/path/to/Figures/Review/figure2/chromosome_location.pdf',height=8,width=8)

ggplot(dist_to_telomere_count,aes(DIST_TELOMERE_BIN/1000000,N))+
    geom_point(alpha=0.6,size=2)+
    geom_vline(xintercept = 5,linetype="dashed",color="red",size=0.3)+
    theme_bw(base_size=20)+
    theme(aspect.ratio=1)+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    labs(x="Telomere Distance (Mb)",y="SVs per 500 Kb")

dev.off()

# SV length histogram - small variants resolution

vars_nr_sm_len <- vars_nr[FILTER=="PASS" & SVTYPE %in% c("INS","DEL") & LEN_LABEL=="< 50bp",.(UNIQUE_ID,SVTYPE,SVLEN,REPEAT_STATUS)]
vars_nr_sm_len[SVTYPE=="DEL",SVLEN:=SVLEN*(-1)]
vars_nr_sm_len$SVTYPE <- factor(vars_nr_sm_len$SVTYPE,levels=c("DEL","INS"))

pdf('/path/to/Figures/Review/figure2/repeat_length_small_variants.pdf',height=6,width=12)

ggplot(vars_nr_sm_len[(SVLEN>=20 & SVLEN<50) | ((SVLEN<-50 & SVLEN<=-20))],aes(SVLEN,fill=REPEAT_STATUS))+
          geom_histogram(bins=80)+
          scale_fill_manual(values=colors)+     
          scale_y_continuous(limits=c(0,8000),labels = function(l) { trans = l / 1000},breaks=c(0,4000,8000))+
          theme_bw(base_size=20)+
          facet_wrap(~SVTYPE,scales="free")+
          theme(aspect.ratio=1/2)+
          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
          theme(legend.position = "bottom")+
          labs(x="SV length (bp)",y="Count (x1,000)",fill="Repeat status")
dev.off()

# SV length histogram - SINE resolution

vars_nr_len <- vars_nr[FILTER=="PASS" & SVTYPE %in% c("INS","DEL") & LEN_LABEL==">= 50bp",.(UNIQUE_ID,SVTYPE,SVLEN,REPEAT_STATUS)]
vars_nr_len[SVTYPE=="DEL",SVLEN:=SVLEN*(-1)]
vars_nr_len$SVTYPE <- factor(vars_nr_len$SVTYPE,levels=c("DEL","INS"))

pdf('/path/to/Figures/Review/figure2/repeat_length_SINE.pdf',height=6,width=12)

ggplot(vars_nr_len[(SVLEN>=50 & SVLEN<1000) | ((SVLEN<=-50 & SVLEN>-1000))],aes(SVLEN,fill=REPEAT_STATUS))+
          geom_histogram(bins=80)+
          scale_fill_manual(values=colors)+     
          scale_y_continuous(limits=c(0,6500),labels = function(l) { trans = l / 1000},breaks=c(0,3000,6000))+
          theme_bw(base_size=20)+
          facet_wrap(~SVTYPE,scales="free")+
          theme(aspect.ratio=1/2)+
          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
          theme(legend.position = "bottom")+
          labs(x="SV length (bp)",y="Count (x1,000)",fill="Repeat status")

dev.off()

pdf('/path/to/Figures/Review/figure2/repeat_length_LINE.pdf',height=6,width=12)

ggplot(vars_nr_len[(SVLEN>=1000 & SVLEN<8000) | ((SVLEN<=-1000 & SVLEN>-8000))],aes(SVLEN,fill=REPEAT_STATUS))+
          geom_histogram(bins=80)+
          scale_fill_manual(values=colors)+     
          scale_y_continuous(limits=c(0,800),labels = function(l) { trans = l / 100},breaks=c(0,400,800))+
          theme_bw(base_size=20)+
          facet_wrap(~SVTYPE,scales="free")+
          theme(aspect.ratio=1/2)+
          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
          theme(legend.position = "bottom")+
          labs(x="SV length (bp)",y="Count (x100)",fill="Repeat status")

dev.off()

# Density repeat length 

vars_nr_len <- vars_nr[FILTER=="PASS" & SVTYPE %in% c("INS","DEL"),.(UNIQUE_ID,SVTYPE,SVLEN,REPEAT_STATUS,INTERSPERSED_COMPLETENESS)]
vars_nr_len$SVTYPE <- factor(vars_nr_len$SVTYPE,levels=c("DEL","INS"))

pdf('/path/to/Figures/Review/figure2/repeat_length_density.pdf',height=8,width=12)

ggplot()+
          geom_density(data=vars_nr_len[INTERSPERSED_COMPLETENESS!="INCOMPLETE",.(SVTYPE,REPEAT_STATUS,INTERSPERSED_COMPLETENESS,SVLEN=ifelse(SVLEN<0,-1*SVLEN,SVLEN))],aes(log2(SVLEN),..scaled..,fill=REPEAT_STATUS),size=0.1)+
          geom_density(data=vars_nr_len[INTERSPERSED_COMPLETENESS=="INCOMPLETE",.(SVTYPE,REPEAT_STATUS,INTERSPERSED_COMPLETENESS,SVLEN=ifelse(SVLEN<0,-1*SVLEN,SVLEN))],aes(log2(SVLEN),..scaled..,fill=REPEAT_STATUS),alpha=0.2,size=0.1)+
          scale_fill_manual(values=colors)+
          scale_y_continuous(breaks=c(0,0.5,1))+
          scale_x_continuous(limits=c(log2(20),log2(8000)),breaks=c(4,8,12))+
          geom_vline(xintercept = log2(50),linetype="dashed",size=0.2)+
          geom_vline(xintercept = log2(330),linetype="dashed",size=0.2)+
          geom_vline(xintercept = log2(6000),linetype="dashed",size=0.2)+
          theme_bw(base_size=20)+
          theme(strip.text.x = element_blank())+
          facet_wrap(~REPEAT_STATUS,ncol=1)+
          guides(color="none")+
          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
          labs(x="SV length (Log2)",y="Density",fill="Repeat status")

dev.off()

# Integrate Gnomad and decode results

vars_novelty <- vars_nr[FILTER=="PASS" & LEN_LABEL==">= 50bp" & SVTYPE %in% c("INS","DEL"),.(LEN_LABEL,REPEAT_STATUS,UCSC_CHAIN_LIFTOVER,UCSC_GNOMAD,UCSC_HGSVC2,UCSC_DECODE,UCSC_ALL_DATABASES)]

vars_novelty[UCSC_CHAIN_LIFTOVER=="Deleted_in_new",UCSC_ALL_DATABASES:="Deleted"]
vars_novelty[UCSC_CHAIN_LIFTOVER=="Partially_deleted_in_new",UCSC_ALL_DATABASES:="Partially deleted"]
vars_novelty[UCSC_CHAIN_LIFTOVER=="Split_in_new",UCSC_ALL_DATABASES:="Split"]

vars_novelty[UCSC_CHAIN_LIFTOVER!="Lifted",UCSC_CHAIN_LIFTOVER:="Not lifted"]

vars_novelty <- vars_novelty[,.N,by=list(LEN_LABEL,REPEAT_STATUS,UCSC_CHAIN_LIFTOVER,UCSC_ALL_DATABASES)]

vars_novelty_total <- vars_novelty[,.(N=sum(N)),by=UCSC_CHAIN_LIFTOVER]
vars_novelty_total <- vars_novelty_total[order(N,decreasing = TRUE)]
vars_novelty_total[,UCSC_ALL_DATABASES:=NA]

vars_novelty$UCSC_CHAIN_LIFTOVER <- factor(vars_novelty$UCSC_CHAIN_LIFTOVER,levels=c('Not lifted','Lifted'))
#vars_novelty$UCSC_BOTH <- factor(vars_novelty$UCSC_ALL_DATABASES,levels=c('Found - both','Found - Gnomad only','Found - DECODE only','Not found','Deleted','Partially deleted','Split'))
vars_novelty$UCSC_ALL_DATABASES <- factor(vars_novelty$UCSC_ALL_DATABASES,levels=c("No overlap","Low overlap","Moderate overlap","High overlap",'Deleted','Partially deleted','Split'))

pdf('/path/to/Figures/Review/figure2/liftover_gnomad_and_decode_intersection.pdf',height=8,width=12)

ggplot(vars_novelty,aes(UCSC_CHAIN_LIFTOVER,N,fill=UCSC_ALL_DATABASES))+
    geom_bar(stat="identity")+
    geom_text(data=vars_novelty_total,aes(x=UCSC_CHAIN_LIFTOVER,y=N+3000,label = N ), vjust = 0, size=5)+
    scale_fill_manual(values=colors_novelty_extended)+
    facet_grid(~LEN_LABEL)+
    theme_bw(base_size=20)+
    theme(aspect.ratio=0.25)+
    theme(legend.position="bottom")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    coord_flip()+
    labs(x="Liftover",y="Count",fill=NULL)

dev.off()

vars_novelty_percent <- vars_novelty[UCSC_CHAIN_LIFTOVER=="Lifted",.(N=sum(N)),by=UCSC_ALL_DATABASES]
vars_novelty_percent[,PERCENT:=N/sum(vars_novelty_percent$N)*100]

vars_novelty_type_lifted <- merge(vars_novelty[UCSC_CHAIN_LIFTOVER=="Lifted"],vars_novelty[UCSC_CHAIN_LIFTOVER=="Lifted",.(TOTAL=sum(N)),by=REPEAT_STATUS])
vars_novelty_type_lifted[,PERCENTAGE:=N/TOTAL*100]
vars_novelty_type[UCSC_CHAIN_LIFTOVER=="Not lifted",UCSC_ALL_DATABASES:="Not lifted"]

pdf('/path/to/Figures/Review/figure2/liftover_gnomad_and_decode_intersection.per_type.lifted.pdf',height=8,width=12)


ggplot(vars_novelty_type_lifted,aes(REPEAT_STATUS,PERCENTAGE,fill=UCSC_ALL_DATABASES))+
    geom_bar(stat="identity")+
    scale_fill_manual(values=colors_novelty)+
    theme_bw(base_size=20)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(aspect.ratio=1)+
    theme(axis.text.x=element_text(angle=45,hjust=1))+    
    labs(x=NULL,y="Percentage (%)",fill="Annotation")

dev.off()

vars_novelty_type_all <- vars_novelty[,.(N=sum(N)),by=list(REPEAT_STATUS,UCSC_CHAIN_LIFTOVER)]
vars_novelty_type_all <- merge(vars_novelty_type_all,vars_novelty_type_all[,.(TOTAL=sum(N)),by=REPEAT_STATUS])
vars_novelty_type_all[,Percentage:=N/TOTAL*100]

pdf('/path/to/Figures/Review/figure2/liftover_gnomad_and_decode_intersection.per_type.all.pdf',height=8,width=12)

ggplot(vars_novelty_type_all,aes(REPEAT_STATUS,Percentage,fill=UCSC_CHAIN_LIFTOVER))+
    geom_bar(stat="identity")+
    #scale_fill_manual(values=colors_novelty)+
    theme_bw(base_size=20)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(aspect.ratio=1)+
    theme(axis.text.x=element_text(angle=45,hjust=1))+    
    labs(x=NULL,y="Percentage (%)",fill="Annotation")

dev.off()

vars_novelty[,Percentage:=N/sum(vars_novelty[,N])*100]

#####################################################

# Statistics

# Not lifted - deleted
vars_novelty[UCSC_BOTH=="Deleted",N]/sum(vars_novelty$N)*100

# Not lifted - partially deleted
vars_novelty[UCSC_BOTH=="Partially deleted",N]/sum(vars_novelty$N)*100

# Number of lifted variants
sum(vars_novelty[UCSC_CHAIN_LIFTOVER=="Lifted",N])

# Total of INS and DEL
sum(vars_novelty[,N])

# Number of known variants
sum(vars_novelty[UCSC_CHAIN_LIFTOVER=="Lifted"][UCSC_BOTH!="Not found",N])

# Proportion of known variants
sum(vars_novelty[UCSC_CHAIN_LIFTOVER=="Lifted"][UCSC_BOTH!="Not found",N])/sum(vars_novelty[,N])*100
