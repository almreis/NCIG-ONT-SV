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
colors_community <- c('NCIGP1'='#F07DD1','NCIGP2'='#986BE3','NCIGP3'='#8DBAF3','NCIGP4'='#6DB3AC','non-NCIG'='#FCB48D')
colors_aus <- c('Community-specific'='#F4C573','Widespread'='#D95050','Private'='#7ECBE6','Ancestral'='#82C493')

# Load master table

vars <- fread('/path/to/Figures/all.filtered.joint_call.v5.master_table.tsv')
vars_individuals <- unique(vars[,.(UNIQUE_ID,GUAID,RunID,Index,Community,Sex,FILTER,SVTYPE,SVLEN,LEN_LABEL,REPEAT_STATUS,REPEAT_LABEL,GLOBAL_RANGE,COMMUNITY_RANGE,COMMUNITY_RANGE_LABEL,SHAREDNESS_FULL,SHAREDNESS_FULL_LABEL,INTERSPERSED_COMPLETENESS,REPEAT_LABEL)])

vars_nr <- unique(vars[,.(UNIQUE_ID,CHROM,POS,ID,QUAL,FILTER,SVTYPE,SVLEN,LEN_LABEL,REPEAT_STATUS,REPEAT_LABEL,GLOBAL_RANGE,COMMUNITY_RANGE,COMMUNITY_RANGE_LABEL,SHAREDNESS_FULL,SHAREDNESS_FULL_LABEL,SHAREDNESS_NCIG,SHAREDNESS_FULL_LABEL,SHAREDNESS_NON_NCIG,SHAREDNESS_NON_NCIG_LABEL,CHROM_SIZE,DIST_TELOMERE,DIST_TELOMERE_BIN,UCSC_CHAIN_LIFTOVER,UCSC_GNOMAD,UCSC_HGSVC2,UCSC_DECODE,UCSC_ALL_DATABASES,CHM13_SPECIFIC_REGION,INTERSPERSED_COMPLETENESS)])

vars_nr$SVTYPE <- factor(vars_nr$SVTYPE,levels=c("INS","DEL","INV","DUP","TRA")) 
vars_nr$REPEAT_STATUS <- factor(vars_nr$REPEAT_STATUS,levels=c("0","HOMO","STR","TR","LTR","LINE","SINE","Retroposon","DNA","NA"))
vars_nr$SHAREDNESS_FULL_LABEL <- factor(vars_nr$SHAREDNESS_FULL_LABEL,levels=c("singleton","polymorphic","major","shared"))


# Degree of sharedness

sharedness_full <- vars_nr[FILTER=="PASS" & LEN_LABEL==">= 50bp",.N,by=list(LEN_LABEL,SHAREDNESS_FULL_LABEL,SHAREDNESS_FULL)]

pdf('/path/to/Figures/Review/figure3/degree_of_sharedness.pdf',height=8,width=8)

ggplot(sharedness_full,aes(SHAREDNESS_FULL,N,fill=SHAREDNESS_FULL_LABEL))+
    geom_bar(stat="identity")+
    scale_fill_manual(values=colors_sharedness)+
    scale_y_continuous(limits=c(0,45000),labels = function(l) { trans = l / 1000},breaks=c(0,20000,40000))+
    theme_bw(base_size=20)+
    theme(aspect.ratio=1)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(legend.position="bottom")+
    labs(x="Degree of sharedness",y="Count (x1,000)",fill=NULL)

dev.off()

# Statistics 
sharedness_full <- vars_nr[FILTER=="PASS" & LEN_LABEL==">= 50bp" & SVTYPE %in% c("INS","DEL"),.N,by=list(LEN_LABEL,SHAREDNESS_FULL_LABEL,SHAREDNESS_FULL)]
sharedness_stats <- sharedness_full[,.(N=sum(N)),by=list(LEN_LABEL,SHAREDNESS_FULL_LABEL)]
sharedness_stats[,TOTAL:=sum(sharedness_stats$N)]
sharedness_stats[,PERCENT:=N/TOTAL*100]

# Repeat status per individual

vars_individuals[REPEAT_LABEL=="Non-repetitive",REPEAT_LABEL:="0"]
vars_individuals[REPEAT_LABEL=="Tandem repeats" & REPEAT_STATUS=="STR",REPEAT_LABEL:="STR"]
vars_individuals[REPEAT_LABEL=="Tandem repeats" & REPEAT_STATUS=="HOMO",REPEAT_LABEL:="STR"]
vars_individuals[REPEAT_LABEL=="Tandem repeats" & REPEAT_STATUS=="TR",REPEAT_LABEL:="TR"]
vars_individuals[REPEAT_LABEL=="Mobile elements" & INTERSPERSED_COMPLETENESS=="INCOMPLETE",REPEAT_LABEL:="ME - fragment"]
vars_individuals[REPEAT_LABEL=="Mobile elements" & INTERSPERSED_COMPLETENESS=="COMPLETE",REPEAT_LABEL:="ME - complete"]
vars_repeat_per_individual <- vars_individuals[FILTER=="PASS" & SVTYPE %in% c("INS","DEL"),.N,by=list(RunID,Index,Community,Sex,LEN_LABEL,REPEAT_LABEL)]


# Statistics discussion
vars_repeat_per_individual_stats <- merge(vars_repeat_per_individual,vars_repeat_per_individual[LEN_LABEL==">= 50bp",.(TOTAL=sum(N)),by=Index])
vars_repeat_per_individual_stats[LEN_LABEL==">= 50bp" & REPEAT_STATUS %in% c("TR","STR")][,.(N=sum(N),TOTAL=TOTAL[1]),by=Index][,.(Mean=mean(N/TOTAL*100),SD=sd(N/TOTAL*100))]

pdf('/path/to/Figures/Review/figure3/variants_per_individual.repeat_status.pdf',height=6,width=12)

ggplot(vars_repeat_per_individual[LEN_LABEL==">= 50bp"],aes(Index,N,fill=REPEAT_LABEL))+
    geom_bar(stat="identity")+
    facet_grid(~Community,scale="free_x",space="free")+
    scale_fill_manual(values=colors_mobile)+
    scale_y_continuous(limits=c(0,22000),labels = function(l) { trans = l / 1000},breaks=c(0,10000,20000))+
    theme_bw(base_size=20)+
    theme(legend.position="bottom")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.text.x = element_blank())+
    theme(axis.ticks.x = element_blank())+
    theme(axis.title.x = element_text(margin=margin(t=20)))+
    labs(x="Individuals",y="Count (x1,000)",fill="Repeat status")

dev.off()

# Repeat classification by repeat status

sharedness_repeat <- vars_nr[FILTER=="PASS" & SVTYPE %in% c("INS","DEL"),.N,by=list(LEN_LABEL,SHAREDNESS_FULL_LABEL,SHAREDNESS_FULL,REPEAT_STATUS)]
sharedness_repeat <- merge(sharedness_repeat,sharedness_repeat[,.(TOTAL=sum(N)),by=list(LEN_LABEL,SHAREDNESS_FULL_LABEL,SHAREDNESS_FULL)],by=c("LEN_LABEL","SHAREDNESS_FULL_LABEL","SHAREDNESS_FULL"))
sharedness_repeat[,PERCENT:=N/TOTAL*100]

pdf('/path/to/Figures/Review/figure3/degree_of_sharedness.repeats.pdf',height=8,width=12)

ggplot(sharedness_repeat[LEN_LABEL==">= 50bp"],aes(SHAREDNESS_FULL,PERCENT,fill=REPEAT_STATUS))+
    geom_bar(stat="identity")+
    scale_fill_manual(values=colors)+
    scale_y_continuous(breaks=c(0,50,100))+
    theme_bw(base_size=20)+
    theme(aspect.ratio=1)+
    theme(legend.position="bottom")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x="Degree of sharedness",y="Percentage (%)",fill="Repeat status")

dev.off()

sharedness_repeat <- vars_nr[FILTER=="PASS" & SVTYPE %in% c("INS","DEL") & LEN_LABEL==">= 50bp",.(LEN_LABEL,SHAREDNESS_FULL_LABEL,SHAREDNESS_FULL,REPEAT_STATUS,REPEAT_LABEL,INTERSPERSED_COMPLETENESS)]
sharedness_repeat[REPEAT_LABEL=="Non-repetitive",REPEAT_LABEL:="0"]
sharedness_repeat[REPEAT_LABEL=="Tandem repeats" & REPEAT_STATUS=="STR",REPEAT_LABEL:="STR"]
sharedness_repeat[REPEAT_LABEL=="Tandem repeats" & REPEAT_STATUS=="HOMO",REPEAT_LABEL:="STR"]
sharedness_repeat[REPEAT_LABEL=="Tandem repeats" & REPEAT_STATUS=="TR",REPEAT_LABEL:="TR"]
sharedness_repeat[REPEAT_LABEL=="Mobile elements" & INTERSPERSED_COMPLETENESS=="INCOMPLETE",REPEAT_LABEL:="ME - fragment"]
sharedness_repeat[REPEAT_LABEL=="Mobile elements" & INTERSPERSED_COMPLETENESS=="COMPLETE",REPEAT_LABEL:="ME - complete"]
sharedness_repeat_count <- sharedness_repeat[!is.na(REPEAT_LABEL),.N,by=list(SHAREDNESS_FULL_LABEL,SHAREDNESS_FULL,REPEAT_LABEL)]
sharedness_repeat_count <- merge(sharedness_repeat_count,sharedness_repeat_count[,.(TOTAL=sum(N)),by=list(SHAREDNESS_FULL_LABEL,SHAREDNESS_FULL)])
sharedness_repeat_count[,PERCENT:=N/TOTAL*100]

sharedness_repeat_count$REPEAT_LABEL <- factor(sharedness_repeat_count$REPEAT_LABEL,levels=c("0","STR","TR","ME - fragment","ME - complete"))

pdf('/path/to/Figures/Review/figure3/degree_of_sharedness.repeats.complete_and_fragment.pdf',height=8,width=12)

ggplot(sharedness_repeat_count,aes(SHAREDNESS_FULL,PERCENT,fill=REPEAT_LABEL))+
    geom_bar(stat="identity")+
    scale_fill_manual(values=colors_mobile)+
    scale_y_continuous(breaks=c(0,50,100))+
    theme_bw(base_size=20)+
    theme(aspect.ratio=1)+
    theme(legend.position="bottom")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x="Degree of sharedness",y="Percentage (%)",fill="Repeat status")

dev.off()

# Statistics

# Proportion of tandem repeats among private and shared variants
sum(sharedness_repeat_count[SHAREDNESS_FULL==min(sharedness_repeat_count$SHAREDNESS_FULL) & REPEAT_LABEL %in% c("TR","STR"),PERCENT])
sum(sharedness_repeat_count[SHAREDNESS_FULL==max(sharedness_repeat_count$SHAREDNESS_FULL) & REPEAT_LABEL %in% c("TR","STR"),PERCENT])

# Proportion of mobile elements in private and shared variants
sum(sharedness_repeat_count[SHAREDNESS_FULL==min(sharedness_repeat_count$SHAREDNESS_FULL) & REPEAT_LABEL %in% c("ME - complete","ME - fragment"),PERCENT])
sum(sharedness_repeat_count[SHAREDNESS_FULL==max(sharedness_repeat_count$SHAREDNESS_FULL) & REPEAT_LABEL %in% c("ME - complete","ME - fragment"),PERCENT])

# Proportion of non-repetitive in private and shared variants
sum(sharedness_repeat_count[SHAREDNESS_FULL==min(sharedness_repeat_count$SHAREDNESS_FULL) & REPEAT_LABEL %in% c("0"),PERCENT])
sum(sharedness_repeat_count[SHAREDNESS_FULL==max(sharedness_repeat_count$SHAREDNESS_FULL) & REPEAT_LABEL %in% c("0"),PERCENT])

# Looking at variant frequency per sharedness

sharedness_global_range <- vars_nr[FILTER=="PASS",.N,by=list(LEN_LABEL,SHAREDNESS_FULL,SHAREDNESS_FULL_LABEL,GLOBAL_RANGE)]
sharedness_global_range <- merge(sharedness_global_range,sharedness_global_range[,.(TOTAL=sum(N)),by=list(LEN_LABEL,SHAREDNESS_FULL)],by=c("LEN_LABEL","SHAREDNESS_FULL"))
sharedness_global_range[,PERCENT:=N/TOTAL*100]

pdf('/path/to/Figures/Review/figure3/degree_of_sharedness.global_range.pdf',height=8,width=12)

ggplot(sharedness_global_range[LEN_LABEL==">= 50bp"],aes(SHAREDNESS_FULL,PERCENT,fill=GLOBAL_RANGE))+
    geom_bar(stat="identity")+
    scale_fill_manual(values=colors_range)+
    scale_y_continuous(breaks=c(0,50,100))+
    theme_bw(base_size=20)+
    theme(aspect.ratio=1)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(legend.position="bottom")+
    labs(x="Degree of sharedness",y="Percentage (%)",fill="Range")

dev.off()

# Range per individual 

range_per_individual <- unique(vars[FILTER=="PASS" & LEN_LABEL==">= 50bp",.N,by=list(RunID,Index,Community_Country_Level,GLOBAL_RANGE)])
range_per_individual <- merge(range_per_individual,range_per_individual[,.(Total=sum(N)),by=list(Community_Country_Level,RunID)],by=c("Community_Country_Level","RunID"))
range_per_individual[,PERCENT:=N/Total*100]

pdf('/path/to/Figures/Review/figure3/range_per_individual.pdf',height=8,width=8)

ggplot(range_per_individual,aes(Index,PERCENT,fill=GLOBAL_RANGE))+
    geom_bar(stat="identity")+
    facet_grid(~Community_Country_Level,space="free",scale="free_x")+
    scale_y_continuous(breaks=c(0,50,100))+
    scale_fill_manual(values=colors_range)+
    theme_bw(base_size=20)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(legend.position="bottom")+
    theme(axis.text.x = element_blank())+
    theme(axis.ticks.x = element_blank())+
    theme(axis.title.x = element_text(margin=margin(t=20)))+
    labs(x="Individuals",y="Percentage (%)",fill="Range")

dev.off()

# Statistics

# Average percentage per individual
range_per_individual <- unique(vars[FILTER=="PASS" & LEN_LABEL==">= 50bp" & SVTYPE %in% c("INS","DEL"),.N,by=list(RunID,Index,Community_Country_Level,GLOBAL_RANGE)])
range_per_individual <- merge(range_per_individual,range_per_individual[,.(Total=sum(N)),by=list(Community_Country_Level,RunID)],by=c("Community_Country_Level","RunID"))
range_per_individual[,PERCENT:=N/Total*100]
range_per_individual[,.(Mean=round(mean(PERCENT),2),SD=sd(PERCENT)),by=GLOBAL_RANGE]

# Range across cohort - proportion of aus-specific
range_stats <- vars_nr[FILTER=="PASS" & LEN_LABEL==">= 50bp" & SVTYPE %in% c("INS","DEL"),.N,by=(GLOBAL_RANGE)]
range_stats[,TOTAL:=sum(range_stats$N)]
range_stats[,PERCENT:=N/TOTAL*100]

# Inspect Aus variants that are Novel

vars_range <- vars[FILTER=="PASS" & LEN_LABEL==">= 50bp" & SVTYPE %in% c("INS","DEL"),.(RunID,Index,ID,LEN_LABEL,REPEAT_STATUS,SHAREDNESS_FULL,SHAREDNESS_FULL_LABEL,SHAREDNESS_NCIG,SHAREDNESS_NCIG_LABEL,GLOBAL_RANGE,UCSC_CHAIN_LIFTOVER,UCSC_GNOMAD,UCSC_HGSVC2,UCSC_DECODE,UCSC_ALL_DATABASES,SVLEN)]
vars_range[UCSC_CHAIN_LIFTOVER!="Lifted" & is.na(UCSC_ALL_DATABASES),UCSC_ALL_DATABASES:="Not lifted"]

vars_range_count <- vars_range[,.N,by=list(Index,GLOBAL_RANGE,UCSC_ALL_DATABASES)]

# Statistics

# Total number of NCIG-only SVs
dim(vars_nr[FILTER=="PASS" & LEN_LABEL==">= 50bp" & GLOBAL_RANGE=="Aus" & SVTYPE %in% c("INS","DEL")])[1]

# Number of SVs per NCIG individual
vars[FILTER=="PASS" & LEN_LABEL==">= 50bp" & Community_Country_Level=="NCIG" & SVTYPE %in% c("INS","DEL"),.N,by=Index][,.(Mean=mean(N),SD=sd(N))]

# Number large indels per NCIG individual
vars[FILTER=="PASS" & LEN_LABEL=="< 50bp" & Community_Country_Level=="NCIG" & SVTYPE %in% c("INS","DEL"),.N,by=Index][,.(Mean=mean(N),SD=sd(N))]

# Size large indels and SV per individual
vars[FILTER=="PASS" & Community_Country_Level=="NCIG" & SVTYPE %in% c("INS","DEL"),.(SIZE=sum(abs(SVLEN))),by=Index][,.(Mean=mean(SIZE),SD=sd(SIZE))]

# Aus specific variants per individual
vars_range_count[GLOBAL_RANGE=="Aus"][,.(N=sum(N)),by=Index][,.(Mean=mean(N),SD=sd(N))]

# Aus specific and novel variants per individual
vars_range_count[UCSC_BOTH=="Not found" & GLOBAL_RANGE=="Aus",.(Mean=mean(N),SD=sd(N))]

# Aus specific variants per individual - cumulative size
vars_range[GLOBAL_RANGE=="Aus",.(SIZE=sum(as.numeric(SVLEN))),by=Index][,.(Mean=mean(as.numeric(SIZE)),SD=sd(as.numeric(SIZE)))]

# Aus specific and novel variants per individual - cumulative size
#vars_range[GLOBAL_RANGE=="Aus" & UCSC_BOTH=="Not found",.(SIZE=sum(as.numeric(SVLEN))),by=Index][,.(Mean=mean(as.numeric(SIZE)),SD=sd(as.numeric(SIZE)))]


vars_range_count <- merge(vars_range_count,vars_range_count[,.(TOTAL=sum(N)),by=list(GLOBAL_RANGE,Index)])
vars_range_count$UCSC_ALL_DATABASES <- factor(vars_range_count$UCSC_ALL_DATABASES,levels=c("No overlap","Low overlap","Moderate overlap","High overlap",'Not lifted'))


pdf('/path/to/Figures/Review/figure3/aus_specific_enriched_novelty.pdf',height=8,width=12)

ggplot(vars_range_count,aes(Index,N/TOTAL*100,fill=UCSC_ALL_DATABASES))+
    geom_bar(stat="identity")+
    facet_grid(~GLOBAL_RANGE,space="free",scale="free")+
    scale_fill_manual(values=colors_novelty)+
    scale_y_continuous(breaks=c(0,50,100))+
    theme_bw(base_size=20)+
    theme(legend.position="bottom")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.text.x=element_blank(),axis.ticks.x = element_blank())+
    theme(axis.title.x = element_text(margin=margin(t=20)))+
    labs(x="Individuals",y="Percentage (%)",fill=NULL)

dev.off()


# Ratio of known and novel variants

vars_range_count[UCSC_ALL_DATABASES=="High overlap",LABEL:="Found"]
vars_range_count[UCSC_ALL_DATABASES=="Moderate overlap",LABEL:="Not found"]
vars_range_count[UCSC_ALL_DATABASES=="Low overlap",LABEL:="Not found"]
vars_range_count[UCSC_ALL_DATABASES=="No overlap",LABEL:="Not found"]

vars_range_ratio <- vars_range_count[!is.na(LABEL),.(N=sum(N)),by=list(Index,GLOBAL_RANGE,LABEL,TOTAL)]
vars_range_ratio[,PERCENT:=N/TOTAL*100]
vars_range_ratio <- data.table(dcast(vars_range_ratio,Index+GLOBAL_RANGE~LABEL,value.var="PERCENT"))
vars_range_ratio[,RATIO:=`Not found`/Found]
                                                          
# Statistics

# Fold-change between NCIG-only ratio and global/NCIG-absent
vars_range_ratio_stats <- vars_range_ratio[,.(Mean=mean(RATIO),Median=median(RATIO)),by=GLOBAL_RANGE]                                                          
mean(vars_range_ratio_stats[GLOBAL_RANGE=="Aus",Median]/vars_range_ratio_stats[GLOBAL_RANGE!="Aus",Median])
                                                          
pdf('/path/to/Figures/Review/figure3/aus_specific_enriched_novelty.ratio.pdf',height=8,width=6)

ggplot(vars_range_ratio,aes(GLOBAL_RANGE,RATIO,fill=GLOBAL_RANGE))+
    geom_boxplot()+
    scale_fill_manual(values=colors_range)+
    scale_y_continuous(limits=c(0.3,1.2))+                                                      
    theme_bw(base_size=20)+
    theme(aspect.ratio=2)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    guides(fill="none")+
    labs(x="Range",y="Not found / Found ratio",fill=NULL)

dev.off()

# Aus specific and novel variants

aus_specific_novel <- unique(vars_range[GLOBAL_RANGE=="Aus",.(ID,LEN_LABEL,REPEAT_STATUS,SHAREDNESS_FULL,SHAREDNESS_FULL_LABEL,SHAREDNESS_NCIG,SHAREDNESS_NCIG_LABEL,UCSC_ALL_DATABASES)])
aus_specific_novel[UCSC_ALL_DATABASES=="High overlap",LABEL:="Found"]
aus_specific_novel[UCSC_ALL_DATABASES=="Moderate overlap",LABEL:="Not found"]
aus_specific_novel[UCSC_ALL_DATABASES=="Low overlap",LABEL:="Not found"]
aus_specific_novel[UCSC_ALL_DATABASES=="No overlap",LABEL:="Not found"]
aus_specific_novel_count <- aus_specific_novel[!is.na(LABEL),.N,by=list(LEN_LABEL,SHAREDNESS_NCIG_LABEL,LABEL,REPEAT_STATUS)]

aus_specific_novel_count$SHAREDNESS_NCIG_LABEL <- factor(aus_specific_novel_count$SHAREDNESS_NCIG_LABEL,levels=c("singleton","polymorphic","major","shared"))


aus_specific_novel_count <- merge(aus_specific_novel_count,aus_specific_novel_count[,.(TOTAL=sum(N)),by=list(SHAREDNESS_NCIG_LABEL,LABEL)])

pdf('/path/to/Figures/Review/figure3/aus_specific_ncig_sharedness.repeats.pdf',height=10,width=8)

ggplot(aus_specific_novel_count,aes(SHAREDNESS_NCIG_LABEL,N,fill=REPEAT_STATUS))+
    geom_bar(stat="identity")+
    facet_grid(~LABEL)+
    geom_text(aes(x=SHAREDNESS_NCIG_LABEL,y=TOTAL+200,label =TOTAL ),check_overlap=TRUE,vjust = 0, size=5)+
    scale_fill_manual(values=colors)+
    scale_y_continuous(limits=c(0,20000),labels = function(l) { trans = l / 1000},breaks=c(0,10000,20000))+
    theme_bw(base_size=20)+
    theme(legend.position="bottom")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.text.x=element_text(angle=45,hjust=1))+
    labs(x=NULL,y="Counts (x1,000)",fill="Repeat status")

dev.off()
