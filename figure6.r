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
library(Biostrings)
library(gtools)
library(gplots)


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

roll_seq <- function(seq,index) {
    letters <- unlist(str_split(seq,''))
    all_combo <- list()
    for (i in letters) {
        instance <- paste0(letters, collapse='')
        letters <- c(letters[2:length(letters)],letters[1])
        all_combo <- append(all_combo, instance)
    }
    all_combo <- unlist(all_combo)
    rev_comp <- sapply(all_combo,function(x) as.character(reverseComplement(DNAString(x))))
    names(rev_comp) <- NULL
    all_combo <- c(all_combo,rev_comp)
    all_combo_dt <- data.table(SEQ=seq,ROLLED=all_combo,INDEX=index)
    return(all_combo_dt[SEQ!=ROLLED])
}

find_rolled_label <- function(len) {                       
    combo <- permutations(4,len,c('A','C','G','T'),repeats.allowed=TRUE)
    combo_seq <- apply(combo, 1, function(x)paste0(x, collapse=''))
    a <- paste0(rep('A',len),collapse='')
    c <- paste0(rep('C',len),collapse='') 
    g <- paste0(rep('G',len),collapse='')
    t <- paste0(rep('T',len),collapse='')                   
    combo_seq <- combo_seq[!(combo_seq %in% c(a,c,g,t))]                           
    rolled_dt <- rbindlist(lapply(seq_along(combo_seq), function(x) roll_seq(combo_seq[x],x)))
    rolling_dict <- rbindlist(lapply(seq_along(combo_seq), function(x) data.table(SEQ=combo_seq[x],INDEX=x,ROLLED_INDEX=min(rolled_dt[ROLLED==combo_seq[x],INDEX]))))
    rolling_dict[,FINAL_INDEX:=ifelse(INDEX<=ROLLED_INDEX,INDEX,ROLLED_INDEX)]                                 
    rolling_dict[,LABEL:=combo_seq[FINAL_INDEX]]
    return(rolling_dict[,.(SEQ,LABEL)])                                 
} 

# Define colors

colors <- c('0'='#8DD3C7','HOMO'='#BEBADA','STR'='#FB8072','TR'='#80B1D3','LTR'='#FDB462','LINE'='#B3DE69','SINE'='#FCCDE5','Retroposon'='#D9D9D9','DNA'='#FFFFB3')
colors_mobile <- c('0'='#8DD3C7','STR'='#FB8072','TR'='#80B1D3','ME - fragment'='#E0BBE4','ME - complete'='#957DAD')
colors_sharedness <- c('singleton'='#E2CFFF','polymorphic'='#7ADEA8','major'='#E995BA','shared'='#3DC1CB')
colors_range <- c('Aus'='#9292D1','Aus-absent'='#D1F0A4','Global'='#C6E8EE')
colors_novelty <- c('Found - both'='#E0D0F5','Found - Gnomad only'='#FCE4BA','Found - DECODE only'='#E1B894','Not found'='#F3AFCC','Not lifted'='#B7CAB4')
colors_community <- c('NCIGP1'='#F07DD1','NCIGP2'='#986BE3','NCIGP3'='#8DBAF3','NCIGP4'='#6DB3AC','non-NCIG'='#FCB48D')
colors_aus <- c('Community-specific'='#F4C573','Widespread'='#D95050','Private'='#7ECBE6','Ancestral'='#82C493')                                             
colors_discovery <- c('nonredundant'='#F7ECB2','shared'='#3DC1CB')
colors_community_extended <- c('NCIGP1'='#F07DD1','NCIGP2'='#986BE3','NCIGP3'='#8DBAF3','NCIGP4'='#6DB3AC','non-NCIG'='#FCB48D','NCIG'='#64369C')
colors_repeats_summary <- c('Non-repetitive'='#8DD3C7','Tandem repeats'='#E88794','Mobile elements'='#AA8FC4')                                       
colors_gene_overlap <- c('CDS'='#9BAB65','UTR/Intron/2kb flanking'='#F0AA4F','Intergenic'='#5F588A')    
colors_aus_extended <- c('Community-specific'='#F4C573','Widespread'='#D95050','Private'='#7ECBE6','Ancestral'='#82C493','Aus-absent'='#D1F0A4','Global'='#C6E8EE')   
colors_str_tr <- c('STR - (2nt)'='#59130F','STR - (3nt)'='#FDD5BB','STR - (4nt)'='#D03E39','STR - (5nt)'='#CE7D74','STR - (6nt)'='#FA6243','STR - (7-12nt)'='#EA9D90','TR'='#80B1D3') 
colors_bases <- c('A'='#47BC1C','C'='#6360AE','G'='#FF8C33','T'='#FF0013') 
                                     
# Load master table

vars <- fread('/path/to/Figures/all.filtered.joint_call.v3.master_table.tsv')
                                     
# Find label for rolled out microsatellites
                                     
all_rolled_labels <- rbindlist(lapply(2:6,function(x) find_rolled_label(x)))                                        
all_rolled_labels <- all_rolled_labels[,.(MAJOR_MOTIF=SEQ,MOTIF_LABEL=LABEL)]
                                      
# non-redundant STR and TR variants 

vars_repeats <- unique(vars[FILTER=="PASS" & SVTYPE %in% c("INS","DEL") & REPEAT_STATUS %in% c("STR","TR"),.(ID,SVTYPE,LEN_LABEL,SVLEN,EXON_INTERSECT_GENE_ID,INTRON_INTERSECT_GENE_ID,TANDEM_REPEAT_STATUS,REPEAT_STATUS,REPEAT_LABEL,MAJOR_MOTIF,TRF_ELEMENT_SIZE,TRF_INTERSECTION,FRAC_SV_COVERED)])
vars_repeats[MAJOR_MOTIF!="100plus",MOTIF_SIZE:=nchar(MAJOR_MOTIF)]
vars_repeats[is.na(MOTIF_SIZE),MOTIF_SIZE:=100]

# Add rolled out motif label
vars_repeats <- merge(vars_repeats,all_rolled_labels,by='MAJOR_MOTIF',all.x=TRUE)                                                                           
# Sub-category of STRs
vars_repeats[MOTIF_SIZE==2,STR_LABEL:="(2nt)"]                                      
vars_repeats[MOTIF_SIZE==3,STR_LABEL:="(3nt)"]                                      
vars_repeats[MOTIF_SIZE==4,STR_LABEL:="(4nt)"]                  
vars_repeats[MOTIF_SIZE==5,STR_LABEL:="(5nt)"]                                        
vars_repeats[MOTIF_SIZE==6,STR_LABEL:="(6nt)"]
vars_repeats[MOTIF_SIZE>6 & MOTIF_SIZE<=12,STR_LABEL:="(7-12nt)"]                           

vars_repeats[REPEAT_STATUS=="STR" & (MOTIF_SIZE>=2 & MOTIF_SIZE<=12),STR_LABEL:=paste(REPEAT_STATUS,STR_LABEL,sep=" - ")]                                      
vars_repeats[is.na(STR_LABEL),STR_LABEL:=REPEAT_STATUS]                                     
                                      
vars_repeats[SVTYPE=="DEL",TRF_INTERSECTION:=TRF_INTERSECTION*(-1)]
vars_repeats$SVTYPE <- factor(vars_repeats$SVTYPE,levels=c("DEL","INS"))

# Cap length at 100bp                                      
                                      
vars_repeats_len <- data.table(vars_repeats)                                      
vars_repeats_len[TRF_INTERSECTION>=100,TRF_INTERSECTION:=100]
vars_repeats_len[TRF_INTERSECTION<=-100,TRF_INTERSECTION:=-100]
 
# Plot variant length distribution                                    
                                      
pdf('/path/to/Figures/main_figures_v2/figure6/repeat_length_STR_and_TR.ins.pdf',height=6,width=12)
                                      
ggplot(vars_repeats_len[(TRF_INTERSECTION>=20 & TRF_INTERSECTION<=101) | ((TRF_INTERSECTION>=-101 & TRF_INTERSECTION<=-20))][SVTYPE=="INS"],aes(TRF_INTERSECTION,fill=STR_LABEL))+
    geom_histogram(bins=80)+
    facet_wrap(~SVTYPE,scales="free")+
    scale_y_continuous(limits=c(0,4000),labels = function(l) { trans = l / 1000},breaks=c(0,2000,4000))+
    scale_fill_manual(values=colors_str_tr)+
    theme(aspect.ratio=1/2)+
    theme_bw(base_size=20)+
    theme(legend.position="bottom")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+                          
    labs(x="Repeat length",y="Count (x1,000)",fill=NULL)
                                      
dev.off()   
                                      
# Statistics
                                      
# Number of STR expansion or contraction
dim(vars_repeats_len[REPEAT_STATUS!="TR" & SVTYPE=="INS"])[1]                                      
dim(vars_repeats_len[REPEAT_STATUS!="TR" & SVTYPE=="DEL"])[1]
                                      
                                      
pdf('/path/to/Figures/main_figures_v2/figure6/repeat_length_STR_and_TR.del.pdf',height=6,width=12)
                                      
ggplot(vars_repeats_len[(TRF_INTERSECTION>=20 & TRF_INTERSECTION<=101) | ((TRF_INTERSECTION>=-101 & TRF_INTERSECTION<=-20))][SVTYPE=="DEL"],aes(TRF_INTERSECTION,fill=STR_LABEL))+
    geom_histogram(bins=80)+
    facet_wrap(~SVTYPE,scales="free")+
    scale_y_continuous(limits=c(0,4000),labels = function(l) { trans = l / 1000},breaks=c(0,2000,4000))+
    scale_fill_manual(values=colors_str_tr)+
    theme(aspect.ratio=1/2)+
    theme_bw(base_size=20)+
    theme(legend.position="bottom")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+                          
    labs(x="Repeat length",y="Count (x1,000)",fill=NULL)
                                      
dev.off()
                                      
################################                                      
                                      
# Plot motif length distribution
#vars_repeats <- unique(vars[FILTER=="PASS" & SVTYPE %in% c("INS","DEL") & REPEAT_STATUS %in% c("STR","TR"),.(ID,SVTYPE,LEN_LABEL,SVLEN,EXON_INTERSECT_GENE_ID,INTRON_INTERSECT_GENE_ID,TANDEM_REPEAT_STATUS,REPEAT_STATUS,REPEAT_LABEL,MAJOR_MOTIF,TRF_ELEMENT_SIZE,TRF_INTERSECTION,FRAC_SV_COVERED)])
#vars_repeats[MAJOR_MOTIF!="100plus",MOTIF_SIZE:=nchar(MAJOR_MOTIF)]
#vars_repeats[is.na(MOTIF_SIZE),MOTIF_SIZE:=100]
                                      
vars_repeats$SVTYPE <- factor(vars_repeats$SVTYPE,levels=c("DEL","INS"))
vars_repeats[EXON_INTERSECT_GENE_ID!="0" & INTRON_INTERSECT_GENE_ID!="0",INTRON_INTERSECT_GENE_ID:="0"]
#vars_repeats[,STR_LABEL:=MOTIF_SIZE] 
                                      
vars_repeats_long <- data.table(melt(vars_repeats[,.(ID,SVTYPE,SVLEN,LEN_LABEL,REPEAT_STATUS,STR_LABEL,EXON_INTERSECT_GENE_ID,INTRON_INTERSECT_GENE_ID)],id.vars=c("ID","LEN_LABEL","SVLEN","SVTYPE","REPEAT_STATUS","STR_LABEL"),variable.name="FEATURE",value.name="ANNOTATION"))
vars_repeats_long[ANNOTATION=="0",OVERLAP:="No"]
vars_repeats_long[ANNOTATION!="0",OVERLAP:="Yes"]

vars_repeats_long[FEATURE=="EXON_INTERSECT_GENE_ID",FEATURE:="CDS"]                                             
vars_repeats_long[FEATURE=="INTRON_INTERSECT_GENE_ID",FEATURE:="UTR/Intron/2kb flanking"]                                      

vars_repeats_long <- rbind(vars_repeats_long[OVERLAP=="Yes"],vars_repeats[ID %in% vars_repeats_long[OVERLAP=="No",.N,by=list(LEN_LABEL,ID)][N==2,ID],.(ID,LEN_LABEL,SVLEN,SVTYPE,REPEAT_STATUS,STR_LABEL,FEATURE="Intergenic",ANNOTATION="0",OVERLAP="Yes")]) 

vars_repeats_count <- vars_repeats_long[,.N,by=list(SVTYPE,REPEAT_STATUS,STR_LABEL)]                                      
vars_repeats_count <- merge(vars_repeats_count,vars_repeats_count[,.(TOTAL=sum(N)),by=list(SVTYPE,REPEAT_STATUS,STR_LABEL)])                                      

# Plot frequency of different STRs                                    
pdf('/path/to/Figures/main_figures_v2/figure6/motif_length_STR_count.absolute_counts.INS.pdf',height=6,width=12)
                                      
ggplot(vars_repeats_count[REPEAT_STATUS=="STR"][SVTYPE=="INS"],aes(STR_LABEL,N,fill=STR_LABEL))+
    geom_bar(stat="identity")+
    scale_fill_manual(values=colors_str_tr)+                                      
    facet_wrap(~SVTYPE)+                                  
    theme_bw(base_size=20)+
    theme(aspect.ratio=1/2)+
    guides(fill="none")+
    scale_y_continuous(limits=c(0,20000),labels = function(l) { trans = l / 1000},breaks=c(0,10000,20000))+
    theme(axis.text.x=element_text(angle=45,hjust=1))+                                  
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x=NULL,y="Counts (x1,000)",fill=NULL)
                                      
dev.off()
                                      
pdf('/path/to/Figures/main_figures_v2/figure6/motif_length_STR_count.absolute_counts.DEL.pdf',height=6,width=12)
                                      
ggplot(vars_repeats_count[REPEAT_STATUS=="STR"][SVTYPE=="DEL"],aes(STR_LABEL,N,fill=STR_LABEL))+
    geom_bar(stat="identity")+
    scale_fill_manual(values=colors_str_tr)+                                      
    facet_wrap(~SVTYPE)+                                  
    theme_bw(base_size=20)+
    theme(aspect.ratio=1/2)+
    guides(fill="none")+
    scale_y_continuous(limits=c(0,20000),labels = function(l) { trans = l / 1000},breaks=c(0,10000,20000))+
    theme(axis.text.x=element_text(angle=45,hjust=1))+                                  
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x=NULL,y="Counts (x1,000)",fill=NULL)
                                      
dev.off() 
                                      
# Plot motif length relative distribution (CDS, UTR/Intron/2kb flanking, Intergenic)
                                      
vars_repeats$SVTYPE <- factor(vars_repeats$SVTYPE,levels=c("DEL","INS"))
vars_repeats[EXON_INTERSECT_GENE_ID!="0" & INTRON_INTERSECT_GENE_ID!="0",INTRON_INTERSECT_GENE_ID:="0"]                                
vars_repeats_long <- data.table(melt(vars_repeats[,.(ID,SVTYPE,SVLEN,LEN_LABEL,REPEAT_STATUS,STR_LABEL,EXON_INTERSECT_GENE_ID,INTRON_INTERSECT_GENE_ID)],id.vars=c("ID","LEN_LABEL","SVLEN","SVTYPE","REPEAT_STATUS","STR_LABEL"),variable.name="FEATURE",value.name="ANNOTATION"))
vars_repeats_long[ANNOTATION=="0",OVERLAP:="No"]
vars_repeats_long[ANNOTATION!="0",OVERLAP:="Yes"]

vars_repeats_long[FEATURE=="EXON_INTERSECT_GENE_ID",FEATURE:="CDS"]                                             
vars_repeats_long[FEATURE=="INTRON_INTERSECT_GENE_ID",FEATURE:="UTR/Intron/2kb flanking"]                                      

vars_repeats_long <- rbind(vars_repeats_long[OVERLAP=="Yes"],vars_repeats[ID %in% vars_repeats_long[OVERLAP=="No",.N,by=list(LEN_LABEL,ID)][N==2,ID],.(ID,LEN_LABEL,SVLEN,SVTYPE,REPEAT_STATUS,STR_LABEL,FEATURE="Intergenic",ANNOTATION="0",OVERLAP="Yes")]) 

vars_repeats_count <- vars_repeats_long[,.N,by=list(SVTYPE,REPEAT_STATUS,STR_LABEL,FEATURE)]                                      
vars_repeats_count <- merge(vars_repeats_count,vars_repeats_count[,.(TOTAL=sum(N)),by=list(SVTYPE,REPEAT_STATUS,STR_LABEL)])                                        
                                      
pdf('/path/to/Figures/main_figures_v2/figure6/motif_length_STR_count.relative_frequency.INS.pdf',height=6,width=12)
                                      
ggplot(vars_repeats_count[REPEAT_STATUS=="STR"][SVTYPE=="INS"],aes(STR_LABEL,N/TOTAL*100,fill=FEATURE))+
    geom_bar(stat="identity")+
    facet_wrap(~SVTYPE)+                                  
    theme_bw(base_size=20)+
    theme(aspect.ratio=1/2)+
    scale_y_continuous(breaks=c(0,50,100))+ 
    scale_fill_manual(values=colors_gene_overlap)+                                  
    theme(axis.text.x=element_text(angle=45,hjust=1))+                                  
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(legend.position="bottom")+
    labs(x=NULL,y="Percentage (%)",fill=NULL)
                                      
dev.off()                                        
  
# Statistics

# Proportion of intergenic STRs                                       
vars_repeats_long_stat <- vars_repeats_long[REPEAT_STATUS=="STR",.N,by=FEATURE]                                      
vars_repeats_long_stat[,TOTAL:=sum(vars_repeats_long_stat[,N])]                                      
vars_repeats_long_stat[,PERCENT:=N/TOTAL*100]
                                      
# Proportion of intronic pentanucleotides
vars_repeats_count[REPEAT_STATUS=="STR"][SVTYPE=="INS" & FEATURE=="UTR/Intron/2kb flanking",.(STR_LABEL,FEATURE,PERCENT=N/TOTAL*100)]    

# Mean and SD for other STRs                                      
vars_repeats_count[REPEAT_STATUS=="STR"][SVTYPE=="INS" & FEATURE=="UTR/Intron/2kb flanking" & STR_LABEL!="STR - (5nt)",.(STR_LABEL,FEATURE,PERCENT=N/TOTAL*100)][,.(Mean=mean(PERCENT),SD=sd(PERCENT))]                              
                                      
pdf('/path/to/Figures/main_figures_v2/figure6/motif_length_STR_count.relative_frequency.DEL.pdf',height=6,width=12)
                                      
ggplot(vars_repeats_count[REPEAT_STATUS=="STR"][SVTYPE=="DEL"],aes(STR_LABEL,N/TOTAL*100,fill=FEATURE))+
    geom_bar(stat="identity")+
    facet_wrap(~SVTYPE)+                                  
    theme_bw(base_size=20)+
    theme(aspect.ratio=1/2)+
    scale_y_continuous(breaks=c(0,50,100))+ 
    scale_fill_manual(values=colors_gene_overlap)+                                  
    theme(axis.text.x=element_text(angle=45,hjust=1))+                                  
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(legend.position="bottom")+
    labs(x=NULL,y="Percentage (%)",fill=NULL)
                                      
dev.off()  
                                      
# Plot motif length relative distribution (CDS)
                                      
vars_repeats_count <- vars_repeats_long[,.N,by=list(SVTYPE,REPEAT_STATUS,STR_LABEL,FEATURE)]                                      
vars_repeats_count <- merge(vars_repeats_count,vars_repeats_count[,.(TOTAL=sum(N)),by=list(SVTYPE,REPEAT_STATUS,STR_LABEL)])
vars_repeats_count <- rbind(vars_repeats_count,
                            data.table(SVTYPE=c("INS","INS"),REPEAT_STATUS=c("STR","STR"),STR_LABEL=c("STR - (4nt)","STR - (4nt)"),FEATURE=c("CDS","CDS"),N=c(0,0),TOTAL=c(1,1)),
                            data.table(SVTYPE=c("INS","INS"),REPEAT_STATUS=c("STR","STR"),STR_LABEL=c("STR - (2nt)","STR - (2nt)"),FEATURE=c("CDS","CDS"),N=c(0,0),TOTAL=c(1,1)),
                            data.table(SVTYPE=c("DEL","DEL"),REPEAT_STATUS=c("STR","STR"),STR_LABEL=c("STR - (4nt)","STR - (4nt)"),FEATURE=c("CDS","CDS"),N=c(0,0),TOTAL=c(1,1)),
                            data.table(SVTYPE=c("DEL","DEL"),REPEAT_STATUS=c("STR","STR"),STR_LABEL=c("STR - (2nt)","STR - (2nt)"),FEATURE=c("CDS","CDS"),N=c(0,0),TOTAL=c(1,1)),
                            data.table(SVTYPE=c("DEL","DEL"),REPEAT_STATUS=c("STR","STR"),STR_LABEL=c("STR - (5nt)","STR - (5nt)"),FEATURE=c("CDS","CDS"),N=c(0,0),TOTAL=c(1,1)))

pdf('/path/to/Figures/main_figures_v2/figure6/motif_length_STR_count.relative_frequency.CDS.INS.pdf',height=6,width=12)
                                      
ggplot(vars_repeats_count[REPEAT_STATUS=="STR" & FEATURE=="CDS"][SVTYPE=="INS"],aes(STR_LABEL,N/TOTAL*100,fill=STR_LABEL))+
    geom_bar(stat="identity")+
    facet_wrap(~SVTYPE)+                                  
    theme_bw(base_size=20)+
    theme(aspect.ratio=1/2)+
    scale_fill_manual(values=colors_str_tr)+   
    guides(fill="none")+

    theme(axis.text.x=element_text(angle=45,hjust=1))+                                  
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x=NULL,y="Percentage (%)",fill=NULL)   

dev.off()

pdf('/path/to/Figures/main_figures_v2/figure6/motif_length_STR_count.relative_frequency.CDS.DEL.pdf',height=6,width=12)                                      
                                      
ggplot(vars_repeats_count[REPEAT_STATUS=="STR" & FEATURE=="CDS"][SVTYPE=="DEL"],aes(STR_LABEL,N/TOTAL*100,fill=STR_LABEL))+
    geom_bar(stat="identity")+
    facet_wrap(~SVTYPE)+                                  
    theme_bw(base_size=20)+
    theme(aspect.ratio=1/2)+
    scale_fill_manual(values=colors_str_tr)+   
    guides(fill="none")+

    theme(axis.text.x=element_text(angle=45,hjust=1))+                                  
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x=NULL,y="Percentage (%)",fill=NULL)  
                                      
dev.off() 
                                      
# Frequency of triplets

triplets <- vars_repeats[MOTIF_SIZE==3]
triplets[EXON_INTERSECT_GENE_ID==0 & INTRON_INTERSECT_GENE_ID==0,FEATURE:="Intergenic"]
triplets[INTRON_INTERSECT_GENE_ID!=0,FEATURE:="UTR/Intron/2kb flanking"]
triplets[EXON_INTERSECT_GENE_ID!=0,FEATURE:="CDS"]
                                      
triplets_counts <- triplets[,.N,by=list(SVTYPE,MOTIF_LABEL)]                                      

pathogenic_motifs <- sort(c('TGC','CAG','GCC','TAG','GAA'))
rolled_pathogenic_motifs <- rbindlist(lapply(seq_along(pathogenic_motifs),function(x) roll_seq(pathogenic_motifs[x],x)))                                      

#triplets_counts[(MOTIF_LABEL %in% rolled_pathogenic_motifs$SEQ) | (MOTIF_LABEL %in% rolled_pathogenic_motifs$ROLLED),STATUS:="pathogenic"]                       
#triplets_counts[is.na(STATUS),STATUS:="non-pathogenic"]                      

pdf('/path/to/Figures/main_figures_v2/figure6/frequency_of_triplets.INS.pdf',height=6,width=12)                                        
                                      
ggplot(triplets_counts[SVTYPE=="INS"],aes(reorder(MOTIF_LABEL,-N),N))+
    geom_bar(stat="identity")+
    #facet_wrap(~FEATURE)+     
    #scale_y_continuous(breaks=c(0,300,600))+                                  
    theme_bw(base_size=20)+
    theme(aspect.ratio=1/2)+
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+                                  
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(legend.position="bottom")+                                         
    labs(x="Trinucleotides",y="Counts",fill=NULL)
                                      
dev.off()                                              
                                             
                                             
#pdf('/path/to/Figures/main_figures_v2/figure6/frequency_of_triplets.non_cds.INS.pdf',height=6,width=12)                                        
                                      
#ggplot(triplets_counts[SVTYPE=="INS" & FEATURE=="UTR/Intron/2kb flanking"],aes(reorder(MOTIF_LABEL,-N),N,fill=STATUS))+
#    geom_bar(stat="identity")+
#    facet_wrap(~FEATURE)+     
#    #scale_y_continuous(breaks=c(0,300,600))+                                  
#    theme_bw(base_size=20)+
#    theme(aspect.ratio=1/2)+
#    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+                                  
#    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#    theme(legend.position="bottom")+                                         
#    labs(x="Trinucleotides",y="Counts",fill=NULL)
                                      
#dev.off() 

#pdf('/path/to/Figures/main_figures_v2/figure6/frequency_of_triplets.intergenic.INS.pdf',height=6,width=12)                                                                                                                                  
#ggplot(triplets_counts[SVTYPE=="INS" & FEATURE=="Intergenic"],aes(reorder(MOTIF_LABEL,-N),N,fill=STATUS))+
#    geom_bar(stat="identity")+
#    facet_wrap(~FEATURE)+     
#    #scale_y_continuous(breaks=c(0,300,600))+                                  
#    theme_bw(base_size=20)+
#    theme(aspect.ratio=1/2)+
#    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+                                  
#    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#    theme(legend.position="bottom")+                                         
#    labs(x="Trinucleotides",y="Counts",fill=NULL)                                             

#dev.off()                                             
                                             
#pdf('/path/to/Figures/main_figures_v2/figure6/frequency_of_triplets.non_cds.DEL.pdf',height=6,width=12)                                        
                                      
#ggplot(triplets_counts[SVTYPE=="DEL" & FEATURE=="UTR/Intron/2kb flanking"],aes(reorder(MOTIF_LABEL,-N),N,fill=STATUS))+
#    geom_bar(stat="identity")+
#    facet_wrap(~FEATURE)+     
#    #scale_y_continuous(breaks=c(0,300,600))+                                  
#    theme_bw(base_size=20)+
#    theme(aspect.ratio=1/2)+
#    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+                                  
#    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#    theme(legend.position="bottom")+                                         
#    labs(x="Trinucleotides",y="Counts",fill=NULL)
                                      
#dev.off() 

#pdf('/path/to/Figures/main_figures_v2/figure6/frequency_of_triplets.intergenic.DEL.pdf',height=6,width=12)                                                                                                                                  
#ggplot(triplets_counts[SVTYPE=="DEL" & FEATURE=="Intergenic"],aes(reorder(MOTIF_LABEL,-N),N,fill=STATUS))+
#    geom_bar(stat="identity")+
#    facet_wrap(~FEATURE)+     
#    #scale_y_continuous(breaks=c(0,300,600))+                                  
#    theme_bw(base_size=20)+
#    theme(aspect.ratio=1/2)+
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+                                  
#    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#    theme(legend.position="bottom")+                                         
#    labs(x="Trinucleotides",y="Counts",fill=NULL)                                             

#dev.off()    
                                             
# Frequency of pentamers
                                             
pentamers <- vars_repeats[MOTIF_SIZE==5]
pentamers[EXON_INTERSECT_GENE_ID==0 & INTRON_INTERSECT_GENE_ID==0,FEATURE:="Intergenic"]
pentamers[INTRON_INTERSECT_GENE_ID!=0,FEATURE:="UTR/Intron/2kb flanking"]
pentamers[EXON_INTERSECT_GENE_ID!=0,FEATURE:="CDS"]
pentamers_counts <- pentamers[MOTIF_SIZE==5,.N,by=list(SVTYPE,MOTIF_LABEL)]                                      

#pathogenic_motifs <- sort(c('ATTTT','ATTTC','TGGAA','AAAAG','AAAGG','AAGGG','ACAGG'))
#rolled_pathogenic_motifs <- rbindlist(lapply(seq_along(pathogenic_motifs),function(x) roll_seq(pathogenic_motifs[x],x)))                                      

#pentamers_counts[(MOTIF_LABEL %in% rolled_pathogenic_motifs$SEQ) | (MOTIF_LABEL %in% rolled_pathogenic_motifs$ROLLED),STATUS:="pathogenic"]                       
#pentamers_counts[is.na(STATUS),STATUS:="non-pathogenic"]
 
pdf('/path/to/Figures/main_figures_v2/figure6/frequency_of_pentamers.INS.pdf',height=6,width=12)                                        
                                      
ggplot(pentamers_counts[SVTYPE=="INS"],aes(reorder(MOTIF_LABEL,-N),log2(N)))+
    geom_bar(stat="identity")+
    #facet_wrap(~FEATURE)+     
    #scale_y_continuous(breaks=c(0,300,600))+                                  
    theme_bw(base_size=20)+
    theme(aspect.ratio=1/2)+
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=7))+                                  
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(legend.position="bottom")+                                         
    labs(x="Pentanucleotides",y="Counts",fill=NULL)
                                      
dev.off()                                              
                                             
#pdf('/path/to/Figures/main_figures_v2/figure6/frequency_of_pentamers.non_cds.INS.pdf',height=6,width=12)                                        
                                      
#ggplot(pentamers_counts[SVTYPE=="INS" & FEATURE=="UTR/Intron/2kb flanking"],aes(reorder(MOTIF_LABEL,-N),N,fill=STATUS))+
#    geom_bar(stat="identity")+
#    facet_wrap(~FEATURE)+     
#    #scale_y_continuous(breaks=c(0,300,600))+                                  
#    theme_bw(base_size=20)+
#    theme(aspect.ratio=1/2)+
#    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=7))+                                  
#    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#    theme(legend.position="bottom")+                                         
#    labs(x="Pentanucleotides",y="Counts",fill=NULL)
                                      
#dev.off() 

#pdf('/path/to/Figures/main_figures_v2/figure6/frequency_of_pentamers.intergenic.INS.pdf',height=6,width=12)                                                                                                                                  
#ggplot(pentamers_counts[SVTYPE=="INS" & FEATURE=="Intergenic"],aes(reorder(MOTIF_LABEL,-N),N,fill=STATUS))+
#    geom_bar(stat="identity")+
#    facet_wrap(~FEATURE)+     
#    #scale_y_continuous(breaks=c(0,300,600))+                                  
#    theme_bw(base_size=20)+
#    theme(aspect.ratio=1/2)+
#    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=7))+                                  
#    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#    theme(legend.position="bottom")+                                         
#    labs(x="Pentanucleotides",y="Counts",fill=NULL)                                             

#dev.off()                                             
                                             
#pdf('/path/to/Figures/main_figures_v2/figure6/frequency_of_pentamers.non_cds.DEL.pdf',height=6,width=12)                                        
                                      
#ggplot(pentamers_counts[SVTYPE=="DEL" & FEATURE=="UTR/Intron/2kb flanking"],aes(reorder(MOTIF_LABEL,-N),N,fill=STATUS))+
#    geom_bar(stat="identity")+
#    facet_wrap(~FEATURE)+     
#    #scale_y_continuous(breaks=c(0,300,600))+                                  
#    theme_bw(base_size=20)+
#    theme(aspect.ratio=1/2)+
#    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=7))+                                  
#    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#    theme(legend.position="bottom")+                                         
#    labs(x="Pentanucleotides",y="Counts",fill=NULL)
                                      
#dev.off() 

#pdf('/path/to/Figures/main_figures_v2/figure6/frequency_of_pentamers.intergenic.DEL.pdf',height=6,width=12)                                                                                                                                  
#ggplot(pentamers_counts[SVTYPE=="DEL" & FEATURE=="Intergenic"],aes(reorder(MOTIF_LABEL,-N),N,fill=STATUS))+
#    geom_bar(stat="identity")+
#    facet_wrap(~FEATURE)+     
#    #scale_y_continuous(breaks=c(0,300,600))+                                  
#    theme_bw(base_size=20)+
#    theme(aspect.ratio=1/2)+
#    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=7))+                                  
#    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#    theme(legend.position="bottom")+                                         
#    labs(x="Pentanucleotides",y="Counts",fill=NULL)                                             

#dev.off()                                             
           
                                             
# Significantly variable STR expansions

alleles <- fread('/path/to/str_analysis_chm13/FULL_PARALLEL/all_allele_sequences.v2.tsv')
alleles[,unique_id:=paste(chrom,start,end,id,sep=";")]
#alleles[str_detect(unique_id,"ATXN3"),sequence:=rev.comp(sequence)]
alleles[,length:=nchar(sequence)-100]
alleles[,RunID:=sample]
communities <- fread('/path/to/Figures/cohort_sex_info.v2.indexed.tab')
alleles <- merge(alleles,communities,by="RunID") 
                                             
alleles_anova <- alleles[,.(pvalue=summary(aov(length~Community))[[1]][[5]][1]),by=unique_id]
alleles_anova[,adjusted_pvalue:=p.adjust(pvalue,method="fdr")]
                                             
alleles_group <- alleles[,.(unique_id,RunID,allele,length,Community)][,.(mean=mean(length),sd=sd(length)),by=list(unique_id,Community)]
alleles_group_significant <- alleles_group[unique_id %in% alleles_anova[pvalue<0.05,unique_id]]
alleles_group_significant <- merge(alleles_group_significant,alleles_group_significant[,.(min=min(sd),max=max(sd)),by=unique_id])
alleles_group_significant[,scaled:=(sd-min)/(max-min)]

alleles_group_wide <- data.table(dcast(alleles_group_significant,unique_id~Community,value.var="scaled"))
alleles_group_wide[,id:=1:.N]
alleles_group_wide_matrix <- as.matrix(alleles_group_wide[,.(NCIGP1,NCIGP2,NCIGP3,NCIGP4,`non-NCIG`)])
row.names(alleles_group_wide_matrix) <- alleles_group_wide$id
                                             
pdf('/path/to/Figures/main_figures_v2/figure6/signficantly_variable_str_sites.pdf',height=12,width=8)

heatmap.2(
    alleles_group_wide_matrix,Colv=NA,col=colorRampPalette(c("blue3", "white", "red3")),
    Rowv=as.dendrogram(hclust(dist(as.matrix(alleles_group_wide[,.(NCIGP1,NCIGP2,NCIGP3,NCIGP4,`non-NCIG`)])))),
    trace="none",key=TRUE,density.info="none",keysize = 1.5,lhei = c(1,8),cexCol=1,cexRow=0.25
 )

dev.off()
                                             
# Statistics

# Number of STR expansions that we investigated                                             
dim(unique(alleles[,.(unique_id,id)]))[1]
                                             
# Number of significantly varied sites
dim(alleles_group_wide_matrix)[1]        
                                             
# Non-ncig more diverse
dim(alleles_group_wide[`non-NCIG`==1])[1]                                             
                                             
# Ncig more diverse
dim(alleles_group_wide_matrix)[1] - dim(alleles_group_wide[`non-NCIG`==1])[1]                                             
