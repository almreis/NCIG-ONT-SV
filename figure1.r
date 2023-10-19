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

# Functions

addline_format <- function(x,...){
    gsub('\\s','\n',x)
}

make_pie_chart <- function(slices,labels,community) {
    dir <- '/path/to/Figures/main_figures_v2/figure1/'
    out <- sprintf('%s/%s_sex_pie_chart.pdf',dir,community)
    n <- sum(slices)
    title <- sprintf('n = %d',n)
    pdf(out,height=6,width=6)
    print(pie(slices, labels = labels, main = title))
    dev.off()
}

# Colors

colors <- c('0'='#8DD3C7','HOMO'='#BEBADA','STR'='#FB8072','TR'='#80B1D3','LTR'='#FDB462','LINE'='#B3DE69','SINE'='#FCCDE5','Retroposon'='#D9D9D9','DNA'='#FFFFB3')
colors_mobile <- c('0'='#8DD3C7','STR'='#FB8072','TR'='#80B1D3','ME - fragment'='#E0BBE4','ME - complete'='#957DAD')
colors_sharedness <- c('singleton'='#E2CFFF','polymorphic'='#7ADEA8','major'='#E995BA','shared'='#3DC1CB')
colors_range <- c('Aus'='#9292D1','Aus-absent'='#D1F0A4','Global'='#C6E8EE')
colors_novelty <- c('Found - both'='#E0D0F5','Found - Gnomad only'='#FCE4BA','Found - DECODE only'='#E1B894','Not found'='#F3AFCC','Not lifted'='#B7CAB4')
colors_novelty_extended <- c('Found - both'='#E0D0F5','Found - Gnomad only'='#FCE4BA','Found - DECODE only'='#E1B894','Not found'='#F3AFCC','Deleted'='#7ECCC5','Partially deleted'='#91F2FF','Split'='#83A1FF')
colors_community <- c('NCIGP1'='#F07DD1','NCIGP2'='#986BE3','NCIGP3'='#8DBAF3','NCIGP4'='#6DB3AC','non-NCIG'='#FCB48D')
colors_aus <- c('Community-specific'='#F4C573','Widespread'='#D95050','Private'='#7ECBE6','Ancestral'='#82C493')                                             
colors_discovery <- c('nonredundant'='#F7ECB2','shared'='#3DC1CB')
colors_community_extended <- c('NCIGP1'='#F07DD1','NCIGP2'='#986BE3','NCIGP3'='#8DBAF3','NCIGP4'='#6DB3AC','non-NCIG'='#FCB48D','NCIG'='#64369C')
colors_repeats_summary <- c('Non-repetitive'='#8DD3C7','Tandem repeats'='#E88794','Mobile elements'='#AA8FC4')                                       
colors_gene_overlap <- c('CDS'='#9BAB65','UTR/Intron/2kb flanking'='#F0AA4F')    
colors_aus_extended <- c('Community-specific'='#F4C573','Widespread'='#D95050','Private'='#7ECBE6','Ancestral'='#82C493','Aus-absent'='#D1F0A4','Global'='#C6E8EE') 
colors_genome <- c('HG38'='#789CBB','CHM13'='#D9ADB6')
colors_svtype <- c('DEL'='#EA8975','INS'='#F7D08D','INV'='#5567AA','DUP'='#A2DEB0','BND'='#4BADC4')

#############################################################

cohort <- fread('/path/to/Figures/cohort_sex_info.v2.tab',sep="\t",header=FALSE)
cohort <- cohort[,.(GUAID=V2,Sex=V3,RunID=V4,Community=V5)]
cohort[Community=="TOB",Community:="non-NCIG"]

# Make pie charts show sex representation for each community

cohort_count <- cohort[,.N,by=list(Community,Sex)]
cohort_count <- cohort_count[order(Community,Sex)]

lapply(unique(cohort_count$Community), function(x) make_pie_chart(cohort_count[Community==x,N],c('XX','XY'),x))

# Relationship between coverage and read length

coverage <- fread('/path/to/Figures/main_figures_v2/figure1/coverage_stats.tab',sep="\t",header=TRUE)
coverage <- merge(coverage,cohort[,.(GUAID,RunID,Community)])
coverage <- coverage[order(Community,L1kb)]
coverage[,Index:=1:.N]

coverage_long <- data.table(melt(coverage[,c("Index","Community",names(coverage)[startsWith(names(coverage),"L")]),with=FALSE],id.vars=c("Index","Community")))
coverage_long[,variable:=str_remove(variable,"kb")]
coverage_long[,variable:=str_remove(variable,"L")]
coverage_long[,variable:=as.numeric(variable)]

pdf('/path/to/Figures/main_figures_v2/figure1/coverage_v_readlength.pdf',height=6,width=12)

ggplot(coverage_long,aes(variable,value,by=as.factor(Index),color=Community))+
    geom_path(alpha=0.6)+
    facet_wrap(~Community,nrow=1)+
    scale_color_manual(values=colors_community)+
    scale_y_continuous(breaks=c(0,25,50))+
    scale_x_continuous(breaks=c(0,15,30))+
    geom_hline(yintercept = 9.5,linetype="dashed",size=0.2)+
    geom_vline(xintercept = 5,linetype="dashed",size=0.2)+
    theme_bw(base_size=20)+
    theme(aspect.ratio=1,legend.position="bottom")+
    theme(strip.text.x = element_blank())+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x="Minimal read length (Kb)",y="Average coverage",fill="Group")

dev.off()

#############################################################

# Average MAPQ

# HG001 + CHM13 + IL
dt1 <- fread('/path/to/data/alignments/bwamem2/2.2.1/chm13/HG001/20221015/HG001_illumina.mapq_fixed_windows.tab')
dt1[,LABEL:="HG001_Illumina_CHM13"]
# HG001 + HG38 + IL
dt2 <- fread('/path/to/data/alignments/bwamem2/2.2.1/hg38/HG001/20221015/HG001_illumina.mapq_fixed_windows.tab')
dt2[,LABEL:="HG001_Illumina_HG38"]

# HG002 + CHM13 + IL
dt3 <- fread('/path/to/data/alignments/bwamem2/2.2.1/chm13/HG002/20221015/HG002_illumina.mapq_fixed_windows.tab')
dt3[,LABEL:="HG002_Illumina_CHM13"]
# HG002 + HG38 + IL
dt4 <- fread('/path/to/data/alignments/bwamem2/2.2.1/hg38/HG002/20221015/HG002_illumina.mapq_fixed_windows.tab')
dt4[,LABEL:="HG002_Illumina_HG38"]

# HG001 + CHM13 + ONT
dt5 <- fread('/path/to/data/alignments/minimap2/v2.22-r1101/<sample_id>/20220428/PLPN239131_pass.mapq_fixed_windows.tab')
dt5[,LABEL:="HG001_ONT_CHM13"]

# HG001 + HG38 + ONT
dt6 <- fread('/path/to/data/alignments/minimap2/v2.22-r1101/hg38/<sample_id>/20221013/PLPN239131_pass.mapq_fixed_windows.tab')
dt6[,LABEL:="HG001_ONT_HG38"]

# HG002 + CHM13 + ONT
dt7 <- fread('/path/to/data/alignments/minimap2/v2.22-r1101/<sample_id>/20220428/<sample_id>_pass.mapq_fixed_windows.tab')
dt7[,LABEL:="HG002_ONT_CHM13"]

# HG002 + HG38 + ONT
dt8 <- fread('/path/to/data/alignments/minimap2/v2.22-r1101/hg38/<sample_id>/20221013/<sample_id>_pass.mapq_fixed_windows.tab')
dt8[,LABEL:="HG002_ONT_HG38"]

all <- rbind(dt1,dt2,dt3,dt4,dt5,dt6,dt7,dt8)

all[as.numeric(V4)<=5 & as.numeric(V4)==0,MAPQ:="0-5"]
all[as.numeric(V4)>5 & as.numeric(V4)!=0,MAPQ:="1-29"]
all[as.numeric(V4)>=30,MAPQ:="30-60"]

mapq_count <- all[!is.na(MAPQ),.N,by=list(LABEL,MAPQ)]

mapq_count <- merge(mapq_count,mapq_count[,.(TOTAL=sum(N)),by=LABEL])

mapq_count[,c("SAMPLE","PLATFORM","GENOME"):=tstrsplit(LABEL,"_",fixed=TRUE)]

mapq_count[,FACET:=paste(SAMPLE,PLATFORM," ")]

pdf('/path/to/Figures/main_figures_v2/figure1/percentage_mapq0.pdf',height=8,width=8)

ggplot(mapq_count[MAPQ=="0-5"],aes(FACET,N/TOTAL*100,fill=GENOME))+
    geom_bar(stat="identity",position=position_dodge())+
    scale_fill_manual(values=colors_genome)+
    theme_bw(base_size=20)+
    theme(aspect.ratio=1)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    scale_y_continuous(breaks=c(0,2,4))+
    scale_x_discrete(breaks=unique(mapq_count$FACET),labels=addline_format(unique(mapq_count$FACET)))+
    theme(legend.position="bottom")+
    labs(x=NULL,y="Percentage (%)",fill="Reference")

dev.off()
       
# Regions of low mappability       
mapq_count[MAPQ=="0-5",.(LABEL,N/TOTAL*100)]

#############################################################

# Percentage of genome with 0 coverage

# HG001 + CHM13 + IL
dt1 <- fread('/path/to/data/alignments/bwamem2/2.2.1/chm13/HG001/20221015/HG001_illumina.bedcov_fixed_windows.tab')
dt1[,LABEL:="HG001_Illumina_CHM13"]
# HG001 + HG38 + IL
dt2 <- fread('/path/to/data/alignments/bwamem2/2.2.1/hg38/HG001/20221015/HG001_illumina.bedcov_fixed_windows.tab')
dt2[,LABEL:="HG001_Illumina_HG38"]

# HG002 + CHM13 + IL
dt3 <- fread('/path/to/data/alignments/bwamem2/2.2.1/chm13/HG002/20221015/HG002_illumina.bedcov_fixed_windows.tab')
dt3[,LABEL:="HG002_Illumina_CHM13"]
# HG002 + HG38 + IL
dt4 <- fread('/path/to/data/alignments/bwamem2/2.2.1/hg38/HG002/20221015/HG002_illumina.bedcov_fixed_windows.tab')
dt4[,LABEL:="HG002_Illumina_HG38"]

# HG001 + CHM13 + ONT
dt5 <- fread('/path/to/data/alignments/minimap2/v2.22-r1101/<sample_id>/20220428/PLPN239131_pass.bedcov_fixed_windows.tab')
dt5[,LABEL:="HG001_ONT_CHM13"]

# HG001 + HG38 + ONT
dt6 <- fread('/path/to/data/alignments/minimap2/v2.22-r1101/hg38/<sample_id>/20221013/PLPN239131_pass.bedcov_fixed_windows.tab')
dt6[,LABEL:="HG001_ONT_HG38"]

# HG002 + CHM13 + ONT
dt7 <- fread('/path/to/data/alignments/minimap2/v2.22-r1101/<sample_id>/20220428/<sample_id>_pass.bedcov_fixed_windows.tab')
dt7[,LABEL:="HG002_ONT_CHM13"]

# HG002 + HG38 + ONT
dt8 <- fread('/path/to/data/alignments/minimap2/v2.22-r1101/hg38/<sample_id>/20221013/<sample_id>_pass.bedcov_fixed_windows.tab')
dt8[,LABEL:="HG002_ONT_HG38"]

all <- rbind(dt1,dt2,dt3,dt4,dt5,dt6,dt7,dt8)

coverage_count <- all[,.(COVERED=sum(V5),UNCOVERED=sum(V6-V5)),by=LABEL]

coverage_count <- data.table(melt(coverage_count,id.vars="LABEL"))

coverage_count <- merge(coverage_count,coverage_count[,.(TOTAL=sum(value)),by=LABEL])

coverage_count[,c("SAMPLE","PLATFORM","GENOME"):=tstrsplit(LABEL,"_",fixed=TRUE)]

coverage_count[,FACET:=paste(SAMPLE,PLATFORM," ")]

pdf('/path/to/Figures/main_figures_v2/figure1/percentage_genome_0_coverage.pdf',height=8,width=8)

ggplot(coverage_count[variable=="UNCOVERED"],aes(FACET,value/TOTAL*100,fill=GENOME))+
    geom_bar(stat="identity",position=position_dodge())+
    scale_fill_manual(values=colors_genome)+
    theme_bw(base_size=20)+
    theme(aspect.ratio=1)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    scale_y_continuous(breaks=c(0,6,12))+
    scale_x_discrete(breaks=unique(coverage_count$FACET),labels=addline_format(unique(coverage_count$FACET)))+
    theme(legend.position="bottom")+
    labs(x=NULL,y="Percentage (%)",fill="Reference")

dev.off()
       
# Percent of regions of 0 coverage
       
coverage_count[variable=="UNCOVERED",.(LABEL,value/TOTAL*100)]
       
# Number of additional Mbases
       
additional_bases <- merge(coverage_count[variable=="UNCOVERED",.(LABEL,ZERO_COV=value,TOTAL)],mapq_count[MAPQ=="0-5",.(LABEL,LOW_MAPQ=N*1000)])
additional_bases[,ALIGNABLE:=(TOTAL-(ZERO_COV+LOW_MAPQ))]
additional_bases_hg001 <- (additional_bases[LABEL=="HG001_ONT_CHM13",ALIGNABLE]-additional_bases[LABEL=="HG001_ONT_HG38",ALIGNABLE])/10^6
additional_bases_hg002 <- (additional_bases[LABEL=="HG002_ONT_CHM13",ALIGNABLE]-additional_bases[LABEL=="HG002_ONT_HG38",ALIGNABLE])/10^6       
#############################################################

# Number and types of variants identified

extract_vcf <- function(dt) {
    dt[,SVTYPE:=str_extract(INFO,"SVTYPE=([A-Z]{3})")]
    dt[,SVTYPE:=str_remove(SVTYPE,"SVTYPE=")]
    dt[,SVLEN:=str_extract(INFO,"SVLEN=([-]?[0-9]+);")]
    dt[,SVLEN:=str_remove(SVLEN,"SVLEN=")]
    dt[,SVLEN:=abs(as.numeric(str_remove(SVLEN,";")))]
    dt <- dt[,.(ID,SVTYPE,SVLEN)]
}

# HG001 + CHM13 + IL
dt1 <- fread('/path/to/data/sv_calls/smoove/0.2.6/chm13/HG001/20221015/HG001-smoove.vcf.gz',skip="#CHROM")
dt1 <- extract_vcf(dt1)
dt1[,LABEL:="HG001 Illumina CHM13"]

# HG001 + HG38 + IL
dt2 <- fread('/path/to/data/sv_calls/smoove/0.2.6/hg38/HG001/20221015/HG001-smoove.vcf.gz',skip="#CHROM")
dt2 <- extract_vcf(dt2)
dt2[,LABEL:="HG001 Illumina HG38"]

# HG002 + CHM13 + IL
dt3 <- fread('/path/to/data/sv_calls/smoove/0.2.6/chm13/HG002/20221015/HG002-smoove.vcf.gz',skip="#CHROM")
dt3 <- extract_vcf(dt3)
dt3[,LABEL:="HG002 Illumina CHM13"]

# HG002 + HG38 + IL
dt4 <- fread('/path/to/data/sv_calls/smoove/0.2.6/hg38/HG002/20221015/HG002-smoove.vcf.gz',skip="#CHROM")
dt4 <- extract_vcf(dt4)
dt4[,LABEL:="HG002 Illumina HG38"]

# HG001 + CHM13 + ONT
dt5 <- fread('/path/to/data/sv_calls/cutesv/v1.0.13/<sample_id>/20220429/PLPN239131_pass.vcf',skip="#CHROM")
dt5 <- extract_vcf(dt5)
dt5[,LABEL:="HG001 ONT CHM13"]

# HG001 + HG38 + ONT
dt6 <- fread('/path/to/data/sv_calls/cutesv/v1.0.13/hg38/<sample_id>/20221013/PLPN239131_pass.vcf',skip="#CHROM")
dt6 <- extract_vcf(dt6)
dt6[,LABEL:="HG001 ONT HG38"]

# HG002 + CHM13 + ONT
dt7 <- fread('/path/to/data/sv_calls/cutesv/v1.0.13/<sample_id>/20220428/<sample_id>_pass.vcf',skip="#CHROM")
dt7 <- extract_vcf(dt7)
dt7[,LABEL:="HG002 ONT CHM13"]

# HG002 + HG38 + ONT
dt8 <- fread('/path/to/data/sv_calls/cutesv/v1.0.13/hg38/<sample_id>/20221013/<sample_id>_pass.vcf',skip="#CHROM")
dt8 <- extract_vcf(dt8)
dt8[,LABEL:="HG002 ONT HG38"]

all <- rbind(dt1,dt2,dt3,dt4,dt5,dt6,dt7,dt8)

sv_count <- all[SVLEN>50 | SVTYPE=="BND" | SVTYPE=="TRA",.N,by=list(LABEL)]

sv_count[,c("SAMPLE","PLATFORM","GENOME"):=tstrsplit(LABEL," ",fixed=TRUE)]

sv_count[,FACET:=paste(SAMPLE,PLATFORM," ")]

pdf('/path/to/Figures/main_figures_v2/figure1/SV_count_by_type.pdf',height=8,width=8)

ggplot(sv_count,aes(FACET,N/1000,fill=GENOME))+
    geom_bar(stat="identity",position=position_dodge())+
    scale_y_continuous(breaks=c(0,10,20))+
    scale_fill_manual(values=colors_genome)+
    theme_bw(base_size=20)+
    theme(aspect.ratio=1)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.text.x=element_text(angle=45,hjust=1))+
    theme(legend.position="bottom")+
    labs(x=NULL,y="Counts (x1,000)",fill="Reference")

dev.off()
       
# INS and DEL imbalance

ins_del_imbalance <- dcast(sv_count[SVTYPE %in% c("INS","DEL") & PLATFORM=="ONT",.(LABEL,SVTYPE,N)],LABEL~SVTYPE,value.var="N")
       
# Statistics
       
female <- 3054832041
male <- 3117292070
       
alignment <- fread('/path/to/Figures/main_figures_v2/figure1/alignment_stats.tab',sep="\t",header=TRUE)
alignment <- merge(alignment,cohort[,.(GUAID,Sex,RunID,Community)])
coverage <- alignment[,.(RunID,Community,Sex,Coverage=ifelse(Sex=="XX",Mapped_Bases/female,Mapped_Bases/male))]
coverage_quantile <- quantile(coverage[,Coverage])
min_coverage <- min(coverage[,Coverage])
max_coverage <- max(coverage[,Coverage])
       
n50 <- fread('/path/to/Figures/main_figures_v2/figure1/coverage_stats.tab',sep="\t",header=TRUE)
n50 <- merge(n50,cohort[,.(GUAID,Sex,RunID,Community)])
n50 <- n50[,.(RunID,Community,Sex,N50)]
n50_quantile <- quantile(n50[,N50])
min_n50 <- min(n50[,N50])
max_n50 <- max(n50[,N50])
       
    
