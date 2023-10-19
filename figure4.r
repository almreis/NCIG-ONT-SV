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

sv_discovery_rate <- function(discovery_rate,sample_names,sampling,label) {
  discovery_rate_sorted <- discovery_rate[order(RunID)]
  discovery_rate_sorted <- merge(discovery_rate_sorted,unique(discovery_rate_sorted[,.(RunID)])[,.(RunID,SAMPLE_ORDER=seq(1:.N))])
  discovery_rate_sorted[,ORDER_APPEARANCE:=seq(1,.N),by=ID]
  discovery_rate_summary <- rbindlist(lapply(unique(discovery_rate_sorted$SAMPLE_ORDER),function(x) discovery_rate_sorted[SAMPLE_ORDER<=x][,.N,by=ID][,.(ID,N,LABEL=ifelse(N==x,"shared","nonredundant"))][,.N,by=LABEL][,.(SAMPLE_ORDER=x,LABEL,N)]))
  discovery_rate_summary[,SAMPLING:=sampling]
  discovery_rate_summary[,VARIABLE:=label]                                          
  discovery_rate_summary <- merge(discovery_rate_summary,unique(discovery_rate_sorted[,.(RunID)])[,.(RunID,SAMPLE_ORDER=seq(1:.N))],by="SAMPLE_ORDER")
  return(discovery_rate_summary)
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
colors_sharedness <- c('singleton'='#E2CFFF','polymorphic'='#7ADEA8','major'='#E995BA','shared'='#3DC1CB')
colors_range <- c('Aus'='#9292D1','Aus-absent'='#D1F0A4','Global'='#C6E8EE')
colors_novelty <- c('Found - both'='#E0D0F5','Found - Gnomad only'='#FCE4BA','Found - DECODE only'='#E1B894','Not found'='#F3AFCC','Not lifted'='#B7CAB4')
colors_community <- c('NCIGP1'='#F07DD1','NCIGP2'='#986BE3','NCIGP3'='#8DBAF3','NCIGP4'='#6DB3AC','non-NCIG'='#FCB48D')
colors_aus <- c('Community-specific'='#F4C573','Widespread'='#D95050','Private'='#7ECBE6','Ancestral'='#82C493')                                             
colors_discovery <- c('nonredundant'='#F7ECB2','shared'='#3DC1CB')
colors_community_extended <- c('NCIGP1'='#F07DD1','NCIGP2'='#986BE3','NCIGP3'='#8DBAF3','NCIGP4'='#6DB3AC','non-NCIG'='#FCB48D','NCIG'='#FD2449')
colors_community_extended2 <- c('NCIGP1'='#F07DD1','NCIGP2'='#986BE3','NCIGP3'='#8DBAF3','NCIGP4'='#6DB3AC','non-NCIG'='#FCB48D','NCIGP1/2/3'='#7fffd4')
colors_repeats_summary <- c('Non-repetitive'='#8DD3C7','Tandem repeats'='#E88794','Mobile elements'='#AA8FC4')                                       
                                             
# Load master table

vars <- fread('/path/to/Figures/all.filtered.joint_call.v5.master_table.tsv')

# PCOA plot

### with all groups

sv_matrix <- data.table(dcast(unique(vars[LEN_LABEL==">= 50bp",.(UNIQUE_ID,RunID,LEN_LABEL,CALL=1)]),UNIQUE_ID+LEN_LABEL~RunID,value.var="CALL"))
sv_matrix[is.na(sv_matrix)] <- 0
sv_matrix[,c("UNIQUE_ID","LEN_LABEL"):=NULL]

sv.dist <- vegdist(t(as.matrix(sv_matrix)), "bray")
res <- pcoa(sv.dist)
res.vectors <- data.table(res$vectors)
res.vectors[,SAMPLE:=rownames(res$vectors)]
                                             

res.vectors <- merge(res.vectors,unique(vars[,.(SAMPLE=RunID,Community)]))

pcoa1.var <- round(res$values[,"Eigenvalues"][1]/sum(res$values[,"Eigenvalues"])*100,2)                                             
pcoa2.var <- round(res$values[,"Eigenvalues"][2]/sum(res$values[,"Eigenvalues"])*100,2)
                                             
pdf('/path/to/Figures/Review/figure4/principal_coordinate_analysis.pdf',height=8,width=8)

ggplot(res.vectors,aes(Axis.1,Axis.2,color=Community))+
    geom_point(alpha=0.6,size=4)+
    scale_color_manual(values=colors_community)+
    theme_bw(base_size=20)+
    theme(aspect.ratio=1)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x=sprintf("PCoA1 (%.2f%%)",pcoa1.var),y=sprintf("PCoA2 (%.2f%%)",pcoa2.var))

dev.off() 
                                             
# PCOA plot - with only complete mobile elements                                             
sv_matrix <- data.table(dcast(unique(vars[INTERSPERSED_COMPLETENESS=="COMPLETE",.(UNIQUE_ID,RunID,LEN_LABEL,CALL=1)]),UNIQUE_ID+LEN_LABEL~RunID,value.var="CALL"))
sv_matrix[is.na(sv_matrix)] <- 0
sv_matrix[,c("UNIQUE_ID","LEN_LABEL"):=NULL]

sv.dist <- vegdist(t(as.matrix(sv_matrix)), "bray")
res <- pcoa(sv.dist)
res.vectors <- data.table(res$vectors)
res.vectors[,SAMPLE:=rownames(res$vectors)]

res.vectors <- merge(res.vectors,unique(vars[,.(SAMPLE=RunID,Community)]))                                             
ggplot(res.vectors,aes(Axis.1,Axis.2,color=Community))+
    geom_point(alpha=0.6,size=4)+
    scale_color_manual(values=colors_community)+
    theme_bw(base_size=20)+
    theme(aspect.ratio=1)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x="PCoA1",y="PCoA2")
                                             
# Variants shared between communities
                                  
aus_specific_nr <- unique(vars[FILTER=="PASS" & GLOBAL_RANGE=="Aus" & SVTYPE %in% c("INS","DEL"),.(UNIQUE_ID,SVTYPE,LEN_LABEL,REPEAT_STATUS,GLOBAL_RANGE,NCIGP1,NCIGP2,NCIGP3,NCIGP4,COMMUNITY_RANGE,COMMUNITY_RANGE_LABEL)])
                                             
                                    
aus_specific_upset <- aus_specific_nr[LEN_LABEL==">= 50bp",.(UNIQUE_ID,LEN_LABEL,NCIGP1,NCIGP2,NCIGP3,NCIGP4,COMMUNITY_RANGE_LABEL)]                                    
set_size(8, 5)

pdf('/path/to/Figures/Review/figure4/community_sharedness_upset.SVs.pdf',height=8,width=12)                                    

                                             
ComplexUpset::upset(
    aus_specific_upset,sort_intersections_by = c("degree","cardinality"),sort_intersections="ascending",
    c("NCIGP1","NCIGP2","NCIGP3","NCIGP4"),
    base_annotations=list(
        'Intersection size'=intersection_size(
            counts=FALSE,
            mapping=aes(fill=COMMUNITY_RANGE_LABEL)
        )+scale_y_continuous(breaks=c(0,5000,10000))+scale_fill_manual(values=colors_aus)+labs(fill=NULL)+theme(legend.position="bottom")
    ),
    width_ratio=0.1,
    themes=upset_default_themes(text=element_text(size=30)),
    set_sizes=FALSE
)   
                                                           
dev.off()  
                                             
aus_specific_upset <- aus_specific_nr[LEN_LABEL=="< 50bp",.(UNIQUE_ID,LEN_LABEL,NCIGP1,NCIGP2,NCIGP3,NCIGP4,COMMUNITY_RANGE_LABEL)] 
aus_specific_upset_percent <- aus_specific_upset[,.N,by=COMMUNITY_RANGE_LABEL][,.(COMMUNITY_RANGE_LABEL,N,PERCENT=N/sum(aus_specific_upset[,.N,by=COMMUNITY_RANGE_LABEL][,N])*100)]
set_size(8, 5)      
                                             
pdf('/path/to/Figures/Review/figure4/community_sharedness_upset.small_variants.pdf',height=8,width=12)                                    

                                             
ComplexUpset::upset(
    aus_specific_upset,sort_intersections_by = c("degree","cardinality"),sort_intersections="ascending",
    c("NCIGP1","NCIGP2","NCIGP3","NCIGP4"),
    base_annotations=list(
        'Intersection size'=intersection_size(
            counts=FALSE,
            mapping=aes(fill=COMMUNITY_RANGE_LABEL)
        )+scale_y_continuous(breaks=c(0,5000,10000))+scale_fill_manual(values=colors_aus)+labs(fill=NULL)+theme(legend.position="bottom")
    ),
    width_ratio=0.1,
    themes=upset_default_themes(text=element_text(size=30)),
    set_sizes=FALSE
)   
                                                           
dev.off()                                               
                                             
                                             
# Statistics

# Proportion of community-specific variants                                             
aus_specific_upset_stats <- aus_specific_upset[LEN_LABEL==">= 50bp",.N,by=COMMUNITY_RANGE_LABEL]                                                                                        
sum(aus_specific_upset_stats[COMMUNITY_RANGE_LABEL %in% c("Private","Community-specific"),N])/sum(aus_specific_upset_stats$N)*100

# Proportion of ancestral variants                                             
sum(aus_specific_upset_stats[COMMUNITY_RANGE_LABEL %in% c("Ancestral"),N])/sum(aus_specific_upset_stats$N)*100                                             
                                             
# Aus-specific variants as boxplot                                    

vars_reach <- vars[FILTER=="PASS",.(UNIQUE_ID,RunID,Index,SVTYPE,SVLEN,LEN_LABEL,REPEAT_STATUS,REPEAT_LABEL,INTERSPERSED_COMPLETENESS,Community,GLOBAL_RANGE,NCIGP1,NCIGP2,NCIGP3,NCIGP4,COMMUNITY_RANGE,COMMUNITY_RANGE_LABEL)]
                                             
aus_specific <- vars_reach[GLOBAL_RANGE=="Aus"]
                                    
aus_specific_indivuals <- aus_specific[,.N,by=list(LEN_LABEL,Index,RunID,Community,COMMUNITY_RANGE_LABEL)]

aus_specific_indivuals <- merge(aus_specific_indivuals,aus_specific_indivuals[,.(TOTAL=sum(N)),by=list(LEN_LABEL,Index)],by=c("LEN_LABEL","Index"))    
                                             
pdf('/path/to/Figures/Review/figure4/aus_specific_individuals.boxplot.pdf',height=8,width=8)                                                                            
                                    
ggplot(aus_specific_indivuals[LEN_LABEL==">= 50bp"],aes(COMMUNITY_RANGE_LABEL,N/TOTAL*100,fill=Community))+
  geom_boxplot()+   
  scale_fill_manual(values=colors_community[names(colors_community)!="non-NCIG"])+
  scale_y_continuous(breaks=c(10,30,50))+
  theme_bw(base_size=20)+
  theme(aspect.ratio=1)+                               
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  theme(axis.text.x=element_text(angle=45,hjust=1))+                                     
  labs(x=NULL,y="Percentage (%)",fill="Community")
                                 
                                    
dev.off()                                               

# Statistics

# Number of community-specific variants per individual                                             
aus_specific_stat <- aus_specific[GLOBAL_RANGE=="Aus" & LEN_LABEL==">= 50bp" & SVTYPE %in% c("INS","DEL")]
aus_specific_stat[COMMUNITY_RANGE_LABEL %in% c("Private","Community-specific"),.N,by=Index][,.(Mean=mean(N),SD=sd(N))]
                                             
# Cumulative size of community-specific variants per individual                                             
aus_specific_stat[COMMUNITY_RANGE_LABEL %in% c("Private","Community-specific"),.(SIZE=sum(SVLEN)),by=Index][,.(Mean=mean(SIZE),SD=sd(SIZE))]  

# Proportion of community-specific variation per individual 
aus_specific_indivuals[LEN_LABEL==">= 50bp" & COMMUNITY_RANGE_LABEL %in% c("Private","Community-specific"),.(N=sum(N),TOTAL=TOTAL[1]),by=Index][,.(Mean=mean(N/TOTAL*100),SD=sd(N/TOTAL*100))]                                             
                                             
# Repeat status for Aus-specific variants
 
aus_specific <- vars_reach[GLOBAL_RANGE=="Aus"]                                    
aus_specific_nr <- unique(aus_specific[SVTYPE %in% c("INS","DEL"),.(UNIQUE_ID,SVTYPE,LEN_LABEL,REPEAT_STATUS,REPEAT_LABEL,INTERSPERSED_COMPLETENESS,GLOBAL_RANGE,NCIGP1,NCIGP2,NCIGP3,NCIGP4,COMMUNITY_RANGE,COMMUNITY_RANGE_LABEL)])   

aus_specific_nr[REPEAT_LABEL=="Non-repetitive",REPEAT_LABEL:="0"]
aus_specific_nr[REPEAT_LABEL=="Tandem repeats" & REPEAT_STATUS=="STR",REPEAT_LABEL:="STR"]
aus_specific_nr[REPEAT_LABEL=="Tandem repeats" & REPEAT_STATUS=="HOMO",REPEAT_LABEL:="STR"]
aus_specific_nr[REPEAT_LABEL=="Tandem repeats" & REPEAT_STATUS=="TR",REPEAT_LABEL:="TR"]
aus_specific_nr[REPEAT_LABEL=="Mobile elements" & INTERSPERSED_COMPLETENESS=="INCOMPLETE",REPEAT_LABEL:="ME - fragment"]
aus_specific_nr[REPEAT_LABEL=="Mobile elements" & INTERSPERSED_COMPLETENESS=="COMPLETE",REPEAT_LABEL:="ME - complete"]     
                                             
aus_specific_nr_repeats <- aus_specific_nr[,.N,by=list(COMMUNITY_RANGE_LABEL,LEN_LABEL,REPEAT_LABEL)]
aus_specific_nr_repeats <- merge(aus_specific_nr_repeats,aus_specific_nr_repeats[,.(TOTAL=sum(N)),by=list(LEN_LABEL,COMMUNITY_RANGE_LABEL)],by=c("LEN_LABEL","COMMUNITY_RANGE_LABEL"))
        
total_label <- unique(aus_specific_nr_repeats[,.(COMMUNITY_RANGE_LABEL,LEN_LABEL,TOTAL)])[,.(LEN_LABEL,COMMUNITY_RANGE_LABEL,REPEAT_STATUS="LINE",y=102,TOTAL)]                                     
       
aus_specific_nr_repeats[,PERCENT:=round(N/TOTAL*100)]   
#percent_label <- aus_specific_nr_repeats[LEN_LABEL==">= 50bp" & REPEAT_STATUS %in% c("SINE","TR","STR")]               
#percent_label <- percent_label[order(COMMUNITY_RANGE_LABEL,REPEAT_STATUS,PERCENT,decreasing=TRUE)]                     
#percent_label[,CUMSUM:=cumsum(PERCENT),by=COMMUNITY_RANGE_LABEL]                                    
#percent_label[,Y:=CUMSUM-((CUMSUM-(CUMSUM-PERCENT))/2)]

                                      
                                             
pdf('/path/to/Figures/Review/figure4/aus_specific_repeat_type.pdf',height=8,width=8)                                                                                                                
                                    
ggplot(aus_specific_nr_repeats[LEN_LABEL==">= 50bp"],aes(COMMUNITY_RANGE_LABEL,N/TOTAL*100,fill=REPEAT_LABEL))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=colors_mobile)+   
  scale_y_continuous(breaks=c(0,50,100))+                                           
  #geom_text(data=total_label[LEN_LABEL==">= 50bp"],aes(x=COMMUNITY_RANGE_LABEL,y=y,label=TOTAL),vjust = 0,size=5)+                           
  #geom_text(data=percent_label[REPEAT_STATUS %in% c("SINE","TR")],aes(x=COMMUNITY_RANGE_LABEL,y=as.numeric(Y),label = paste0(PERCENT,"%")),size = 5,color="white")+                                                                  
  theme_bw(base_size=20)+
  theme(aspect.ratio=1)+                            
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x=element_text(angle=45,hjust=1))+                                                                     
  labs(x=NULL,y="Percentage (%)",fill=NULL)   
                                    
dev.off()                                             
                                                                                   
                                             
# SV discovery rate
                                             
iterations <- 2                                                    

discovery_rate <- vars[FILTER=="PASS" & SVTYPE %in% c("INS","DEL") & LEN_LABEL==">= 50bp" & Community_Country_Level=="NCIG",.(RunID,ID,SVTYPE,SVLEN,CALL,REPEAT_STATUS,SHAREDNESS_FULL_LABEL,LEN_LABEL)]
sample_names <- unique(discovery_rate$RunID)
                                             
## get curves SVs                                             

VAR <- "LEN_LABEL" 
                                                                                        
discovery_curve <- rbindlist(lapply(unique(discovery_rate[,get(VAR)]), function(y)                                            
                        rbindlist(lapply(1:iterations, function(x) sv_discovery_rate(discovery_rate[get(VAR)==y],sample(sample_names),x,y)))))                                                                                       
discovery_curve_mean <- discovery_curve[,.(N=mean(N)),by=list(SAMPLE_ORDER,LABEL,VARIABLE)]
discovery_curve_mean[,TOTAL:=max(discovery_curve_mean[,.(N=sum(as.numeric(N))),by=list(VARIABLE,SAMPLE_ORDER)][,N])]                                         
                                         
pdf('/path/to/Figures/Review/figure4/variant_discovery_saturation_curves.pdf',height=8,width=8)                                            
                                            
ggplot(discovery_curve_mean,aes(SAMPLE_ORDER,N/TOTAL,fill=LABEL))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=colors_discovery)+  
  scale_y_continuous(breaks=c(0,0.5,1))+                                       
  theme_bw(base_size=20)+
  theme(aspect.ratio=1)+
  theme(legend.position="bottom")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x="Samples",y="Discovery set proportion",fill=NULL) 
                                            
dev.off()   

# SV discovery rate curve per community
                                         
discovery_rate <- vars[FILTER=="PASS" & SVTYPE %in% c("INS","DEL") & LEN_LABEL==">= 50bp",.(RunID,ID,SVTYPE,SVLEN,CALL,REPEAT_STATUS,LEN_LABEL,Community_Country_Level,Community)]
sample_names <- unique(discovery_rate$RunID)

VAR <- "Community" 
                                                    
discovery_curve_community <- rbindlist(lapply(unique(discovery_rate[,get(VAR)]), function(y)                                            
                        rbindlist(lapply(1:iterations, function(x) sv_discovery_rate(discovery_rate[get(VAR)==y],sample(sample_names),x,y)))))
                                         
discovery_curve_community_mean <- discovery_curve_community[,.(N=mean(N)),by=list(SAMPLE_ORDER,LABEL,VARIABLE)]
discovery_curve_community_mean[,TOTAL:=max(discovery_curve_community_mean[,.(N=sum(as.numeric(N))),by=list(VARIABLE,SAMPLE_ORDER)][,N])]
                                         
discovery_curve_community_sample <-  discovery_curve_community_mean[,.(TOTAL=sum(as.numeric(N))),by=list(SAMPLE_ORDER,VARIABLE)]
   
labels_community <- unique(discovery_curve_community_sample[,VARIABLE])                                      
                                         
x1 <- seq(1,70,length.out = 1000)  
                                         
regression_models_community <- lapply(labels_community,function(x) lm(TOTAL ~ log(SAMPLE_ORDER),data=discovery_curve_community_sample[VARIABLE==x]))                                                
prediction_data_community <- rbindlist(lapply(1:length(labels_community),function(y) data.table(x=x1)[,.(x,y=predict(regression_models_community[[y]],newdata = data.table(SAMPLE_ORDER=x1)),community=labels_community[y])]))                      

# Regression for NCIG samples combined
                                              
discovery_rate <- vars[FILTER=="PASS" & SVTYPE %in% c("INS","DEL") & LEN_LABEL==">= 50bp" & Community_Country_Level=="NCIG",.(RunID,ID,SVTYPE,SVLEN,CALL,REPEAT_STATUS,LEN_LABEL,Community_Country_Level,Community)]
sample_names <- unique(discovery_rate$RunID)

VAR <- "Community_Country_Level" 
                                                    
discovery_curve_ncig <- rbindlist(lapply(unique(discovery_rate[,get(VAR)]), function(y)                                            
                        rbindlist(lapply(1:iterations, function(x) sv_discovery_rate(discovery_rate[get(VAR)==y],sample(sample_names),x,y)))))
                                         
discovery_curve_ncig_mean <- discovery_curve_ncig[,.(N=mean(N)),by=list(SAMPLE_ORDER,LABEL,VARIABLE)]
discovery_curve_ncig_mean[,TOTAL:=max(discovery_curve_ncig_mean[,.(N=sum(as.numeric(N))),by=list(VARIABLE,SAMPLE_ORDER)][,N])]
                                         
discovery_curve_ncig_sample <-  discovery_curve_ncig_mean[,.(TOTAL=sum(as.numeric(N))),by=list(SAMPLE_ORDER,VARIABLE)]
   
labels_ncig <- unique(discovery_curve_ncig_sample[,VARIABLE])                                      
                                         
x1 <- seq(1,70,length.out = 1000)  
                                         
regression_models_ncig <- lapply(labels_ncig,function(x) lm(TOTAL ~ log(SAMPLE_ORDER),data=discovery_curve_ncig_sample[VARIABLE==x]))                                                
prediction_data_ncig <- rbindlist(lapply(1:length(labels_ncig),function(y) data.table(x=x1)[,.(x,y=predict(regression_models_ncig[[y]],newdata = data.table(SAMPLE_ORDER=x1)),community=labels_ncig[y])]))                     
# Regression for NCIGP1,2,3 samples combined                                         
discovery_rate <- vars[FILTER=="PASS" & SVTYPE %in% c("INS","DEL") & LEN_LABEL==">= 50bp" & Community %in% c("NCIGP1","NCIGP2","NCIGP3"),.(RunID,ID,SVTYPE,SVLEN,CALL,REPEAT_STATUS,LEN_LABEL,Community_Country_Level,Community)]
discovery_rate[,Community_Country_Level:="NCIGP1/2/3"]                                         
sample_names <- unique(discovery_rate$RunID)

VAR <- "Community_Country_Level" 
                                                    
discovery_curve_ncig123 <- rbindlist(lapply(unique(discovery_rate[,get(VAR)]), function(y)                                            
                        rbindlist(lapply(1:iterations, function(x) sv_discovery_rate(discovery_rate[get(VAR)==y],sample(sample_names),x,y)))))
                                         
discovery_curve_ncig123_mean <- discovery_curve_ncig123[,.(N=mean(N)),by=list(SAMPLE_ORDER,LABEL,VARIABLE)]
discovery_curve_ncig123_mean[,TOTAL:=max(discovery_curve_ncig123_mean[,.(N=sum(as.numeric(N))),by=list(VARIABLE,SAMPLE_ORDER)][,N])]
                                         
discovery_curve_ncig123_sample <-  discovery_curve_ncig123_mean[,.(TOTAL=sum(as.numeric(N))),by=list(SAMPLE_ORDER,VARIABLE)]
   
labels_ncig123 <- unique(discovery_curve_ncig123_sample[,VARIABLE])                                      
                                         
x1 <- seq(1,70,length.out = 1000)  
                                         
regression_models_ncig123 <- lapply(labels_ncig123,function(x) lm(TOTAL ~ log(SAMPLE_ORDER),data=discovery_curve_ncig123_sample[VARIABLE==x]))                                                
prediction_data_ncig123 <- rbindlist(lapply(1:length(labels_ncig),function(y) data.table(x=x1)[,.(x,y=predict(regression_models_ncig123[[y]],newdata = data.table(SAMPLE_ORDER=x1)),community=labels_ncig123[y])]))    
                                         
                                         
                                         
          
prediction_data <- rbind(prediction_data_community,prediction_data_ncig,prediction_data_ncig123)

pdf('/path/to/Figures/Review/figure4/variant_discovery_saturation_curves.per_community.pdf',height=8,width=8)                                                                                                                          
ggplot(prediction_data[community!="NCIG"],aes(x,y,color=community))+
  geom_line(size=2,alpha=0.6)+
  scale_color_manual(values=colors_community_extended2)+ 
  scale_y_continuous(labels = function(l) { trans = l / 1000},breaks=c(0,50000,100000))+  
  theme_bw(base_size=20)+
  theme(aspect.ratio=1)+                                        
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x="Samples",y="Non-redundant variants (x1,000)",color=NULL)                                  
                                         
dev.off()
                                         
# Regression for Aus-specific x Aus-absent variants      
                                         
discovery_rate <- vars[FILTER=="PASS" & SVTYPE %in% c("INS","DEL") & LEN_LABEL==">= 50bp",.(RunID,ID,SVTYPE,SVLEN,CALL,REPEAT_STATUS,LEN_LABEL,Community_Country_Level,GLOBAL_RANGE)]
sample_names <- unique(discovery_rate[,.(RunID,Community_Country_Level)])
                                         
equal_sampling <- function(dt) {
    non_ncig <- dt[Community_Country_Level=="non-NCIG",RunID]
    ncig <- sample(dt[Community_Country_Level=="NCIG",RunID],length(non_ncig))
    return(c(non_ncig,ncig))
}

VAR <- "GLOBAL_RANGE" 
                                                    
discovery_curve <- rbindlist(lapply(unique(discovery_rate[,get(VAR)]), function(y)                                            
                        rbindlist(lapply(1:iterations, function(x) sv_discovery_rate(discovery_rate[get(VAR)==y],equal_sampling(sample_names),x,y)))))
discovery_curve[,LEN_LABEL:=">= 50bp"]
                                         
discovery_curve_mean <- discovery_curve[,.(N=mean(N)),by=list(SAMPLE_ORDER,LABEL,VARIABLE,LEN_LABEL)]
discovery_curve_mean <- merge(discovery_curve_mean,discovery_curve_mean[,.(TOTAL=max(N)),by=LEN_LABEL],by=c("LEN_LABEL"))  

# Log models                                        
                                         
discovery_curve_sample <-  discovery_curve_mean[,.(TOTAL=sum(as.numeric(N))),by=list(LEN_LABEL,SAMPLE_ORDER,VARIABLE)]
   
labels <- unique(discovery_curve_sample[,VARIABLE])
                                         
x1 <- seq(1,70,length.out = 1000)    
                                         
regression_models <- lapply(labels,function(x) lm(TOTAL ~ log(SAMPLE_ORDER),data=discovery_curve_sample[LEN_LABEL==">= 50bp" & VARIABLE==x]))                                                    
prediction_data<- rbindlist(lapply(1:length(labels),function(y) data.table(x=x1)[,.(x,y=predict(regression_models[[y]],newdata = data.table(SAMPLE_ORDER=x1)),community=labels[y])]))   
prediction_data[,LEN_LABEL:=">= 50bp"]     

pdf('/path/to/Figures/Review/figure4/variant_discovery_saturation_curves.aus_v_aus_absent_regression.pdf',height=8,width=8)  
                                   
ggplot(prediction_data,aes(x,y,color=community))+
  geom_line(size=2)+                              
  scale_y_continuous(limits=c(0,70000),labels = function(l) { trans = l / 1000},breaks=c(0,25000,50000))+                                 
  scale_color_manual(values=colors_range)+
  theme_bw(base_size=20)+
  theme(aspect.ratio=1)+                                        
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x="Samples",y="Non-redundant variants (x1,000)",color=NULL)  
                                   
dev.off()           
                                   
##################################################                                         
   
# SV discovery curves parsed by (i) non-repetitive (ii) interspersed repeat (iii) tandem repeat                               
                                     
discovery_rate <- vars[FILTER=="PASS" & SVTYPE %in% c("INS","DEL") & LEN_LABEL==">= 50bp",.(RunID,ID,SVTYPE,SVLEN,CALL,REPEAT_STATUS,LEN_LABEL,REPEAT_LABEL)]
sample_names <- unique(discovery_rate$RunID)
                                     
VAR <- "REPEAT_LABEL" 
                                                    
discovery_curve <- rbindlist(lapply(unique(discovery_rate[,get(VAR)]), function(y)                                            
                        rbindlist(lapply(1:iterations, function(x) sv_discovery_rate(discovery_rate[get(VAR)==y],sample(sample_names),x,y)))))
                            
discovery_curve_mean <- discovery_curve[,.(N=mean(N)),by=list(SAMPLE_ORDER,LABEL,VARIABLE)]
discovery_curve_mean[,TOTAL:=max(discovery_curve_mean[,.(N=sum(as.numeric(N))),by=list(VARIABLE,SAMPLE_ORDER)][,N])]
                                         
discovery_curve_mean$VARIABLE <- factor(discovery_curve_mean$VARIABLE,levels=c("Non-repetitive","Mobile elements","Tandem repeats"))                                                                                                              
# Log models                                        
                                         
discovery_curve_sample <-  discovery_curve_mean[,.(TOTAL=sum(as.numeric(N))),by=list(SAMPLE_ORDER,VARIABLE)]
   
labels <- unique(discovery_curve_sample[,VARIABLE])
                                         
x1 <- seq(1,70,length.out = 1000)
                                         
# SVs
                                     
regression_models <- lapply(labels,function(x) lm(TOTAL ~ log(SAMPLE_ORDER),data=discovery_curve_sample[VARIABLE==x]))                        
                            
prediction_data <- rbindlist(lapply(1:length(labels),function(y) data.table(x=x1)[,.(x,y=predict(regression_models[[y]],newdata = data.table(SAMPLE_ORDER=x1)),community=labels[y])]))                       
                                    
pdf('/path/to/Figures/Review/figure4/variant_discovery_saturation_curves.repeat_labels_regression.pdf',height=8,width=12)                                                                        
ggplot(prediction_data,aes(x,y,color=community))+
  geom_line(size=2)+
  scale_color_manual(values = colors_repeats_summary)+
  scale_y_continuous(labels = function(l) { trans = l / 1000},breaks=c(0,40000,80000))+                                 
  theme_bw(base_size=20)+
  theme(aspect.ratio=1)+                                        
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x="Samples",y="Non-redundant variants (x1,000)",color=NULL)                                        
                                    
dev.off()                                
