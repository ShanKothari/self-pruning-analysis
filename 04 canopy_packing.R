setwd("C:/Users/querc/Dropbox/PostdocProjects/SelfPruning")

self_pruning<-read.csv("SelfPruningData/self_pruning_processed.csv")
## get indicators of species composition
self_pruning_list<-split(self_pruning,f = self_pruning$Plot)
self_pruning_comp<-unlist(lapply(self_pruning_list,
                                 function(plot) paste(unique(plot$Species),collapse="|")))

C0B<-0.196
C1B<-0.511

species.list<-c("ABBA","ACRU","ACSA","BEAL","BEPA","LALA",
           "PIGL","PIRE","PIRU","PIST","QURU","THOC")
Tj<-c(0.278,0.536,0.560,0.592,0.435,0.308,
      0.278,0.287,0.324,0.417,0.538,0.251)

crown_shape<-data.frame(species=species.list,
                        Tj=Tj)
crown_shape$Bj<-(1-crown_shape$Tj)*C0B+crown_shape$Tj*C1B
self_pruning$Bj<-crown_shape$Bj[match(self_pruning$Species,crown_shape$species)]

crown_vol<-function(CD,CR,beta){(pi*CR^2*CD)/(2*beta+1)}

self_pruning$crown_vol<-sapply(1:nrow(self_pruning),
                               function(i) {
                                 beta_i<-self_pruning$Bj[i]
                                 crown_i<-crown_vol(CD=self_pruning$CrownDepth[i]/100,
                                                    CR=self_pruning$CR_average[i]/100,
                                                    beta=beta_i)
                                 return(crown_i)
                               })

self_pruning$crown_vol[is.na(self_pruning$crown_vol)]<-0

crown_vol_agg<-aggregate(self_pruning$crown_vol,
                         by=list(self_pruning$Block,
                                 self_pruning$Plot,
                                 self_pruning$Species,
                                 self_pruning$nbsp),
                         FUN=mean,na.rm=T)

colnames(crown_vol_agg)<-c("Block","Plot","Species","Richness","crown_vol")
crown_vol_agg$sp_comp<-self_pruning_comp[match(crown_vol_agg$Plot,names(self_pruning_comp))]

crown_vol_agg$mono.means<-apply(crown_vol_agg,1,
                                function(x) {
                                  crown_vol_agg$crown_vol[crown_vol_agg$sp_comp==x["Species"] & crown_vol_agg$Block==x["Block"]]
                                })
crown_vol_agg$OY<-crown_vol_agg$crown_vol-crown_vol_agg$mono.means

crown_OY_plot<-aggregate(crown_vol_agg$OY,
                         by=list(crown_vol_agg$Plot,
                                 crown_vol_agg$Richness),
                         FUN=sum)

## to do:
## add measures of functional diversity and
## heterogeneity in shade tolerance
## (non-abundance weighted?)
## test how much overyielding you get when holding
## self-shading behavior constant (i.e. no plasticity)

######################################
## canopy complementarity sandbox

CRmax_1<-50
CB_1<-500
CD_1<-200
Bj_1<-0.3

CRmax_2<-100
CB_2<-600
CD_2<-50
Bj_2<-0.3

h_intersect<-function(h) CRmax_1*((CD_1+CB_1-h)/CD_1)^Bj_1-CRmax_2*((CD_2+CB_2-h)/CD_2)^Bj_2
uniroot(h_intersect,lower=max(CB_1,CB_2),upper=min(CB_1+CD_1,CB_2+CD_2))

####################################
## calculate canopy complementarity
## volume of overlap (min(radius) across height) divided by
## sum of two individual columes
