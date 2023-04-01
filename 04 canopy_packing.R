setwd("C:/Users/querc/Dropbox/PostdocProjects/SelfPruning")

library(mosaic)

########################################
## data on crown shape and size

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

######################################
## calculate canopy packing

crown_area<-function(CD,CR,beta){(CR*CD)/(beta+1)}
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

## for the purpose of checking,
## note that Bj = 0.5 is parabolic (V = 1/2*pi*r^2*h)
## and Bj = 1 is conical (V = 1/3*pi*r^2*h)

## long, narrow crown
tree1<-list(CRmax=50,
            CB=500,
            CD=200,
            Bj=0.5)

## shallow, wide crown
tree2<-list(CRmax=100,
            CB=600,
            CD=50,
            Bj=0.5)

tree3<-list(CRmax=30,
            CB=400,
            CD=200,
            Bj=0.3)

tree4<-list(CRmax=30,
            CB=500,
            CD=100,
            Bj=0.5)

## to see the effects
x1<-with(tree3,CB:(CB+CD))
y1<-with(tree3,CRmax*((CD+CB-x1)/CD)^Bj)
plot(y1~x1,ylim=c(0,tree3$CRmax))

x2<-with(tree4,CB:(CB+CD))
y2<-with(tree4,CRmax*((CD+CB-x2)/CD)^Bj)
points(y2~x2,col="red")

calculate_CCI<-function(two_trees){
  area1<-crown_area(CD=two_trees[[1]]$CD,CR=two_trees[[1]]$CRmax,beta=two_trees[[1]]$Bj)
  area2<-crown_area(CD=two_trees[[2]]$CD,CR=two_trees[[2]]$CRmax,beta=two_trees[[2]]$Bj)
  overlap_2D<-calculate_2D_overlap(two_trees)
  
  vol1<-crown_vol(CD=two_trees[[1]]$CD,CR=two_trees[[1]]$CRmax,beta=two_trees[[1]]$Bj)
  vol2<-crown_vol(CD=two_trees[[2]]$CD,CR=two_trees[[2]]$CRmax,beta=two_trees[[2]]$Bj)
  overlap_3D<-calculate_3D_overlap(two_trees)
  
  perc_overlap_2D<-overlap_2D/(area1+area2)
  CCI_2D<-1-2*perc_overlap_2D
  
  perc_overlap_3D<-overlap_3D/(vol1+vol2)
  CCI_3D<-1-2*perc_overlap_3D
  CCI_list<-list(area1=area1,
                 area2=area2,
                 overlap_2D=overlap_2D,
                 CCI_2D=CCI_2D,
                 vol1=vol1,
                 vol2=vol2,
                 overlap_3D=overlap_3D,
                 CCI_3D=CCI_3D)
  return(CCI_list)
}

self_pruning_alive<-self_pruning[-which(is.na(self_pruning$CrownDepth)),]

tree_samp_1<-self_pruning_alive[sample(1:nrow(self_pruning_alive),1),]
tree_samp_2<-self_pruning_alive[sample(1:nrow(self_pruning_alive),1),]

tree_list<-list(list(CRmax=tree_samp_1$CR_average,
                     CD=tree_samp_1$CrownDepth,
                     CB=tree_samp_1$HeightBase,
                     Bj=tree_samp_1$Bj),
                list(CRmax=tree_samp_2$CR_average,
                     CD=tree_samp_2$CrownDepth,
                     CB=tree_samp_2$HeightBase,
                     Bj=tree_samp_2$Bj))

calculate_CCI(tree_list)
