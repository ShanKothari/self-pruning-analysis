setwd("C:/Users/querc/Dropbox/PostdocProjects/SelfPruning")

library(mosaic)

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

## long, narrow crown
tree1<-list(CRmax=50,
            CB=500,
            CD=200,
            Bj=0.3)

## shallow, wide crown
tree2<-list(CRmax=100,
            CB=600,
            CD=50,
            Bj=0.3)

tree3<-list(CRmax=50,
            CB=400,
            CD=100,
            Bj=0.4)

tree4<-list(CRmax=75,
            CB=600,
            CD=100,
            Bj=0.2)

calculate_overlap<-function(treeA,treeB){
  CRmax_A<-treeA$CRmax
  CB_A<-treeA$CB
  CD_A<-treeA$CD
  Bj_A<-treeA$Bj
  H_A<-CB_A+CD_A
  
  CRmax_B<-treeB$CRmax
  CB_B<-treeB$CB
  CD_B<-treeB$CD
  Bj_B<-treeB$Bj
  H_B<-CB_B+CD_B
  
  max_base<-max(CB_A,CB_B)
  min_height<-min(H_A,H_B)
  
  ## if one tree's crown is entirely below the other
  ## return 0
  if(H_A<=CB_B | H_B<=CB_A){
    return(0)
  }
  
  ## otherwise we have to solve for whether their
  ## crown traces intersect
  trace_int<-as.numeric(findZeros(CRmax_A*((CD_A+CB_A-h)/CD_A)^Bj_A-CRmax_B*((CD_B+CB_B-h)/CD_B)^Bj_B~h,
                                 xlim=c(max_base,min_height)))
  
  if(length(trace_int)==0){
    ## integrate the tree with the shorter height from
    ## its top to the higher of the two crown bases
  }
  
  if(length(trace_int)==1){
    if(isTRUE(all.equal(H_A,H_B)) &&
       isTRUE(all.equal(H_A,trace_int))){
      ## integrate the tree with the smaller CRmax from the
      ## higher of the two crown bases to the top
    }
    else{
      ## integrate the shorter tree from min_height to the
      ## intersection point
      ## then the taller tree from the intersection point
      ## to max_base
    }
  }
  
  if(length(trace_int)>1){
    stop("multiple intersections of crown trace")
  }

}

####################################
## calculate canopy complementarity
## volume of overlap (min(radius) across height) divided by
## sum of two individual columes
