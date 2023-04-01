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

calculate_overlap<-function(two_trees){
  
  treeA<-two_trees[[1]]
  treeB<-two_trees[[2]]
  
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
  if(min_height < max_base){
    return(0)
  }
  
  ## otherwise we have to solve for whether (and at what height)
  ## their crown traces intersect
  trace_int<-unlist(findZeros(CRmax_A*((CD_A+CB_A-h)/CD_A)^Bj_A-CRmax_B*((CD_B+CB_B-h)/CD_B)^Bj_B~h,
                              xlim=c(max_base,min_height)))
  
  ## if they don't intersect...
  if(length(trace_int)==0){
    ## integrate the tree with the shorter height from
    ## its top to the higher of the two crown bases
    
    shorter_tree<-two_trees[[which.min(c(H_A,H_B))]]
    CRmax_s<-shorter_tree$CRmax
    CB_s<-shorter_tree$CB
    CD_s<-shorter_tree$CD
    Bj_s<-shorter_tree$Bj
    
    constant_frac<- (-pi*CRmax_s^2*CD_s)/(2*Bj_s+1)
    integral<-((CB_s+CD_s-min_height)/CD_s)^(2*Bj_s+1)-((CB_s+CD_s-max_base)/CD_s)^(2*Bj_s+1)
    
    full_integral<-constant_frac*integral
    return(full_integral)
  }
  
  if(length(trace_int)==1){
    
    ## if the one intersection point is at the top
    ## due to equal heights...
    
    if(isTRUE(all.equal(min_height,trace_int))){
      
      ## the narrower tree is the one with the smaller CRmax;
      ## if the CRmax values are also equal, pick the one with
      ## the larger Bj
      
      if(isTRUE(all.equal(CRmax_A,CRmax_B))){
        narrower_tree<-two_trees[[which.max(c(Bj_A,Bj_B))]]
      } else {
        narrower_tree<-two_trees[[which.min(c(CRmax_A,CRmax_B))]]
      }
      
      ## integrate the narrower tree from the
      ## higher of the two crown bases to the top
      
      CRmax_n<-narrower_tree$CRmax
      CB_n<-narrower_tree$CB
      CD_n<-narrower_tree$CD
      Bj_n<-narrower_tree$Bj
      
      constant_frac<- (-pi*CRmax_n^2*CD_n)/(2*Bj_n+1)
      integral<-((CB_n+CD_n-min_height)/CD_n)^(2*Bj_n+1)-((CB_n+CD_n-max_base)/CD_n)^(2*Bj_n+1)
      
      full_integral<-constant_frac*integral
      return(full_integral)
    }
    
    ## if the one intersection point is not at the top...
    
    else{
    
      ## integrate the shorter tree from min_height to the
      ## intersection point
      
      shorter_tree<-two_trees[[which.min(c(H_A,H_B))]]
      taller_tree<-two_trees[[which.max(c(H_A,H_B))]]
      
      CRmax_s<-shorter_tree$CRmax
      CB_s<-shorter_tree$CB
      CD_s<-shorter_tree$CD
      Bj_s<-shorter_tree$Bj
      
      constant_frac_s<- (-pi*CRmax_s^2*CD_s)/(2*Bj_s+1)
      integral_s<-((CB_s+CD_s-min_height)/CD_s)^(2*Bj_s+1)-((CB_s+CD_s-trace_int)/CD_s)^(2*Bj_s+1)
      
      full_integral_s<-constant_frac_s*integral_s
      
      ## then integrate the taller tree from the
      ## intersection point to max_base
      
      CRmax_t<-taller_tree$CRmax
      CB_t<-taller_tree$CB
      CD_t<-taller_tree$CD
      Bj_t<-taller_tree$Bj
      
      constant_frac_t<- (-pi*CRmax_t^2*CD_t)/(2*Bj_t+1)
      integral_t<-((CB_t+CD_t-trace_int)/CD_t)^(2*Bj_t+1)-((CB_t+CD_t-max_base)/CD_t)^(2*Bj_t+1)
      
      full_integral_t<-constant_frac_t*integral_t
      
      full_integral<-full_integral_t+full_integral_s
      
      return(full_integral)
    }
  }
  
  ## I think it should be impossible for the two traces
  ## to have more than two points of intersection
  ## so we just deal with the case where there are exactly two
  if(length(trace_int)==2){
    
    ## if the trees have the same height
    ## so one of the points of intersection is min_height...
    if(any(sapply(trace_int,function(x) isTRUE(all.equal(x,min_height))))){

      ## figure out which crown is "narrower" in each of
      ## the upper and lower portions
      
      hu<-runif(1,min = min(trace_int),max=min_height)
      r1A<-CRmax_A*((CD_A+CB_A-hu)/CD_A)^Bj_A
      r1B<-CRmax_B*((CD_B+CB_B-hu)/CD_B)^Bj_B
      
      narrow_upper_tree<-two_trees[[which.min(c(r1A,r1B))]]
      narrow_lower_tree<-two_trees[[which.max(c(r1A,r1B))]]
      
      CRmax_nu<-narrow_upper_tree$CRmax
      CB_nu<-narrow_upper_tree$CB
      CD_nu<-narrow_upper_tree$CD
      Bj_nu<-narrow_upper_tree$Bj
      
      constant_frac_nu<- (-pi*CRmax_nu^2*CD_nu)/(2*Bj_nu+1)
      integral_nu<-((CB_nu+CD_nu-min_height)/CD_nu)^(2*Bj_nu+1)-((CB_nu+CD_nu-min(trace_int))/CD_nu)^(2*Bj_nu+1)
      
      full_integral_nu<-constant_frac_nu*integral_nu
      
      ## then integrate the other tree from the
      ## intersection point to max_base
      
      CRmax_nl<-narrow_lower_tree$CRmax
      CB_nl<-narrow_lower_tree$CB
      CD_nl<-narrow_lower_tree$CD
      Bj_nl<-narrow_lower_tree$Bj
      
      constant_frac_nl<- (-pi*CRmax_nl^2*CD_nl)/(2*Bj_nl+1)
      integral_nl<-((CB_nl+CD_nl-min(trace_int))/CD_nl)^(2*Bj_nl+1)-((CB_nl+CD_nl-max_base)/CD_nl)^(2*Bj_nl+1)
      
      full_integral_nl<-constant_frac_nl*integral_nl
      
      full_integral<-full_integral_nu+full_integral_nl
      return(full_integral)
      
    } else {
      
      ## if neither of the points of intersection is
      ## at min_height...
      
      ## integrate the shorter tree from
      ## first intersection to min_height
      ## then the taller between the intersections
      ## then the shorter again from the
      ## second intersection to max base
      
      shorter_tree<-two_trees[[which.min(c(H_A,H_B))]]
      taller_tree<-two_trees[[which.max(c(H_A,H_B))]]
      
      CRmax_s<-shorter_tree$CRmax
      CB_s<-shorter_tree$CB
      CD_s<-shorter_tree$CD
      Bj_s<-shorter_tree$Bj
      
      constant_frac_s<- (-pi*CRmax_s^2*CD_s)/(2*Bj_s+1)
      integral_s_upper<-((CB_s+CD_s-min_height)/CD_s)^(2*Bj_s+1)-((CB_s+CD_s-max(trace_int))/CD_s)^(2*Bj_s+1)
      integral_s_lower<-((CB_s+CD_s-min(trace_int))/CD_s)^(2*Bj_s+1)-((CB_s+CD_s-max_base)/CD_s)^(2*Bj_s+1)
      
      full_integral_s<-constant_frac_s*(integral_s_upper+integral_s_lower)
      
      ## then integrate the taller tree from the
      ## intersection point to max_base
      
      CRmax_t<-taller_tree$CRmax
      CB_t<-taller_tree$CB
      CD_t<-taller_tree$CD
      Bj_t<-taller_tree$Bj
      
      constant_frac_t<- (-pi*CRmax_t^2*CD_t)/(2*Bj_t+1)
      integral_t<-((CB_t+CD_t-max(trace_int))/CD_t)^(2*Bj_t+1)-((CB_t+CD_t-min(trace_int))/CD_t)^(2*Bj_t+1)
      
      full_integral_t<-constant_frac_t*integral_t
      
      full_integral<-full_integral_s+full_integral_t
      return(full_integral)
    }
    
  }
  
  if(length(trace_int)>2){
    stop("more than two intersections between crown traces???")
  }

  return("no output")
  
}

calculate_CCI<-function(two_trees){
  overlap<-calculate_overlap(two_trees)
  vol1<-crown_vol(CD=two_trees[[1]]$CD,CR=two_trees[[1]]$CR,beta=two_trees[[1]]$Bj)
  vol2<-crown_vol(CD=two_trees[[2]]$CD,CR=two_trees[[2]]$CR,beta=two_trees[[2]]$Bj)
  
  perc_overlap<-overlap/(vol1+vol2)
  CCI<-1-2*perc_overlap
  CCI_list<-list(vol1=vol1,
                 vol2=vol2,
                 overlap=overlap,
                 CCI=CCI)
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

calculate_overlap(tree_list)
calculate_CCI(tree_list)
