setwd("C:/Users/querc/Dropbox/PostdocProjects/SelfPruning")

library(mosaic)

## to dos:
## calculate functional diversity/heterogeneity in shade tolerance
## at the plot scale for complementarity analyses

########################################
## data on crown shape and size

self_pruning<-read.csv("SelfPruningData/self_pruning_processed.csv")
## get indicators of species composition
self_pruning_list<-split(self_pruning,f = self_pruning$Plot)
self_pruning_comp<-unlist(lapply(self_pruning_list,
                                 function(plot) paste(unique(plot$Species),collapse="|")))

## read in mortality data
DB_mortality<-read.csv("IDENTMontrealData/mortality_2018.csv")
DB_mortality$plot_sp<-paste(DB_mortality$Block,
                             DB_mortality$Plot,
                             DB_mortality$CodeSp,
                             sep="_")

## parameters from Purves et al. paper
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

## formulas for area/volume of crown
crown_area<-function(CD,CR,beta){(CR*CD)/(beta+1)}
crown_vol<-function(CD,CR,beta){(pi*CR^2*CD)/(2*beta+1)}

## through rows of the self-pruning data
self_pruning$crown_vol<-sapply(1:nrow(self_pruning),
                               function(i) {
                                 crown_i<-crown_vol(CD=self_pruning$CrownDepth[i]/100,
                                                    CR=self_pruning$CR_average[i]/100,
                                                    beta=self_pruning$Bj[i])
                                 return(crown_i)
                               })


## NOTE: does this aggregation account for unequal numbers of planted individuals
## of various species within a plot?

crown_vol_agg<-aggregate(self_pruning$crown_vol,
                         by=list(self_pruning$Block,
                                 self_pruning$Plot,
                                 self_pruning$Species,
                                 self_pruning$nbsp),
                         FUN=mean,na.rm=T)

colnames(crown_vol_agg)<-c("Block","Plot","Species","Richness","crown_vol")
crown_vol_agg$sp_comp<-self_pruning_comp[match(crown_vol_agg$Plot,names(self_pruning_comp))]

## add mortality rate column and multiply volume by mortality rate
crown_vol_agg$plot_sp<-paste(crown_vol_agg$Block,
                             crown_vol_agg$Plot,
                             crown_vol_agg$Species,
                             sep="_")

crown_vol_agg$prop_alive<-DB_mortality$Alive[match(crown_vol_agg$plot_sp,
                                                   DB_mortality$plot_sp)]
crown_vol_agg$crown_vol_adj<-crown_vol_agg$crown_vol*crown_vol_agg$prop_alive

crown_vol_agg$mono.means<-apply(crown_vol_agg,1,
                                function(x) {
                                  crown_vol_agg$crown_vol_adj[crown_vol_agg$sp_comp==x["Species"] & crown_vol_agg$Block==x["Block"]]
                                })
crown_vol_agg$OY<-crown_vol_agg$crown_vol_adj-crown_vol_agg$mono.means

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

######################################
## function to calculate complementarity

calculate_CCI<-function(two_trees){
  area1<-crown_area(CD=two_trees[[1]]$CD,CR=two_trees[[1]]$CRmax,beta=two_trees[[1]]$Bj)
  area2<-crown_area(CD=two_trees[[2]]$CD,CR=two_trees[[2]]$CRmax,beta=two_trees[[2]]$Bj)
  overlap_2D<-calculate_2D_overlap(two_trees)
  
  vol1<-crown_vol(CD=two_trees[[1]]$CD,CR=two_trees[[1]]$CRmax,beta=two_trees[[1]]$Bj)
  vol2<-crown_vol(CD=two_trees[[2]]$CD,CR=two_trees[[2]]$CRmax,beta=two_trees[[2]]$Bj)
  overlap_3D<-calculate_3D_overlap(two_trees)
  
  CCI_2D<-1-2*overlap_2D/(area1+area2)
  CCI_min_2D<-1-overlap_2D/min(c(area1,area2))
  
  CCI_3D<-1-2*overlap_3D/(vol1+vol2)
  CCI_min_3D<-1-overlap_3D/min(c(vol1,vol2))
  
  ## here I make a crucial assumption:
  ## if both trees are dead, return NA for all CCI metric
  ## if one tree is dead, return 1 (perfect complementarity)
  
  ## how many trees are dead?
  dead_trees<-sum(c(isTRUE(all.equal(area1,0)),
                    isTRUE(all.equal(area2,0))))
  
  if(dead_trees==2){
    CCI_2D<-NA
    CCI_min_2D<-NA
    CCI_3D<-NA
    CCI_min_3D<-NA
  }
  
  if(dead_trees==1){
    CCI_min_2D<-1
    CCI_min_3D<-1
  }
  
  CCI_list<-list(dead_trees=dead_trees,
                 area1=area1,
                 area2=area2,
                 overlap_2D=overlap_2D,
                 CCI_2D=CCI_2D,
                 CCI_min_2D=CCI_min_2D,
                 vol1=vol1,
                 vol2=vol2,
                 overlap_3D=overlap_3D,
                 CCI_3D=CCI_3D,
                 CCI_min_3D=CCI_min_3D)
  return(unlist(CCI_list))
}

tree_samp_1<-self_pruning[sample(1:nrow(self_pruning),1),]
tree_samp_2<-self_pruning[sample(1:nrow(self_pruning),1),]

tree_list<-list(list(CRmax=tree_samp_1$CR_average,
                     CD=tree_samp_1$CrownDepth,
                     CB=tree_samp_1$HeightBase,
                     Bj=tree_samp_1$Bj),
                list(CRmax=tree_samp_2$CR_average,
                     CD=tree_samp_2$CrownDepth,
                     CB=tree_samp_2$HeightBase,
                     Bj=tree_samp_2$Bj))

calculate_CCI(tree_list)

#####################################
## applying the complementarity function
## to calculate the complementarity of each plot

self_pruning_split<-split(self_pruning,f = ~Block+Plot)

## this procedure needs to be updated because mortalities
## in the self-pruning data are not "true" mortalities
## so we need a sampling procedure to determine whether
## to select dead trees based on actual mortality rates
## in the full community survey

plot_CCI<-function(plot) {
  
  ## get every pair of trees sampled in the plot
  plot_combos<-t(combn(1:nrow(plot),2))
  
  tree_pair_CCI<-list()
  ## loop through pairs and calculate CCI for each one
  for(i in 1:nrow(plot_combos)){
    tree_a<-list(CRmax=plot$CR_average[plot_combos[i,1]],
                 CB=plot$HeightBase[plot_combos[i,1]],
                 CD=plot$CrownDepth[plot_combos[i,1]],
                 Bj=plot$Bj[plot_combos[i,1]])
    tree_b<-list(CRmax=plot$CR_average[plot_combos[i,2]],
                 CB=plot$HeightBase[plot_combos[i,2]],
                 CD=plot$CrownDepth[plot_combos[i,2]],
                 Bj=plot$Bj[plot_combos[i,2]])
    tree_list<-list(tree_a,tree_b)
    tree_pair_CCI[[i]]<-calculate_CCI(tree_list)
  }
  
  tree_pair_df<-do.call(rbind.data.frame, tree_pair_CCI)
  colnames(tree_pair_df)<-names(tree_pair_CCI[[i]])
  return(tree_pair_df)
}

plot_CCI_all<-lapply(self_pruning_split,plot_CCI)
plot_CCI_means<-lapply(plot_CCI_all,colMeans,na.rm=T)
plot_CCI_mean_df<-data.frame(do.call(rbind,plot_CCI_means))
plot_CCI_mean_df$plot<-rownames(plot_CCI_mean_df)

############################################
## applying the complementarity function
## to calculate the simulated complementarity
## of each plot mixture plot by drawing
## monoculture trees

## this procedure needs to be updated because mortalities
## in the self-pruning data are not "true" mortalities
## so we need a sampling procedure to determine whether
## to select dead trees based on actual mortality rates
## in the full community survey

plot_simulator<-function(plot,full_df) {

  ## simulate a new plot by combining the monocultures
  ## of the constituent species in the same block
  
  ## extract features of the mixture plot
  plot_comp<-unique(plot$Species)
  plot_block<-plot$Block[1]
  
  ## this works since monoculture plots have a name
  ## that is simply the one species they contain
  sim_plot<-self_pruning[which(full_df$Block==plot_block & full_df$Plot %in% plot_comp),]
  return(sim_plot)
  
}

split_nbsp<-unlist(lapply(self_pruning_split,function(x) x$nbsp[1]))
self_pruning_mix<-self_pruning_split[-which(split_nbsp==1)]

self_pruning_sim<-lapply(self_pruning_mix,
                         plot_simulator,
                         full_df=self_pruning)
plot_CCI_sim<-lapply(self_pruning_sim,plot_CCI)
