setwd("C:/Users/querc/Dropbox/PostdocProjects/SelfPruning")

library(mosaic)

########################################
## data on crown shape and size

## this self-pruning data includes random sample of living trees
## but has no information on mortality rates, which we will have
## to account for later
self_pruning<-read.csv("SelfPruningData/self_pruning_processed.csv")

## in these three plots, there are species that were living
## in the plot in 2018 that were accidentally not measured 
## in the self-pruning survey. in A_4N8, THOC is missing;
## in D_4N3, PIST is missing; in A_2N3, BEPA (misplanted)
## is missing. these issues would otherwise threaten to
## bias crown packing and complementarity calculations
self_pruning<-subset(self_pruning,!(unique_plot %in% c("A_2N3","A_4N8","D_4N3")))

## get indicators of species composition
self_pruning_list<-split(self_pruning,f = self_pruning$Plot)
self_pruning_comp<-unlist(lapply(self_pruning_list,
                                 function(plot) paste(unique(plot$Species),collapse="|")))
self_pruning$sp_comp<-self_pruning_comp[match(self_pruning$Plot,names(self_pruning_comp))]

## parameters from Purves et al. paper
C0B<-0.196
C1B<-0.511

species.list<-c("ABBA","ACRU","ACSA","BEAL","BEPA","LALA",
                "PIGL","PIRE","PIRU","PIST","QURU","THOC")
Tj<-c(0.278,0.536,0.560,0.592,0.435,0.308,
      0.278,0.287,0.324,0.417,0.538,0.251)

crown_shape<-data.frame(species=species.list,Tj=Tj)
crown_shape$Bj<-(1-crown_shape$Tj)*C0B+crown_shape$Tj*C1B
self_pruning$Bj<-crown_shape$Bj[match(self_pruning$Species,crown_shape$species)]

## formulas for area/volume of crown
crown_area<-function(CD,CR,beta){(CR*CD)/(beta+1)}
crown_vol<-function(CD,CR,beta){(pi*CR^2*CD)/(2*beta+1)}

## calculate for each row of the self-pruning data
self_pruning$crown_vol<-crown_vol(CD=self_pruning$CrownDepth/100,
                                  CR=self_pruning$CR_average/100,
                                  beta=self_pruning$Bj)

## read in mortality data
mortality_2018<-read.csv("IDENTMontrealData/mortality_2018.csv")
mortality_2018$unique_plot<-paste(mortality_2018$Block,
                                 mortality_2018$Plot,
                                 sep="_")

mortality_2018$plot_sp<-paste(mortality_2018$unique_plot,
                              mortality_2018$CodeSp,
                              sep="_")

## proportion of species planted within the inner 36 and
## proportion of planted individuals of the species still alive
mortality_2018$prop_planted<-mortality_2018$num_planted/36
mortality_2018$prop_alive<-mortality_2018$num_alive/mortality_2018$num_planted

######################################
## calculate canopy packing

## take the mean crown volume for each species in each plot
crown_vol_agg<-aggregate(crown_vol~Block+Plot+Species+nbsp,
                         data=self_pruning,FUN=mean)

colnames(crown_vol_agg)<-c("Block","Plot","Species","Richness","crown_vol")

## columns to match against for adding mortality and monoculture data
crown_vol_agg$block_plot<-paste(crown_vol_agg$Block,
                                crown_vol_agg$Plot,
                                sep="_")

crown_vol_agg$block_sp<-paste(crown_vol_agg$Block,
                              crown_vol_agg$Species,
                              sep="_")

crown_vol_agg$plot_sp<-paste(crown_vol_agg$block_plot,
                             crown_vol_agg$Species,
                             sep="_")

## attach numbers alive and planted in mixture
mortality_match<-match(crown_vol_agg$plot_sp,mortality_2018$plot_sp)
crown_vol_agg$num_alive<-mortality_2018$num_alive[mortality_match]
crown_vol_agg$num_planted<-mortality_2018$num_planted[mortality_match]
crown_vol_agg$prop_alive<-mortality_2018$prop_alive[mortality_match]

## attach mean crown volume and mortality rates in monoculture
## here we take advantage of the fact that monocultures have a
## plot name that is the same as their species code
mono_match<-match(crown_vol_agg$block_sp,crown_vol_agg$block_plot)
crown_vol_agg$crown_vol_mono<-crown_vol_agg$crown_vol[mono_match]
crown_vol_agg$prop_alive_mono<-crown_vol_agg$prop_alive[mono_match]

## multiply by proportions alive to account for mortality
crown_vol_agg$crown_vol_mono_live<-crown_vol_agg$crown_vol_mono*crown_vol_agg$prop_alive_mono
crown_vol_agg$crown_vol_live<-crown_vol_agg$crown_vol*crown_vol_agg$prop_alive
crown_vol_agg$total_sp_crown_vol<-crown_vol_agg$crown_vol_live*crown_vol_agg$num_planted

##########################################
## simulations of canopy packing holding crown depth constant

## number of crowns to simulate for each species in each plot
n_crowns<-250

## first, we simulate trees where the crown radius is sampled from
## species in the given plot and crown depth is sampled from
## the corresponding monocultures
sim_crown_vols<-list()
for(i in crown_vol_agg$plot_sp){

  plot_sp_split<-strsplit(i,"_")[[1]]
  block<-plot_sp_split[1]
  plot<-plot_sp_split[2]
  species<-plot_sp_split[3]
  
  ## selecting mixtures & monocultures
  mix_ids<-which(self_pruning$Block==block & self_pruning$Plot==plot & self_pruning$Species==species)
  ## here we take advantage of the fact that monoculture plot names are species names
  mono_ids<-which(self_pruning$Block==block & self_pruning$Plot==species & self_pruning$Species==species)
  
  ## draw crown radius from mixtures
  ## but crown depth from monocultures
  sim_plot_sp<-data.frame(CR_average=sample(self_pruning$CR_average[mix_ids],
                                            size=n_crowns,replace=T),
                          CrownDepth=sample(self_pruning$CrownDepth[mono_ids],
                                            size=n_crowns,replace=T),
                          Bj=sample(self_pruning$Bj[mono_ids],
                                    size=n_crowns,replace=T))
  
  sim_crown_vols[i]<-mean(crown_vol(CD=sim_plot_sp$CrownDepth/100,
                                    CR=sim_plot_sp$CR_average/100,
                                    beta=sim_plot_sp$Bj))
  
}

## next, we simulate trees where the both the crown radius
## and the crown depth are sampled from species in the given plot

## the purpose of this simulation is to test sensitivity to
## the potential bias caused by random matching of crown radius and crown depth
null_sim_crown_vols<-list()
for(i in crown_vol_agg$plot_sp){
  
  plot_sp_split<-strsplit(i,"_")[[1]]
  block<-plot_sp_split[1]
  plot<-plot_sp_split[2]
  species<-plot_sp_split[3]
  
  mix_ids<-with(self_pruning,which(Block==block & Plot==plot & Species==species))
  mono_ids<-with(self_pruning,which(Block==block & Plot==species & Species==species))
  
  ## draw crown radius and crown depth from mixtures
  null_sim_plot_sp<-data.frame(CR_average=sample(self_pruning$CR_average[mix_ids],
                                            size=n_crowns,replace=T),
                          CrownDepth=sample(self_pruning$CrownDepth[mix_ids],
                                            size=n_crowns,replace=T),
                          Bj=sample(self_pruning$Bj[mono_ids],
                                    size=n_crowns,replace=T))
  
  null_sim_crown_vols[i]<-mean(crown_vol(CD=null_sim_plot_sp$CrownDepth/100,
                                    CR=null_sim_plot_sp$CR_average/100,
                                    beta=null_sim_plot_sp$Bj))
  
}

## multiply mean volume by mixture mortality rate
crown_vol_agg$sim_vol<-unlist(sim_crown_vols)
crown_vol_agg$sim_vol_live<-crown_vol_agg$sim_vol*crown_vol_agg$prop_alive
crown_vol_agg$total_sp_sim_vol<-crown_vol_agg$sim_vol_live*crown_vol_agg$num_planted

crown_vol_agg$null_sim_vol<-unlist(null_sim_crown_vols)
crown_vol_agg$null_sim_vol_live<-crown_vol_agg$null_sim_vol*crown_vol_agg$prop_alive
crown_vol_agg$total_sp_null_sim_vol<-crown_vol_agg$null_sim_vol_live*crown_vol_agg$num_planted

######################################
## aggregating to the plot scale

## aggregating overyielding to the species x plot level
## rather than just individuals
crown_vol_agg$OY_actual<-with(crown_vol_agg,
                              num_planted*(crown_vol_live-crown_vol_mono_live))
crown_vol_agg$OY_sim<-with(crown_vol_agg,
                           num_planted*(sim_vol_live-crown_vol_mono_live))
crown_vol_agg$OY_null_sim<-with(crown_vol_agg,
                           num_planted*(null_sim_vol_live-crown_vol_mono_live))

## and aggregating to the whole plot (inner 6 x 6 trees)
## since planting numbers / mortalities are only calculated
## within those inner 6 x 6
crown_vol_agg_sub<-crown_vol_agg[,c("Block","Plot","Richness",
                                    "OY_actual","OY_sim","OY_null_sim",
                                    "total_sp_crown_vol",
                                    "total_sp_sim_vol",
                                    "total_sp_null_sim_vol")]

crown_vol_plot<-aggregate(.~Block+Plot+Richness,
                          data=crown_vol_agg_sub,
                          FUN=sum)

crown_vol_plot$OY_actual[which(crown_vol_plot$Richness==1)]<-NA
crown_vol_plot$OY_sim[which(crown_vol_plot$Richness==1)]<-NA
crown_vol_plot$OY_null_sim[which(crown_vol_plot$Richness==1)]<-NA

## all values should be divided by 9 for analysis per m^2
## since upscaling was done to the 6 x 6 inner plot (no edge)

## read in plot-level heterogeneity measures
plot_vars<-read.csv("IDENTMontrealData/plot_vars.csv")

crown_vol_plot$unique_plot<-paste(crown_vol_plot$Block,crown_vol_plot$Plot,sep="_")
crown_vol_plot$FDis<-plot_vars$FDis[match(crown_vol_plot$unique_plot,plot_vars$unique_plot)]
crown_vol_plot$FTD<-plot_vars$FTD[match(crown_vol_plot$unique_plot,plot_vars$unique_plot)]

##########################################
## plotting canopy packing results



######################################
## function to calculate complementarity
## for a pair of trees

source("Scripts/self-pruning-analysis/complementarity_functions.R")

calculate_CCI<-function(tree1,tree2){
  area1<-crown_area(CD=tree1$CD,CR=tree1$CRmax,beta=tree1$Bj)
  area2<-crown_area(CD=tree2$CD,CR=tree2$CRmax,beta=tree2$Bj)
  overlap_2D<-calculate_2D_overlap(list(tree1,tree2))
  
  vol1<-crown_vol(CD=tree1$CD,CR=tree1$CRmax,beta=tree1$Bj)
  vol2<-crown_vol(CD=tree2$CD,CR=tree2$CRmax,beta=tree2$Bj)
  overlap_3D<-calculate_3D_overlap(list(tree1,tree2))
  
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
    
    ## no need to define CCI_2D and CCI_3D explicitly
    ## because they can be calculated
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

tree_list_1<-list(CRmax=tree_samp_1$CR_average,
                  CD=tree_samp_1$CrownDepth,
                  CB=tree_samp_1$HeightBase,
                  Bj=tree_samp_1$Bj)

tree_list_2<-list(CRmax=tree_samp_2$CR_average,
                  CD=tree_samp_2$CrownDepth,
                  CB=tree_samp_2$HeightBase,
                  Bj=tree_samp_2$Bj)

calculate_CCI(tree_list_1,tree_list_2)

#####################################
## applying the complementarity function
## to calculate the complementarity of each plot

## split self-pruning data by unique plot
self_pruning_split<-split(self_pruning,f = ~unique_plot)

## split mortality data by unique plot, and
## ensure that retained plots (and their order)
## match the split self-pruning data
mortality_split<-split(mortality_2018,f = ~unique_plot)[names(self_pruning_split)]

## check that plots are in the same order; should be 69
sum(names(self_pruning_split)==names(mortality_split))

## species is a vector of species in the plot
## nums_alive is a vector of numbers of surviving trees of each species in the plot
## nums_planted is a vector of numbers of those species planted within the plot
## species, props_planted, and props_alive MUST be in the same order

## (NOTE: we use raw numbers planted rather than proportions because we may
## want to create simulated plots by combining monocultures
## so if you just combine those plots via rbind, numbers may exceed 1)

plot_combos<-function(species,nums_alive,nums_planted,n_pairs=100){
  
  ## calculate proportion of living trees of each
  ## species out of total trees of all species 
  props_alive_sp<-nums_alive/sum(nums_planted)
  
  ## proportion of dead trees (of any species)
  ## out of those planted
  prop_dead<-1-sum(props_alive_sp)
  
  ## vector of probabilities for outcomes
  probs_vec<-c(props_alive_sp,prop_dead)
  sample_vec<-c(species,NA)
  
  ## sample from sample_vec according to probs_vec
  ## returns species label for living tree of given species
  ## returns NA for dead tree of either species
  sp_sample<-as.factor(sample(x=sample_vec,
                              size=n_pairs*2,
                              replace=T,
                              prob=probs_vec))
  
  sp_sample_df<-data.frame(ind1=sp_sample[1:n_pairs],
                           ind2=sp_sample[(n_pairs+1):(2*n_pairs)])
  
  return(sp_sample_df)
}

## sp_plot is used for drawing information related to crown shape
## and would usually be in the same format as self_pruning

## mortality_plot is used for drawing information related to mortality
## it must have one row per species with columns 'CodeSp',
## 'num_alive', and 'num_planted'

plot_CCI<-function(sp_plot,mortality_plot,n_pairs=100){
  
  ## generate 100 pairs of two species sampled
  ## based on real planting numbers and mortality rates
  sp_pairs<-plot_combos(species=mortality_plot$CodeSp,
                        nums_alive=mortality_plot$num_alive,
                        nums_planted=mortality_plot$num_planted,
                        n_pairs=n_pairs)
  
  ## output object with outcomes
  tree_pair_list<-list()
  
  for(i in 1:nrow(sp_pairs)){
    
    ## extract the two species from this row
    sp1<-sp_pairs[i,1]
    sp2<-sp_pairs[i,2]
    
    ## if both sampled trees are dead (NA)
    if(is.na(sp1) & is.na(sp2)){
      
      tree_a<-list(CRmax=0,CB=0,CD=0,Bj=0)
      tree_b<-list(CRmax=0,CB=0,CD=0,Bj=0)
      tree_pair_list[[i]]<-calculate_CCI(tree_a,tree_b)
      
    }
    
    ## if the first tree is dead,
    ## we sample the second tree from
    ## measured individuals of the species
    ## in the plot
    if(is.na(sp1) & !is.na(sp2)){
      
      tree_a<-list(CRmax=0,CB=0,CD=0,Bj=0)
      
      sp_plot_sub2<-sp_plot[sp_plot$Species==sp2,]
      sample2<-sample(1:nrow(sp_plot_sub2),size=1)
      tree_b<-list(CRmax=sp_plot_sub2$CR_average[sample2],
                   CB=sp_plot_sub2$HeightBase[sample2],
                   CD=sp_plot_sub2$CrownDepth[sample2],
                   Bj=sp_plot_sub2$Bj[sample2])
      
      tree_pair_list[[i]]<-calculate_CCI(tree_a,tree_b)
      
    }
    
    ## if the second tree is dead
    ## we sample the first tree from
    ## measured individuals of the species
    ## in the plot
    if(!is.na(sp1) & is.na(sp2)){
      
      sp_plot_sub1<-sp_plot[sp_plot$Species==sp1,]
      sample1<-sample(1:nrow(sp_plot_sub1),size=1)
      tree_a<-list(CRmax=sp_plot_sub1$CR_average[sample1],
                   CB=sp_plot_sub1$HeightBase[sample1],
                   CD=sp_plot_sub1$CrownDepth[sample1],
                   Bj=sp_plot_sub1$Bj[sample1])
      
      tree_b<-list(CRmax=0,CB=0,CD=0,Bj=0)
      
      tree_pair_list[[i]]<-calculate_CCI(tree_a,tree_b)
      
    }
    
    ## if neither tree is dead, we sample both trees
    if(!is.na(sp1) & !is.na(sp2)){
      
      sp_plot_sub1<-sp_plot[sp_plot$Species==sp1,]
      sp_plot_sub2<-sp_plot[sp_plot$Species==sp2,]
      sample1<-sample(1:nrow(sp_plot_sub1),size=1)
      sample2<-sample(1:nrow(sp_plot_sub2),size=1)
      
      ## if the same individual is sampled twice above
      ## redo the sampling
      while(sp1==sp2 & sample1==sample2){
        sample1<-sample(1:nrow(sp_plot_sub1),size=1)
        sample2<-sample(1:nrow(sp_plot_sub2),size=1)
      }
      
      tree_a<-list(CRmax=sp_plot_sub1$CR_average[sample1],
                   CB=sp_plot_sub1$HeightBase[sample1],
                   CD=sp_plot_sub1$CrownDepth[sample1],
                   Bj=sp_plot_sub1$Bj[sample1])
      
      tree_b<-list(CRmax=sp_plot_sub2$CR_average[sample2],
                   CB=sp_plot_sub2$HeightBase[sample2],
                   CD=sp_plot_sub2$CrownDepth[sample2],
                   Bj=sp_plot_sub2$Bj[sample2])
      
      tree_pair_list[[i]]<-calculate_CCI(tree_a,tree_b)
      
    }
  }
  
  ## turn to data frame for output
  tree_pair_df<-do.call(rbind.data.frame, tree_pair_list)
  colnames(tree_pair_df)<-names(tree_pair_list[[1]])
  return(tree_pair_df)
  
}

plot_CCI_all<-lapply(1:length(self_pruning_split),function(i){
  print(i)
  return(plot_CCI(self_pruning_split[[i]],mortality_split[[i]]))
})

plot_CCI_means<-lapply(plot_CCI_all,colMeans,na.rm=T)
plot_CCI_df<-data.frame(do.call(rbind,plot_CCI_means))
plot_CCI_df$unique_plot<-names(self_pruning_split)

############################################
## applying the complementarity function
## to calculate the simulated complementarity
## of each mixture plot by drawing
## monoculture trees

## simulate a new plot by combining the monocultures
## of the constituent species in the same block
plot_simulator<-function(plot,sp_df,mortality_df) {
  
  ## extract species composition and block of plot
  plot_comp<-unique(plot$Species)
  plot_block<-plot$Block[1]
  plot_name<-plot$Plot[1]
  
  ## this works since monoculture plots have a name
  ## that is simply the one species they contain
  sim_sp_match<-which(sp_df$Block==plot_block & sp_df$Plot %in% plot_comp)
  sim_sp<-sp_df[sim_sp_match,]
  
  sim_mortality_match<-with(mortality_df,which(Block==plot_block & Plot %in% plot_comp))
  sim_mortality<-mortality_df[sim_mortality_match,]
  
  orig_sp_match<-which(sp_df$Block==plot_block & sp_df$Plot==plot_name)
  orig_sp<-sp_df[orig_sp_match,]
  
  orig_mortality_match<-with(mortality_df,which(Block==plot_block & Plot==plot_name))
  orig_mortality<-mortality_df[orig_mortality_match,]
  
  sim_all<-list(sim_sp=sim_sp,
                sim_mortality=sim_mortality,
                orig_sp=orig_sp,
                orig_mortality=orig_mortality)
  
  return(sim_all)
  
}

## optionally can run on just mixture plots
# split_nbsp<-unlist(lapply(self_pruning_split,function(x) x$nbsp[1]))
# self_pruning_mix<-self_pruning_split[-which(split_nbsp==1)]

self_pruning_sim<-lapply(self_pruning_split,
                         plot_simulator,
                         sp_df=self_pruning,
                         mortality_df=mortality_2018)

plot_CCI_sim<-lapply(self_pruning_sim,function(x) plot_CCI(sp_plot=x$sp_plot,
                                                           mortality_plot = x$mortality_plot))
plot_CCI_sim_means<-lapply(plot_CCI_sim,colMeans,na.rm=T)
plot_CCI_sim_df<-data.frame(do.call(rbind,plot_CCI_sim_means))
plot_CCI_sim_df$unique_plot<-names(self_pruning_sim)

# sim_match<-match(plot_CCI_df$unique_plot,plot_CCI_sim_df$unique_plot)
# plot_CCI_df$sim_CCI_2D<-plot_CCI_sim_df$CCI_2D[sim_match]
# plot_CCI_df$sim_CCI_3D<-plot_CCI_sim_df$CCI_3D[sim_match]
# plot_CCI_df$sim_CCI_min_2D<-plot_CCI_sim_df$CCI_min_2D[sim_match]
# plot_CCI_df$sim_CCI_min_3D<-plot_CCI_sim_df$CCI_min_3D[sim_match]

################################################
## applying the complementarity function
## to calculate the simulated complementarity
## of each mixture plot by drawing mixture trees
## but with crown depth from monoculture

## here we need a function to pull both monoculture and mixture plots

## then sample 