## species is a vector of species in the plot
## nums_planted is a vector of numbers of those species planted within the plot
## (we use raw numbers planted rather than proportions because we may
## want to create simulated plots by combining monocultures
## so if you just combine those plots via rbind, numbers may exceed 1)
## props_alive is a vector of proportions of planted trees of the species that survive
## species, props_planted, and props_alive MUST be in the same order

## it's not strictly necessary to sample trees like this as pairs
## rather than pairing them up later, but I find that it makes things
## more intuitive, even if a bit more complicated later on

plot_combos<-function(species,nums_planted,props_alive,n_pairs=100){
  
  props_planted<-nums_planted/sum(nums_planted)
  
  ## vector of probabilities for outcomes
  ## each species planted * prob of being still alive
  probs_vec<-c(props_planted*props_alive,props_planted*(1-props_alive))
  sample_vec<-c(species,rep(NA,times=length(species)))
  
  ## sample from sample_vec according to probs_vec
  ## returns species label for living tree of given species
  ## returns NA for dead tree of either species
  sp_sample<-as.factor(sample(sample_vec,size=n_pairs*2,
                              replace=T,prob=probs_vec))
  
  sp_sample_df<-data.frame(ind1=sp_sample[1:n_pairs],
                           ind2=sp_sample[(n_pairs+1):(2*n_pairs)])
  return(sp_sample_df)
  
}


##################################
## for a given plot, calculate CCI

## for testing
unique_plot<-self_pruning_split$D.4N2
mortality_plot<-mortality_2018[mortality_2018$Block=="D" & mortality_2018$Plot=="4N2",]

## sp_plot is used for drawing information related to crown shape
## and would usually be in the same format as self_pruning

## mortality_plot is used for drawing information related to mortality
## it must have one row per species with columns 'CodeSp',
## 'prop_planted', and 'prop_alive'

plot_CCI<-function(sp_plot,mortality_plot){
  
  ## generate 100 pairs of two species sampled
  ## based on real planting numbers and mortality rates
  sp_pairs<-plot_combos(species=mortality_plot$CodeSp,
                        nums_planted=mortality_plot$num_planted,
                        props_alive=mortality_plot$prop_alive)
  
  ## output object with outcomes
  tree_pair_list<-list()
  
  for(i in 1:nrow(sp_pairs)){
    
    print(i)
    
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
      
      tree_pair_list[[i]]<-NA

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
  
  tree_pair_df<-do.call(rbind.data.frame, tree_pair_list)
  colnames(tree_pair_df)<-names(tree_pair_list[[1]])
  return(tree_pair_df)
}
