## species is a vector of species in the plot
## props_planted is a vector of proportions of those species planted within the plot
## props_alive is a vector of proportions of planted trees of the species that survive
## species, props_planted, and props_alive MUST be in the same order

plot_combos<-function(species,props_planted,props_alive,nsamp=100){
  
  if(sum(props_planted)!=1){
    stop("proportions of planted species do not sum to 1")
  }
  
  ## vector of probabilities for outcomes
  ## each species planted * prob of being still alive
  probs_vec<-c(props_planted*props_alive,props_planted*(1-props_alive))
  sample_vec<-c(species,rep(NA,times=length(species)))
  
  ## sample from sample_vec according to probs_vec
  ## twice independently to generate two columns
  ## returns species label for living tree of given species
  ## returns NA for dead tree of either species
  sp_sample1<-as.factor(sample(sample_vec,size=nsamp,replace=T,prob=probs_vec))
  sp_sample2<-as.factor(sample(sample_vec,size=nsamp,replace=T,prob=probs_vec))
  
  sp_sample_df<-data.frame(ind1=sp_sample1,
                           ind2=sp_sample2)
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
  
  sp_pairs<-plot_combos(species=mortality_plot$CodeSp,
                        props_planted=mortality_plot$prop_planted,
                        props_alive=mortality_plot$prop_alive)
  
  for(i in 1:nrow(sp_pairs)){
    
    sp1<-sp_pairs[i,1]
    sp2<-sp_pairs[i,2]
    sp_plot_sub1<-sp_plot[sp_plot$Species==sp1,]
    sp_plot_sub2<-sp_plot[sp_plot$Species==sp2,]
    
    sample1<-sample(1:nrow(sp_plot_sub1),sp_pairs[i,1])
    sample2<-sample(1:nrow(sp_plot_sub2),sp_pairs[i,2])
  }
  
   
}
