setwd("C:/Users/Shan Kothari/Dropbox/PostdocProjects/SelfPruning/")
source("../../DecomposingFD/R/AlphaFD.R")

library(stringr)
library(dplyr)
library(labdsv)
library(FD)
library(ggplot2)

Inventory2018<-read.csv("IDENTMontrealData/Inventory2018_cleaned.csv")

## extract row and column from position indicator
Inventory2018$Col<-str_sub(Inventory2018$Pos,1,1)
Inventory2018$Row<-as.numeric(str_sub(Inventory2018$Pos,2,2))

#############################
## finding each tree's neighbors using
## X_Pos and Y_Pos columns to set an explicit radius

neighbor.finder<-function(dat,outcome.var,sp.var,radius){
  outcome.list<-list()
  
  for(i in 1:nrow(dat)){
    xpos<-dat$X_Pos[i]
    ypos<-dat$Y_Pos[i]

    ## get IDs of neighbors (within radius)
    distances<-sqrt((dat$X_Pos-xpos)^2+(dat$Y_Pos-ypos[1])^2)
    neighbors<-dat[which(distances <= radius),
                   c(outcome.var,sp.var,"unique_tree")]
    colnames(neighbors)<-c("outcome","species","unique_tree")
    neighbors$distance<-distances[which(distances <= radius)]
    
    ## replace NAs in the outcome with 0
    ## if there was a planted tree
    neighbor.na<-which(is.na(neighbors$outcome) & !is.na(neighbors$species))
    neighbors$outcome[neighbor.na]<-0
    
    ## add to list
    outcome.list[[i]]<-neighbors
  }
  
  names(outcome.list)<-dat$unique_tree
  return(outcome.list)
}

#####################################
## calculate NCI

## grab basal areas of neighbors
## within radius of 1.1 m
## this is the 90th percentile of measured crown radii
## (concatenating CR1, CR2, CR3, CR4 from self pruning data)
neighbor.area.NCI<-neighbor.finder(Inventory2018,"BasalArea","CodeSp",radius=1.1)

## delete focal trees in plots we don't care about
## based on plot composition identifiers
del_species<-c("ACPL","LADE","TICO","PIAB","PISY","QURO","PIOM","BELE","BEPE","PIMA")
del_plots<-unique(Inventory2018$Plot[Inventory2018$CodeSp %in% del_species])
del_plot_pattern<-paste(paste("_",del_plots,"_",sep=""),collapse="|")
neighbor.area.NCI<-neighbor.area.NCI[-which(grepl(del_plot_pattern,names(neighbor.area.NCI)))]

## a function to calculate a simplified NCI
## (without species-specific parameters)
simp.NCI<-function(neighborhood,num_exp=1,den_exp=1){
  
  ## ignore the focal tree itself
  cl<-neighborhood[-which(neighborhood$distance==0),]
  comp.index<-sum(cl$outcome^num_exp/cl$distance^den_exp,na.rm=T)
  return(comp.index)
}

## calculate the simplified NCI for each neighborhood
neighbor.comp<-unlist(lapply(neighbor.area.NCI,simp.NCI))

###########################################
## calculate functional diversity and CWMs

## we divide this into a separate section in case we
## want to use different radii and outcome variables
## to calculate NCI and neighborhood function

neighbor.area.FD<-neighbor.finder(Inventory2018,"BasalArea","CodeSp",radius=1.1)

## get rid of the same plots we dropped for NCI calculations
neighbor.area.FD<-neighbor.area.FD[-which(grepl(del_plot_pattern,names(neighbor.area.FD)))]

## row bind all the neighborhoods, with a new column
## that designates the ID of the focal tree
## 'cl' objects omit the focal tree itself (centerless)
neighbor.join<-bind_rows(neighbor.area.FD, .id = "unique_tree")
neighbor.join.cl<-neighbor.join[-which(neighbor.join$distance==0),]

## for analyses based on numbers of
## planted trees rather than sizes
neighbor.join$planted<-1
neighbor.join.cl$planted<-1

## summing trees of the same species in the same neighborhood
## because otherwise matrify won't work properly
neighbor.agg.num<-aggregate(neighbor.join$outcome,
                            by=list(neighbor.join$unique_tree,
                                    neighbor.join$species),
                            FUN=sum)

# neighbor.agg.cl<-aggregate(neighbor.join.cl$outcome,
#                            by=list(neighbor.join.cl$unique_tree,
#                                    neighbor.join.cl$species),
#                            FUN=sum)

neighbor.agg.clnum<-aggregate(neighbor.join$outcome,
                              by=list(neighbor.join$unique_tree,
                                      neighbor.join$species),
                              FUN=sum)

## note that this community matrix includes the focal
## tree of each neighborhood!
neighbor.df.num<-matrify(neighbor.agg.num)
# neighbor.df.cl<-matrify(neighbor.agg.cl)
neighbor.df.clnum<-matrify(neighbor.agg.clnum)

## make sure columns are in the same order
neighbor.df.clnum<-neighbor.df.clnum[,colnames(neighbor.df.num)]

## the number of species planted in the immediate neighborhood
neighbor.richness<-rowSums(neighbor.df.num>0)

#########################################
## calculate functional diversity

trait_summary<-read.csv("TraitData/trait_summary.csv")
trait_summary<-trait_summary[match(colnames(neighbor.df.num),trait_summary$SpeciesCode),]

core_traits<-trait_summary[,c("LDMC","N","SLA","SRL","WD")]
rownames(core_traits)<-trait_summary$SpeciesCode

## calculate CWM of first trait PCA component, omitting the central tree
## this gives as the functional identity of the neighborhood
neighbor.prop<-neighbor.df.clnum/rowSums(neighbor.df.clnum)
neighbor.CWM1<-as.matrix(neighbor.prop) %*% trait_summary$PC1

## for trees with no neighbors, assign a neighborhood functional ID of 0
## generally not needed if neighborhood is greater than just immediate neighbors
# neighbor.CWM1[is.na(neighbor.CWM1)]<-0

## calculate FDis using FD package, including the central tree
## note that this will throw up an error if any tree has no 
## living trees in its neighborhood
# neighbor.FDis<-fdisp(dist(core_traits),
#                      a = as.matrix(neighbor.df.num))

## calculate Scheiner's metrics
neighbor.FTD<-FTD.comm(tdmat=dist(core_traits),
                       spmat = as.matrix(neighbor.df.num),
                       abund=T)

neighbor.data<-data.frame(unique_tree=names(neighbor.area.NCI),
                          neighbor.richness=neighbor.richness,
                          NCI=neighbor.comp,
                          neighbor.FI1=neighbor.CWM1,
                          # FDis=neighbor.FDis$FDis,
                          qDTM=neighbor.FTD$com.FTD$qDTM)

## write data
# write.csv(neighbor.data,"IDENTMontrealData/neighborhood_vars.csv")
