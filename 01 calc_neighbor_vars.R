setwd("C:/Users/querc/Dropbox/PostdocProjects/SelfPruning/")
source("../../DecomposingFD/R/AlphaFD.R")

library(stringr)
library(dplyr)
library(labdsv)
library(FD)
library(ggplot2)

Inventory2018<-read.csv("IDENTMontrealData/Inventory2018_cleaned.csv")
Inventory2018$Col<-str_sub(Inventory2018$Pos,1,1)
Inventory2018$Row<-as.numeric(str_sub(Inventory2018$Pos,2,2))
Inventory2018$UniqueTreeID<-paste(Inventory2018$Block,Inventory2018$Plot,Inventory2018$Pos,sep="_")

## calculate basal area in cm^2 from basal diameter
Inventory2018$BasalArea<-(Inventory2018$BasalDiam_cleaned/20)^2*pi

#############################
## finding each tree's neighbors using
## X_Pos and Y_Pos columns to set an explicit radius

neighbor.finder<-function(dat,outcome.var,sp.var,radius){
  outcome.list<-list()
  
  for(i in 1:nrow(dat)){
    xpos<-dat$X_Pos[i]
    ypos<-dat$Y_Pos[i]

    ## get neighbor IDs
    distances<-sqrt((dat$X_Pos-xpos)^2+(dat$Y_Pos-ypos[1])^2)
    neighbors<-dat[which(distances <= radius),
                   c(outcome.var,sp.var,"UniqueTreeID")]
    colnames(neighbors)<-c("outcome","species","UniqueTreeID")
    neighbors$distance<-distances[which(distances <= radius)]
    
    ## replace NAs in the outcome with 0
    ## if there was a planted tree
    neighbor.na<-which(is.na(neighbors$outcome) & !is.na(neighbors$species))
    neighbors$outcome[neighbor.na]<-0
    
    ## add to list
    outcome.list[[i]]<-neighbors
  }
  
  names(outcome.list)<-dat$UniqueTreeID
  return(outcome.list)
}

## grab basal areas of neighbors
## within radius of 2 m
neighbor.area<-neighbor.finder(Inventory2018,"BasalArea","CodeSp",radius=2)

## delete focal trees in plots we don't care about
## based on plot composition identifiers
del_species<-c("ACPL","LADE","TICO","PIAB","PISY","QURO","PIOM","BELE","BEPE","PIMA")
del_plots<-unique(Inventory2018$Plot[Inventory2018$CodeSp %in% del_species])
del_plot_pattern<-paste(paste("_",del_plots,"_",sep=""),collapse="|")
neighbor.area<-neighbor.area[-which(grepl(del_plot_pattern,names(neighbor.area)))]

## delete neighborhoods with no living trees
neighbor.total<-unlist(lapply(neighbor.area,
                              function(neighborhood) sum(neighborhood$outcome,na.rm=T)))
if(sum(neighbor.total==0)>0){
  neighbor.area<-neighbor.area[-which(neighbor.total==0)]
}

########################################
## calculate various competition indices
## including a distance-weighted sum
## and Hegyi (also incorporates focal tree size)

neighbor.comp<-unlist(lapply(neighbor.area,
                             function(neighborhood){
                               centerless<-neighborhood[-which(neighborhood$distance==0),]
                               comp.index<-sum(centerless$outcome/centerless$distance,na.rm=T)
                               return(comp.index)
                             }))

###########################################
## get community data into the right format(s)
## for calculating functional diversity

## row bind all the neighborhoods, with a new column
## that designates the ID of the focal tree
## 'centerless' objects omit the focal tree itself
neighbor.join<-bind_rows(neighbor.area, .id = "UniqueTreeID")
neighbor.join.centerless<-neighbor.join[-which(neighbor.join$distance==0),]
neighbor.join$distance<-NULL
neighbor.join.centerless$distance<-NULL

## summing trees of the same species in the same neighborhood
## because otherwise matrify won't work properly
neighbor.agg<-aggregate(neighbor.join$outcome,
                        by=list(neighbor.join$UniqueTreeID,
                                neighbor.join$species),
                        FUN=sum)

neighbor.agg.centerless<-aggregate(neighbor.join.centerless$outcome,
                                   by=list(neighbor.join.centerless$UniqueTreeID,
                                           neighbor.join.centerless$species),
                                   FUN=sum)

## note that this community matrix includes the focal
## tree of each neighborhood!
neighbor.df<-matrify(neighbor.agg)
neighbor.df.centerless<-matrify(neighbor.agg.centerless)
## make sure columns are in the same order
neighbor.df.centerless<-neighbor.df.centerless[,colnames(neighbor.df)]
## the number of species in the immediate neighborhood
neighbor.richness<-rowSums(neighbor.df>0)

#########################################
## calculate functional diversity

trait_summary<-read.csv("TraitData/trait_summary.csv")
trait_summary<-trait_summary[match(colnames(neighbor.df),trait_summary$SpeciesCode),]

core_traits<-trait_summary[,c("LDMC","N","SLA","SRL","WD")]
rownames(core_traits)<-trait_summary$SpeciesCode

## calculate CWM of first trait PCA component, omitting the central tree
## this gives as the functional identity of the neighborhood
neighbor.prop<-neighbor.df.centerless/rowSums(neighbor.df.centerless)
neighbor.CWM1<-as.matrix(neighbor.prop) %*% trait_summary$PC1

## for trees with no neighbors, assign a neighborhood functional ID of 0
## generally not needed if neighborhood is greater than just immediate neighbors
# neighbor.CWM1[is.na(neighbor.CWM1)]<-0

## calculate FDis using FD package
## including the central tree
neighbor.FDis<-fdisp(dist(core_traits),
                     a = as.matrix(neighbor.df))

## calculate Scheiner's metrics
neighbor.FTD<-FTD.comm(tdmat=dist(core_traits),
                       spmat = as.matrix(neighbor.df),
                       abund=T)

neighbor.data<-data.frame(UniqueTreeID=names(neighbor.area),
                          neighbor.richness=neighbor.richness,
                          NCI=neighbor.comp,
                          neighbor.FI1=neighbor.CWM1,
                          FDis=neighbor.FDis$FDis,
                          qDTM=neighbor.FTD$com.FTD$qDTM)

## to do:
## calculate CWMs of major traits

## write data
write.csv(neighbor.data,"IDENTMontrealData/neighborhood_vars.csv")

#################################
## output mortalities

## need to resolve these issues
## and maybe do some data cleaning (comparisons to 2017 and 2019?)
Inventory2018[which(Inventory2018$StateDesc!="Dead" & is.na(Inventory2018$BasalArea)),]
Inventory2018[which(Inventory2018$StateDesc=="Dead" & !is.na(Inventory2018$BasalArea)),]

## ignoring that for now
## count a tree as dead if its status is "dead"
Inventory2018$Alive<-ifelse(Inventory2018$StateDesc!="Dead",yes=1,no=0)
DB_mortality<-aggregate(Alive~Block+Plot+CodeSp,
                        data=Inventory2018,
                        FUN=mean)

write.csv(DB_mortality,"IDENTMontrealData/mortality_2018.csv")
