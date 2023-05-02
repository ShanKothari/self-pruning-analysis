setwd("C:/Users/querc/Dropbox/PostdocProjects/SelfPruning/")
source("../../DecomposingFD/R/AlphaFD.R")

library(stringr)
library(dplyr)
library(labdsv)
library(FD)
library(ggplot2)

## to dos:
## increase NCI radius (to 2 m or so?)
## use a radius-based approach rather than a grid-based approach

DB_community<-read.csv("IDENTMontrealData/Inventory2018.csv")
DB_community$Col<-str_sub(DB_community$Pos,1,1)
DB_community$Row<-as.numeric(str_sub(DB_community$Pos,2,2))
DB_community$UniqueTreeID<-paste(DB_community$Block,DB_community$Plot,DB_community$Pos,sep="_")

## calculate basal area from basal diameter
DB_community$BasalArea<-(DB_community$BasalDiam/2)^2*pi

## plot composition identifiers for those that include
## species we don't care about for the self-pruning paper
## none of these are neighbors of focal species
del_species<-c("ACPL","LADE","TICO","PIAB","PISY","QURO","PIOM","BELE","BEPE","PIMA")
del_plots<-unique(DB_community$Plot[DB_community$CodeSp %in% del_species])
DB_community<-DB_community[-which(DB_community$Plot %in% del_plots),]

let2num<-function(let) return(match(let,LETTERS[1:26]))
num2let<-function(num) return(ifelse(num<=0,NA,LETTERS[num]))

#############################
## finding each tree's immediate neighbors
## by looking for those a certain number of
## letters and/or numbers away

## an alternate approach would be to use the
## X_Pos and Y_Pos columns to set an explicit radius

neighbor.id<-function(row,column,radius=1){
  
  rows_out<-rep((row-radius):(row+radius),each=radius*2+1)
  col_nums<-rep(let2num(column)+((-1*radius):radius),radius*2+1)
  cols_out<-sapply(col_nums,num2let)
  
  row_diff<-rows_out-row
  col_diff<-col_nums-let2num(column)
    
  return(data.frame(rows=rows_out,
                    cols=cols_out,
                    row_diff=row_diff,
                    col_diff=col_diff))
}

## this function uses the neighbor.id function to
## pull the outcome variable of interest for the trees
## a certain radius around the focal tree
neighbor.finder<-function(dat,outcome.var,sp.var,radius=1,grid.size=0.5){
  outcome.list<-list()
  
  for(i in 1:nrow(dat)){
    block<-dat$Block[i]
    plot<-dat$Plot[i]
    row<-dat$Row[i]
    col<-dat$Col[i]

    ## get neighbor IDs
    neighbors<-neighbor.id(row=row,column=col,radius=radius)
    neighbors$cols<-factor(neighbors$cols,levels=LETTERS[1:8])
    neighbors$blocks<-block
    neighbors$plots<-plot
    
    neighbors$position<-paste0(neighbors$cols,neighbors$rows)
    neighbors$UniqueTreeID<-paste(neighbors$blocks,
                                  neighbors$plots,
                                  neighbors$position,
                                  sep="_")
    
    ## extract outcome variables and species
    out<-sapply(1:nrow(neighbors),
                function(x) which(dat$UniqueTreeID==neighbors$UniqueTreeID[x]))
    out[which(out=="integer(0)")]<-NA
    out<-unlist(out)
    
    outcome.vector<-dat[out,outcome.var]
    outcome.species<-as.factor(dat[out,sp.var])
    
    ## replace NAs in the outcome with 0
    ## if there was a planted tree
    outcome.vector[is.na(outcome.vector) & !is.na(outcome.species)]<-0

    ## calculate distances to each neighbor
    outcome.dist<-grid.size*sqrt(neighbors$row_diff^2+neighbors$col_diff^2)
    outcome.list[[i]]<-data.frame(distance=outcome.dist,
                                  species=outcome.species,
                                  outcome=outcome.vector)
  }
  
  names(outcome.list)<-dat$UniqueTreeID
  return(outcome.list)
}

## grab basal areas of neighbors
neighbor.area<-neighbor.finder(DB_community,"BasalArea","CodeSp")

## delete neighborhoods with no living trees
neighbor.total<-unlist(lapply(neighbor.area,
                              function(neighborhood) sum(neighborhood$outcome,na.rm=T)))
neighbor.area<-neighbor.area[-which(neighbor.total==0)]

########################################
## calculate various competition indices
## including a distance-weighted sum
## and Hegyi (also incorporates focal tree size)

neighbor.comp<-lapply(neighbor.area,
                      function(neighborhood){
                        center<-neighborhood[which(neighborhood$distance==0),]
                        centerless<-neighborhood[-which(neighborhood$distance==0),]
                        comp.index<-sum(centerless$outcome/centerless$distance,na.rm=T)
                        Hegyi.index<-comp.index/center$outcome
                        return(list(comp.index=comp.index,
                                    Hegyi.index=Hegyi.index))
                      })

neighbor.comp.df<-data.frame(comp.index=unlist(lapply(neighbor.comp, function(x) x$comp.index)),
                             Hegyi.index=unlist(lapply(neighbor.comp, function(x) x$Hegyi.index)))

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

#######################################
## calculate index of functional identity from trait data

## read in trait data and subset by relevant species and traits
traits<-read.csv("TraitData/IDENT_TRAIT_DATABASE_2020-10-20.csv")
traits<-traits[-which(traits$SpeciesCode %in% del_species),]

traits<-traits[,c("SpeciesCode","LL","LDMC","Leaf_N_mass","SLA..all.include.",
                  "SRL..fine.root.","SSD...WD")]

## match trait row order to columns of community matrix
traits<-traits[match(colnames(neighbor.df),traits$SpeciesCode),]
rownames(traits)<-as.character(traits$SpeciesCode)
traits$SpeciesCode<-NULL

## just extract whether species are deciduous or evergreen
## for plotting purposes (but not for the PCA)
leaf_habit<-ifelse(traits$LL>12,
                   yes = "Evergreen",
                   no = "Deciduous")
traits$LL<-NULL

## do PCA of traits
trait.pca<-prcomp(traits,scale. = T)
trait.pca.scores<-data.frame(species = rownames(trait.pca$x),
                             trait.pca$x,
                             leaf_habit=leaf_habit)
trait.pca.loadings<-data.frame(variables = rownames(trait.pca$rotation),
                               trait.pca$rotation)
trait.pca.loadings$variables<-c("LDMC","%N","SLA","SRL","WD")
trait.pca.perc<-trait.pca$sdev^2/sum(trait.pca$sdev^2)*100

leaf_habit_cols<-c("Deciduous"="#a44f30",
                   "Evergreen"="#60941a")

trait.pca.plot<-ggplot(trait.pca.scores, 
                       aes(x = -PC1, y = PC2)) +
  geom_text(size = 3.5,label = trait.pca.scores$species,
            aes(color = trait.pca.scores$leaf_habit)) +
  geom_segment(data = trait.pca.loadings,
               aes(x = 0, y = 0, xend = -PC1*5, yend = PC2*5),
               arrow = arrow(length = unit(1/2, "picas")),
               size=1,color="#3366FF") +
  annotate("text",
           x = -trait.pca.loadings$PC1*5.5,
           y = trait.pca.loadings$PC2*5.5,
           label = trait.pca.loadings$variables,
           color="#3366FF")+
  theme_bw()+theme(text=element_text(size=15),
                   panel.background = element_rect(fill='transparent'), #transparent panel bg
                   plot.background = element_rect(fill='transparent', color=NA))+
  coord_fixed(ratio=trait.pca.perc[2]/trait.pca.perc[1])+
  guides(color="none")+
  scale_color_manual(values = leaf_habit_cols)+
  labs(x=paste("CP1 (",round(trait.pca.perc[1],1),"% variance)",sep=""),
       y=paste("CP2 (",round(trait.pca.perc[2],1),"% variance)",sep=""))

ggsave("Images/trait_pca_plot.png", trait.pca.plot, bg='transparent',
       dpi=600,width=7,height=4)

write.csv(trait.pca.scores,"TraitData/trait_pca_scores.csv")

#########################################
## calculate functional diversity

## calculate CWM of first trait PCA component, omitting the central tree
## this gives as the functional identity of the neighborhood
neighbor.prop<-neighbor.df.centerless/rowSums(neighbor.df.centerless)
neighbor.CWM1<-data.frame(as.matrix(neighbor.prop) %*% trait.pca.scores)$PC1
## for trees with no neighbors, assign a neighborhood functional ID of 0
neighbor.CWM1[is.na(neighbor.CWM1)]<-0

## calculate FDis using FD package
## including the central tree
neighbor.FDis<-fdisp(dist(traits),
                     a = as.matrix(neighbor.df))

## calculate Scheiner's metrics
neighbor.FTD<-FTD.comm(tdmat=dist(traits),
                       spmat = as.matrix(neighbor.df),
                       abund=T)

neighbor.data<-data.frame(UniqueTreeID=names(neighbor.area),
                          neighbor.richness=neighbor.richness,
                          comp.index=neighbor.comp.df$comp.index,
                          Hegyi.index=neighbor.comp.df$Hegyi.index,
                          neighbor.FI1=neighbor.CWM1,
                          FDis=neighbor.FDis$FDis,
                          qDTM=neighbor.FTD$com.FTD$qDTM)

## to do:
## calculate CWMs of major traits

## write data
write.csv(neighbor.data,"IDENTMontrealData/neighborhood_vars.csv")
