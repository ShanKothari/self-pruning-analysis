setwd("C:/Users/querc/Dropbox/PostdocProjects/SelfPruning")

library(ggplot2)

## read and clean data
self_pruning<-read.csv("SelfPruningData/Self_Pruning_DATA.csv")
self_pruning<-self_pruning[-which(self_pruning$Species==""),]
self_pruning$Plot<-gsub(" ","",self_pruning$Plot)
self_pruning$UniqueTreeID<-paste(self_pruning$Block,
                                 self_pruning$Plot,
                                 self_pruning$TreeID,
                                 sep="_")

## produce new self-pruning and predictor variables
## 'pseudo-LAI' style vars have better statistical properties
self_pruning$CD<-self_pruning$HeightTop-self_pruning$HeightBase
self_pruning$pseudoLAI_top<- -log(self_pruning$LightTop../100)
self_pruning$pseudoLAI_base<- -log(self_pruning$LightBase./100)

## attach neighborhood-level variables
neighbor.data<-read.csv("IDENTMontrealData/neighborhood_vars.csv")
self_pruning$neighbor.richness<-neighbor.data$neighbor.richness[match(self_pruning$UniqueTreeID,
                                                                      neighbor.data$UniqueTreeID)]
self_pruning$FDis<-neighbor.data$FDis[match(self_pruning$UniqueTreeID,
                                            neighbor.data$UniqueTreeID)]
self_pruning$qDTM<-neighbor.data$qDTM[match(self_pruning$UniqueTreeID,
                                            neighbor.data$UniqueTreeID)]
self_pruning$neighborID<-neighbor.data$neighbor.FI1[match(self_pruning$UniqueTreeID,
                                                          neighbor.data$UniqueTreeID)]
self_pruning$neighbor.comp<-neighbor.data$comp.index[match(self_pruning$UniqueTreeID,
                                                           neighbor.data$UniqueTreeID)]

## one BEPA has a neighborhood competition index more than twice the others
## so we can eliminate it as a potential outlier
neighbor_outlier<-which(self_pruning$neighbor.comp>60000)
self_pruning<-self_pruning[-neighbor_outlier,]

## read in traits to get shade tolerance
traits<-read.csv("TraitData/IDENT_TRAIT_DATABASE_2020-10-20.csv")
self_pruning$ShadeTol<-traits$Shade.tolerance[match(self_pruning$Species,
                                                    traits$SpeciesCode)]

## read in trait PCA to get focal sp functional ID
## focal tree identity of PC1
trait.pca.scores<-read.csv("TraitData/trait_pca_scores.csv")
self_pruning$focalID<-trait.pca.scores$PC1[match(self_pruning$Species,
                                                 trait.pca.scores$X)]

write.csv(self_pruning,"SelfPruningData/self_pruning_processed.csv",row.names=F)

#######################################
## working with species means

## could look only at trees in monoculture plots
## although I currently do not
self_pruning_mono<-self_pruning[self_pruning$nbsp==1,]

## species means of crown base pseudo-LAI
pseudoLAI_table<-aggregate(self_pruning$pseudoLAI_base,
                           by=list(self_pruning$Species),
                           FUN=mean,na.rm=T)
colnames(pseudoLAI_table)<-c("Species","PseudoLAI")
pseudoLAI_table$ShadeTol<-traits$Shade.tolerance[match(pseudoLAI_table$Species,
                                                       traits$SpeciesCode)]
pseudoLAI_table$focalID<-trait.pca.scores$PC1[match(pseudoLAI_table$Species,
                                                    trait.pca.scores$X)]

summary(lm(PseudoLAI~ShadeTol+focalID,data=pseudoLAI_table))

ggplot(data=pseudoLAI_table,
       aes(x=ShadeTol,y=PseudoLAI,label=Species))+
  geom_smooth(method="lm")+geom_text()+
  theme_bw()+
  theme(text=element_text(size=15))+
  labs(x="Shade tolerance",
       y="Pseudo-LAI above crown base")

ggplot(data=pseudoLAI_table,
       aes(x=focalID,y=PseudoLAI,label=Species))+
  geom_smooth(method="lm")+geom_text()+
  theme_bw()+
  theme(text=element_text(size=15))+
  labs(x="Functional identity",
       y="Pseudo-LAI above crown base")

########################################
## examining simple relationships within the bivariate data

## as neighbor competition increases
## height and crown depth both decrease
## but the living fraction of the crown increases
## and the light at the crown base increases a bit?
ggplot(self_pruning,
       aes(x=neighbor.comp,
           y=pseudoLAI_base,
           color=Species))+
  geom_point()+geom_smooth(method="lm",se=F)

## this plot is kind of odd and visually striking...
ggplot(self_pruning,
       aes(x=HeightTop,
           y=pseudoLAI_base,
           color=Species))+
  geom_point()+geom_smooth(method="lm",se=F)

self_pruning_sp<-split(self_pruning,
                       f = self_pruning$Species)
light_height_slopes<-unlist(lapply(self_pruning_sp,
                                   function(x) {
                                     reg<-lm(pseudoLAI_base~HeightTop,data=x)
                                     return(reg$coefficients[2])
                                     }))
pseudoLAI_table$light_height_slope<-light_height_slopes[match(pseudoLAI_table$Species,
                                                              names(self_pruning_sp))]
ggplot(data=pseudoLAI_table,
       aes(x=focalID,y=light_height_slope,label=Species))+
  geom_smooth(method="lm")+geom_text()+
  theme_bw()+
  theme(text=element_text(size=15))+
  labs(x="Functional identity",
       y="Change in pseudo-LAI at base with top height")


##########
## to do: check that neighbor comp is calculated correctly