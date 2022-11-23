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

## read in traits to get shade tolerance
traits<-read.csv("TraitData/IDENT_TRAIT_DATABASE_2020-10-20.csv")
self_pruning$ShadeTol<-traits$Shade.tolerance[match(self_pruning$Species,
                                                    traits$SpeciesCode)]

## read in trait PCA to get focal sp functional ID
## focal tree identity of PC1
trait.pca.scores<-read.csv("TraitData/trait_pca_scores.csv")
self_pruning$focalID<-trait.pca.scores$PC1[match(self_pruning$Species,
                                                 trait.pca.scores$X)]

## species means of crown base pseudo-LAI
pseudoLAI_table<-aggregate(self_pruning$pseudoLAI_base,
                           by=list(self_pruning$Species),
                           FUN=mean,na.rm=T)
colnames(pseudoLAI_table)<-c("Species","PseudoLAI")
pseudoLAI_table$ShadeTol<-traits$Shade.tolerance[match(pseudoLAI_table$Species,
                                                       traits$SpeciesCode)]
pseudoLAI_table$FocalID<-trait.pca.scores$PC1[match(pseudoLAI_table$Species,
                                                    trait.pca.scores$X)]

summary(lm(PseudoLAI~ShadeTol,data=pseudoLAI_table))

ggplot(data=self_pruning,aes(x=Species,y=pseudoLAI_base))+
  geom_violin()+geom_point()

ggplot(data=pseudoLAI_table,
       aes(x=ShadeTol,y=PseudoLAI,label=Species))+
  geom_smooth(method="lm")+geom_text()+
  theme_bw()+
  theme(text=element_text(size=15))+
  labs(x="Shade tolerance",
       y="Pseudo-LAI above crown base")

ggplot(data=pseudoLAI_table,
       aes(x=FocalID,y=PseudoLAI,label=Species))+
  geom_smooth(method="lm")+geom_text()+
  theme_bw()+
  theme(text=element_text(size=15))+
  labs(x="Functional identity",
       y="Pseudo-LAI above crown base")
