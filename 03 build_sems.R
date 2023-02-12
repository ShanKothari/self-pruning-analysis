setwd("C:/Users/querc/Dropbox/PostdocProjects/SelfPruning")

library(piecewiseSEM)
library(tidySEM)
library(lavaan) ## for testing
library(ggplot2)

self_pruning<-read.csv("SelfPruningData/self_pruning_processed.csv")

## look at simple pairwise relationships

## as neighbor competition increases
## height and crown depth both decrease
## but the living fraction of the crown increases
## and the light at the crown base increases a bit?
ggplot(self_pruning,
       aes(x=neighbor.comp,
           y=pseudoLAI_base,
           color=Species))+
  geom_point()+geom_smooth(method="lm",se=F)

## z-standardize important variables
self_pruning_standard<-self_pruning
standard_cols<-c("neighbor.comp","FDis","qDTM",
                 "neighbor.richness","HeightTop",
                 "HeightBase","CR_average","ShadeTol",
                 "pseudoLAI_base","pseudoLAI_top",
                 "focalID","neighborID")
self_pruning_standard[,standard_cols]<-scale(self_pruning_standard[,standard_cols])

## drop dead trees
self_pruning_standard<-self_pruning_standard[-which(toupper(self_pruning_standard$TreeID)=="DEAD"),]

## try some exploratory multivariate analyses I guess
summary(lm(HeightBase~HeightTop+neighbor.comp+qDTM+neighborID+ShadeTol+focalID,
           data=self_pruning_standard))
summary(lm(pseudoLAI_base~HeightTop+neighbor.comp+qDTM+neighborID+ShadeTol+focalID,
           data=self_pruning_standard))
summary(lm(AliveCrown....total.height.~HeightTop+neighbor.comp+qDTM+neighborID+ShadeTol+focalID,
           data=self_pruning_standard))

summary(lm(pseudoLAI_base~HeightTop+neighbor.comp+qDTM+neighborID+ShadeTol+focalID,
           data=self_pruning_standard))

## specify SEMs
## should I include qDTM directly in predicting HeightTop for neighborhood models?

m_base_light_lavaan<-'
neighbor.comp~1
HeightTop~1+ShadeTol+focalID
CR_average~1+ShadeTol+focalID
pseudoLAI_top~1+HeightTop
pseudoLAI_base~1+HeightTop+pseudoLAI_top+CR_average+ShadeTol+focalID
'
fit_base_light <- sem(m_base_light_lavaan,
                      data=self_pruning_standard)

m_base_light_neighbor_lavaan<-'
neighbor.comp~1+qDTM+neighborID
HeightTop~1+neighbor.comp+ShadeTol+focalID
CR_average~1+neighbor.comp+ShadeTol+focalID
pseudoLAI_top~1+HeightTop+neighbor.comp+qDTM+neighborID
pseudoLAI_base~1+HeightTop+pseudoLAI_top+CR_average+ShadeTol+focalID+neighbor.comp+qDTM+neighborID
'
fit_base_light_neighbor <- sem(m_base_light_neighbor_lavaan,
                               data=self_pruning_standard)

m_base_height_lavaan<-'
neighbor.comp~1
HeightTop~1+ShadeTol+focalID
CR_average~1+ShadeTol+focalID
pseudoLAI_top~1+HeightTop
HeightBase~1+HeightTop+pseudoLAI_top+CR_average+ShadeTol+focalID
'
fit_base_height <- sem(m_base_height_lavaan,
                       data=self_pruning_standard)

m_base_height_neighbor_lavaan<-'
neighbor.comp~1+qDTM+neighborID
HeightTop~1+neighbor.comp+ShadeTol+focalID
CR_average~1+neighbor.comp+ShadeTol+focalID
pseudoLAI_top~1+HeightTop+neighbor.comp+qDTM+neighborID
HeightBase~1+HeightTop+pseudoLAI_top+CR_average+ShadeTol+focalID+neighbor.comp+qDTM+neighborID
'
fit_base_height_neighbor <- sem(m_base_height_neighbor_lavaan,
                                data=self_pruning_standard)

m_base_height_light_lavaan<-'
neighbor.comp~1
HeightTop~ShadeTol+focalID
CR_average~ShadeTol+focalID
pseudoLAI_top~HeightTop
pseudoLAI_base~HeightTop+pseudoLAI_top+CR_average+ShadeTol+focalID
HeightBase~pseudoLAI_base+HeightTop+CR_average+ShadeTol+focalID
'
fit_base_height_light <- sem(m_base_height_light_lavaan,
                             data=self_pruning_standard)

m_base_height_light_neighbor_lavaan<-'
neighbor.comp~qDTM+neighborID
HeightTop~neighbor.comp+ShadeTol+focalID
CR_average~neighbor.comp+ShadeTol+focalID
pseudoLAI_top~HeightTop+neighbor.comp+qDTM+neighborID
pseudoLAI_base~HeightTop+pseudoLAI_top+CR_average+ShadeTol+focalID+neighbor.comp+qDTM+neighborID
HeightBase~pseudoLAI_base+HeightTop+CR_average+ShadeTol+focalID+neighbor.comp+qDTM+neighborID
'
fit_base_height_light_neighbor <- sem(m_base_height_light_neighbor_lavaan,
                                      data=self_pruning_standard)

########################################
## species-specific fits

## these fits remove the focal ID and shade tolerance
## variables since these are constant within species

m_base_light_lavaan_sp<-'
neighbor.comp~1
HeightTop~1
CR_average~1
pseudoLAI_top~1+HeightTop
pseudoLAI_base~1+HeightTop+pseudoLAI_top+CR_average
'

m_base_light_neighbor_lavaan_sp<-'
neighbor.comp~1+qDTM+neighborID
HeightTop~1+neighbor.comp
CR_average~1+neighbor.comp
pseudoLAI_top~1+HeightTop+neighbor.comp+qDTM+neighborID
pseudoLAI_base~1+HeightTop+pseudoLAI_top+CR_average+neighbor.comp+qDTM+neighborID
'

m_base_height_lavaan_sp<-'
neighbor.comp~1
HeightTop~1
CR_average~1
pseudoLAI_top~1+HeightTop
HeightBase~1+HeightTop+pseudoLAI_top+CR_average
'

m_base_height_neighbor_lavaan_sp<-'
neighbor.comp~1+qDTM+neighborID
HeightTop~1+neighbor.comp
CR_average~1+neighbor.comp
pseudoLAI_top~1+HeightTop+neighbor.comp+qDTM+neighborID
HeightBase~1+HeightTop+pseudoLAI_top+CR_average+neighbor.comp+qDTM+neighborID
'

## split up data
self_pruning_standard_sp<-split(self_pruning_standard,
                                f = self_pruning_standard$Species)

fit_base_height_sp_list<-lapply(self_pruning_standard_sp,
                                function(sp_df) sem(m_base_height_lavaan_sp,data=sp_df))
fit_base_height_neighbor_sp_list<-lapply(self_pruning_standard_sp,
                                         function(sp_df) sem(m_base_height_neighbor_lavaan_sp,data=sp_df))


########################################
## sandbox for simulations to predict on alternate data
## to evaluate whether functional variation predicts shade tolerance

sp_standard_nodiv<-self_pruning_standard
sp_standard_nodiv$qDTM<-0

tmp<-predict(m_base_height_light_neighbor,newdata=sp_standard_nodiv)

#####################################
## local estimation with piecewise SEM

m_base_light<-psem(
  lm(HeightTop~ShadeTol+focalID,data=self_pruning_standard),
  lm(CR_average~ShadeTol+focalID,data=self_pruning_standard),
  lm(pseudoLAI_top~HeightTop,data=self_pruning_standard),
  lm(pseudoLAI_base~HeightTop+pseudoLAI_top+CR_average+ShadeTol+focalID,
     data=self_pruning_standard)
)

m_base_light_neighbor<-psem(
  lm(neighbor.comp~qDTM+neighborID,data=self_pruning_standard),
  lm(HeightTop~neighbor.comp+ShadeTol+focalID,data=self_pruning_standard),
  lm(CR_average~neighbor.comp+ShadeTol+focalID,data=self_pruning_standard),
  lm(pseudoLAI_top~HeightTop+neighbor.comp+qDTM+neighborID,data=self_pruning_standard),
  lm(pseudoLAI_base~HeightTop+pseudoLAI_top+CR_average+ShadeTol+focalID+neighbor.comp+qDTM+neighborID,
     data=self_pruning_standard)
)

m_base_height<-psem(
  lm(HeightTop~ShadeTol+focalID,data=self_pruning_standard),
  lm(CR_average~ShadeTol+focalID,data=self_pruning_standard),
  lm(pseudoLAI_top~HeightTop,data=self_pruning_standard),
  lm(HeightBase~HeightTop+pseudoLAI_top+CR_average+ShadeTol+focalID,data=self_pruning_standard)
)

m_base_height_neighbor<-psem(
  lm(neighbor.comp~qDTM+neighborID,data=self_pruning_standard),
  lm(HeightTop~neighbor.comp+ShadeTol+focalID,data=self_pruning_standard),
  lm(CR_average~neighbor.comp+ShadeTol+focalID,data=self_pruning_standard),
  lm(pseudoLAI_top~HeightTop+neighbor.comp+qDTM+neighborID,data=self_pruning_standard),
  lm(HeightBase~HeightTop+pseudoLAI_top+CR_average+ShadeTol+focalID+neighbor.comp+qDTM+neighborID,
     data=self_pruning_standard)
)

m_base_height_light<-psem(
  lm(HeightTop~ShadeTol+focalID,data=self_pruning_standard),
  lm(CR_average~ShadeTol+focalID,data=self_pruning_standard),
  lm(pseudoLAI_top~HeightTop,data=self_pruning_standard),
  lm(pseudoLAI_base~HeightTop+pseudoLAI_top+CR_average+ShadeTol+focalID,
     data=self_pruning_standard),
  lm(HeightBase~pseudoLAI_base+HeightTop+CR_average+ShadeTol+focalID,
     data=self_pruning_standard)
)

m_base_height_light_neighbor<-psem(
  lm(neighbor.comp~qDTM+neighborID,data=self_pruning_standard),
  lm(HeightTop~neighbor.comp+ShadeTol+focalID,data=self_pruning_standard),
  lm(CR_average~neighbor.comp+ShadeTol+focalID,data=self_pruning_standard),
  lm(pseudoLAI_top~HeightTop+neighbor.comp+qDTM+neighborID,data=self_pruning_standard),
  lm(pseudoLAI_base~HeightTop+pseudoLAI_top+CR_average+ShadeTol+focalID+neighbor.comp+qDTM+neighborID,
     data=self_pruning_standard),
  lm(HeightBase~pseudoLAI_base+HeightTop+CR_average+ShadeTol+focalID+neighbor.comp+qDTM+neighborID,
     data=self_pruning_standard)
)
