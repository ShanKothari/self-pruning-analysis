setwd("C:/Users/querc/Dropbox/PostdocProjects/SelfPruning")

library(piecewiseSEM)
library(tidySEM)

## z-standardize important variables
self_pruning_standard<-self_pruning
standard_cols<-c("neighbor.comp","FDis","qDTM",
                 "neighbor.richness","HeightTop",
                 "HeightBase","CR_average","ShadeTol",
                 "pseudoLAI_base","pseudoLAI_top",
                 "focalID","neighborID")
self_pruning_standard[,standard_cols]<-scale(self_pruning_standard[,standard_cols])

## drop rows with NAs
## right now there are four living trees with missing base height
## that are also being dropped
self_pruning_standard<-self_pruning_standard[-unique(which(is.na(self_pruning_standard),arr.ind=T)[,1]),]

## try some exploratory analyses I guess
summary(lm(HeightBase~HeightTop+neighbor.comp+qDTM+neighborID+ShadeTol+focalID,
           data=self_pruning_standard))
summary(lm(pseudoLAI_base~HeightTop+neighbor.comp+qDTM+neighborID+ShadeTol+focalID,
           data=self_pruning_standard))

## specify SEMs
## should I include qDTM directly in predicting HeightTop for neighborhood models?

m_base_light<-psem(
  lm(HeightTop~1+ShadeTol+focalID,data=self_pruning_standard),
  lm(pseudoLAI_top~1+HeightTop,data=self_pruning_standard),
  lm(pseudoLAI_base~1+HeightTop+CR_average+ShadeTol+focalID+pseudoLAI_top,
     data=self_pruning_standard)
)

m_base_light_neighbor<-psem(
  lm(neighbor.comp~1+qDTM+neighborID,data=self_pruning_standard),
  lm(HeightTop~1+neighbor.comp+ShadeTol+focalID,data=self_pruning_standard),
  lm(pseudoLAI_top~1+HeightTop+neighbor.comp+qDTM+neighborID,data=self_pruning_standard),
  lm(pseudoLAI_base~1+HeightTop+pseudoLAI_top+neighbor.comp+qDTM+neighborID+CR_average+ShadeTol+focalID,
     data=self_pruning_standard)
)

m_base_height<-psem(
  lm(HeightTop~1+ShadeTol+focalID,data=self_pruning_standard),
  lm(pseudoLAI_top~1+HeightTop,data=self_pruning_standard),
  lm(HeightBase~1+HeightTop+CR_average+ShadeTol+focalID+pseudoLAI_top,data=self_pruning_standard)
)

m_base_height_neighbor<-psem(
  lm(neighbor.comp~1+qDTM+neighborID,data=self_pruning_standard),
  lm(HeightTop~1+neighbor.comp+ShadeTol+focalID,data=self_pruning_standard),
  lm(pseudoLAI_top~1+HeightTop+neighbor.comp+qDTM+neighborID,data=self_pruning_standard),
  lm(HeightBase~1+HeightTop+neighbor.comp+qDTM+neighborID+CR_average+ShadeTol+focalID+pseudoLAI_top,
     data=self_pruning_standard)
)

m_base_height_light<-psem(
  lm(HeightTop~1+ShadeTol+focalID,data=self_pruning_standard),
  lm(pseudoLAI_top~1+HeightTop,data=self_pruning_standard),
  lm(pseudoLAI_base~1+HeightTop+CR_average+ShadeTol+focalID+pseudoLAI_top,
     data=self_pruning_standard),
  lm(HeightBase~1+pseudoLAI_base+HeightTop+CR_average+ShadeTol+focalID,
     data=self_pruning_standard),
)

m_base_height_light_neighbor<-psem(
  lm(neighbor.comp~1+qDTM+neighborID,data=self_pruning_standard),
  lm(HeightTop~1+neighbor.comp+ShadeTol+focalID,data=self_pruning_standard),
  lm(pseudoLAI_top~1+HeightTop+neighbor.comp+qDTM+neighborID,data=self_pruning_standard),
  lm(pseudoLAI_base~1+HeightTop+pseudoLAI_top+neighbor.comp+qDTM+neighborID+CR_average+ShadeTol+focalID,
     data=self_pruning_standard),
  lm(HeightBase~1+pseudoLAI_base+HeightTop+neighbor.comp+qDTM+neighborID+CR_average+ShadeTol+focalID,
     data=self_pruning_standard)
)

## predict on alternate data
sp_standard_nodiv<-self_pruning_standard
sp_standard_nodiv$qDTM<-0

tmp<-predict(fit_base_height_neighbor,newdata=sp_standard_nodiv)
