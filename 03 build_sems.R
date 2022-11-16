library(lavaan)
library(tidySEM)

## z-standardize important variables
self_pruning_standard<-self_pruning
standard_cols<-c("neighbor.comp","FDis","qDTM",
                 "neighbor.richness","HeightTop",
                 "HeightBase","CR_average","ShadeTol",
                 "pseudoLAI_base","pseudoLAI_top",
                 "focalID","neighborID")
self_pruning_standard[,standard_cols]<-scale(self_pruning_standard[,standard_cols])

## try some exploratory analyses I guess
summary(lm(HeightBase~HeightTop+neighbor.comp+qDTM+neighborID+ShadeTol+focalID,
           data=self_pruning_standard))
summary(lm(pseudoLAI_base~HeightTop+neighbor.comp+qDTM+neighborID+ShadeTol+focalID,
           data=self_pruning_standard))

## specify SEMs
## should I include qDTM directly in predicting HeightTop for neighborhood models?

m_base_light<-'
neighbor.comp~1
HeightTop~1+ShadeTol+focalID
pseudoLAI_top~1+HeightTop
pseudoLAI_base~1+HeightTop+CR_average+ShadeTol+focalID+pseudoLAI_top
'

m_base_light_neighbor<-'
neighbor.comp~1+qDTM+neighborID
HeightTop~1+neighbor.comp+ShadeTol+focalID
pseudoLAI_top~1+HeightTop+neighbor.comp+qDTM+neighborID
pseudoLAI_base~1+HeightTop+pseudoLAI_top+neighbor.comp+qDTM+neighborID+CR_average+ShadeTol+focalID
'

m_base_height<-'
neighbor.comp~1
HeightTop~1+ShadeTol+focalID
pseudoLAI_top~1+HeightTop
HeightBase~1+HeightTop+CR_average+ShadeTol+focalID+pseudoLAI_top
'

m_base_height_neighbor<-'
neighbor.comp~1+qDTM+neighborID
HeightTop~1+neighbor.comp+ShadeTol+focalID
pseudoLAI_top~1+HeightTop+neighbor.comp+qDTM+neighborID
HeightBase~1+HeightTop+neighbor.comp+qDTM+neighborID+CR_average+ShadeTol+focalID+pseudoLAI_top
'

m_base_height_light<-'
neighbor.comp~1
HeightTop~1+ShadeTol+focalID
pseudoLAI_top~1+HeightTop
pseudoLAI_base~1+HeightTop+CR_average+ShadeTol+focalID+pseudoLAI_top
HeightBase~1+pseudoLAI_base+HeightTop+CR_average+ShadeTol+focalID
'

m_base_height_light_neighbor<-'
neighbor.comp~1+qDTM+neighborID
HeightTop~1+neighbor.comp+ShadeTol+focalID
pseudoLAI_top~1+HeightTop+neighbor.comp+qDTM+neighborID
pseudoLAI_base~1+HeightTop+pseudoLAI_top+neighbor.comp+qDTM+neighborID+CR_average+ShadeTol+focalID
HeightBase~1+pseudoLAI_base+HeightTop+neighbor.comp+qDTM+neighborID+CR_average+ShadeTol+focalID
'

## fit SEMs
fit_base_height <- sem(m_base_height,data=self_pruning_standard)
fit_base_height_neighbor <- sem(m_base_height_neighbor,data=self_pruning_standard)

graph_sem(fit_base_height_neighbor)

## predict on alternate data
sp_standard_nodiv<-self_pruning_standard
sp_standard_nodiv$qDTM<-0

tmp<-predict(fit_base_height_neighbor,newdata=sp_standard_nodiv)
