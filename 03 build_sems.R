library(lavaan)
library(tidySEM)

## create dummy variables for each species
# self_pruning$ACRU<-ifelse(self_pruning$Species=="ACRU",1,0)
# self_pruning$ACSA<-ifelse(self_pruning$Species=="ACSA",1,0)
# self_pruning$BEAL<-ifelse(self_pruning$Species=="BEAL",1,0)
# self_pruning$BEPA<-ifelse(self_pruning$Species=="BEPA",1,0)
# self_pruning$LALA<-ifelse(self_pruning$Species=="LALA",1,0)
# self_pruning$PIGL<-ifelse(self_pruning$Species=="PIGL",1,0)
# self_pruning$PIRE<-ifelse(self_pruning$Species=="PIRE",1,0)
# self_pruning$PIRU<-ifelse(self_pruning$Species=="PIRU",1,0)
# self_pruning$PIST<-ifelse(self_pruning$Species=="PIST",1,0)
# self_pruning$QURU<-ifelse(self_pruning$Species=="QURU",1,0)
# self_pruning$THOC<-ifelse(self_pruning$Species=="THOC",1,0)

self_pruning_standard<-self_pruning
standard_cols<-c("neighbor.comp","FDis","qDTM",
                 "neighbor.richness","HeightTop",
                 "HeightBase","CR_average","ShadeTol",
                 "pseudoLAI_base","pseudoLAI_top",
                 "focalID","neighborID")
self_pruning_standard[,standard_cols]<-scale(self_pruning_standard[,standard_cols])

# m_baseheight_sp<-'
# neighbor.comp~1+FDis+neighbor.richness
# HeightTop~1+neighbor.comp+FDis+neighbor.richness+ACRU+ACSA+BEAL+BEPA+LALA+PIGL+PIRE+PIRU+PIST+QURU+THOC
# HeightBase~1+HeightTop+neighbor.comp+FDis+neighbor.richness+CR_average+ACRU+ACSA+BEAL+BEPA+LALA+PIGL+PIRE+PIRU+PIST+QURU+THOC
# '

# m_baseheight<-'
# neighbor.comp~1
# HeightTop~1+ShadeTol
# pseudoLAI_top~1+HeightTop
# pseudoLAI_base~1+HeightTop+CR_average+ShadeTol+pseudoLAI_top
# HeightBase~1+pseudoLAI_base+HeightTop+CR_average+ShadeTol
# '

m_baseheight<-'
neighbor.comp~1
HeightTop~1+ShadeTol+focalID
pseudoLAI_top~1+HeightTop
pseudoLAI_base~1+HeightTop+CR_average+ShadeTol+focalID+pseudoLAI_top
HeightBase~1+pseudoLAI_base+HeightTop+CR_average+ShadeTol+focalID
'

m_baseheight_shade_tol_orig<-'
neighbor.comp~1+FDis+neighbor.richness
HeightTop~1+neighbor.comp+FDis+neighbor.richness+ShadeTol
pseudoLAI_top~1+HeightTop+neighbor.comp+FDis+neighbor.richness
pseudoLAI_base~1+HeightTop+neighbor.comp+FDis+neighbor.richness+CR_average+ShadeTol+pseudoLAI_top
HeightBase~1+pseudoLAI_base+HeightTop+neighbor.comp+FDis+neighbor.richness+CR_average+ShadeTol
'

m_baseheight_shade_tol<-'
neighbor.comp~1+qDTM+neighborID
HeightTop~1+neighbor.comp+qDTM+ShadeTol+focalID
pseudoLAI_top~1+HeightTop+neighbor.comp+qDTM+neighborID
pseudoLAI_base~1+HeightTop+pseudoLAI_top+neighbor.comp+qDTM+neighborID+CR_average+ShadeTol+focalID
HeightBase~1+pseudoLAI_base+HeightTop+neighbor.comp+qDTM+neighborID+CR_average+ShadeTol+focalID
'

fit_baseheight <- sem(m_baseheight,data=self_pruning_standard)
fit_baseheight_shade_tol <- sem(m_baseheight_shade_tol,data=self_pruning_standard)

graph_sem(fit_baseheight_shade_tol)

summary(lm(HeightBase~HeightTop+neighbor.comp+FDis+neighbor.richness+Species,
           data=self_pruning))
summary(lm(pseudoLAI_base~HeightTop+neighbor.comp+FDis+neighbor.richness+Species,
           data=self_pruning))
