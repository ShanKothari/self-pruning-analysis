setwd("C:/Users/querc/Dropbox/PostdocProjects/SelfPruning")

library(piecewiseSEM)
library(tidySEM)
library(lavaan)
library(ggplot2)
library(lme4)
library(lmerTest)

self_pruning<-read.csv("SelfPruningData/self_pruning_processed.csv")

## drop dead trees
self_pruning<-self_pruning[-which(toupper(self_pruning$TreeID)=="DEAD"),]

## z-standardize important variables
self_pruning_standard<-self_pruning
standard_cols<-c("neighbor_comp","FDis","qDTM",
                 "neighbor_richness","HeightTop",
                 "HeightBase","CR_average","shade_tol",
                 "logLightBase","logLightTop",
                 "focal_acq","neighbor_acq")
self_pruning_standard[,standard_cols]<-scale(self_pruning_standard[,standard_cols])

## try some exploratory multivariate analyses I guess
summary(lm(HeightBase~HeightTop+neighbor_comp+qDTM+neighbor_acq+shade_tol+focal_acq,
           data=self_pruning_standard))
summary(lm(logLightBase~HeightTop+neighbor_comp+qDTM+neighbor_acq+shade_tol+focal_acq,
           data=self_pruning_standard))
summary(lm(AliveCrown....total.height.~HeightTop+neighbor_comp+qDTM+neighbor_acq+shade_tol+focal_acq,
           data=self_pruning_standard))
summary(lm(CrownDepth~HeightTop+neighbor_comp+qDTM+neighbor_acq+shade_tol+focal_acq,
           data=self_pruning_standard))

summary(lmer(logLightBase~HeightTop+neighbor_comp+qDTM+neighbor_acq+(1|Species),
           data=self_pruning_standard))

########################################
## local estimation of SEMs

## specify SEMs
## should I include qDTM directly in predicting HeightTop for neighborhood models?

m_base_light_lavaan<-'
neighbor_comp~1
HeightTop~1+shade_tol+focal_acq
CR_average~1+shade_tol+focal_acq
logLightTop~1+HeightTop
logLightBase~1+HeightTop+logLightTop+CR_average+shade_tol+focal_acq
'
fit_base_light <- sem(m_base_light_lavaan,
                      data=self_pruning_standard)

m_base_light_neighbor_lavaan<-'
neighbor_comp~1+qDTM+neighbor_acq
HeightTop~1+neighbor_comp+shade_tol+focal_acq
CR_average~1+neighbor_comp+shade_tol+focal_acq
logLightTop~1+HeightTop+neighbor_comp+qDTM+neighbor_acq
logLightBase~1+HeightTop+logLightTop+CR_average+shade_tol+focal_acq+neighbor_comp+qDTM+neighbor_acq
'
fit_base_light_neighbor <- sem(m_base_light_neighbor_lavaan,
                               data=self_pruning_standard)

m_base_height_lavaan<-'
neighbor_comp~1
HeightTop~1+shade_tol+focal_acq
CR_average~1+shade_tol+focal_acq
logLightTop~1+HeightTop
HeightBase~1+HeightTop+logLightTop+CR_average+shade_tol+focal_acq
'
fit_base_height <- sem(m_base_height_lavaan,
                       data=self_pruning_standard)

m_base_height_neighbor_lavaan<-'
neighbor_comp~1+qDTM+neighbor_acq
HeightTop~1+neighbor_comp+shade_tol+focal_acq
CR_average~1+neighbor_comp+shade_tol+focal_acq
logLightTop~1+HeightTop+neighbor_comp+qDTM+neighbor_acq
HeightBase~1+HeightTop+logLightTop+CR_average+shade_tol+focal_acq+neighbor_comp+qDTM+neighbor_acq
'
fit_base_height_neighbor <- sem(m_base_height_neighbor_lavaan,
                                data=self_pruning_standard)

m_base_height_light_lavaan<-'
neighbor_comp~1
HeightTop~shade_tol+focal_acq
CR_average~shade_tol+focal_acq
logLightTop~HeightTop
logLightBase~HeightTop+logLightTop+CR_average+shade_tol+focal_acq
HeightBase~logLightBase+HeightTop+CR_average+shade_tol+focal_acq
'
fit_base_height_light <- sem(m_base_height_light_lavaan,
                             data=self_pruning_standard)

m_base_height_light_neighbor_lavaan<-'
neighbor_comp~qDTM+neighbor_acq
HeightTop~neighbor_comp+shade_tol+focal_acq
CR_average~neighbor_comp+shade_tol+focal_acq
logLightTop~HeightTop+neighbor_comp+qDTM+neighbor_acq
logLightBase~HeightTop+logLightTop+CR_average+shade_tol+focal_acq+neighbor_comp+qDTM+neighbor_acq
HeightBase~logLightBase+HeightTop+CR_average+shade_tol+focal_acq+neighbor_comp+qDTM+neighbor_acq
'
fit_base_height_light_neighbor <- sem(m_base_height_light_neighbor_lavaan,
                                      data=self_pruning_standard)

########################################
## species-specific fits

## these fits remove the focal ID and shade tolerance
## variables since these are constant within species

m_base_light_lavaan_sp<-'
neighbor_comp~1
HeightTop~1
CR_average~1
logLightTop~1+HeightTop
logLightBase~1+HeightTop+logLightTop+CR_average
'

m_base_light_neighbor_lavaan_sp<-'
neighbor_comp~1+qDTM+neighbor_acq
HeightTop~1+neighbor_comp
CR_average~1+neighbor_comp
logLightTop~1+HeightTop+neighbor_comp+qDTM+neighbor_acq
logLightBase~1+HeightTop+logLightTop+CR_average+neighbor_comp+qDTM+neighbor_acq
'

m_base_height_lavaan_sp<-'
neighbor_comp~1
HeightTop~1
CR_average~1
logLightTop~1+HeightTop
HeightBase~1+HeightTop+logLightTop+CR_average
'

m_base_height_neighbor_lavaan_sp<-'
neighbor_comp~1+qDTM+neighbor_acq
HeightTop~1+neighbor_comp
CR_average~1+neighbor_comp
logLightTop~1+HeightTop+neighbor_comp+qDTM+neighbor_acq
HeightBase~1+HeightTop+logLightTop+CR_average+neighbor_comp+qDTM+neighbor_acq
'

## split up data
self_pruning_standard_sp<-split(self_pruning_standard,
                                f = self_pruning_standard$Species)

fit_base_height_sp_list<-lapply(self_pruning_standard_sp,
                                function(sp_df) sem(m_base_height_lavaan_sp,data=sp_df))
fit_base_height_neighbor_sp_list<-lapply(self_pruning_standard_sp,
                                         function(sp_df) sem(m_base_height_neighbor_lavaan_sp,data=sp_df))
data.frame(no.neighbor=unlist(lapply(fit_base_height_sp_list,AIC)),
           neighbor=unlist(lapply(fit_base_height_neighbor_sp_list,AIC)))

fit_base_light_sp_list<-lapply(self_pruning_standard_sp,
                               function(sp_df) sem(m_base_light_lavaan_sp,data=sp_df))
fit_base_light_neighbor_sp_list<-lapply(self_pruning_standard_sp,
                                        function(sp_df) sem(m_base_light_neighbor_lavaan_sp,data=sp_df))
data.frame(no.neighbor=unlist(lapply(fit_base_light_sp_list,AIC)),
           neighbor=unlist(lapply(fit_base_light_neighbor_sp_list,AIC)))

########################################
## sandbox for simulations to predict on alternate data
## to evaluate whether functional variation predicts shade tolerance

sp_standard_nodiv<-self_pruning_standard
sp_standard_nodiv$qDTM<-0

tmp<-predict(m_base_height_light_neighbor,newdata=sp_standard_nodiv)

#####################################
## local estimation with piecewise SEM

m_base_light<-psem(
  lm(HeightTop~shade_tol+focal_acq,data=self_pruning_standard),
  lm(CR_average~shade_tol+focal_acq,data=self_pruning_standard),
  lm(logLightTop~HeightTop,data=self_pruning_standard),
  lm(logLightBase~HeightTop+logLightTop+CR_average+shade_tol+focal_acq,
     data=self_pruning_standard)
)

m_base_light_neighbor<-psem(
  lm(neighbor_comp~qDTM+neighbor_acq,data=self_pruning_standard),
  lm(HeightTop~neighbor_comp+shade_tol+focal_acq,data=self_pruning_standard),
  lm(CR_average~neighbor_comp+shade_tol+focal_acq,data=self_pruning_standard),
  lm(logLightTop~HeightTop+neighbor_comp+qDTM+neighbor_acq,data=self_pruning_standard),
  lm(logLightBase~HeightTop+logLightTop+CR_average+shade_tol+focal_acq+neighbor_comp+qDTM+neighbor_acq,
     data=self_pruning_standard)
)

m_base_height<-psem(
  lm(HeightTop~shade_tol+focal_acq,data=self_pruning_standard),
  lm(CR_average~shade_tol+focal_acq,data=self_pruning_standard),
  lm(logLightTop~HeightTop,data=self_pruning_standard),
  lm(HeightBase~HeightTop+logLightTop+CR_average+shade_tol+focal_acq,data=self_pruning_standard)
)

m_base_height_neighbor<-psem(
  lm(neighbor_comp~qDTM+neighbor_acq,data=self_pruning_standard),
  lm(HeightTop~neighbor_comp+shade_tol+focal_acq,data=self_pruning_standard),
  lm(CR_average~neighbor_comp+shade_tol+focal_acq,data=self_pruning_standard),
  lm(logLightTop~HeightTop+neighbor_comp+qDTM+neighbor_acq,data=self_pruning_standard),
  lm(HeightBase~HeightTop+logLightTop+CR_average+shade_tol+focal_acq+neighbor_comp+qDTM+neighbor_acq,
     data=self_pruning_standard)
)

m_base_height_light<-psem(
  lm(HeightTop~shade_tol+focal_acq,data=self_pruning_standard),
  lm(CR_average~shade_tol+focal_acq,data=self_pruning_standard),
  lm(logLightTop~HeightTop,data=self_pruning_standard),
  lm(logLightBase~HeightTop+logLightTop+CR_average+shade_tol+focal_acq,
     data=self_pruning_standard),
  lm(HeightBase~logLightBase+HeightTop+CR_average+shade_tol+focal_acq,
     data=self_pruning_standard)
)

m_base_height_light_neighbor<-psem(
  lm(neighbor_comp~qDTM+neighbor_acq,data=self_pruning_standard),
  lm(HeightTop~neighbor_comp+shade_tol+focal_acq,data=self_pruning_standard),
  lm(CR_average~neighbor_comp+shade_tol+focal_acq,data=self_pruning_standard),
  lm(logLightTop~HeightTop+neighbor_comp+qDTM+neighbor_acq,data=self_pruning_standard),
  lm(logLightBase~HeightTop+logLightTop+CR_average+shade_tol+focal_acq+neighbor_comp+qDTM+neighbor_acq,
     data=self_pruning_standard),
  lm(HeightBase~logLightBase+HeightTop+CR_average+shade_tol+focal_acq+neighbor_comp+qDTM+neighbor_acq,
     data=self_pruning_standard)
)
