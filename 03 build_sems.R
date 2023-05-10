setwd("C:/Users/querc/Dropbox/PostdocProjects/SelfPruning")

library(piecewiseSEM)
library(tidySEM)
library(lavaan)
library(ggplot2)
library(lme4)
library(lmerTest)

####################################
## to dos

## try adding functional distance of neighbors from focal individual to SEMs
## consider non-abundance weighted metrics?

## maybe use crown depth as the main "position" variable in SEMs?

###################################

self_pruning<-read.csv("SelfPruningData/self_pruning_processed.csv")

## drop dead trees
self_pruning<-self_pruning[-which(toupper(self_pruning$TreeID)=="DEAD"),]

## z-standardize important variables
self_pruning_standard<-self_pruning
standard_cols<-c("neighbor_comp","FDis","qDTM",
                 "neighbor_richness","HeightTop",
                 "HeightBase","CrownDepth","CR_average",
                 "logLightBase","logLightTop",
                 "shade_tol","focal_acq","neighbor_acq",
                 "acq_dist","acq_dist_abs")
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

## decide between ML and REML
summary(lmer(logLightBase~HeightTop+neighbor_comp+qDTM+neighbor_acq+(1|Species),
           data=self_pruning_standard,REML=F))

########################################
## local estimation of SEMs

## specify SEMs
## should I include qDTM directly in predicting HeightTop for neighborhood models?

m_base_light_lavaan<-'
neighbor_comp~1
HeightTop~1+shade_tol+focal_acq
logLightTop~1+HeightTop
logLightBase~1+HeightTop+logLightTop+shade_tol+focal_acq
'
fit_base_light <- sem(m_base_light_lavaan,
                      data=self_pruning_standard)

m_base_light_neighbor_lavaan<-'
neighbor_comp~1+qDTM+neighbor_acq
HeightTop~1+neighbor_comp+shade_tol+focal_acq
logLightTop~1+HeightTop+neighbor_comp+qDTM+neighbor_acq
logLightBase~1+HeightTop+logLightTop+shade_tol+focal_acq+neighbor_comp+qDTM+neighbor_acq
'
fit_base_light_neighbor <- sem(m_base_light_neighbor_lavaan,
                               data=self_pruning_standard)

m_base_height_lavaan<-'
neighbor_comp~1
HeightTop~1+shade_tol+focal_acq
logLightTop~1+HeightTop
HeightBase~1+HeightTop+logLightTop+shade_tol+focal_acq
'
fit_base_height <- sem(m_base_height_lavaan,
                       data=self_pruning_standard)

m_base_height_neighbor_lavaan<-'
neighbor_comp~1+qDTM+neighbor_acq
HeightTop~1+neighbor_comp+shade_tol+focal_acq
logLightTop~1+HeightTop+neighbor_comp+qDTM+neighbor_acq
HeightBase~1+HeightTop+logLightTop+shade_tol+focal_acq+neighbor_comp+qDTM+neighbor_acq
'
fit_base_height_neighbor <- sem(m_base_height_neighbor_lavaan,
                                data=self_pruning_standard)

########################################
## species-specific fits

## these fits remove the focal ID and shade tolerance
## variables since these are constant within species

m_base_light_lavaan_sp<-'
neighbor_comp~1
HeightTop~1
logLightTop~1+HeightTop
logLightBase~1+HeightTop+logLightTop
'

m_base_light_neighbor_lavaan_sp<-'
neighbor_comp~1+qDTM+neighbor_acq
HeightTop~1+neighbor_comp
logLightTop~1+HeightTop+neighbor_comp+qDTM+neighbor_acq
logLightBase~1+HeightTop+logLightTop+neighbor_comp+qDTM+neighbor_acq
'

m_base_height_lavaan_sp<-'
neighbor_comp~1
HeightTop~1
logLightTop~1+HeightTop
HeightBase~1+HeightTop+logLightTop
'

m_base_height_neighbor_lavaan_sp<-'
neighbor_comp~1+qDTM+neighbor_acq
HeightTop~1+neighbor_comp
logLightTop~1+HeightTop+neighbor_comp+qDTM+neighbor_acq
HeightBase~1+HeightTop+logLightTop+neighbor_comp+qDTM+neighbor_acq
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

## using multigroup analysis instead
fit_base_light_multi<-sem(m_base_light_lavaan_sp,
                          data=self_pruning_standard,
                          group="Species")

fit_base_light_neighbor_multi<-sem(m_base_light_neighbor_lavaan_sp,
                          data=self_pruning_standard,
                          group="Species")

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
  lm(logLightTop~HeightTop,data=self_pruning_standard),
  lm(logLightBase~HeightTop+logLightTop+shade_tol+focal_acq,
     data=self_pruning_standard)
)

m_base_light_neighbor<-psem(
  lmer(neighbor_comp~qDTM+neighbor_acq+(1|UniquePlot),data=self_pruning_standard),
  lmer(HeightTop~neighbor_comp+shade_tol+focal_acq+(1|UniquePlot),data=self_pruning_standard),
  lmer(logLightTop~HeightTop+neighbor_comp+qDTM+neighbor_acq+(1|UniquePlot),data=self_pruning_standard),
  lmer(logLightBase~HeightTop+logLightTop+shade_tol+focal_acq+neighbor_comp+qDTM+neighbor_acq+(1|UniquePlot),
     data=self_pruning_standard)
)

m_base_height<-psem(
  lm(HeightTop~shade_tol+focal_acq,data=self_pruning_standard),
  lm(logLightTop~HeightTop,data=self_pruning_standard),
  lm(HeightBase~HeightTop+logLightTop+shade_tol+focal_acq,data=self_pruning_standard)
)

m_base_height_neighbor<-psem(
  lm(neighbor_comp~qDTM+neighbor_acq,data=self_pruning_standard),
  lm(HeightTop~neighbor_comp+shade_tol+focal_acq,data=self_pruning_standard),
  lm(logLightTop~HeightTop+neighbor_comp+qDTM+neighbor_acq,data=self_pruning_standard),
  lm(HeightBase~HeightTop+logLightTop+shade_tol+focal_acq+neighbor_comp+qDTM+neighbor_acq,
     data=self_pruning_standard)
)

## try multigroup SEM?
m_base_light_neighbor_multi<-psem(
  lm(neighbor_comp~qDTM+neighbor_acq,data=self_pruning_standard),
  lm(HeightTop~neighbor_comp,data=self_pruning_standard),
  lm(logLightTop~HeightTop+neighbor_comp+qDTM+neighbor_acq,data=self_pruning_standard),
  lm(logLightBase~HeightTop+logLightTop+neighbor_comp+qDTM+neighbor_acq,
       data=self_pruning_standard)
)

summary(m_base_light_neighbor_multi)
tmp<-multigroup(m_base_light_neighbor_multi,group = "Species")
