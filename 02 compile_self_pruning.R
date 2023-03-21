setwd("C:/Users/querc/Dropbox/PostdocProjects/SelfPruning")

library(ggplot2)

## to do:
## also calculate non-abundance-weighted values
## of neighbor function/diversity?

## read and clean data
self_pruning<-read.csv("SelfPruningData/Self_Pruning_DATA.csv")
self_pruning<-self_pruning[-which(self_pruning$Species==""),]
self_pruning$Plot<-gsub(" ","",self_pruning$Plot)
self_pruning$UniqueTreeID<-paste(self_pruning$Block,
                                 self_pruning$Plot,
                                 self_pruning$TreeID,
                                 sep="_")

####################################
## produce new self-pruning and predictor variables
self_pruning$CrownDepth<-self_pruning$HeightTop-self_pruning$HeightBase

## basal diameter is missing from one tree, which presumably died
## between Jon's measurements and the fall survey
self_pruning$BasalDiam<-as.numeric(self_pruning$BasalDiam)
self_pruning$BasalArea<-(self_pruning$BasalDiam/2)^2*pi

## log-transforming light shows better statistical properties
## and is justifiable via the Beer-Lambert idealization of
## canopy light transmission
self_pruning$logLightTop<-log(self_pruning$LightTop../100)
self_pruning$logLightBase<-log(self_pruning$LightBase./100)

## attach neighborhood-level variables
neighbor.data<-read.csv("IDENTMontrealData/neighborhood_vars.csv")
self_pruning$neighbor_richness<-neighbor.data$neighbor.richness[match(self_pruning$UniqueTreeID,
                                                                      neighbor.data$UniqueTreeID)]
self_pruning$FDis<-neighbor.data$FDis[match(self_pruning$UniqueTreeID,
                                            neighbor.data$UniqueTreeID)]
self_pruning$qDTM<-neighbor.data$qDTM[match(self_pruning$UniqueTreeID,
                                            neighbor.data$UniqueTreeID)]
self_pruning$neighbor_acq<- -neighbor.data$neighbor.FI1[match(self_pruning$UniqueTreeID,
                                                          neighbor.data$UniqueTreeID)]
self_pruning$neighbor_comp<-neighbor.data$comp.index[match(self_pruning$UniqueTreeID,
                                                           neighbor.data$UniqueTreeID)]

## one BEPA has a neighborhood competition index more than twice
## the others so we can eliminate it as a potential outlier
## NOTE: this code will probably need changing based on
## recalculation of NCI
neighbor_outlier<-which(self_pruning$neighbor_comp>60000)
self_pruning<-self_pruning[-neighbor_outlier,]

## read in traits to get shade tolerance
traits<-read.csv("TraitData/IDENT_TRAIT_DATABASE_2020-10-20.csv")
self_pruning$shade_tol<-traits$Shade.tolerance[match(self_pruning$Species,
                                                    traits$SpeciesCode)]

## read in trait PCA to get focal sp functional ID
## focal tree identity of PC1
trait.pca.scores<-read.csv("TraitData/trait_pca_scores.csv")
self_pruning$focal_acq<- -trait.pca.scores$PC1[match(self_pruning$Species,
                                                     trait.pca.scores$X)]
self_pruning$focal_fundist<- self_pruning$focal_acq-self_pruning$neighbor_acq
self_pruning$focal_fundist_abs<- abs(self_pruning$focal_fundist)

# write.csv(self_pruning,"SelfPruningData/self_pruning_processed.csv",row.names=F)

#######################################
## working with species means

## could look only at trees in monoculture plots
## although I currently do not
self_pruning_mono<-self_pruning[self_pruning$nbsp==1,]

## species means of crown base light
self_pruning_sub<-self_pruning[,c("logLightBase","Species","HeightTop","BasalArea")]
species_means<-aggregate(.~Species,data = self_pruning_sub,
                         FUN = mean,na.rm = T)
species_means$shade_tol<-traits$Shade.tolerance[match(species_means$Species,
                                                     traits$SpeciesCode)]
species_means$focal_acq<- -trait.pca.scores$PC1[match(species_means$Species,
                                                  trait.pca.scores$X)]

summary(lm(logLightBase~shade_tol+focal_acq,data=species_means))

# png("Images/lfbase_st_sp.png",width=5,height=5,units="in",res=150)
# ggplot(data=species_means,
#        aes(x=shade_tol,y=logLightBase,label=Species))+
#   geom_smooth(method="lm")+geom_text()+
#   theme_bw()+
#   theme(text=element_text(size=15))+
#   labs(x="Shade tolerance",
#        y="log(light fraction) at crown base")
# dev.off()

# png("Images/lfbase_acq_sp.png",width=5,height=5,units="in",res=150)
# ggplot(data=species_means,
#        aes(x=focal_acq,y=logLightBase,label=Species))+
#   geom_smooth(method="lm")+geom_text()+
#   theme_bw()+
#   theme(text=element_text(size=15))+
#   labs(x="Focal tree acquisitiveness",
#        y="log(light fraction) at crown base")
# dev.off()

########################################
## examining simple relationships within the bivariate data

## as neighbor competition increases
## height and crown depth both decrease
## but the living fraction of the crown increases
## and the light at the crown base increases a bit?

ggplot(self_pruning,
       aes(x=FDis,
           y=logLightBase,
           color=Species))+
  geom_point()+geom_smooth(method="lm",se=F)

ggplot(self_pruning,
       aes(x=neighbor_comp,
           y=logLightBase,
           color=Species))+
  geom_point()+geom_smooth(method="lm",se=F)

png("Images/lfbase_nacq_ind.png",width=7,height=5,units="in",res=150)
ggplot(self_pruning,
       aes(x=neighbor_acq,
           y=logLightBase,
           color=Species))+
  geom_point()+geom_smooth(method="lm",se=F)+
  theme_bw()+theme(text=element_text(size=20))+
  labs(x="Neighbor acquisitiveness",
       y="log(light fraction) at crown base")
dev.off()

ggplot(self_pruning,
       aes(x=focal_fundist_abs,
           y=logLightBase,
           color=Species))+
  geom_point()+geom_smooth(method="lm",se=F)+
  theme_bw()+theme(text=element_text(size=20))+
  labs(x="Functional distance from neighbors",
       y="log(light fraction) at crown base")

## test of correlative inhibition:
## we should expect a positive slope here
ggplot(self_pruning,
       aes(x=logLightTop,
           y=logLightBase,
           color=Species))+
  geom_point()+geom_smooth(method="lm",se=F)+
  theme_bw()+theme(text=element_text(size=20))+
  labs(x="log(light fraction) at crown top",
       y="log(light fraction) at crown base")

## this plot is kind of odd and visually striking...
ggplot(self_pruning,
       aes(x=HeightTop,
           y=logLightBase,
           color=Species))+
  geom_point()+geom_smooth(method="lm",se=F)

## pull out the species-specific slopes
self_pruning_sp<-split(self_pruning,
                       f = self_pruning$Species)

light_height_slopes<-unlist(lapply(self_pruning_sp,
                                   function(x) {
                                     reg<-lmer(logLightBase~HeightTop+(1|Plot),data=x)
                                     return(fixef(reg)[2])
                                   }))
species_means$light_height_slope<-light_height_slopes[match(species_means$Species,
                                                            names(self_pruning_sp))]

light_neighbor_acq_slopes<-unlist(lapply(self_pruning_sp,
                                         function(x) {
                                           reg<-lmer(logLightBase~neighbor_acq+(1|Plot),data=x)
                                           return(fixef(reg)[2])
                                         }))
species_means$light_neighbor_acq_slope<-light_neighbor_acq_slopes[match(species_means$Species,
                                                                        names(self_pruning_sp))]

## which species may show some evidence of correlative inhibition?
light_toplight_slopes<-unlist(lapply(self_pruning_sp,
                                     function(x) {
                                       reg<-lmer(logLightBase~logLightTop+(1|Plot),data=x)
                                       return(fixef(reg)[2])
                                     }))
species_means$light_toplight_slope<-light_toplight_slopes[match(species_means$Species,
                                                                names(self_pruning_sp))]

ggplot(data=species_means,
       aes(x=focal_acq,y=light_height_slope,label=Species))+
  geom_smooth(method="lm")+geom_text()+
  theme_bw()+
  theme(text=element_text(size=15))+
  labs(x="Functional identity",
       y="Change in log(LF) at base with top height")

ggplot(data=species_means,
       aes(x=focal_acq,y=light_neighbor_acq_slope,label=Species))+
  geom_smooth(method="lm")+geom_text()+
  theme_bw()+
  theme(text=element_text(size=15))+
  labs(x="Functional identity",
       y="Change in log(LF) at base with neighbor acquisitiveness")

ggplot(data=species_means,
       aes(x=focal_acq,y=light_toplight_slope,label=Species))+
  geom_smooth(method="lm")+geom_text()+
  theme_bw()+
  theme(text=element_text(size=15))+
  labs(x="Functional identity",
       y="Change in log(LF) at base with log(LF) at top")

##########
## to do: check that neighbor comp is calculated correctly