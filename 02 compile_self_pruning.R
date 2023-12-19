setwd("C:/Users/querc/Dropbox/PostdocProjects/SelfPruning")

library(ggplot2)
library(ggpubr)
library(lme4)
library(chron)
library(fishmethods)

## manual changes to data files sent by Jon:
## corrected "DEAd" to "DEAD" for A	4N8	4	THOC
## fixed typo in base measurement time for D	2N7	2	LALA
## added correct values for ABBA block D 2NR7A

#######################################
## read and clean data

self_pruning<-read.csv("SelfPruningData/Self_Pruning_DATA_TimeINFO.csv")
self_pruning$Plot<-gsub(" ","",self_pruning$Plot)
self_pruning$unique_plot<-paste(self_pruning$Block,self_pruning$Plot,sep="_")
self_pruning$unique_tree<-paste(self_pruning$unique_plot,
                                self_pruning$TreeID,
                                sep="_")

## drop dead trees
## mortalities in the self-pruning dataset are uninformative
## about true mortality rates, so they convey no real information
self_pruning<-self_pruning[-which(toupper(self_pruning$TreeID)=="DEAD"),]

####################################
## calculate solar zenith angle from time of day

## IDENT-Montreal coordinates
latitude<-45.425
longitude<- -73.939

## time sampled at the base
base_sampling<-strsplit(self_pruning$Time.Sampled,split = "-")
base_start<-unlist(lapply(base_sampling,function(x) ifelse(length(x) > 1, x[1], NA)))
base_start<-sapply(base_start,function(x) ifelse(!is.na(x),paste(x,":00",sep=""),NA))

base_end<-unlist(lapply(base_sampling,function(x) ifelse(length(x) > 1, x[2], NA)))
base_end<-sapply(base_end,function(x) ifelse(!is.na(x),paste(x,":00",sep=""),NA))

base_times<-data.frame(base_start=chron(times=base_start),
                       base_end=chron(times=base_end),
                       date_sampled=self_pruning$Day.Sampled)
date_split<-strsplit(base_times$date_sampled,split="/")
base_times$day<-unlist(lapply(date_split,function(x) as.numeric(ifelse(length(x) > 0, x[1], NA))))
base_times$month<-unlist(lapply(date_split,function(x) as.numeric(ifelse(length(x) > 0, x[2], NA))))
base_times$year<-unlist(lapply(date_split,function(x) as.numeric(ifelse(length(x) > 0, x[3], NA))))

base_times$mean_decimal<-rowMeans(base_times[,c("base_start","base_end")])
base_times$mean_time<-times(base_times$mean_decimal)
base_times$mean_hours<-base_times$mean_decimal*24

base_times$date_time<-paste(base_times$date_sampled,base_times$mean_time,sep=" ")
base_times$date_time_POSIX<-as.POSIXct(base_times$date_time,format="%d/%m/%Y %H:%M:%S")

base_times_sub<-base_times[which(!is.na(base_times$day)),]
## Montreal is always GMT -4 during the months of July and August
base_times_angle<-astrocalc4r(day=base_times_sub$day,
                              month=base_times_sub$month,
                              year=base_times_sub$year,
                              hour=base_times_sub$mean_hours,
                              timezone = rep(-4,times=nrow(base_times_sub)),
                              lat = rep(latitude,times=nrow(base_times_sub)),
                              lon = rep(longitude,times=nrow(base_times_sub)),
                              seaorland = "continental")

## clunky!
self_pruning$zenith<-NULL
self_pruning$zenith[which(!is.na(base_times$day))]<-base_times_angle$zenith

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
neighbor.match<-match(self_pruning$unique_tree,neighbor.data$unique_tree)
self_pruning$neighbor_richness<-neighbor.data$neighbor.richness[neighbor.match]
self_pruning$FDis<-neighbor.data$FDis[neighbor.match]
self_pruning$qDTM<-neighbor.data$qDTM[neighbor.match]
self_pruning$neighbor_acq<- -neighbor.data$neighbor.FI1[neighbor.match]
self_pruning$neighbor_comp<-neighbor.data$NCI[neighbor.match]

## read in trait data to get shade tolerance
trait_summary<-read.csv("TraitData/trait_summary.csv")
trait.match<-match(self_pruning$Species,trait_summary$SpeciesCode)
self_pruning$shade_tol<-trait_summary$shade_tol[trait.match]

## and get focal species functional ID (PC1)
self_pruning$focal_acq<- -trait_summary$PC1[trait.match]
self_pruning$acq_dist<- self_pruning$focal_acq-self_pruning$neighbor_acq
self_pruning$acq_dist_abs<- abs(self_pruning$acq_dist)

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
species_means$shade_tol<-trait_summary$shade_tol[match(species_means$Species,
                                                       trait_summary$SpeciesCode)]
species_means$focal_acq<- -trait_summary$PC1[match(species_means$Species,
                                                   trait_summary$SpeciesCode)]
species_means$PC2<- -trait_summary$PC2[match(species_means$Species,
                                             trait_summary$SpeciesCode)]
species_means$leaf_lifespan<- trait_summary$LL[match(species_means$Species,
                                                     trait_summary$SpeciesCode)]
species_means$leaf_habit<- trait_summary$leaf_habit[match(species_means$Species,
                                                          trait_summary$SpeciesCode)]

# summary(lm(logLightBase~shade_tol+focal_acq,data=species_means))

leaf_habit_cols<-c("Deciduous"="chocolate4",
                   "Evergreen"="#60941a")

st_acq_sp<-ggplot(data=species_means,
                     aes(x=shade_tol,
                         y=focal_acq,
                         label=Species))+
  geom_smooth(method="lm")+
  geom_text(size=5,aes(color=leaf_habit))+
  theme_bw()+
  theme(text=element_text(size=20))+
  guides(color="none")+
  scale_color_manual(values = leaf_habit_cols)+
  labs(x="Shade tolerance",
       y="Acquisitiveness")
# ggsave(filename = "Images/st_acq_sp.png",plot = st_acq_sp,
#        width=6,height=5,units="in",dpi=600)

lfbase_ll_sp<-ggplot(data=species_means,
                     aes(x=log(leaf_lifespan),
                         y=logLightBase,
                         label=Species))+
  geom_smooth(method="lm")+
  geom_text(size=5,aes(color=leaf_habit))+
  theme_bw()+
  theme(text=element_text(size=20))+
  coord_cartesian(ylim=c(-6.5,-2.5))+
  guides(color="none")+
  scale_color_manual(values = leaf_habit_cols)+
  labs(x="log(leaf lifespan [months])",
       y=expression(paste("log (",italic(L[base]),")")))
# ggsave(filename = "Images/lfbase_ll_sp.png",plot = lfbase_ll_sp,
#        width=6,height=5,units="in",dpi=600)

lfbase_st_sp<-ggplot(data=species_means,
                     aes(x=shade_tol,
                         y=logLightBase,
                         label=Species))+
  geom_smooth(method="lm")+
  geom_text(size=5,aes(color=leaf_habit))+
  theme_bw()+
  theme(text=element_text(size=20))+
  coord_cartesian(ylim=c(-6.5,-2.5))+
  guides(color="none")+
  scale_color_manual(values = leaf_habit_cols)+
  # labs(x="Tol\u00e9rance \u00e0 l'ombre",
  #      y=expression(paste("log (",italic(L[base]),")")))
  labs(x="Shade tolerance",
       y=expression(paste("log (",italic(L[base]),")")))
# ggsave(filename = "Images/lfbase_st_sp_FR.png",plot = lfbase_st_sp,
#        width=6,height=5,units="in",dpi=600)

lfbase_acq_sp<-ggplot(data=species_means,
                      aes(x=focal_acq,
                          y=logLightBase,
                          label=Species))+
  geom_smooth(method="lm")+
  geom_text(size=5,aes(color=leaf_habit))+
  theme_bw()+
  theme(text=element_text(size=20))+
  coord_cartesian(ylim=c(-6.5,-2.5))+
  guides(color="none")+
  scale_color_manual(values = leaf_habit_cols)+
  # labs(x="Tendance acquisitive (CP1)",
  #      y=expression(paste("log (",italic(L[base]),")")))
  labs(x="Focal tree acquisitiveness",
       y=expression(paste("log (",italic(L[base]),")")))
# ggsave(filename = "Images/lfbase_acq_sp_FR.png", plot = lfbase_acq_sp,
#        width=6,height=5,units="in",dpi=600)

## PC2 from the trait data (mostly LDMC)
## is actually correlated with shade tolerance!
## but PC2 is a poor predictor of light at the crown base
# ggplot(data=species_means,
#        aes(x=shade_tol,y=PC2,label=Species))+
#   geom_smooth(method="lm")+geom_text()+
#   theme_bw()+
#   theme(text=element_text(size=15))+
#   labs(x="Shade tolerance",
#        y="PC2")

########################################
## examining simple relationships within the bivariate data

## as neighbor competition increases
## height and crown depth both decrease
## but the living fraction of the crown increases
## and the light at the crown base increases a bit?

qDTM_plastic<-ggplot(self_pruning,
                     aes(x=qDTM,
                         y=logLightBase,
                         color=Species))+
  geom_point()+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=15))

NCI_plastic<-ggplot(self_pruning,
                    aes(x=neighbor_comp,
                        y=logLightBase,
                        color=Species))+
  geom_point()+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=20))+
  labs(x="Indice de comp\u00e9tition (NCI)",
       y=expression(paste("log (",italic(L[base]),")")))+
  # labs(x="NCI",
  #      y="log(light fraction) at crown base")+
  guides(color="none")
# ggsave(filename="Images/NCI_plastic.png",NCI_plastic,
#        width=6,height=5,dpi=600)

neighbor_acq_plastic<-ggplot(self_pruning,
                             aes(x=neighbor_acq,
                                 y=logLightBase,
                                 color=Species))+
  geom_point()+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=15))+
  labs(x="Neighbor acquisitiveness",
       y=expression(paste("log (",italic(L[base]),")")))+
  guides(color="none")

# png("Images/lfbase_nacq_ind.png",width=7,height=5,units="in",res=150)
# neighbor_acq_plastic
# dev.off()

acq_dist_plastic<-ggplot(self_pruning,
                      aes(x=acq_dist_abs,
                          y=logLightBase,
                          color=Species))+
  geom_point()+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=15))+
  labs(x="Functional distance from neighbors",
       y="log(light fraction) at crown base")+
  guides(color="none")

## test of correlative inhibition:
## we should expect a positive slope here
light_top_plastic<-ggplot(self_pruning,
                          aes(x=logLightTop,
                              y=logLightBase,
                              color=Species))+
  geom_point()+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=15))+
  labs(x="log(light fraction) at crown top",
       y="log(light fraction) at crown base")

## this plot is kind of odd and visually striking...
height_plastic<-ggplot(self_pruning,
                       aes(x=HeightTop,
                           y=logLightBase,
                           color=Species))+
  geom_point()+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=15))+
  labs(x="Tree height",
       y="log(light fraction) at crown base")+
  guides(color="none")

## relationships with crown depth

NCI_plastic_CD<-ggplot(self_pruning,
                    aes(x=neighbor_comp,
                        y=CrownDepth,
                        color=Species))+
  geom_point()+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=15))+
  labs(x="NCI",
       y="Crown depth (cm)")+
  guides(color="none")

neighbor_acq_plastic_CD<-ggplot(self_pruning,
                                aes(x=neighbor_acq,
                                    y=CrownDepth,
                                    color=Species))+
  geom_point()+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=15))+
  labs(x="Neighbor acquisitiveness",
       y="Crown depth (cm)")+
  guides(color="none")

light_top_plastic_CD<-ggplot(self_pruning,
                          aes(x=logLightTop,
                              y=CrownDepth,
                              color=Species))+
  geom_point()+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=15))+
  labs(x="log(light fraction) at crown top",
       y="Crown depth (cm)")

height_plastic_CD<-ggplot(self_pruning,
                       aes(x=HeightTop,
                           y=CrownDepth,
                           color=Species))+
  geom_point()+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=15))+
  labs(x="Tree height",
       y="Crown depth (cm)")+
  guides(color="none")

################################################
## pull out the species-specific slopes from mixed-effects models

light_neighbor_acq_sp<-lmer(logLightBase~neighbor_acq*Species+(1|unique_plot),data=self_pruning)
species_means$light_neighbor_acq_slope<-rep(fixef(light_neighbor_acq_sp)[2],12)+c(0,fixef(light_neighbor_acq_sp)[14:24])
anova(light_neighbor_acq_sp, type="III")

light_NCI_sp<-lmer(logLightBase~neighbor_comp*Species+(1|unique_plot),data=self_pruning)
species_means$light_NCI_slope<-rep(fixef(light_NCI_sp)[2],12)+c(0,fixef(light_NCI_sp)[14:24])
anova(light_NCI_sp, type="III")

light_height_sp<-lmer(logLightBase~HeightTop*Species+(1|unique_plot),data=self_pruning)
species_means$light_height_slope<-rep(fixef(light_height_sp)[2],12)+c(0,fixef(light_height_sp)[14:24])
anova(light_height_sp, type="III")

## which species may show some evidence of correlative inhibition?
light_toplight_sp<-lmer(logLightBase~logLightTop*Species+(1|unique_plot),data=self_pruning)
species_means$light_toplight_slope<-rep(fixef(light_toplight_sp)[2],12)+c(0,fixef(light_toplight_sp)[15:25])
anova(light_toplight_sp, type="III")

neighbor_acq_slopes<-ggplot(data=species_means,
                            aes(x=focal_acq,
                                y=light_neighbor_acq_slope,
                                label=Species))+
  geom_smooth(method="lm")+
  geom_text(size=5,aes(color=leaf_habit))+
  scale_color_manual(values = leaf_habit_cols)+
  theme_bw()+
  theme(text=element_text(size=20))+
  labs(x="Focal acquisitiveness",
       y=expression(paste("Slopes: ",italic(L[base])," ~ neighbor acquisitiveness")))
# ggsave(filename = "Images/neighbor_acq_slopes.png",neighbor_acq_slopes,
#        dpi=600,width = 6,height=5)

NCI_slopes<-ggplot(data=species_means,
                   aes(x=focal_acq,
                       y=light_NCI_slope,
                       label=Species))+
  geom_smooth(method="lm")+
  geom_text(size=5,aes(color=leaf_habit))+
  theme_bw()+
  theme(text=element_text(size=20))+
  guides(color="none")+
  scale_color_manual(values = leaf_habit_cols)+
  labs(x="Focal acquisitiveness",
       y=expression(paste("Slopes: ",italic(L[base])," ~ NCI")))
  # labs(x="Tendance acquisitive (CP1)",
  #      y=expression(paste("Pentes : ",italic(L[base])," ~ NCI")))
# ggsave(filename = "Images/NCI_slopes_FR.png",NCI_slopes,
#        dpi=600,width = 6,height=5)

height_slopes<-ggplot(data=species_means,
                      aes(x=focal_acq,
                          y=light_height_slope,
                          label=Species))+
  geom_smooth(method="lm")+
  geom_text(size=5,aes(color=leaf_habit))+
  scale_color_manual(values = leaf_habit_cols)+
  theme_bw()+
  theme(text=element_text(size=20))+
  labs(x="Focal acquisitiveness",
       y=expression(paste("Slopes: ",italic(L[base])," ~ top height")))

light_top_slopes<-ggplot(data=species_means,
                         aes(x=focal_acq,
                             y=light_toplight_slope,
                             label=Species))+
  geom_smooth(method="lm")+
  geom_text(size=5,aes(color=leaf_habit))+
  scale_color_manual(values = leaf_habit_cols)+
  theme_bw()+
  theme(text=element_text(size=20))+
  labs(x="Focal acquisitiveness",
       y=expression(paste("Slopes: ",italic(L[base])," ~ ",italic(L[top]))))

plot_compile<-ggpubr::ggarrange(neighbor_acq_plastic,neighbor_acq_slopes,
                                NCI_plastic,NCI_slopes,
                                height_plastic,height_slopes,
                                light_top_plastic,light_top_slopes,
                                ncol=2,nrow=4,
                                common.legend=T,legend="bottom")

pdf("Images/Fig2XX.pdf",height=16,width=8)
plot_compile
dev.off()