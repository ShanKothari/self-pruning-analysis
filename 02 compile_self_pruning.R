setwd("C:/Users/Shan Kothari/Dropbox/PostdocProjects/SelfPruning")

library(ggplot2)
library(chron)
library(fishmethods)
library(plotrix)
library(patchwork)

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
neighbor_match<-match(self_pruning$unique_tree,neighbor.data$unique_tree)
self_pruning$neighbor_richness<-neighbor.data$neighbor.richness[neighbor_match]
self_pruning$FDis<-neighbor.data$FDis[neighbor_match]
self_pruning$qDTM<-neighbor.data$qDTM[neighbor_match]
self_pruning$neighbor_acq<- -neighbor.data$neighbor.FI1[neighbor_match]
self_pruning$neighbor_comp<-neighbor.data$NCI[neighbor_match]

## read in trait data to get shade tolerance
## and focal species functional ID (PC1)
trait_summary<-read.csv("TraitData/trait_summary.csv")
trait_match<-match(self_pruning$Species,trait_summary$SpeciesCode)
self_pruning$shade_tol<-trait_summary$shade_tol[trait_match]
self_pruning$focal_acq<- -trait_summary$PC1[trait_match]
self_pruning$acq_dist<- self_pruning$focal_acq-self_pruning$neighbor_acq
self_pruning$acq_dist_abs<- abs(self_pruning$acq_dist)

# write.csv(self_pruning,"SelfPruningData/self_pruning_processed.csv",row.names=F)

#######################################
## working with species means

## species means of crown base light
species_means<-aggregate(logLightBase~Species,data = self_pruning,
                         FUN = mean,na.rm = T)

## standard errors of crown base light
species_se<-aggregate(logLightBase~Species,data = self_pruning,
                         FUN = std.error,na.rm = T)
species_means$logLightBase_se<-species_se$logLightBase[match(species_means$Species,
                                                               species_se$Species)]

## species means of crown base light
## from only monoculture plots
self_pruning_mono<-self_pruning[self_pruning$nbsp==1,]
mono_means<-aggregate(logLightBase~Species,data = self_pruning_mono,
                      FUN = mean,na.rm = T)

species_means$logLightBase_mono<-mono_means$logLightBase[match(species_means$Species,
                                                               mono_means$Species)]

## attach traits to species means
trait_match_means<-match(species_means$Species,trait_summary$SpeciesCode)
species_means$shade_tol<-trait_summary$shade_tol[trait_match_means]
species_means$focal_acq<- -trait_summary$PC1[trait_match_means]
species_means$PC2<- -trait_summary$PC2[trait_match_means]
species_means$leaf_lifespan<- trait_summary$LL[trait_match_means]
species_means$leaf_habit<- trait_summary$leaf_habit[trait_match_means]

# write.csv(species_means,"SelfPruningData/species_means.csv",row.names=F)

# summary(lm(logLightBase~shade_tol+focal_acq,data=species_means))

leaf_habit_cols<-c("Deciduous"="chocolate4",
                   "Evergreen"="#60941a")

st_acq_sp<-ggplot(data=species_means,
                     aes(x=shade_tol,
                         y=focal_acq,
                         label=Species))+
  geom_smooth(method="lm")+
  geom_point(size=1.5)+
  geom_text(size=5,aes(color=leaf_habit))+
  theme_bw()+
  theme(text=element_text(size=20))+
  guides(color="none")+
  scale_color_manual(values = leaf_habit_cols)+
  labs(x="Shade tolerance",
       y="Acquisitiveness")
# ggsave(filename = "Images/st_acq_sp.png",plot = st_acq_sp,
#        width=6,height=5,units="in",dpi=600)

lfbase_st_sp<-ggplot(data=species_means,
                     aes(x=shade_tol,
                         y=logLightBase,
                         label=Species))+
  geom_smooth(method="lm")+
  geom_point(size=1.5)+
  geom_text(size=5,aes(color=leaf_habit))+
  theme_bw()+
  theme(text=element_text(size=20),
        panel.grid.minor = element_blank())+
  coord_cartesian(ylim=c(-6.5,-2.5),
                  xlim=c(0.8,5.2))+
  guides(color="none")+
  scale_color_manual(values = leaf_habit_cols)+
  # labs(x="Tol\u00e9rance \u00e0 l'ombre",
  #      y=expression(italic(L[base])))
  labs(x="Shade tolerance",
       y=expression(italic(L[base])))
# ggsave(filename = "Images/lfbase_st_sp.png",plot = lfbase_st_sp,
#        width=6,height=5,units="in",dpi=600)

lfbase_acq_sp<-ggplot(data=species_means,
                      aes(x=focal_acq,
                          y=logLightBase,
                          label=Species))+
  geom_smooth(method="lm")+
  geom_point(size=1.5)+
  geom_text(size=5,aes(color=leaf_habit))+
  theme_bw()+
  theme(text=element_text(size=20),
        panel.grid.minor = element_blank())+
  coord_cartesian(ylim=c(-6.5,-2.5),
                  xlim=c(-2.7,2.5))+
  guides(color="none")+
  scale_color_manual(values = leaf_habit_cols)+
  # labs(x="Tendance acquisitive (CP1)",
  #      y=expression(italic(L[base])))
  labs(x="Focal tree acquisitiveness",
       y=expression(italic(L[base])))
# ggsave(filename = "Images/lfbase_acq_sp.png", plot = lfbase_acq_sp,
#        width=6,height=5,units="in",dpi=600)

lfbase_ll_sp<-ggplot(data=species_means,
                     aes(x=log(leaf_lifespan),
                         y=logLightBase,
                         label=Species))+
  geom_smooth(method="lm")+
  geom_smooth(data=species_means[species_means$leaf_habit=="Evergreen",],
              aes(x=log(leaf_lifespan),y=logLightBase),
              method="lm",color="red",se=F)+
  geom_text(size=5,aes(color=leaf_habit))+
  theme_bw()+
  theme(text=element_text(size=20),
        legend.position = c(0.7, 0.8),
        legend.background = element_rect(fill = NA),
        panel.grid.minor = element_blank())+
  coord_cartesian(ylim=c(-6.5,-2.5),
                  xlim=c(1.2,4.5))+
  scale_color_manual(values = leaf_habit_cols)+
  labs(x="log(leaf lifespan [months])",
       y=expression(italic(L[base])),
       color="Leaf habit")
# ggsave(filename = "Images/lfbase_ll_sp.png",plot = lfbase_ll_sp,
#        width=6,height=5,units="in",dpi=600)

mono_comparison<-ggplot(data=species_means,
                        aes(x=logLightBase,
                            y=logLightBase_mono,
                            label=Species))+
  geom_abline(slope=1,intercept=0,linewidth=2,linetype="dashed")+
  geom_smooth(method="lm")+
  geom_smooth(data=species_means[species_means$leaf_habit=="Evergreen",],
              aes(x=log(leaf_lifespan),y=logLightBase),
              method="lm",color="red",se=F)+
  geom_text(size=5,aes(color=leaf_habit))+
  theme_bw()+
  theme(text=element_text(size=20),
        legend.position = c(0.7, 0.3),
        legend.background = element_rect(fill = NA),
        panel.grid.minor = element_blank())+
  coord_cartesian(ylim=c(-6.5,-2.7),
                  xlim=c(-6.5,-2.7))+
  scale_color_manual(values = leaf_habit_cols)+
  labs(y=expression(paste(italic(L[base])," in monocultures")),
       x=expression(italic(L[base])),
       color="Leaf habit")

pdf("Images/FigS_mono_comp.pdf",width = 5,height = 5)
mono_comparison
dev.off()

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

# pdf("Images/Fig2.pdf",width = 5,height=13)
# lfbase_st_sp/lfbase_acq_sp/lfbase_ll_sp
# dev.off()
