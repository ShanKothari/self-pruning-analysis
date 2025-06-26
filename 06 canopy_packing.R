setwd("C:/Users/Shan Kothari/Dropbox/PostdocProjects/SelfPruning")

library(mosaic)
library(reshape2)
library(patchwork)

########################################
## data on crown shape and size

## this self-pruning data includes random sample of living trees
## but has no information on mortality rates, which we will have
## to account for later
self_pruning<-read.csv("SelfPruningData/self_pruning_processed.csv")

## in these three plots, there are species that were living
## in the plot in 2018 that were accidentally not measured 
## in the self-pruning survey. in A_4N8, THOC is missing;
## in D_4N3, PIST is missing; in A_2N3, BEPA (misplanted)
## is missing. these issues would otherwise threaten to
## bias crown packing calculations
self_pruning<-subset(self_pruning,!(unique_plot %in% c("A_2N3","A_4N8","D_4N3")))

## get indicators of species composition
self_pruning_list<-split(self_pruning,f = self_pruning$Plot)
self_pruning_comp<-unlist(lapply(self_pruning_list,
                                 function(plot) paste(unique(plot$Species),collapse="|")))
self_pruning$sp_comp<-self_pruning_comp[match(self_pruning$Plot,names(self_pruning_comp))]

## parameters from Purves et al. paper
C0B<-0.196
C1B<-0.511

species.list<-c("ABBA","ACRU","ACSA","BEAL","BEPA","LALA",
                "PIGL","PIRE","PIRU","PIST","QURU","THOC")
Tj<-c(0.278,0.536,0.560,0.592,0.435,0.308,
      0.278,0.287,0.324,0.417,0.538,0.251)

crown_shape<-data.frame(species=species.list,Tj=Tj)
crown_shape$Bj<-(1-crown_shape$Tj)*C0B+crown_shape$Tj*C1B
self_pruning$Bj<-crown_shape$Bj[match(self_pruning$Species,crown_shape$species)]

## formulas for area/volume of crown
crown_area<-function(CD,CR,beta){(CR*CD)/(beta+1)}
crown_vol<-function(CD,CR,beta){(pi*CR^2*CD)/(2*beta+1)}

## calculate for each row of the self-pruning data
self_pruning$crown_vol<-crown_vol(CD=self_pruning$CrownDepth/100,
                                  CR=self_pruning$CR_average/100,
                                  beta=self_pruning$Bj)

## read in mortality data
mortality_2018<-read.csv("IDENTMontrealData/mortality_2018.csv")
mortality_2018$unique_plot<-paste(mortality_2018$Block,
                                 mortality_2018$Plot,
                                 sep="_")

mortality_2018$plot_sp<-paste(mortality_2018$unique_plot,
                              mortality_2018$CodeSp,
                              sep="_")

## proportion of species planted within the inner 36 and
## proportion of planted individuals of the species still alive
mortality_2018$prop_planted<-mortality_2018$num_planted/36
mortality_2018$prop_alive<-mortality_2018$num_alive/mortality_2018$num_planted

######################################
## calculate canopy packing

## take the mean crown volume for each species in each plot
crown_vol_agg<-aggregate(crown_vol~Block+Plot+Species+nbsp,
                         data=self_pruning,FUN=mean)

colnames(crown_vol_agg)<-c("Block","Plot","Species","Richness","crown_vol")

## columns to match against for adding mortality and monoculture data
crown_vol_agg$block_plot<-paste(crown_vol_agg$Block,
                                crown_vol_agg$Plot,
                                sep="_")

crown_vol_agg$block_sp<-paste(crown_vol_agg$Block,
                              crown_vol_agg$Species,
                              sep="_")

crown_vol_agg$plot_sp<-paste(crown_vol_agg$block_plot,
                             crown_vol_agg$Species,
                             sep="_")

## attach numbers alive and planted in mixture
mortality_match<-match(crown_vol_agg$plot_sp,mortality_2018$plot_sp)
crown_vol_agg$num_alive<-mortality_2018$num_alive[mortality_match]
crown_vol_agg$num_planted<-mortality_2018$num_planted[mortality_match]
crown_vol_agg$prop_alive<-mortality_2018$prop_alive[mortality_match]

## attach mean crown volume and mortality rates in monoculture
## here we take advantage of the fact that monocultures have a
## plot name that is the same as their species code
mono_match<-match(crown_vol_agg$block_sp,crown_vol_agg$block_plot)
crown_vol_agg$crown_vol_mono<-crown_vol_agg$crown_vol[mono_match]
crown_vol_agg$prop_alive_mono<-crown_vol_agg$prop_alive[mono_match]

## multiply by proportions alive to account for mortality
crown_vol_agg$crown_vol_mono_live<-crown_vol_agg$crown_vol_mono*crown_vol_agg$prop_alive_mono
crown_vol_agg$crown_vol_live<-crown_vol_agg$crown_vol*crown_vol_agg$prop_alive
crown_vol_agg$total_sp_crown_vol<-crown_vol_agg$crown_vol_live*crown_vol_agg$num_planted

##########################################
## simulations of canopy packing holding crown depth constant

## number of crowns to simulate for each species in each plot
n_crowns<-500

## first, we simulate trees where the crown radius is sampled from
## species in the given plot and crown depth is sampled from
## the corresponding monocultures
sim_crown_vols<-list()
for(i in crown_vol_agg$plot_sp){

  plot_sp_split<-strsplit(i,"_")[[1]]
  block<-plot_sp_split[1]
  plot<-plot_sp_split[2]
  species<-plot_sp_split[3]
  
  ## selecting mixtures & monocultures
  mix_ids<-which(self_pruning$Block==block & self_pruning$Plot==plot & self_pruning$Species==species)
  ## here we take advantage of the fact that monoculture plot names are species names
  mono_ids<-which(self_pruning$Block==block & self_pruning$Plot==species & self_pruning$Species==species)
  
  ## draw crown radius from mixtures
  ## but crown depth from monocultures
  sim_plot_sp<-data.frame(CR_average=sample(self_pruning$CR_average[mix_ids],
                                            size=n_crowns,replace=T),
                          CrownDepth=sample(self_pruning$CrownDepth[mono_ids],
                                            size=n_crowns,replace=T),
                          Bj=sample(self_pruning$Bj[mono_ids],
                                    size=n_crowns,replace=T))
  
  sim_crown_vols[i]<-mean(crown_vol(CD=sim_plot_sp$CrownDepth/100,
                                    CR=sim_plot_sp$CR_average/100,
                                    beta=sim_plot_sp$Bj))
  
}

## next, we simulate trees where the both the crown radius
## and the crown depth are sampled from species in the given plot

## the purpose of this simulation is to test sensitivity to
## the potential bias caused by random matching of crown radius and crown depth
null_sim_crown_vols<-list()
for(i in crown_vol_agg$plot_sp){
  
  plot_sp_split<-strsplit(i,"_")[[1]]
  block<-plot_sp_split[1]
  plot<-plot_sp_split[2]
  species<-plot_sp_split[3]
  
  mix_ids<-with(self_pruning,which(Block==block & Plot==plot & Species==species))
  mono_ids<-with(self_pruning,which(Block==block & Plot==species & Species==species))
  
  ## draw crown radius and crown depth from mixtures
  null_sim_plot_sp<-data.frame(CR_average=sample(self_pruning$CR_average[mix_ids],
                                            size=n_crowns,replace=T),
                          CrownDepth=sample(self_pruning$CrownDepth[mix_ids],
                                            size=n_crowns,replace=T),
                          Bj=sample(self_pruning$Bj[mono_ids],
                                    size=n_crowns,replace=T))
  
  null_sim_crown_vols[i]<-mean(crown_vol(CD=null_sim_plot_sp$CrownDepth/100,
                                    CR=null_sim_plot_sp$CR_average/100,
                                    beta=null_sim_plot_sp$Bj))
  
}

## multiply mean volume by mixture mortality rate
crown_vol_agg$sim_vol<-unlist(sim_crown_vols)
crown_vol_agg$sim_vol_live<-crown_vol_agg$sim_vol*crown_vol_agg$prop_alive
crown_vol_agg$total_sp_sim_vol<-crown_vol_agg$sim_vol_live*crown_vol_agg$num_planted

crown_vol_agg$null_sim_vol<-unlist(null_sim_crown_vols)
crown_vol_agg$null_sim_vol_live<-crown_vol_agg$null_sim_vol*crown_vol_agg$prop_alive
crown_vol_agg$total_sp_null_sim_vol<-crown_vol_agg$null_sim_vol_live*crown_vol_agg$num_planted

## draw monocultures for simulation
crown_vol_agg$crown_vol_mono_sim<-crown_vol_agg$sim_vol[mono_match]
crown_vol_agg$crown_vol_mono_sim_live<-crown_vol_agg$crown_vol_mono_sim*crown_vol_agg$prop_alive_mono

crown_vol_agg$crown_vol_mono_null_sim<-crown_vol_agg$null_sim_vol[mono_match]
crown_vol_agg$crown_vol_mono_null_sim_live<-crown_vol_agg$crown_vol_mono_null_sim*crown_vol_agg$prop_alive_mono

######################################
## aggregating to the plot scale

## aggregating overyielding to the species x plot level
## rather than just individuals
crown_vol_agg$NBE_actual<-with(crown_vol_agg,
                              num_planted*(crown_vol_live-crown_vol_mono_live))
crown_vol_agg$NBE_sim<-with(crown_vol_agg,
                           num_planted*(sim_vol_live-crown_vol_mono_sim_live))
crown_vol_agg$NBE_null_sim<-with(crown_vol_agg,
                           num_planted*(null_sim_vol_live-crown_vol_mono_null_sim_live))

## and aggregating to the whole plot (inner 6 x 6 trees)
## since planting numbers / mortalities are only calculated
## within those inner 6 x 6
crown_vol_agg_sub<-crown_vol_agg[,c("Block","Plot","Richness",
                                    "NBE_actual","NBE_sim","NBE_null_sim",
                                    "total_sp_crown_vol",
                                    "total_sp_sim_vol",
                                    "total_sp_null_sim_vol")]

crown_vol_plot<-aggregate(.~Block+Plot+Richness,
                          data=crown_vol_agg_sub,
                          FUN=sum)

crown_vol_plot$NBE_actual[which(crown_vol_plot$Richness==1)]<-NA
crown_vol_plot$NBE_sim[which(crown_vol_plot$Richness==1)]<-NA
crown_vol_plot$NBE_null_sim[which(crown_vol_plot$Richness==1)]<-NA

## all values should be divided by 9 for analysis per m^2
## since upscaling was done to the 6 x 6 inner plot (no edge)

## read in plot-level heterogeneity measures
plot_vars<-read.csv("IDENTMontrealData/plot_vars.csv")

crown_vol_plot$unique_plot<-paste(crown_vol_plot$Block,crown_vol_plot$Plot,sep="_")
crown_vol_plot$BasalArea<-plot_vars$BasalArea[match(crown_vol_plot$unique_plot,plot_vars$unique_plot)]

crown_vol_plot$FDis<-plot_vars$FDis[match(crown_vol_plot$unique_plot,plot_vars$unique_plot)]
crown_vol_plot$FTD<-plot_vars$FTD[match(crown_vol_plot$unique_plot,plot_vars$unique_plot)]
crown_vol_plot$FTD_LH<-plot_vars$FTD_LH[match(crown_vol_plot$unique_plot,plot_vars$unique_plot)]

##########################################
## plotting canopy packing results

colorBlind  <- c("#E69F00","#009E73","#56B4E9","#F0E442",
                 "#0072B2","#CC79A7","#D55E00","#999999")

## this should be very tightly clustered around the 1:1 line
## because the 'null simulation' should be close to real values
null_vs_real<-ggplot(crown_vol_plot,aes(x=total_sp_crown_vol/9,
                                        y=total_sp_null_sim_vol/9,
                                        color=as.factor(Richness)))+
  geom_point(size=2)+
  geom_abline(slope=1,intercept=0,linewidth=2,linetype="dashed")+
  theme_bw()+
  theme(text=element_text(size=15),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  coord_cartesian(xlim=c(0,16),ylim=c(0,16))+
  labs(x="Actual crowns approach (m)",
       y="Independent draws with plasticity (m)",
       color="Richness")+
  guides(color="none")+
  scale_color_manual(values=colorBlind)

## this reveals that if you draw crown depth from monocultures,
## you generally get less total crown volume
sim_vs_real<-ggplot(crown_vol_plot,
                    aes(x=total_sp_crown_vol/9,
                        y=total_sp_sim_vol/9,
                        color=as.factor(Richness)))+
  geom_point(size=2)+
  geom_abline(slope=1,intercept=0,linewidth=2,linetype="dashed")+
  theme_bw()+
  theme(text=element_text(size=15))+
  coord_cartesian(xlim=c(0,16),ylim=c(0,16))+
  labs(x="Actual crowns approach (m)",
       y="Independent draws without plasticity (m)",
       color="Richness")+
  scale_color_manual(values=colorBlind)

# pdf("Images/FigS2.pdf",height=9,width=5)
# null_vs_real/sim_vs_real +
#   plot_layout(guides = "collect") &
#   theme(legend.position = 'bottom')
# dev.off()

## wide to long
crown_vol_plot_sub<-crown_vol_plot[,c("unique_plot","Richness","FTD","FTD_LH","NBE_actual","NBE_sim","NBE_null_sim",
                                      "total_sp_crown_vol","total_sp_sim_vol","total_sp_null_sim_vol")]
crown_vol_plot_long<-melt(crown_vol_plot_sub,id.vars=c("unique_plot","Richness","FTD","FTD_LH"))

crown_vol_plot_total<-crown_vol_plot_long[which(crown_vol_plot_long$variable
                                                %in% c("total_sp_null_sim_vol","total_sp_sim_vol")),]
crown_vol_plot_NBE<-crown_vol_plot_long[which(crown_vol_plot_long$variable
                                                %in% c("NBE_null_sim","NBE_sim")),]

ggplot(crown_vol_plot_long,
       aes(x=log(FTD_LH),
           y=value/9,
           color=variable))+
  geom_smooth(method="lm",se=F)+geom_point()+
  theme_bw()+
  labs(x=expression("log("^q*"D(TM)) of "*italic("L"*""[base])),
       y="Crown volume or NBE per ground area (m)")

plasticity_comparison_FT<-ggplot(crown_vol_plot_total,
                                 aes(x=log(FTD),
                                     y=value/9,
                                     color=variable))+
  geom_smooth(method="lm",se=F,linewidth=2)+
  geom_point(size=2)+
  theme_bw()+
  theme(text=element_text(size=20),
        legend.position = "bottom",
        plot.margin = unit(c(0,0.5,0,0), "cm"))+
  scale_color_hue(labels = c("Without plasticity", "With plasticity"))+
  labs(x=expression("log("^q*"D(TM)) of functional traits"),
       y="Crown volume per ground area (m)",
       color="")

plasticity_comparison_Lbase<-ggplot(crown_vol_plot_total,
                              aes(x=log(FTD_LH),
                                  y=value/9,
                                  color=variable))+
  geom_smooth(method="lm",se=F,linewidth=2)+
  geom_point(size=2)+
  theme_bw()+
  theme(text=element_text(size=20),
        legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  scale_color_hue(labels = c("Without plasticity", "With plasticity"))+
  labs(x=expression("log("^q*"D(TM)) of "*italic("L"*""[base])),
       y="Crown volume per ground area (m)",
       color="")

# pdf("Images/Fig6.pdf",height=6,width=10)
# plasticity_comparison_FT + plasticity_comparison_Lbase +
#   plot_layout(guides = "collect") & theme(legend.position = "bottom")
# dev.off()

NBE_FT<-ggplot(crown_vol_plot_NBE,
                 aes(x=log(FTD),
                     y=value/9,
                     color=variable))+
  geom_smooth(method="lm",se=F,linewidth=2)+
  geom_point(size=2)+
  theme_bw()+
  theme(text=element_text(size=20),
        legend.position = "bottom",
        plot.margin = unit(c(0,0.5,0,0), "cm"))+
  scale_color_hue(labels = c("Without plasticity", "With plasticity"))+
  labs(x=expression("log("^q*"D(TM)) of functional traits"),
       y="NBE on crown volume (m)",
       color="")

NBE_Lbase<-ggplot(crown_vol_plot_NBE,
                 aes(x=log(FTD_LH),
                     y=value/9,
                     color=variable))+
  geom_smooth(method="lm",se=F,linewidth=2)+
  geom_point(size=2)+
  theme_bw()+
  theme(text=element_text(size=20),
        legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  scale_color_hue(labels = c("Without plasticity", "With plasticity"))+
  labs(x=expression("log("^q*"D(TM)) of "*italic("L"*""[base])),
       y="NBE on crown volume (m)",
       color="")

# pdf("Images/FigS5.pdf",height=6,width=10)
# NBE_FT + NBE_Lbase +
#   plot_layout(guides = "collect") & theme(legend.position = "bottom")
# dev.off()

summary(lm(total_sp_null_sim_vol~log(FTD_LH),data=crown_vol_plot))
summary(lm(total_sp_sim_vol~log(FTD_LH),data=crown_vol_plot))

summary(lm(total_sp_null_sim_vol~log(FTD),data=crown_vol_plot))
summary(lm(total_sp_sim_vol~log(FTD),data=crown_vol_plot))

t.test(crown_vol_plot$total_sp_null_sim_vol[crown_vol_plot$Richness!=1]/9,
       crown_vol_plot$total_sp_sim_vol[crown_vol_plot$Richness!=1]/9,
       paired=T)

t.test(crown_vol_plot$NBE_null_sim/9)
t.test(crown_vol_plot$NBE_sim/9)