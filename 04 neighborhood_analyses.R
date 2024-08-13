setwd("C:/Users/querc/Dropbox/PostdocProjects/SelfPruning")

library(ggplot2)
library(ggpubr)
library(GGally)
library(lme4)
library(lmerTest)
library(emmeans)

self_pruning<-read.csv("SelfPruningData/self_pruning_processed.csv")

########################################
## relationships among individual-level variables

ind_level_vars<-c("neighbor_acq","neighbor_comp","HeightTop","logLightTop")

pair_plot<-ggpairs(self_pruning,
                   columns=ind_level_vars,
                   aes(color=Species),
                   lower = list(continuous = wrap("smooth",
                                                  method = "lm",
                                                  se=F)))

# pdf("Images/FigS3.pdf",height=16,width=16)
# pair_plot
# dev.off()

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
  theme(text=element_text(size=20))

neighbor_acq_plastic<-ggplot(self_pruning,
                             aes(x=neighbor_acq,
                                 y=logLightBase,
                                 color=Species))+
  geom_point()+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=20))+
  labs(x="Neighbor acquisitiveness",
       y=expression(italic(L[base])))+
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
  theme(text=element_text(size=20))+
  labs(x="Functional distance from neighbors",
       y=expression(italic(L[base])))+
  guides(color="none")

NCI_plastic<-ggplot(self_pruning,
                    aes(x=neighbor_comp,
                        y=logLightBase,
                        color=Species))+
  geom_point()+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=20))+
  # labs(x="Indice de comp\u00e9tition (NCI)",
  #      y=expression(italic(L[base])))+
  labs(x="NCI",
       y=expression(italic(L[base])))+
  guides(color="none")
# ggsave(filename="Images/NCI_plastic.png",NCI_plastic,
#        width=6,height=5,dpi=600)

## test of correlative inhibition:
## we should expect a positive slope here
light_top_plastic<-ggplot(self_pruning,
                          aes(x=logLightTop,
                              y=logLightBase,
                              color=Species))+
  geom_point()+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=20))+
  labs(x=expression(italic(L[top])),
       y=expression(italic(L[base])))

## this plot is kind of odd and visually striking...
height_plastic<-ggplot(self_pruning,
                       aes(x=HeightTop,
                           y=logLightBase,
                           color=Species))+
  geom_point()+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=20))+
  labs(x="Tree height",
       y=expression(italic(L[base])))+
  guides(color="none")

############################################
## relationships with crown depth

neighbor_acq_plastic_CD<-ggplot(self_pruning,
                                aes(x=neighbor_acq,
                                    y=CrownDepth,
                                    color=Species))+
  geom_point()+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=20))+
  labs(x="Neighbor acquisitiveness",
       y="Crown depth (cm)")+
  guides(color="none")

NCI_plastic_CD<-ggplot(self_pruning,
                       aes(x=neighbor_comp,
                           y=CrownDepth,
                           color=Species))+
  geom_point()+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=20))+
  labs(x="NCI",
       y="Crown depth (cm)")+
  guides(color="none")

light_top_plastic_CD<-ggplot(self_pruning,
                             aes(x=logLightTop,
                                 y=CrownDepth,
                                 color=Species))+
  geom_point()+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=20))+
  labs(x=expression(italic(L[top])),
       y="Crown depth (cm)")+
  guides(color="none")

height_plastic_CD<-ggplot(self_pruning,
                          aes(x=HeightTop,
                              y=CrownDepth,
                              color=Species))+
  geom_point()+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=20))+
  labs(x="Tree height",
       y="Crown depth (cm)")

plot_compile_CD<-ggpubr::ggarrange(neighbor_acq_plastic_CD,NCI_plastic_CD,
                                   height_plastic_CD,light_top_plastic_CD,
                                   ncol=1,nrow=4,
                                   common.legend=T,legend="bottom")

# pdf("Images/Fig5.pdf",height=24,width=7)
# plot_compile_CD
# dev.off()

neighbor_acq_plastic_CB<-ggplot(self_pruning,
                                aes(x=neighbor_acq,
                                    y=HeightBase,
                                    color=Species))+
  geom_point()+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=20),
        legend.position = "bottom")+
  labs(x="Neighbor acquisitiveness",
       y="Crown base height (cm)")

# pdf("Images/FigS4.pdf",height=8,width=7)
# neighbor_acq_plastic_CB
# dev.off()

################################################
## pull out the species-specific slopes from mixed-effects models

species_means<-read.csv("SelfPruningData/species_means.csv")

## no comparisons; REML = T by default
light_neighbor_acq_sp<-lmer(logLightBase~neighbor_acq*Species+(1|unique_plot),data=self_pruning)
light_neighbor_acq_trend<-emtrends(light_neighbor_acq_sp,"Species",var="neighbor_acq")
species_means$light_neighbor_acq_slope<-summary(light_neighbor_acq_trend)$neighbor_acq.trend
anova(light_neighbor_acq_sp, type="III")
## very similar results
# Anova(light_neighbor_acq_sp, type="III",test.statistic = "F")

light_NCI_sp<-lmer(logLightBase~neighbor_comp*Species+(1|unique_plot),data=self_pruning)
light_NCI_trend<-emtrends(light_NCI_sp,"Species",var="neighbor_comp")
species_means$light_NCI_slope<-summary(light_NCI_trend)$neighbor_comp.trend
anova(light_NCI_sp, type="III")

light_height_sp<-lmer(logLightBase~HeightTop*Species+(1|unique_plot),data=self_pruning)
light_height_trend<-emtrends(light_height_sp,"Species",var="HeightTop")
species_means$light_height_slope<-summary(light_height_trend)$HeightTop.trend
anova(light_height_sp, type="III")

## which species may show some evidence of correlative inhibition?
light_toplight_sp<-lmer(logLightBase~logLightTop*Species+(1|unique_plot),data=self_pruning)
light_toplight_trend<-emtrends(light_toplight_sp,"Species",var="logLightTop")
species_means$light_toplight_slope<-summary(light_toplight_trend)$logLightTop.trend
anova(light_toplight_sp, type="III")

leaf_habit_cols<-c("Deciduous"="chocolate4",
                   "Evergreen"="#60941a")

neighbor_acq_slopes<-ggplot(data=species_means,
                            aes(x=focal_acq,
                                y=light_neighbor_acq_slope,
                                label=Species))+
  geom_smooth(method="lm")+
  geom_text(size=5,aes(color=leaf_habit))+
  scale_color_manual(values = leaf_habit_cols)+
  theme_bw()+
  theme(text=element_text(size=20))+
  guides(color="none")+
  labs(x="Focal acquisitiveness",
       y=expression(paste("Slopes: ",
                          italic(L[base]),
                          " ~ neighbor acquisitiveness")))
# ggsave(filename = "Images/neighbor_acq_slopes.png",neighbor_acq_slopes,
#        dpi=600,width = 6,height=5)

NCI_slopes<-ggplot(data=species_means,
                   aes(x=focal_acq,
                       y=light_NCI_slope*100,
                       label=Species))+
  geom_smooth(method="lm")+
  geom_text(size=5,aes(color=leaf_habit))+
  theme_bw()+
  theme(text=element_text(size=20))+
  guides(color="none")+
  scale_color_manual(values = leaf_habit_cols)+
  labs(x="Focal acquisitiveness",
       y=expression(paste("Slopes (" %*% "100): ",
                          italic(L[base]),
                          " ~ NCI")))
# labs(x="Tendance acquisitive (CP1)",
#      y=expression(paste("Pentes : ",italic(L[base])," ~ NCI")))
# ggsave(filename = "Images/NCI_slopes_FR.png",NCI_slopes,
#        dpi=600,width = 6,height=5)

height_slopes<-ggplot(data=species_means,
                      aes(x=focal_acq,
                          y=light_height_slope*1000,
                          label=Species))+
  geom_smooth(method="lm")+
  geom_text(size=5,aes(color=leaf_habit))+
  scale_color_manual(values = leaf_habit_cols)+
  theme_bw()+
  theme(text=element_text(size=20))+
  guides(color="none")+
  labs(x="Focal acquisitiveness",
       y=expression(paste("Slopes (" %*% "1000): ",
                          italic(L[base]),
                          " ~ top height")))

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
       y=expression(paste("Slopes: ",
                          italic(L[base]),
                          " ~ ",italic(L[top]))),
       color="Leaf habit")

plot_compile_Lbase<-ggpubr::ggarrange(neighbor_acq_plastic,neighbor_acq_slopes,
                                      NCI_plastic,NCI_slopes,
                                      height_plastic,height_slopes,
                                      light_top_plastic,light_top_slopes,
                                      ncol=2,nrow=4,
                                      common.legend=T,legend="bottom")

# pdf("Images/Fig4.pdf",height=24,width=12)
# plot_compile_Lbase
# dev.off()

CD_neighbor_acq_sp<-lmer(CrownDepth~neighbor_acq*Species+(1|unique_plot),data=self_pruning)
species_means$CD_neighbor_acq_slope<-rep(fixef(CD_neighbor_acq_sp)[2],12)+c(0,fixef(CD_neighbor_acq_sp)[14:24])
anova(CD_neighbor_acq_sp, type="III")

CD_NCI_sp<-lmer(CrownDepth~neighbor_comp*Species+(1|unique_plot),data=self_pruning)
species_means$CD_NCI_slope<-rep(fixef(CD_NCI_sp)[2],12)+c(0,fixef(CD_NCI_sp)[14:24])
anova(CD_NCI_sp, type="III")

CD_height_sp<-lmer(CrownDepth~HeightTop*Species+(1|unique_plot),data=self_pruning)
species_means$CD_height_slope<-rep(fixef(CD_height_sp)[2],12)+c(0,fixef(CD_height_sp)[14:24])
anova(CD_height_sp, type="III")

CD_toplight_sp<-lmer(CrownDepth~logLightTop*Species+(1|unique_plot),data=self_pruning)
species_means$CD_toplight_slope<-rep(fixef(CD_toplight_sp)[2],12)+c(0,fixef(CD_toplight_sp)[14:24])
anova(CD_toplight_sp, type="III")

