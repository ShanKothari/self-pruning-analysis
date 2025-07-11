setwd("C:/Users/Shan Kothari/Dropbox/PostdocProjects/SelfPruning")

library(ggplot2)
library(ggpubr)
library(GGally)
library(lme4)
library(lmerTest)
library(emmeans)
library(patchwork)

self_pruning<-read.csv("SelfPruningData/self_pruning_processed.csv")
species_means<-read.csv("SelfPruningData/species_means.csv")
self_pruning$leaf_habit<-species_means$leaf_habit[match(self_pruning$Species,species_means$Species)]

self_pruning$Species<-factor(self_pruning$Species,
                             levels=c("ACRU","ACSA","BEAL","BEPA","LALA","QURU",
                                      "ABBA","PIGL","PIRE","PIRU","PIST","THOC"))

trial_scale<-c("#D17711","#FFFF00","#F2681F",
               "#FDA600","#FFC58F","chocolate4",
               "#003C96","#11A69C","#924AF7",
               "#0081FE","#60941A","#A4F9AC")

########################################
## relationships among individual-level variables

ind_level_vars<-c("neighbor_acq","neighbor_comp","HeightTop","logLightTop")

pair_plot<-ggpairs(self_pruning,
                   columns=ind_level_vars,
                   aes(color=Species),
                   lower = list(continuous = wrap("smooth",
                                                  method = "lm",
                                                  se=F)))+
  theme_bw()+
  theme(text=element_text(size=25))+
  scale_color_manual(values=trial_scale)

# pdf("Images/FigS3.pdf",height=14,width=15)
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
  geom_point(alpha=0.5)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=20))

neighbor_acq_plastic<-ggplot(self_pruning,
                             aes(x=neighbor_acq,
                                 y=logLightBase,
                                 color=Species))+
  geom_point(alpha=0.5)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=20))+
  labs(x="Neighbor acquisitiveness",
       y=expression(italic(L[base])),
       title="Individuals")+
  guides(color="none")+
  scale_color_manual(values=trial_scale)

# png("Images/lfbase_nacq_ind.png",width=7,height=5,units="in",res=150)
# neighbor_acq_plastic
# dev.off()

acq_dist_plastic<-ggplot(self_pruning,
                         aes(x=acq_dist_abs,
                             y=logLightBase,
                             color=Species))+
  geom_point(alpha=0.5)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=20))+
  labs(x="Functional distance from neighbors",
       y=expression(italic(L[base])))+
  guides(color="none")+
  scale_color_manual(values=trial_scale)

NCI_plastic<-ggplot(self_pruning,
                    aes(x=neighbor_comp,
                        y=logLightBase,
                        color=Species))+
  geom_point(alpha=0.5)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=20))+
  # labs(x="Indice de comp\u00e9tition (NCI)",
  #      y=expression(italic(L[base])))+
  labs(x="NCI",
       y=expression(italic(L[base])))+
  guides(color="none")+
  scale_color_manual(values=trial_scale)

# ggsave(filename="Images/NCI_plastic.png",NCI_plastic,
#        width=6,height=5,dpi=600)

## test of correlative inhibition:
## we should expect a positive slope here
light_top_plastic<-ggplot(self_pruning,
                          aes(x=logLightTop,
                              y=logLightBase,
                              color=Species))+
  geom_point(alpha=0.5)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=20))+
  labs(x=expression(italic(L[top])),
       y=expression(italic(L[base])))+
  guides(color = guide_legend(position = "bottom"))+
  scale_color_manual(values=trial_scale)

## this plot is kind of odd and visually striking...
height_plastic<-ggplot(self_pruning,
                       aes(x=HeightTop,
                           y=logLightBase,
                           color=Species))+
  geom_point(alpha=0.5)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=20))+
  labs(x="Tree height (cm)",
       y=expression(italic(L[base])))+
  guides(color="none")+
  scale_color_manual(values=trial_scale)

depth_plastic<-ggplot(self_pruning,
                       aes(y=CrownDepth,
                           x=logLightBase,
                           color=Species))+
  geom_point(alpha=0.5)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=20),
        legend.position = "bottom")+
  labs(y="Crown depth (cm)",
       x=expression(italic(L[base])))+
  scale_color_manual(values=trial_scale)

pdf("Images/depth_plastic.pdf",height=8,width=7)
depth_plastic
dev.off()

############################################
## relationships with crown depth

neighbor_acq_plastic_CD<-ggplot(self_pruning,
                                aes(x=neighbor_acq,
                                    y=CrownDepth,
                                    color=Species))+
  geom_point(alpha=0.5)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=20))+
  labs(x="Neighbor acquisitiveness",
       y="Crown depth (cm)")+
  guides(color="none")+
  scale_color_manual(values=trial_scale)

NCI_plastic_CD<-ggplot(self_pruning,
                       aes(x=neighbor_comp,
                           y=CrownDepth,
                           color=Species))+
  geom_point(alpha=0.5)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=20),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())+
  labs(x="NCI",
       y="Crown depth (cm)")+
  guides(color="none")+
  scale_color_manual(values=trial_scale)

height_plastic_CD<-ggplot(self_pruning,
                          aes(x=HeightTop,
                              y=CrownDepth,
                              color=Species))+
  geom_point(alpha=0.5)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=20))+
  labs(x="Tree height (cm)",
       y="Crown depth (cm)")+
  scale_color_manual(values=trial_scale)

light_top_plastic_CD<-ggplot(self_pruning,
                             aes(x=logLightTop,
                                 y=CrownDepth,
                                 color=Species))+
  geom_point(alpha=0.5)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=20),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())+
  labs(x=expression(italic(L[top])),
       y="Crown depth (cm)")+
  guides(color="none")+
  scale_color_manual(values=trial_scale)

# pdf("Images/Fig4.pdf",height=11,width=10)
# (neighbor_acq_plastic_CD + NCI_plastic_CD)/
#   (height_plastic_CD + light_top_plastic_CD) +
#   plot_layout(guides = "collect") & theme(legend.position = "bottom")
# dev.off()

neighbor_acq_plastic_CB<-ggplot(self_pruning,
                                aes(x=neighbor_acq,
                                    y=HeightBase,
                                    color=Species))+
  geom_point(alpha=0.5)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=20),
        legend.position = "bottom")+
  labs(x="Neighbor acquisitiveness",
       y="Crown base height (cm)")+
  scale_color_manual(values=trial_scale)

# pdf("Images/FigS6.pdf",height=8,width=7)
# neighbor_acq_plastic_CB
# dev.off()

## radius ~ depth allometry
radius_depth<-ggplot(self_pruning,
                     aes(y=CR_average,
                         x=CrownDepth,
                         color=Species))+
  geom_point(alpha=0.5)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(text=element_text(size=20),
        legend.position = "bottom")+
  labs(y="Average crown radius (cm)",
       x="Crown depth (cm)")+
  scale_color_manual(values=trial_scale)

# pdf("Images/FigS_radius_depth.pdf",height=8,width=7)
# radius_depth
# dev.off()

################################################
## pull out the species-specific slopes from mixed-effects models

sp_match<-match(species_means$Species,levels(self_pruning$Species))

## no comparisons; REML = T by default
light_neighbor_acq_sp<-lmer(logLightBase~neighbor_acq*Species+(1|unique_plot),data=self_pruning)
light_neighbor_acq_trend<-emtrends(light_neighbor_acq_sp,"Species",var="neighbor_acq")
species_means$light_neighbor_acq_slope<-summary(light_neighbor_acq_trend)$neighbor_acq.trend[sp_match]
anova(light_neighbor_acq_sp, type="III")
## very similar results
# Anova(light_neighbor_acq_sp, type="III",test.statistic = "F")

light_NCI_sp<-lmer(logLightBase~neighbor_comp*Species+(1|unique_plot),data=self_pruning)
light_NCI_trend<-emtrends(light_NCI_sp,"Species",var="neighbor_comp")
species_means$light_NCI_slope<-summary(light_NCI_trend)$neighbor_comp.trend[sp_match]
anova(light_NCI_sp, type="III")

light_height_sp<-lmer(logLightBase~HeightTop*Species+(1|unique_plot),data=self_pruning)
light_height_trend<-emtrends(light_height_sp,"Species",var="HeightTop")
species_means$light_height_slope<-summary(light_height_trend)$HeightTop.trend[sp_match]
anova(light_height_sp, type="III")

## which species may show some evidence of correlative inhibition?
light_toplight_sp<-lmer(logLightBase~logLightTop*Species+(1|unique_plot),data=self_pruning)
light_toplight_trend<-emtrends(light_toplight_sp,"Species",var="logLightTop")
species_means$light_toplight_slope<-summary(light_toplight_trend)$logLightTop.trend[sp_match]
anova(light_toplight_sp, type="III")

leaf_habit_cols<-c("Deciduous"="chocolate4",
                   "Evergreen"="#60941a")

neighbor_acq_slopes<-ggplot(data=species_means,
                            aes(x=focal_acq,
                                y=light_neighbor_acq_slope,
                                label=Species))+
  geom_smooth(method="lm")+
  geom_point(size=1.5)+
  geom_text(size=5,aes(color=leaf_habit))+
  scale_color_manual(values = leaf_habit_cols)+
  theme_bw()+
  theme(text=element_text(size=20))+
  coord_cartesian(xlim = c(-2.6,2.4))+
  guides(color="none")+
  labs(x="Focal acquisitiveness",
       y=expression(paste(italic(L[base]),
                          " ~ neighbor acquisitiveness")),
       title="Species slopes")
# ggsave(filename = "Images/neighbor_acq_slopes.png",neighbor_acq_slopes,
#        dpi=600,width = 6,height=5)

NCI_slopes<-ggplot(data=species_means,
                   aes(x=focal_acq,
                       y=light_NCI_slope*100,
                       label=Species))+
  geom_smooth(method="lm")+
  geom_point(size=1.5)+
  geom_text(size=5,aes(color=leaf_habit))+
  theme_bw()+
  theme(text=element_text(size=20))+
  coord_cartesian(xlim = c(-2.6,2.4))+
  guides(color="none")+
  scale_color_manual(values = leaf_habit_cols)+
  labs(x="Focal acquisitiveness",
       y=expression(paste(italic(L[base]),
                          " ~ NCI (" %*% "100)")))
# labs(x="Tendance acquisitive (CP1)",
#      y=expression(paste("Pentes : ",italic(L[base])," ~ NCI")))
# ggsave(filename = "Images/NCI_slopes_FR.png",NCI_slopes,
#        dpi=600,width = 6,height=5)

height_slopes<-ggplot(data=species_means,
                      aes(x=focal_acq,
                          y=light_height_slope*1000,
                          label=Species))+
  geom_smooth(method="lm")+
  geom_point(size=1.5)+
  geom_text(size=5,aes(color=leaf_habit))+
  scale_color_manual(values = leaf_habit_cols)+
  theme_bw()+
  theme(text=element_text(size=20))+
  coord_cartesian(xlim = c(-2.6,2.4))+
  guides(color="none")+
  labs(x="Focal acquisitiveness",
       y=expression(paste(italic(L[base]),
                          " ~ top height (" %*% "1000)")))

light_top_slopes<-ggplot(data=species_means,
                         aes(x=focal_acq,
                             y=light_toplight_slope,
                             label=Species))+
  geom_smooth(method="lm")+
  geom_point(size=1.5)+
  geom_text(size=5,aes(color=leaf_habit))+
  scale_color_manual(values = leaf_habit_cols)+
  theme_bw()+
  theme(text=element_text(size=20))+
  coord_cartesian(xlim = c(-2.6,2.4))+
  labs(x="Focal acquisitiveness",
       y=expression(paste(italic(L[base]),
                          " ~ ",italic(L[top]))),
       color="Leaf habit")+
  guides(colour = guide_legend(position = "bottom"))

# pdf("Images/Fig3.pdf",height=22,width=11)
# (neighbor_acq_plastic + neighbor_acq_slopes)/
#   (NCI_plastic + NCI_slopes)/
#   (height_plastic + height_slopes)/
#   (light_top_plastic + light_top_slopes) +
#   plot_layout(guides = "collect") & theme(legend.position = "bottom")
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

