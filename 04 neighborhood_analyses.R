setwd("C:/Users/querc/Dropbox/PostdocProjects/SelfPruning")

library(ggplot2)
library(ggpubr)
library(lme4)

self_pruning<-read.csv("SelfPruningData/self_pruning_processed.csv")

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