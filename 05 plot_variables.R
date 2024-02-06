setwd("C:/Users/querc/Dropbox/PostdocProjects/SelfPruning/")
source("../../DecomposingFD/R/AlphaFD.R")

library(stringr)
library(dplyr)
library(labdsv)
library(FD)
library(ggplot2)

Inventory2018<-read.csv("IDENTMontrealData/Inventory2018_cleaned.csv")

## extract row and column from position indicator
Inventory2018$Col<-str_sub(Inventory2018$Pos,1,1)
Inventory2018$Row<-as.numeric(str_sub(Inventory2018$Pos,2,2))

#################################
## calculate plot-level overyielding

## keep only plots with just native species
del_species<-c("ACPL","LADE","TICO","PIAB","PISY","QURO","PIOM","BELE","BEPE","PIMA")
del_plots<-unique(Inventory2018$Plot[Inventory2018$CodeSp %in% del_species])
del_plot_pattern<-paste(paste("_",del_plots,"_",sep=""),collapse="|")
Inventory2018_sub<-subset(Inventory2018,
                          !grepl(del_plot_pattern,unique_tree))

## we won't remove 12-species plots or blocks B and C
## although we could to correspond with the self-pruning survey
# Inventory2018_sub<-subset(Inventory2018_sub,PlotRichness!=12)
# Inventory2018_sub<-subset(Inventory2018,Block %in% c("A","D"))

## remove edge trees
edge.trees<-c(paste(LETTERS[1:8],"1",sep=""),
              paste(LETTERS[1:8],"8",sep=""),
              paste("A",1:8,sep=""),
              paste("H",1:8,sep=""))
edge_pattern<-paste(edge.trees,collapse="|")
Inventory2018_sub<-Inventory2018_sub[-which(grepl(edge_pattern,Inventory2018_sub$Pos)),]

## create a dummy variable with a value of 1 for all planted trees
Inventory2018_sub$num_planted<-1

## because we start overyielding calculations at the
## individual level, we want to consider trees to have
## 'underperformed' if they died. to this end, we create
## a new basal area column which replaces NAs (for dead
## trees) with 0s, to use in overyielding calculations
Inventory2018_sub$BasalArea_0<-Inventory2018_sub$BasalArea
Inventory2018_sub$BasalArea_0[which(is.na(Inventory2018_sub$BasalArea_0))]<-0

## sum basal area and number of planted individuals by plot x species
plot_sp_ba<-aggregate(cbind(BasalArea_0,num_planted)~Block+Plot+CodeSp+PlotRichness,
                      data=Inventory2018_sub,
                      FUN=sum,na.rm=T)

## calculate planted frequencies (36 individuals per plot)
plot_sp_ba$planted_freq<-plot_sp_ba$num_planted/36

## attach monoculture plot basal area from the same block and species
plot_sp_ba_mono<-subset(plot_sp_ba,PlotRichness==1)
plot_sp_ba$mono_ba<-apply(plot_sp_ba,1,
                          function(x) {
                            plot_sp_ba_mono$BasalArea_0[plot_sp_ba_mono$CodeSp==x["CodeSp"] & 
                                                          plot_sp_ba_mono$Block==x["Block"]]
                          })

## calculate monoculture expectations based on
## monoculture basal area and planted frequency in mixture
plot_sp_ba$mono_exp<-plot_sp_ba$mono_ba*plot_sp_ba$planted_freq
## convert from cm^2 per plot (9 m^2) to m^2 per hectare
plot_sp_ba$spNBE<-(plot_sp_ba$BasalArea_0-plot_sp_ba$mono_exp)/9
plot_sp_ba$BasalArea<-plot_sp_ba$BasalArea_0/9

plot_ba<-aggregate(cbind(BasalArea,spNBE)~Block+Plot+PlotRichness,
                   data=plot_sp_ba,
                   FUN=sum,na.rm=T)
plot_ba$spNBE[plot_ba$PlotRichness==1]<-NA

plot_ba$unique_plot<-paste(plot_ba$Block,
                           plot_ba$Plot,
                           sep="_")

#########################################
## output mean basal area in monoculture in 2018
## as well as mortalities at the species level
## (all plots and only monocultures)

## count a tree as alive if its status is not "dead"
Inventory2018_sub$num_alive<-ifelse(Inventory2018_sub$StateDesc!="Dead",yes=1,no=0)

## mortality by species (across all plots)
mortality_agg<-aggregate(cbind(num_alive,num_planted)~CodeSp,
                         data=Inventory2018_sub,
                         FUN=sum)
mortality_agg$mortality<-with(mortality_agg,1-(num_alive/num_planted))

## mortality by species (monoculture only)
Inventory2018_sub_mono<-subset(Inventory2018_sub,PlotRichness==1)
mortality_mono_agg<-aggregate(cbind(num_alive,num_planted)~CodeSp,
                              data=Inventory2018_sub_mono,
                              FUN=sum)
mortality_mono_agg$mortality<-with(mortality_mono_agg,1-(num_alive/num_planted))

## basal area by species (monoculture only)
BA_mono_agg<-aggregate(BasalArea~CodeSp,
                       data=Inventory2018_sub_mono,
                       FUN=mean,na.rm=T)

## mortality per species per unique plot (all plots)
mortality_2018<-aggregate(cbind(num_alive,num_planted)~Block+Plot+CodeSp,
                          data=Inventory2018_sub,
                          FUN=sum)

# write.csv(mortality_2018,"IDENTMontrealData/mortality_2018.csv",row.names=F)

#########################################
## calculate plot-level functional diversity
## and heterogeneity in shade-tolerance
## based on inner 6 x 6 trees

mortality_2018$unique_plot<-paste(mortality_2018$Block,
                                  mortality_2018$Plot,
                                  sep="_")
## weights are by numbers of planted individuals
planted_comm<-matrify(mortality_2018[,c("unique_plot","CodeSp","num_planted")])

## functional dispersion in traits
## using Laliberte's FDis and Scheiner's qDTM
## the latter is not richness-independent
trait_summary<-read.csv("TraitData/trait_summary.csv")
core_traits<-trait_summary[,c("LDMC","N","SLA","SRL","WD")]
rownames(core_traits)<-trait_summary$SpeciesCode
core_traits_sub<-core_traits[colnames(planted_comm),]

plot_fdis<-fdisp(dist(core_traits_sub),a = as.matrix(planted_comm))$FDis
plot.FTD<-FTD.comm(tdmat=dist(core_traits_sub),
                   spmat = as.matrix(planted_comm),
                   abund=T)

## and heterogeneity in shade tolerance
shade_tol<-setNames(trait_summary$shade_tol,
                    trait_summary$SpeciesCode)[colnames(planted_comm)]
plot_STH<-fdisp(dist(shade_tol),a = as.matrix(planted_comm))$FDis
plot.FTD_STH<-FTD.comm(tdmat=dist(shade_tol),
                       spmat = as.matrix(planted_comm),
                       abund=T)

## and heterogeneity in Lbase
species_means<-read.csv("SelfPruningData/species_means.csv")
Lbase<-setNames(species_means$logLightBase,
                species_means$Species)[colnames(planted_comm)]
plot_LH<-fdisp(dist(Lbase),a = as.matrix(planted_comm))$FDis
plot.FTD_LH<-FTD.comm(tdmat=dist(Lbase),
                      spmat = as.matrix(planted_comm),
                      abund=T)

plot_vars<-data.frame(unique_plot=rownames(planted_comm),
                      FDis=plot_fdis,
                      STH=plot_STH,
                      LH=plot_LH,
                      FTD=plot.FTD$com.FTD$qDTM,
                      FTD_STH=plot.FTD_STH$com.FTD$qDTM,
                      FTD_LH=plot.FTD_LH$com.FTD$qDTM,
                      richness=plot.FTD$com.FTD$nsp)

plot_vars$BasalArea<-plot_ba$BasalArea[match(plot_vars$unique_plot,plot_ba$unique_plot)]
plot_vars$NBE<-plot_ba$spNBE[match(plot_vars$unique_plot,plot_ba$unique_plot)]

## to try block random effects, but always singular
plot_vars$block<-unlist(lapply(strsplit(plot_vars$unique_plot,split = "_"),function(x) x[[1]]))

# write.csv(plot_vars,"IDENTMontrealData/plot_vars.csv",row.names=F)

BA_FTD<-ggplot(plot_vars,aes(x=log(FTD),y=BasalArea))+
  geom_point(size=2)+
  geom_smooth(method="lm")+
  theme_bw()+
  theme(text=element_text(size=20))+
  labs(x=expression("log("^q*"D(TM)) of functional traits"),
       y=expression("Basal area ( "*m^2*" "*ha^{-1}*")"))

BA_FTD_STH<-ggplot(plot_vars,aes(x=log(FTD_STH),y=BasalArea))+
  geom_point(size=2)+
  geom_smooth(method="lm")+
  theme_bw()+
  theme(text=element_text(size=20),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(x=expression("log("^q*"D(TM)) of shade tolerance"),
       y=expression("Basal area ( "*m^2*" "*ha^{-1}*")"))

BA_FTD_LH<-ggplot(plot_vars,aes(x=log(FTD_LH),y=BasalArea))+
  geom_point(size=2)+
  geom_smooth(method="lm")+
  theme_bw()+
  theme(text=element_text(size=20),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(x=expression("log("^q*"D(TM)) of "*italic("L"*""[base])),
       y=expression("Basal area ( "*m^2*" "*ha^{-1}*")"))
  # labs(x=expression("Diversity of "*italic("L"*""[base])),
  #      y=expression("Basal area ( "*m^2*" "*ha^{-1}*")"))
# ggsave(filename = "Images/BA_FTD_LF.png", plot = BA_FTD_LH,
#        width=6,height=5,units="in",dpi=600)

pdf("Images/Fig3.pdf",width=13,height=5)
BA_FTD+BA_FTD_STH+BA_FTD_LH
dev.off()

NBE_FTD<-ggplot(plot_vars,aes(x=log(FTD),y=NBE))+
  geom_point(size=2)+
  geom_smooth(method="lm")+
  theme_bw()+
  theme(text=element_text(size=20))+
  labs(x=expression("log("^q*"D(TM)) of functional traits"),
       y=expression("NBE ( "*m^2*" "*ha^{-1}*")"))

NBE_FTD_STH<-ggplot(plot_vars,aes(x=log(FTD_STH),y=NBE))+
  geom_point(size=2)+
  geom_smooth(method="lm")+
  theme_bw()+
  theme(text=element_text(size=20),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(x=expression("log("^q*"D(TM)) of shade tolerance"),
       y=expression("NBE ( "*m^2*" "*ha^{-1}*")"))

NBE_FTD_LH<-ggplot(plot_vars,aes(x=log(FTD_LH),y=NBE))+
  geom_point(size=2)+
  geom_smooth(method="lm")+
  theme_bw()+
  theme(text=element_text(size=20),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(x=expression("log("^q*"D(TM)) of "*italic("L"*""[base])),
       y=expression("NBE ( "*m^2*" "*ha^{-1}*")"))

pdf("Images/FigS2.pdf",width=13,height=5)
NBE_FTD+NBE_FTD_STH+NBE_FTD_LH
dev.off()
