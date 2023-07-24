setwd("C:/Users/querc/Dropbox/PostdocProjects/SelfPruning/")

library(ggplot2)
library(mice)

## read in trait data and subset by relevant species and traits
traits<-read.csv("TraitData/IDENT_TRAIT_DATABASE_2020-10-20.csv")
del_species<-c("BELE","BEPE","PIMA")
traits<-traits[-which(traits$SpeciesCode %in% del_species),]
traits<-traits[,c("SpeciesCode","LL","LDMC","Leaf_N_mass","SLA..all.include.",
                  "SRL..fine.root.","SSD...WD","Shade.tolerance")]
colnames(traits)<-c("SpeciesCode","LL","LDMC","N","SLA","SRL","WD","shade_tol")

## more TRY data (numbers are dataset IDs)
## 227: Pierce S., Brusa G., Vagge I., Cerabolini B.E.L. (2013) Functional Ecology, 27(4): 1002-1010
## 72: Cornelissen, J. H. C., B. Cerabolini, et al. 2003. Journal of Vegetation Science 14:311-322.
## 281: Hattermann D, et al https://gepris.dfg.de/gepris/projekt/260851423?language=en&selectedSubTab=2
traits$LDMC[traits$SpeciesCode=="TICO"]<-344.177

## just extract whether species are deciduous or evergreen
## for plotting purposes (but not for the PCA)
leaf_habit<-ifelse(traits$LL>12,
                   yes = "Evergreen",
                   no = "Deciduous")

## nativeness for plotting purposes
## add below (so we can de-emphasize non-native species)

## subset traits (no LL or shade tolerance) in preparation for PCA
traits_sub<-traits[,-which(colnames(traits) %in% c("SpeciesCode","LL","shade_tol"))]
rownames(traits_sub)<-traits$SpeciesCode

## linear regression interpolation of missing data
## add interpolated values to summary dataset
traits_sub_comp<-complete(mice(traits_sub,method = "norm.predict",maxit=1))
traits$SRL[traits$SpeciesCode=="PIOM"]<-traits_sub_comp$SRL[rownames(traits_sub_comp)=="PIOM"]
traits$LDMC[traits$SpeciesCode=="PIOM"]<-traits_sub_comp$LDMC[rownames(traits_sub_comp)=="PIOM"]

## do scaled PCA of traits, extract scores and loadings
trait_pca<-prcomp(traits_sub_comp,scale. = T)
trait_pca_scores<-data.frame(species = rownames(trait_pca$x),
                             trait_pca$x,
                             leaf_habit=leaf_habit)
trait_pca_loadings<-data.frame(variables = rownames(trait_pca$rotation),
                               trait_pca$rotation)
trait_pca_loadings$variables<-c("LDMC","%N","SLA","SRL","WD")
trait_pca_perc<-trait_pca$sdev^2/sum(trait_pca$sdev^2)*100

leaf_habit_cols<-c("Deciduous"="#a44f30",
                   "Evergreen"="#60941a")

trait_pca_plot<-ggplot(trait_pca_scores, 
                       aes(x = -PC1, y = PC2)) +
  geom_text(size = 3.5,label = trait_pca_scores$species,
            aes(color = trait_pca_scores$leaf_habit)) +
  geom_segment(data = trait_pca_loadings,
               aes(x = 0, y = 0, xend = -PC1*5, yend = PC2*5),
               arrow = arrow(length = unit(1/2, "picas")),
               size=1,color="#3366FF") +
  annotate("text",
           x = -trait_pca_loadings$PC1*5.5,
           y = trait_pca_loadings$PC2*5.5,
           label = trait_pca_loadings$variables,
           color="#3366FF")+
  theme_bw()+theme(text=element_text(size=15),
                   panel.background = element_rect(fill='transparent'), #transparent panel bg
                   plot.background = element_rect(fill='transparent', color=NA))+
  coord_fixed(ratio=trait_pca_perc[2]/trait_pca_perc[1])+
  guides(color="none")+
  scale_color_manual(values = leaf_habit_cols)+
  labs(x=paste("CP1 (",round(trait_pca_perc[1],1),"% variance)",sep=""),
       y=paste("CP2 (",round(trait_pca_perc[2],1),"% variance)",sep=""))

ggsave("Images/trait_pca_plot.png", trait_pca_plot, bg='transparent',
       dpi=600,width=7,height=4)

trait_summary<-cbind(traits,trait_pca_scores)
trait_summary$species<-NULL

write.csv(trait_summary,"TraitData/trait_summary.csv",row.names = F)
