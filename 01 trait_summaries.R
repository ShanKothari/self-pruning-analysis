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

## linear regression interpolation of missing data
## via core set of five traits
## add interpolated values to summary dataset
traits_comp<-complete(mice(traits[,c("LDMC","N","SRL","SLA","WD")],
                           method = "norm.predict",maxit=1))
traits$SRL[traits$SpeciesCode=="PIOM"]<-traits_comp$SRL[traits$SpeciesCode=="PIOM"]
traits$LDMC[traits$SpeciesCode=="PIOM"]<-traits_comp$LDMC[traits$SpeciesCode=="PIOM"]

## z-standardize the core set of five traits
traits$LDMC<-(traits$LDMC-mean(traits$LDMC))/sd(traits$LDMC)
traits$N<-(traits$N-mean(traits$N))/sd(traits$N)
traits$SRL<-(traits$SRL-mean(traits$SRL))/sd(traits$SRL)
traits$SLA<-(traits$SLA-mean(traits$SLA))/sd(traits$SLA)
traits$WD<-(traits$WD-mean(traits$WD))/sd(traits$WD)

## just extract whether species are deciduous or evergreen
## for plotting purposes (but not for the PCA)
leaf_habit<-ifelse(traits$LL>12,
                   yes = "Evergreen",
                   no = "Deciduous")

## nativeness for plotting purposes
## add below (so we can de-emphasize non-native species)
non_native<-ifelse(traits$SpeciesCode %in% c("ACPL","LADE","PIOM","PIRU",
                                             "PISY","QURO","TICO"),
                   yes="non_native",no="native")

## subset traits (no LL or shade tolerance) in preparation for PCA
traits_sub<-traits[,-which(colnames(traits) %in% c("SpeciesCode","LL","shade_tol"))]
rownames(traits_sub)<-traits$SpeciesCode

## do PCA of traits, extract scores and loadings
## no need to use scaling in the function
## because traits are already scaled
trait_pca<-prcomp(traits_sub)
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
  geom_text(size = 4,label = trait_pca_scores$species,
            aes(color = trait_pca_scores$leaf_habit)) +
  geom_segment(data = trait_pca_loadings,
               aes(x = 0, y = 0, xend = -PC1*5, yend = PC2*5),
               arrow = arrow(length = unit(1/2, "picas")),
               linewidth = 1,color = "#3366FF") +
  annotate("text",
           x = -trait_pca_loadings$PC1*5.5,
           y = trait_pca_loadings$PC2*5.5,
           label = trait_pca_loadings$variables,
           color="#3366FF")+
  theme_bw()+
  theme(text=element_text(size=15),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        legend.position = c(0.1,0.17),
        legend.background = element_rect(fill = NA))+
  coord_fixed(ratio=trait_pca_perc[2]/trait_pca_perc[1])+
  scale_color_manual(values = leaf_habit_cols)+
  labs(x=paste("PC1 (",round(trait_pca_perc[1],1),"% variance)",sep=""),
       y=paste("PC2 (",round(trait_pca_perc[2],1),"% variance)",sep=""),
       color="Leaf habit")

pdf("Images/FigS1.pdf",width=7,height=4)
trait_pca_plot
dev.off()

trait_summary<-cbind(traits,trait_pca_scores)
trait_summary$species<-NULL

write.csv(trait_summary,"TraitData/trait_summary.csv",row.names = F)
