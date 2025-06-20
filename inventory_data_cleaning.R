setwd("C:/Users/querc/Dropbox/PostdocProjects/SelfPruning/")

#######################
## cleaning script for inventory data

## read in data
Inventory2016<-read.csv("IDENTMontrealData/Inventory2016.csv")
Inventory2017<-read.csv("IDENTMontrealData/Inventory2017.csv")
Inventory2018<-read.csv("IDENTMontrealData/Inventory2018.csv")
Inventory2019<-read.csv("IDENTMontrealData/Inventory2019.csv")
Inventory2020<-read.csv("IDENTMontrealData/Inventory2020.csv")

Inventory2018$unique_tree<-with(Inventory2018,paste(Block,Plot,Pos,sep="_"))

## add basal diams from other years to 2018 data
Inventory2018$BasalDiam_2016<-Inventory2016$BasalDiam
Inventory2018$BasalDiam_2017<-Inventory2017$BasalDiam
Inventory2018$BasalDiam_2019<-Inventory2019$BasalDiam
Inventory2018$BasalDiam_2020<-Inventory2020$BasalDiam

## duplicate column for cleaning
Inventory2018$BasalDiam_cleaned<-Inventory2018$BasalDiam

## note to self: we don't take out unused plots, or blocks B and C,
## because some trees in these blocks could be in the radii of focal trees.
## however, this script does not implement cleaning of those plots
## because no annotation was done about what sort of replacement
## to use. I also haven't checked the data sheets.

## Jon's modeled values
## filtered to only include year 2018
## I don't actually use this, since I confirm that
## results are similar to the geometric mean method
## and because it wasn't calculated for plot edges
model_JU<-read.csv("IDENTMontrealData/DB_Corrected_Shan.csv")
model_JU<-model_JU[model_JU$Inv==10,]
model_JU$unique_tree<-with(model_JU,paste(Block,Plot,Pos,sep="_"))

## my own data cleaning notes
## in ToDo column, "Mean" means replace with interpolated value
## "Replace" means replace with specific value
## I input "replace" when I felt like, upon checking the sheets,
## I had a good reason to believe that a specific mistake
## was made in data entry or recording
cleaning_notes<-read.csv("IDENTMontrealData/IDENT_MTL_data_cleaning_sorted.csv")
cleaning_notes$ToDo<-as.factor(cleaning_notes$ToDo)
cleaning_notes<-cleaning_notes[which(cleaning_notes$ToDo!=""),]

for(i in 1:nrow(cleaning_notes)){
  print(i)
  inventory_row<-match(cleaning_notes$UniqueTreeID[i],
                       Inventory2018$unique_tree)
  JU_row<-match(cleaning_notes$UniqueTreeID[i],
                model_JU$unique_tree)
  
  if(cleaning_notes$ToDo[i]=="Replace"){
    Inventory2018$BasalDiam_cleaned[inventory_row]<-cleaning_notes$Replacement[i]
  }
  
  if(cleaning_notes$ToDo[i]=="Mean"){
    ## geometric mean of two years on either side
    Inventory2018$BasalDiam_cleaned[inventory_row]<-exp(mean(log(c(Inventory2018$BasalDiam_2017[inventory_row],
                                                                    Inventory2018$BasalDiam_2019[inventory_row])),
                                                              na.rm=T))
  }
}

## these trees are missing basal diameters so I classify them as dead
## it's possible they are not actually dead (e.g. some have 2019 measurements)
## just for consistency and accurate tallying of mortality
Inventory2018$StateDesc[Inventory2018$unique_tree=="A_2NR7B_H1"]<-"Dead"
Inventory2018$StateDesc[Inventory2018$unique_tree=="A_4N8_G8"]<-"Dead"
Inventory2018$StateDesc[Inventory2018$unique_tree=="D_PIRE_H5"]<-"Dead"

#### these trees had basal diameters but were classified as dead originally
Inventory2018$BasalDiam_cleaned[Inventory2018$unique_tree=="A_4N7_G7"]<-NA
Inventory2018$BasalDiam_cleaned[Inventory2018$unique_tree=="A_PIST_D3"]<-NA
Inventory2018$BasalDiam_cleaned[Inventory2018$unique_tree=="A_PIRU_F8"]<-NA

## found in both 2016 and 2020 with similar basal diameter
Inventory2018$StateDesc[Inventory2018$unique_tree=="D_2N4_H2"]<-"Alive"

## basal diameter probably linked to the wrong tree
Inventory2018$BasalDiam_cleaned[Inventory2018$unique_tree=="D_2N5_E4"]<-43.12
Inventory2018$BasalDiam_cleaned[Inventory2018$unique_tree=="D_2N5_F4"]<-NA

Inventory2018$StateDesc[Inventory2018$unique_tree=="D_2NR5_E1"]<-"Alive"

## probably H6 measured twice
Inventory2018$BasalDiam_cleaned[Inventory2018$unique_tree=="D_4N5_G6"]<-NA

## probably switched with either A4 or A6
## see if I can confirm?
Inventory2018$BasalDiam_cleaned[Inventory2018$unique_tree=="D_PIRE_A4"]<-35.7
Inventory2018$StateDesc[Inventory2018$unique_tree=="D_PIRE_A4"]<-"Alive"
Inventory2018$BasalDiam_cleaned[Inventory2018$unique_tree=="D_PIRE_A5"]<-NA

Inventory2018$BasalDiam_cleaned[Inventory2018$unique_tree=="D_PIRE_G7"]<-NA

###############################
## height estimation

## calculate basal area in cm^2
Inventory2018$BasalArea<-(Inventory2018$BasalDiam_cleaned/20)^2*pi

## building models based only on cleaned data
exotic_species<-c("ACPL","LADE","PIAB","PIOM","PISY","QURO","TICO")
exotic_plots<-c(paste("M",1:8,sep=""),"MM6A","MM6B",exotic_species)
Inventory2018_sub<-subset(Inventory2018,Block %in% c("A","D"))
Inventory2018_sub<-subset(Inventory2018_sub,!(Plot %in% exotic_plots))

## removing outlier for model-building only
Inventory2018_sub$Height[Inventory2018_sub$unique_tree=="A_4NR7_G8"]<-NA

gen.ah.models<-function(df,sp.col,area.col,height.col){
  model.list<-list()
  df.sp.list<-split(df,df[[sp.col]])
  model.list<-lapply(df.sp.list,function(x){
    area<-x[,area.col]
    height<-x[,height.col]
    return(lm(log(height)~log(area)))
  })
  return(model.list)
}

height.models<-gen.ah.models(Inventory2018_sub,"CodeSp","BasalArea","Height")

Inventory2018$Height_est<-NA
for(i in 1:nrow(Inventory2018)){

  if(!is.na(Inventory2018$Height[i])){
    Inventory2018$Height_est[i]<-Inventory2018$Height[i]
  } else if(Inventory2018$CodeSp[i] %in% exotic_species){
    Inventory2018$Height_est[i]<-NA
  } else {
    sp<-as.character(Inventory2018$CodeSp[i])
    model.pred<-height.models[[match(sp,names(height.models))]]
    pred.d<-predict(model.pred,
                    newdata=data.frame(area=Inventory2018$BasalArea[i]))
    Inventory2018$Height_est[i]<-exp(pred.d)
  }
}

## calculate parabolic volume
Inventory2018$PV<-with(Inventory2018,1/2*BasalArea*Height_est)

write.csv(Inventory2018,"IDENTMontrealData/Inventory2018_cleaned.csv",row.names = F)
