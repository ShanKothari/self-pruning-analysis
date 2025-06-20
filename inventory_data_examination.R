## PIGL C_PIGL_E5 high Height?

setwd("C:/Users/querc/Dropbox/PostdocProjects/SelfPruning/")

Inventory2016<-read.csv("IDENTMontrealData/Inventory2016.csv")
Inventory2017<-read.csv("IDENTMontrealData/Inventory2017.csv")
Inventory2018<-read.csv("IDENTMontrealData/Inventory2018.csv")
Inventory2019<-read.csv("IDENTMontrealData/Inventory2019.csv")
Inventory2020<-read.csv("IDENTMontrealData/Inventory2020.csv")

## to check that rows are in the same order
Inventory2016$UniqueTreeID<-with(Inventory2016,paste(Block,Plot,Pos,sep="_"))
Inventory2017$UniqueTreeID<-with(Inventory2017,paste(Block,Plot,Pos,sep="_"))
Inventory2018$UniqueTreeID<-with(Inventory2018,paste(Block,Plot,Pos,sep="_"))
Inventory2019$UniqueTreeID<-with(Inventory2019,paste(Block,Plot,Pos,sep="_"))
Inventory2020$UniqueTreeID<-with(Inventory2020,paste(Block,Plot,Pos,sep="_"))

## add basal diams from 2017 and 2019 to 2018 data
Inventory2018$BasalDiam_2016<-Inventory2016$BasalDiam
Inventory2018$BasalDiam_2017<-Inventory2017$BasalDiam
Inventory2018$BasalDiam_2018<-Inventory2018$BasalDiam
Inventory2018$BasalDiam_2019<-Inventory2019$BasalDiam
Inventory2018$BasalDiam_2020<-Inventory2020$BasalDiam
Inventory2018$BasalDiam<-NULL

## ABBA
plot(BasalDiam_2018~BasalDiam_2017,data=Inventory2018[Inventory2018$CodeSp=="ABBA",])
with(Inventory2018,which(CodeSp=="ABBA" & BasalDiam_2017>60 & BasalDiam_2018<40))
with(Inventory2018,which(CodeSp=="ABBA" & BasalDiam_2018>100))
with(Inventory2018,which(CodeSp=="ABBA" & BasalDiam_2018>50 & BasalDiam_2017<25))

plot(BasalDiam_2019~BasalDiam_2018,data=Inventory2018[Inventory2018$CodeSp=="ABBA",])
with(Inventory2018,which(CodeSp=="ABBA" & BasalDiam_2018<35 & BasalDiam_2019>45))

plot(Height~BasalDiam,data=Inventory2018[Inventory2018$CodeSp=="ABBA",])

## ACRU
Inventory2018_ACRU<-Inventory2018[Inventory2018$CodeSp=="ACRU",]

plot(BasalDiam_2018~BasalDiam_2017,data=Inventory2018_ACRU)
with(Inventory2018,which(CodeSp=="ACRU" & BasalDiam_2017<40 & BasalDiam_2018>60))

tmp<-lm(BasalDiam_2018~BasalDiam_2017,data=Inventory2018_ACRU)
Inventory2018[names(tmp$residuals[which(abs(tmp$residuals)>15)]),
              c("BasalDiam_2017","BasalDiam_2018","BasalDiam_2019","Height","UniqueTreeID")]

tmp<-lm(BasalDiam_2019~BasalDiam,data=Inventory2018_ACRU)
Inventory2018[names(tmp$residuals[which(abs(tmp$residuals)>15)]),
              c("BasalDiam_2017","BasalDiam_2018","BasalDiam_2019","Height","UniqueTreeID")]

tmp<-lm(BasalDiam_2018~Height,data=Inventory2018_ACRU)
Inventory2018[names(tmp$residuals[which(abs(tmp$residuals)>15)]),
              c("BasalDiam_2017","BasalDiam_2018","BasalDiam_2019","Height","UniqueTreeID")]

## ACSA
Inventory2018_ACSA<-Inventory2018[Inventory2018$CodeSp=="ACSA",]

with(Inventory2018_ACSA,plot(BasalDiam_2018~BasalDiam_2017))
#with(Inventory2018,which(CodeSp=="ACSA" & BasalDiam_2017<37 & BasalDiam_2018>50))

tmp<-lm(BasalDiam_2018~BasalDiam_2017,data=Inventory2018_ACSA)
Inventory2018[names(tmp$residuals[which(abs(tmp$residuals)>12)]),
              c("BasalDiam_2017","BasalDiam_2018","BasalDiam_2019","Height","UniqueTreeID")]

tmp<-lm(BasalDiam_2019~BasalDiam,data=Inventory2018_ACSA)
Inventory2018[names(tmp$residuals[which(abs(tmp$residuals)>12)]),
              c("BasalDiam_2017","BasalDiam_2018","BasalDiam_2019","Height","UniqueTreeID")]

plot(Height~BasalDiam,data=Inventory2018_ACSA)

## BEAL
Inventory2018_BEAL<-Inventory2018[Inventory2018$CodeSp=="BEAL",]
plot(BasalDiam_2018~BasalDiam_2017,data=Inventory2018_BEAL)

tmp<-lm(BasalDiam_2018~BasalDiam_2017,data=Inventory2018_BEAL)
Inventory2018[names(tmp$residuals[which(abs(tmp$residuals)>14)]),
              c("BasalDiam_2017","BasalDiam_2018","BasalDiam_2019","Height","UniqueTreeID")]

tmp<-lm(BasalDiam_2019~BasalDiam,data=Inventory2018_BEAL)
Inventory2018[names(tmp$residuals[which(abs(tmp$residuals)>12)]),
              c("BasalDiam_2017","BasalDiam_2018","BasalDiam_2019","Height","UniqueTreeID")]

## BEPA
Inventory2018_BEPA<-Inventory2018[Inventory2018$CodeSp=="BEPA",]
plot(BasalDiam_2018~BasalDiam_2017,data=Inventory2018_BEPA)

tmp<-lm(BasalDiam_2018~BasalDiam_2017,data=Inventory2018_BEPA)
Inventory2018[names(tmp$residuals[which(abs(tmp$residuals)>15)]),
              c("BasalDiam_2017","BasalDiam_2018","BasalDiam_2019","Height","UniqueTreeID")]

tmp<-lm(BasalDiam_2019~BasalDiam,data=Inventory2018_BEPA)
Inventory2018[names(tmp$residuals[which(abs(tmp$residuals)>15)]),
              c("BasalDiam_2017","BasalDiam_2018","BasalDiam_2019","Height","UniqueTreeID")]

tmp<-lm(Height~BasalDiam,data=Inventory2018_BEPA)
Inventory2018[names(tmp$residuals[which(abs(tmp$residuals)>15)]),
              c("BasalDiam_2017","BasalDiam_2018","BasalDiam_2019","Height","UniqueTreeID")]

## LALA
Inventory2018_LALA<-Inventory2018[Inventory2018$CodeSp=="LALA",]
plot(BasalDiam_2018~BasalDiam_2017,data=Inventory2018_LALA)

tmp<-lm(BasalDiam_2018~BasalDiam_2017,data=Inventory2018_LALA)
Inventory2018[names(tmp$residuals[which(abs(tmp$residuals)>12)]),
              c("BasalDiam_2017","BasalDiam_2018","BasalDiam_2019","Height","UniqueTreeID")]

tmp<-lm(BasalDiam_2019~BasalDiam,data=Inventory2018_LALA)
Inventory2018[names(tmp$residuals[which(abs(tmp$residuals)>12)]),
              c("BasalDiam_2017","BasalDiam_2018","BasalDiam_2019","Height","UniqueTreeID")]

## PIGL
Inventory2018_PIGL<-Inventory2018[Inventory2018$CodeSp=="PIGL",]
plot(BasalDiam_2018~BasalDiam_2017,data=Inventory2018_PIGL)

tmp<-lm(BasalDiam_2018~BasalDiam_2017,data=Inventory2018_PIGL)
Inventory2018[names(tmp$residuals[which(abs(tmp$residuals)>10)]),
              c("BasalDiam_2017","BasalDiam_2018","BasalDiam_2019","Height","UniqueTreeID")]

tmp<-lm(BasalDiam_2019~BasalDiam,data=Inventory2018_PIGL)
Inventory2018[names(tmp$residuals[which(abs(tmp$residuals)>10)]),
              c("BasalDiam_2017","BasalDiam_2018","BasalDiam_2019","Height","UniqueTreeID")]


## PIRE
Inventory2018_PIRE<-Inventory2018[Inventory2018$CodeSp=="PIRE",]
plot(BasalDiam_2018~BasalDiam_2017,data=Inventory2018_PIRE)

tmp<-lm(BasalDiam_2018~BasalDiam_2017,data=Inventory2018_PIRE)
Inventory2018[names(tmp$residuals[which(abs(tmp$residuals)>10)]),
              c("BasalDiam_2017","BasalDiam_2018","BasalDiam_2019","Height","UniqueTreeID")]

## PIRU
Inventory2018_PIRU<-Inventory2018[Inventory2018$CodeSp=="PIRU",]
plot(BasalDiam_2018~BasalDiam_2017,data=Inventory2018_PIRU)

tmp<-lm(BasalDiam_2018~BasalDiam_2017,data=Inventory2018_PIRU)
Inventory2018[names(tmp$residuals[which(abs(tmp$residuals)>10)]),
              c("BasalDiam_2017","BasalDiam_2018","BasalDiam_2019","Height","UniqueTreeID")]

## PIST
Inventory2018_PIST<-Inventory2018[Inventory2018$CodeSp=="PIST",]
plot(BasalDiam_2018~BasalDiam_2017,data=Inventory2018_PIST)

tmp<-lm(BasalDiam_2018~BasalDiam_2017,data=Inventory2018_PIST)
Inventory2018[names(tmp$residuals[which(abs(tmp$residuals)>12)]),
              c("BasalDiam_2017","BasalDiam_2018","BasalDiam_2019","Height","UniqueTreeID")]

## QURU
Inventory2018_QURU<-Inventory2018[Inventory2018$CodeSp=="QURU",]
plot(BasalDiam_2018~BasalDiam_2017,data=Inventory2018_QURU)

tmp<-lm(BasalDiam_2018~BasalDiam_2017,data=Inventory2018_QURU)
Inventory2018[names(tmp$residuals[which(abs(tmp$residuals)>15)]),
              c("BasalDiam_2017","BasalDiam_2018","BasalDiam_2019","Height","UniqueTreeID")]

## THOC
Inventory2018_THOC<-Inventory2018[Inventory2018$CodeSp=="THOC",]
plot(BasalDiam_2018~BasalDiam_2017,data=Inventory2018_THOC)

tmp<-lm(BasalDiam_2018~BasalDiam_2017,data=Inventory2018_THOC)
Inventory2018[names(tmp$residuals[which(abs(tmp$residuals)>13)]),
              c("BasalDiam_2017","BasalDiam_2018","BasalDiam_2019","Height","UniqueTreeID")]

write.csv(Inventory2018,"IDENTMontrealData/combined_2017-2019.csv",row.names = F)

## issues related to mortality...
## here, which trees' diameter is inconsistent with status
## could also check for zombie trees
Inventory2018[which(Inventory2018$StateDesc!="Dead" & is.na(Inventory2018$BasalArea)),]
Inventory2018[which(Inventory2018$StateDesc=="Dead" & !is.na(Inventory2018$BasalArea)),]
