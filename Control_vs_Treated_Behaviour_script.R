# Charli Davies
# Meerkat Preg behaviour analysis for C vs F - over 2 pregancy stages (LP, PP) -  all locations, burrow then forage

# last modified last modified 6/1/21

# clear R of all objects
rm(list=ls())

# Load libraries####
library("car")
library("bbmle")
library("glmmADMB")
library("lsmeans")
library("multcomp")
library("predictmeans")
library("ggplot2")
library("broom")
library("tidyverse")

# First set directory as downloaded zipped folder ####

# Load data####

library(readr)
Allbehaviour <- read_csv("Behaviour_dataset.csv", col_types = cols(AM.PM = col_factor(levels = c("AM",
                                                                        "PM")), Age = col_number(),  Fcomp.I = col_number(),
                                          Fcomp.R = col_number(), Focal.Length.Hour = col_number(),
                                          Group.Size = col_number(), HIA.I.All = col_number(),
                                          HIA.R.All = col_number(), Location = col_factor(levels = c("Forage",
                                                                                                     "Burrow")), Olfactory = col_number(),
                                          Preg.stage = col_factor(levels = c("MP",
                                                                             "LP", "PP")), Prosocial.I = col_number(),
                                          Prosocial.R = col_number(), Rank = col_factor(levels = c("D",
                                                                                                   "S")), Sub.I.All = col_number(),
                                          Sub.R.All = col_number(), Total.Monthly.Rainfall = col_number(),
                                          Treatment.4way = col_factor(levels = c("Control",
                                                                                 "Flutamide", "Subordinate", "SubordinateFlut"))),
                         na = "NA")

View(Allbehaviour)
str(Allbehaviour)
head(Allbehaviour)
print.data.frame(Allbehaviour)

################################################# All locations ####

# Create a datset for Control Dominant versus Treated Dominant- all locations #### 
CF_All <- Allbehaviour %>% filter(Rank == "D")%>% droplevels()%>% filter(Preg.stage =! "MP")%>% droplevels()
View(CF_All)

CF_All$KPM.ID=as.factor(CF_All$KMP.ID)
CF_All$Litter.Code=as.factor(CF_All$Litter.Code)
CF_All$Treatment=as.factor(CF_All$Treatment.4way)
CF_All$Location=as.factor(CF_All$Location)
CF_All$AM.PM=as.factor(CF_All$AM.PM)
CF_All$Rank=as.factor(CF_All$Rank)

# Start with HIA.I all locations####
#check data is zero-inflated
z=summary(CF_All$Treatment[CF_All$HIA.I.All==0])
nz=summary(CF_All$Treatment[CF_All$HIA.I.All>0])
nz/z
dotchart(CF_All$HIA.I.All/CF_All$Focal.Length.Hour)
boxplot(CF_All$HIA.I.All/CF_All$Focal.Length.Hour~CF_All$Treatment)

# start models - first check which distribution is best
b2<-glmmadmb(HIA.I.All~Treatment+Location+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data = CF_All )
p2<-glmmadmb(HIA.I.All~Treatment+Location+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = CF_All )
n2<-glmmadmb(HIA.I.All~Treatment+Location+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data = CF_All).

# n1 and p1 will not run - too many variables...... will try without KMP.ID 
p2<-glmmadmb(HIA.I.All~Treatment+Location+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson" , data = CF_All)
n2<-glmmadmb(HIA.I.All~Treatment+Location+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = CF_All)
b2<-glmmadmb(HIA.I.All~Treatment+Location+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data = CF_All )
AICtab(b2,n2, p2, b1)            

# b1/b2 only model that works - makes no difference if litter is removed as a random effect
Anova(b1)
b3<-glmmadmb(HIA.I.All~Treatment+AM.PM+Location+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data = CF_All )
Anova(b3)
b4<-glmmadmb(HIA.I.All~Treatment+AM.PM+Location+Group.Size+Age+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1" , data = CF_All)
Anova(b4)
b5<-glmmadmb(HIA.I.All~Treatment+AM.PM+Location+Group.Size+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1" , data = CF_All)
Anova(b5)
b6<-glmmadmb(HIA.I.All~AM.PM+Location+Group.Size+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data = CF_All )
Anova(b6)
AICtab(b1,b2,b3,b4,b5, b6)
Anova(b6)

# b6 best model = Cf.HIA.I.All####
# confirm model by adding dropped terms back in
CF.HIA.I.All<-glmmadmb(HIA.I.All~AM.PM+Location+Group.Size+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1" )
Anova(CF.HIA.I.All)
coef(summary(CF.HIA.I.All))
write.csv((coef(summary(CF.HIA.I.All))), "CF.HIA.I.All.csv")
summary(CF.HIA.I.All)

# Start with HIA.R all locations####

#check data is zero-inflated
z=summary(CF_All$Preg.stage[CF_All$HIA.R.All==0])
nz=summary(CF_All$Preg.stage[CF_All$HIA.R.All>0])
nz/z
dotchart(HIA.R.All)
dotchart(CF_All$HIA.R.All/CF_All$Focal.Length.Hour)
boxplot(CF_All$HIA.R.All/CF_All$Focal.Length.Hour~CF_All$Treatment)
# one big outlier in burrow control

# start models - first check which distribution is best

p1<-glmmadmb(HIA.R.All~Treatment*Location+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson" , data = CF_All)
n1<-glmmadmb(HIA.R.All~Treatment*Location+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = CF_All )
b1<-glmmadmb(HIA.R.All~Treatment*Location+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data = CF_All )
AICtab(p1,n1,b1)
Anova(p1)
summary(p1)

# only p1 will run - try without interaction? - also not much varaince explained by random factors
p2<-glmmadmb(HIA.R.All~Treatment+Locationf+AM.PMf+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|FKPM.ID/FLitter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson" , data = CF_All)
n2<-glmmadmb(HIA.R.All~Treatment+Locationf+AM.PMf+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|FKPM.ID/FLitter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = CF_All )
b2<-glmmadmb(HIA.R.All~Treatment+Locationf+AM.PMf+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|FKPM.ID/FLitter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1" , data = CF_All)
AICtab(p1,n1,b1)

# model only runs for p1- too many zeroes! - remove age
p3<-glmmadmb(HIA.R.All~Treatment+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson" , data = CF_All)
Anova(p3)
p4<-glmmadmb(HIA.R.All~Treatment+Group.Size+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = CF_All )
Anova(p4)
p5<-glmmadmb(HIA.R.All~Group.Size+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = CF_All )
Anova(p5)
AICtab(p1,p2,p3,p4,p5)
Anova(p4)
summary(p4)

# p4  is best model= CF.HIA.R.All ####
# confirm model by adding dropped terms back in CF.HIA.R.All<-glmmadmb(HIA.R.All~Treatment+Group.Size+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = CF_All )
Anova(CF.HIA.R.All)
coef(summary(CF.HIA.R.All))
write.csv((coef(summary(CF.HIA.R.All))), "CF.HIA.R.All.csv")
summary(CF.HIA.R.All)
ph <- lsmeans(CF.HIA.R.All, ~ Treatment)
summary(as.glht(pairs(ph), by = NULL))

# Start with Olfactory.All locations####

#check data is zero-inflated
z=summary(CF_All$Preg.stage[CF_All$Olfactory.All==0])
nz=summary(CF_All$Preg.stage[CF_All$Olfactory.All>0])
nz/z
dotchart(CF_All$Olfactory.All)
dotchart(CF_All$Olfactory.All/CF_All$Focal.Length.Hour)
boxplot(CF_All$Olfactory.All/CF_All$Focal.Length.Hour~CF_All$Treatment)
# ONE BIG OOUTLIER IN CONTROL

# start models - first check which distribution is best
p1<-glmmadmb(Olfactory.All~Treatment*Location+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson" , data = CF_All)
n1<-glmmadmb(Olfactory.All~Treatment*Location+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data = CF_All)
b1<-glmmadmb(Olfactory.All~Treatment*Location+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1" , data = CF_All)
AICtab(p1,b1, n1)
# try without interaction
p2<-glmmadmb(Olfactory.All~Treatment+Location+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = CF_All )
n2<-glmmadmb(Olfactory.All~Treatment+Location+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = CF_All )
b2<-glmmadmb(Olfactory.All~Treatment+Location+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data = CF_All )
AICtab(n1,n2,b1,b2,p1)
Anova(n2)
Anova(b2)
# n2 best so far
n3<-glmmadmb(Olfactory.All~Treatment+Location+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data = CF_All)
# doesn't work when age removed - try simplifying further
n4<-glmmadmb(Olfactory.All~Treatment+Location+AM.PM+Group.Size+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data = CF_All)
Anova(n4)
AICtab(n1,n2,n4,b1,b2)
summary(n4)

# n4  is best model = Cf.Olfactory.All####
# confirm model by adding dropped terms back in 
CF.Olfactory.All<-glmmadmb(Olfactory.All~Treatment+Location+AM.PM+Group.Size+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data = CF_All)
Anova(CF.Olfactory.All)
coef(summary(CF.Olfactory.All))
write.csv((coef(summary(CF.Olfactory.All))), "CF.Olfactory.All.csv")
summary(CF.Olfactory.All)


################################################## Forage focals only ####

# seperate out focals taken while foraging only

# Create a datset for C vs F - forage only ####
CF_Forage <- CF_All %>% filter(Location != "Burrow")
View(CF_Forage)
str(CF_Forage)

# check for overall colinearity
scatterplotMatrix(~Age+Total.Monthly.Rainfall+Treatment+Preg.stage+AM.PM+Group.Size, data=CF_Forage)

# Start with Fcomp.R ####

#check data is zero-inflated
z=summary(CF_Forage$Treatment[CF_Forage$Fcomp.R==0])
nz=summary(CF_Forage$Treatment[CF_Forage$Fcomp.R>0])
nz/z
dotchart(CF_Forage$Fcomp.R/CF_Forage$Focal.Length.Hour)
boxplot(CF_Forage$Fcomp.R/CF_Forage$Focal.Length.Hour~CF_Forage$Treatment)

# start models - first check which distribution is best
p1<-glmmadmb(Fcomp.R~Treatment+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = CF_Forage)
n1<-glmmadmb(Fcomp.R~Treatment+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = CF_Forage )
b1<-glmmadmb(Fcomp.R~Treatment+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1" , data = CF_Forage)
AICtab(p1,n1, b1)
#b1 looks best 

Anova(b1)
b2<-glmmadmb(Fcomp.R~AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1" , data = CF_Forage)
#  doesnt work - try without KMP.ID as variance was low
p2<-glmmadmb(Fcomp.R~Treatment+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = CF_Forage )
n2<-glmmadmb(Fcomp.R~Treatment+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data = CF_Forage)
b2<-glmmadmb(Fcomp.R~Treatment+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1" , data = CF_Forage)
AICtab(p1,n1,b1, b2,n2,p2)
#n2 looks best 

Anova(n2)
n3<-glmmadmb(Fcomp.R~AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data = CF_Forage)
Anova(n3)
n4<-glmmadmb(Fcomp.R~AM.PM+Group.Size+Age+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data = CF_Forage)
Anova(n4)
n5<-glmmadmb(Fcomp.R~AM.PM+Group.Size+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = CF_Forage )
Anova(n5)
n6<-glmmadmb(Fcomp.R~Group.Size+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = CF_Forage )
Anova(n6)
AICtab(p1,n1,b1, b2,n2,p2, n3, n4, n5, n6)

# n6  is best model  = Cf.Fcomp.R ####
# confirm model by adding dropped terms back in 
CF.Fcomp.R<-glmmadmb(Fcomp.R~Group.Size+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = CF_Forage )
Anova(CF.Fcomp.R)
coef(summary(CF.Fcomp.R))
write.csv((coef(summary(CF.Fcomp.R))), "CF.Fcomp.R.csv")
summary(CF.Fcomp.R)


# Start with Fcomp.I ####

#check data is zero-inflated
z=summary(CF_Forage$Treatment[CF_Forage$Fcomp.I==0])
nz=summary(CF_Forage$Treatment[CF_Forage$Fcomp.I>0])
nz/z
dotchart(CF_Forage$Fcomp.I/CF_Forage$Focal.Length.Hour)
boxplot(CF_Forage$Fcomp.I/CF_Forage$Focal.Length.Hour~CF_Forage$Treatment)

# start models - first check which distribution is best
p1<-glmmadmb(Fcomp.I~Treatment+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson" , data = CF_Forage )
n1<-glmmadmb(Fcomp.I~Treatment+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data = CF_Forage )
b1<-glmmadmb(Fcomp.I~Treatment+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1" , data = CF_Forage )
AICtab(p1,n1, b1)

#none work try reducing random 
p1<-glmmadmb(Fcomp.I~Treatment+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = CF_Forage  )
n1<-glmmadmb(Fcomp.I~Treatment+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = CF_Forage  )
b1<-glmmadmb(Fcomp.I~Treatment+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data = CF_Forage  )
AICtab(p1,n1)
# looks like n1 is best - b won't work

Anova(n1)
n2<-glmmadmb(Fcomp.I~Treatment+AM.PM+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = CF_Forage  )
Anova(n2)
n3<-glmmadmb(Fcomp.I~Treatment+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1FLitter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = CF_Forage  )
Anova(n3)
AICtab(n1,n2,n3,p1)

# n3  is best model  = Cf.Fcomp.I ####
# confirm model by adding dropped terms back in 
CF.Fcomp.I<-glmmadmb(Fcomp.I~Treatment+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data = CF_Forage)
Anova(CF.Fcomp.I)
coef(summary(CF.Fcomp.I))
write.csv((coef(summary(CF.Fcomp.I))), "CF.Fcomp.I.csv")
summary(CF.Fcomp.I)


################################################## Burrow focals only ####

# seperate out focals taken at the burrow only

# Create a datset for C vs F - burrow only ####
CF_Burrow <- CF_All %>% filter(Location != "Forage")
View(CF_Burrow)
str(CF_Burrow)

# check for overall colinearity
scatterplotMatrix(~Age+Total.Monthly.Rainfall+Rank+Preg.stage+AM.PM+Group.Size, data=CF_Burrow)

# Start with SUB.R all####

#check data is zero-inflated
z=summary(CF_Burrow$Treatment[CF_Burrow$Sub.R.All==0])
nz=summary(CF_Burrow$Treatment[CF_Burrow$Sub.R.All>0])
nz/z
dotchart(CF_Burrow$Sub.R.All/CF_Burrow$Focal.Length.Hour)
boxplot(CF_Burrow$Sub.R.All/CF_Burrow$Focal.Length.Hour~CF_Burrow$Treatment)

# start models - first check which distribution is best
p1<-glmmadmb(Sub.R.All~Treatment+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data=CF_Burrow )
n1<-glmmadmb(Sub.R.All~Treatment+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data=CF_Burrow )
b1<-glmmadmb(Sub.R.All~Treatment+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data=CF_Burrow )
AICtab(p1,n1,b1)
# b1 looks best

Anova(b1)
b2<-glmmadmb(Sub.R.All~Treatment+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data=CF_Burrow )
Anova(b2)
b3<-glmmadmb(Sub.R.All~Treatment+AM.PM+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data=CF_Burrow )
Anova(b3)
b4<-glmmadmb(Sub.R.All~Treatment+AM.PM+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1" , data=CF_Burrow)
Anova(b4)
AICtab(n1,p1,b1,b2,b3,b4)

# b4  is best model  - = Cf.Sub.R.Burrow ####
# confirm model by adding dropped terms back in 
CF.Sub.R.Burrow<-glmmadmb(Sub.R.All~Treatment+AM.PM+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data=CF_Burrow )
Anova(CF.Sub.R.Burrow)
coef(summary(CF.Sub.R.Burrow))
write.csv((coef(summary(CF.Sub.R.Burrow))), "CF.Sub.R.Burrow.csv")
summary(CF.Sub.R.Burrow)

# Start with Prosocial.I.All locations####

#check data is zero-inflated
z=summary(CF_Burrow$Treatment[CF_Burrow$Prosocial.I.All==0])
nz=summary(CF_Burrow$Treatment[CF_Burrow$Prosocial.I.All>0])
nz/z
dotchart(CF_Burrow$Prosocial.I.All/CF_Burrow$Focal.Length.Hour)
boxplot(CF_Burrow$Prosocial.I.All/CF_Burrow$Focal.Length.Hour~CF_Burrow$Treatment)

# start models - first check which distribution is best
p1<-glmmadmb(Prosocial.I.All~Treatment+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = CF_Burrow)
n1<-glmmadmb(Prosocial.I.All~Treatment+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = CF_Burrow )
b1<-glmmadmb(Prosocial.I.All~Treatment+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data = CF_Burrow )
AICtab(p1,n1,b1)
# b1 doesnt work - n1 looks best

Anova(n1)
n2<-glmmadmb(Prosocial.I.All~Treatment+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data = CF_Burrow)
Anova(n2)
n3<-glmmadmb(Prosocial.I.All~Treatment+Group.Size+Age+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = CF_Burrow )
Anova(n3)
n4<-glmmadmb(Prosocial.I.All~Treatment+Group.Size+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = CF_Burrow )
Anova(n4)
n5<-glmmadmb(Prosocial.I.All~Group.Size+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = CF_Burrow )
Anova(n5)
AICtab(p1,n1,n2,n3,n4,n5)

# n4  is best model  - = Cf.Prosocial.I.All####
# confirm model by adding dropped terms back in 
CF.Prosocial.I.Burrow<-glmmadmb(Prosocial.I.All~Treatment+Group.Size+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data = CF_Burrow)
Anova(CF.Prosocial.I.Burrow)
coef(summary(CF.Prosocial.I.Burrow))
write.csv((coef(summary(CF.Prosocial.I.Burrow))), "CF.Prosocial.I.Burrow.csv")
summary(CF.Prosocial.I.Burrow)

# Start with Prosocial.R.All locations####

#check data is zero-inflated
z=summary(CF_Burrow$Preg.stage[CF_Burrow$Prosocial.R.All==0])
nz=summary(CF_Burrow$Preg.stage[CF_Burrow$Prosocial.R.All>0])
nz/z
dotchart(Prosocial.R.All)
dotchart(CF_Burrow$Prosocial.R.All/CF_Burrow$Focal.Length.Hour)
boxplot(CF_Burrow$Prosocial.R.All/CF_Burrow$Focal.Length.Hour~CF_Burrow$Treatment)
# ONE BIG OUTLIER IN CONTROL

# start models - first check which distribution is best
p1<-glmmadmb(Prosocial.R.All~Treatment+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson" , data = CF_Burrow)
n1<-glmmadmb(Prosocial.R.All~Treatment+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data = CF_Burrow)
b1<-glmmadmb(Prosocial.R.All~Treatment+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data = CF_Burrow )
AICtab(p1,n1, b1)
# b1 is best

Anova(b1)
b2<-glmmadmb(Prosocial.R.All~Treatment+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data = CF_Burrow )
Anova(b2)
b3<-glmmadmb(Prosocial.R.All~Treatment+Group.Size+Age+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data = CF_Burrow )
Anova(b3)
b4<-glmmadmb(Prosocial.R.All~Treatment+Group.Size+(1|Preg.stage)+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data = CF_Burrow )
Anova(b4)
AICtab(b1,b2,b3,b4)
# b3 and b4 don't converge - try removing KMPID as a random factor
summary(b2)
#remove KMP.ID?
b3<-glmmadmb(Prosocial.R.All~Treatment+Group.Size+Age+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1" , data = CF_Burrow)
Anova(b3)
b4<-glmmadmb(Prosocial.R.All~Treatment+Group.Size+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data = CF_Burrow )
Anova(b4)
AICtab(b1,b2,b3,b4)

# b4  is best model = Cf.Prosocial.R.All ####
CF.Prosocial.R.Burrow<-glmmadmb(Prosocial.R.All~Treatment+Group.Size+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1" , data = CF_Burrow)
Anova(CF.Prosocial.R.Burrow)
coef(summary(CF.Prosocial.R.Burrow))
write.csv((coef(summary(CF.Prosocial.R.Burrow))), "CF.Prosocial.R.Burrow.csv")
summary(CF.Prosocial.R.Burrow)
