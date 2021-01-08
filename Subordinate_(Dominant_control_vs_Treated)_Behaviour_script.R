# Charli Davies
# Meerkat Preg behaviour analysis for subordinate dams when their dominant dam is either C vs F - over 3 pregancy stages (MP, LP, PP) -  all locations, burrow then forage

# last modified 6/1/21

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

# Create a datset for Subordinate dams when their dominant is either a control or Treated - all locations #### 
SCSF_All <- Allbehaviour %>% filter(Rank == "S")%>% droplevels()
View(SCSF_All)

SCSF_All$KPM.ID=as.factor(CF_All$KMP.ID)
SCSF_All$Litter.Code=as.factor(CF_All$Litter.Code)
SCSF_All$Treatment=as.factor(CF_All$Treatment.4way)
SCSF_All$Location=as.factor(CF_All$Location)
SCSF_All$AM.PM=as.factor(CF_All$AM.PM)
SCSF_All$Preg.Stage=as.factor(CF_All$Preg.stage)

# check for overall colinearity
scatterplotMatrix(~Age+Total.Monthly.Rainfall+Treatment+Preg.stage+AM.PM+Group.Size+Location, data=SCSF_All)

# Start with HIA.I all locations####

#check data is zero-inflated
z=summary(SCSF_All$Treatment[SCSF_All$HIA.I.All==0])
nz=summary(SCSF_All$Treatment[SCSF_All$HIA.I.All>0])
nz/z
dotchart(SCSF_All$HIA.I.All/SCSF_All$Focal.Length.Hour)
boxplot(SCSF_All$HIA.I.All/SCSF_All$Focal.Length.Hour~CF_All$Treatment)

# start models - first check which distribution is best
p1<-glmmadmb(HIA.I.All~Treatment*Location+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson" , data = SCSF_All)
n1<-glmmadmb(HIA.I.All~Treatment*Location+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data = SCSF_All)
b1<-glmmadmb(HIA.I.All~Treatment*Location+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1" , data = SCSF_All)
AICtab(p1,n1,b1)
Anova(b1)
# none will run - too many variables for model - remove interaction
p2<-glmmadmb(HIA.I.All~Treatment+Location+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson" , data = SCSF_All)
n2<-glmmadmb(HIA.I.All~Treatment+Location+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data = SCSF_All)
b2<-glmmadmb(HIA.I.All~Treatment+Location+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data = SCSF_All )
AICtab(p2,n2)
#p2 is best

Anova (p2)
p3<-glmmadmb(HIA.I.All~Treatment+Location+AM.PM+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = SCSF_All )
Anova(p3)
p4<-glmmadmb(HIA.I.All~Treatment+Location+AM.PM+Age+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = SCSF_All )
Anova(p4)
p5<-glmmadmb(HIA.I.All~Treatment+Location+AM.PM+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = SCSF_All )
Anova(p5)
p6<-glmmadmb(HIA.I.All~Treatment+AM.PM+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson" , data = SCSF_All)
Anova(p6)
p7<-glmmadmb(HIA.I.All~Treatment+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = SCSF_All )
Anova(p7)
AICtab(p2,p3,p4,p5,p6,p7,n2)


# best model is p5 = S.HIA.I.All####
# confirm model by adding dropped terms back in
S.HIA.I.All<-glmmadmb(HIA.I.All~Treatment+Location+AM.PM+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = SCSF_All )
Anova(S.HIA.I.All)
summary(S.HIA.I.All)
coef(summary(S.HIA.I.All))
write.csv((coef(summary(S.HIA.I.All))), "S.HIA.I.All.csv")

# Start with HIA.R all locations####

#check data is zero-inflated
z=summary(SCSF_All$Treatment[SCSF_All$HIA.R.All==0])
nz=summary(SCSF_All$Treatment[SCSF_All$HIA.R.All>0])
nz/z
dotchart(SCSF_All$HIA.R.All/SCSF_All$Focal.Length.Hour)
boxplot(SCSF_All$HIA.R.All/SCSF_All$Focal.Length.Hour~CF_All$Treatment)
# one big outlier in SF treatments - tried models without - did not make any impact on final result 

# start models - first check which distribution is best
p1<-glmmadmb(HIA.R.All~Treatment*Location+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = SCSF_All )
n1<-glmmadmb(HIA.R.All~Treatment*Location+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data = SCSF_All)
b1<-glmmadmb(HIA.R.All~Treatment*Location+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1" , data = SCSF_All)
AICtab(p1,n1)
# only p1 and n1 will run - remove interaction - also low variance explained by Litter.Code

p2<-glmmadmb(HIA.R.All~Treatment+Location+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = SCSF_All )
n2<-glmmadmb(HIA.R.All~Treatment+Location+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data = SCSF_All)
b2<-glmmadmb(HIA.R.All~Treatment+Location+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1" , data = SCSF_All)
AICtab(p1,n1, p2, n2)
# n2 best - also remove KMPID

p2<-glmmadmb(HIA.R.All~Treatment+Location+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = SCSF_All )
n2<-glmmadmb(HIA.R.All~Treatment+Location+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = SCSF_All )
b2<-glmmadmb(HIA.R.All~Treatment+Location+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1" , data = SCSF_All)
AICtab(p1,n1,n2,b2,p2)
#n2 is best

Anova(n2)
n3<-glmmadmb(HIA.R.All~Treatment+Location+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|Preg.stage)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = SCSF_All )
Anova(n3)
n4<-glmmadmb(HIA.R.All~Treatment+Location+AM.PM+Group.Size+(1|Preg.stage)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = SCSF_All )
Anova(n4)
n5<-glmmadmb(HIA.R.All~Treatment+Location+Group.Size+(1|Preg.stage)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = SCSF_All  )
Anova(n5)
n6<-glmmadmb(HIA.R.All~Location+Group.Size+(1|Preg.stage)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = SCSF_All  )
AICtab(n1,n2,n3,n4,n5,n6,p2,b2)
Anova(n6)
summary(n6)

# best model is n6 = S.HIA.R.All ####
S.HIA.R.All<-glmmadmb(HIA.R.All~Location+Group.Size+(1|Preg.stage)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data = SCSF_All )
Anova(S.HIA.R.All)
summary(S.HIA.R.All)
coef(summary(S.HIA.R.All))
write.csv((coef(summary(S.HIA.R.All))), "S.HIA.R.All.csv")


################################################## Forage focals only ####

# seperate out focals taken while foraging only

# Create a datset for C vs F - forage only ####
SCSF_Forage <- SCSF_All %>% filter(Location != "Burrow")
View(SCSF_Forage)
str(SCSF_Forage)

# check for overall colinearity
scatterplotMatrix(~Age+Total.Monthly.Rainfall+Treatment+Preg.stage+AM.PM+Group.Size, data=SCSF_Forage)

# Start with Fcomp.I ####

#check data is zero-inflated
z=summary(SCSF_Forage$Treatment[SCSF_Forage$Fcomp.I==0])
nz=summary(SCSF_Forage$Treatment[SCSF_Forage$Fcomp.I>0])
nz/z
dotchart(SCSF_Forage$Fcomp.I/SCSF_Forage$Focal.Length.Hour)
# very few non-zeros
boxplot(SCSF_Forage$Fcomp.I/SCSF_Forage$Focal.Length.Hour~SCSF_Forage$Treatment)

# start models - first check which distribution is best
p1<-glmmadmb(Fcomp.I~Treatment+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = SCSF_Forage)
n1<-glmmadmb(Fcomp.I~Treatment+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = SCSF_Forage )
b1<-glmmadmb(Fcomp.I~Treatment+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1" , data = SCSF_Forage)
AICtab(p1,n1,b1)
# p1 is best

Anova(p1)
summary(p1)
p2<-glmmadmb(Fcomp.I~Treatment+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = SCSF_Forage )
Anova(p2)
p3<-glmmadmb(Fcomp.I~Treatment+Group.Size+AM.PM+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = SCSF_Forage )
Anova(p3)
p4<-glmmadmb(Fcomp.I~Treatment+Group.Size+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = SCSF_Forage )
Anova(p4)
p5<-glmmadmb(Fcomp.I~Treatment+AM.PM+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = SCSF_Forage )
Anova(p6)
p6<-glmmadmb(Fcomp.I~Treatment+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = SCSF_Forage )
AICtab(p1,p2,p3,p4,p5,p6)
Anova(p5)

# best model is p5 = S.Fcomp.I####
S.Fcomp.I<-glmmadmb(Fcomp.I~Treatment+AM.PM+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = SCSF_Forage )
Anova(S.Fcomp.I)
summary(S.Fcomp.I)
coef(summary(S.Fcomp.I))
write.csv((coef(summary(S.Fcomp.I))), "S.Fcomp.I.csv")

# Start with Fcomp.R ####

#check data is zero-inflated
z=summary(SCSF_Forage$Treatment[SCSF_Forage$Fcomp.R==0])
nz=summary(SCSF_Forage$Treatment[SCSF_Forage$Fcomp.R>0])
nz/z
dotchart(SCSF_Forage$Fcomp.R/SCSF_Forage$Focal.Length.Hour)
boxplot(SCSF_Forage$Fcomp.R/SCSF_Forage$Focal.Length.Hour~SCSF_Forage$Treatment)

# start models - first check which distribution is best
p1<-glmmadmb(Fcomp.R~Treatment+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = SCSF_Forage )
n1<-glmmadmb(Fcomp.R~Treatment+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = SCSF_Forage )
b1<-glmmadmb(Fcomp.R~Treatment+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data = SCSF_Forage )
#only p1 works

Anova(p1)
p2<-glmmadmb(Fcomp.R~Treatment+AM.PMf+Group.Size+Total.Monthly.Rainfall+(1|Preg.stage)+(1|FLitter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = SCSF_Forage )
Anova(p2)
p3<-glmmadmb(Fcomp.R~Treatment+AM.PMf+Group.Size+(1|Preg.stage)+(1|FLitter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = SCSF_Forage )
Anova(p3)
p4<-glmmadmb(Fcomp.R~Treatment+Group.Size+(1|Preg.stage)+(1|FLitter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson" , data = SCSF_Forage)
Anova(p4)
p5<-glmmadmb(Fcomp.R~Group.Size+(1|Preg.stage)+(1|FLitter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = SCSF_Forage )
Anova(p5)
AICtab(p1,p2,p3,p4,p5)

# best model is p5 = S.Fcomp.R####
S.Fcomp.R<-glmmadmb(Fcomp.R~Group.Size+(1|Preg.stage)+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson" , data = SCSF_Forage)
Anova(S.Fcomp.R)
summary(S.Fcomp.R)
coef(summary(S.Fcomp.R))
write.csv((coef(summary(S.Fcomp.R))), "S.Fcomp.R.csv")
