# Charli Davies
# Meerkat Preg behaviour analysis for C vs S - over 3 pregancy stages (MP, LP, PP) -  all locations, burrow then forage

# last modified 29/05/2019

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
Allbehaviour <- read_csv("Behaviour_dataset.csv",
col_types = cols(AM.PM = col_factor(levels = c("AM",
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

# Create a datset for Control Dominant versus subordinate dams- all locations #### 
CS_All <- Allbehaviour %>% filter(DF.Treatment == "Control")%>% droplevels()
View(CS_All)

CS_All$Preg.stage = factor(CS_All$Preg.stage, levels(CS_All$Preg.stage)[c(1,2,3)])

CS_All$KPM.ID=as.factor(CS_All$KMP.ID)
CS_All$Litter.Code=as.factor(CS_All$Litter.Code)
CS_All$Treatment=as.factor(CS_All$Treatment.4way)
CS_All$Location=as.factor(CS_All$Location)
CS_All$AM.PM=as.factor(CS_All$AM.PM)
CS_All$Rank=as.factor(CS_All$Rank)

# check for overall colinearity
library("psych")
library(car)
scatterplotMatrix(~Age+Total.Monthly.Rainfall+Rank+Preg.stage+AM.PM+Location+Group.Size, data=CS_All)
#age and rank might be correlated
age_rank<-lm(Age~Rank, data=CS_All)
Anova(age_rank)
# age and rank are correlated

# Start with HIA.I all locations ####

#check data is zero-inflated
z=summary(CS_All$Rank[CS_All$HIA.I.All==0])
nz=summary(CS_All$Rank[CS_All$HIA.I.All>0])
nz/z
dotchart(CS_All$HIA.I.All/CS_All$Focal.Length.Hour)
boxplot(CS_All$HIA.I.All/CS_All$Focal.Length.Hour~CS_All$Rank)
boxplot(CS_All$HIA.I.All/CS_All$Focal.Length.Hour~CS_All$Preg.stage)

# start models - first check which distribution is best

I_HIA_p1<-glmmadmb(HIA.I.All~Rank*Preg.stage+Location+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson" , data = CS_All)
I_HIA_n1<-glmmadmb(HIA.I.All~Rank+Preg.stage+Location+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data=CS_All )
I_HIA_b1<-glmmadmb(HIA.I.All~Rank+Preg.stage+Location+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data=CS_All )
AICtab(I_HIA_p1,I_HIA_n1, I_HIA_b1)
Anova(I_HIA_b1)
# I_HIA_b1 has lowest AIC
I_HIA_b2<-glmmadmb(HIA.I.All~Rank+Preg.stage+Location+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1",data=CS_All )
Anova(I_HIA_b2)
I_HIA_b3<-glmmadmb(HIA.I.All~Rank+Preg.stage+Location+AM.PM+Group.Size+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1",data=CS_All )
Anova(I_HIA_b3)
I_HIA_b4<-glmmadmb(HIA.I.All~Rank+Location+AM.PM+Group.Size+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1",data=CS_All )
Anova(I_HIA_b4)
AICtab(I_HIA_b1,I_HIA_b2,I_HIA_b3,I_HIA_b4,I_HIA_b5,I_HIA_p1)

# Best model is b4 = CS.HIA.I.All ####
# confirm model by adding dropped terms back in 
CS.HIA.I.All<-glmmadmb(HIA.I.All~Rank+Location+AM.PM+Group.Size+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data=CS_All)
Anova(CS.HIA.I.All)
summary(CS.HIA.I.All)
coef(summary(CS.HIA.I.All))
write.csv((coef(summary(CS.HIA.I.All))), "CS.HIA.I.All.csv")
plot(cld(summary(glht(CS.HIA.I.All, linfct=mcp(Rank="Tukey")))))
plot(cld(summary(glht(CS.HIA.I.All, linfct=mcp(AM.PM="Tukey")))))
plot(cld(summary(glht(CS.HIA.I.All, linfct=mcp(Location="Tukey")))))



# Start with HIA.R all locations ####

#check data is zero-inflated
z=summary(CS_All$Rank[CS_All$HIA.R.All==0])
nz=summary(CS_All$Rank[CS_All$HIA.R.All>0])
nz/z
dotchart(CS_All$HIA.R.All/CS_All$Focal.Length.Hour)
boxplot(CS_All$HIA.R.All/CS_All$Focal.Length.Hour~CS_All$Rank)
boxplot(CS_All$HIA.R.All/CS_All$Focal.Length.Hour~CS_All$Preg.stage)

# start models - first check which distribution is best

p1<-glmmadmb(HIA.R.All~Rank*Preg.stage*Location+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data=CS_All )
n1<-glmmadmb(HIA.R.All~Rank*Preg.stage*Location+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data=CS_All )
b1<-glmmadmb(HIA.R.All~Rank*Preg.stage*Location+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data=CS_All )
AICtab(p1,n1,b1)
Anova(p1)
summary(p1)
# p1 looks best

p1<-glmmadmb(HIA.R.All~Rank*Preg.stage+Rank*Location+Preg.stage*Location+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson" , data=CS_All)
Anova(p1)
p2<-glmmadmb(HIA.R.All~Rank*Preg.stage+Rank*Location+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data=CS_All )
Anova(p2)
p3<-glmmadmb(HIA.R.All~Rank*Preg.stage+Rank*Location+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data=CS_All )
Anova(p3)
AICtab(p1,p2,p3)

# Best model is p3= CS.HIA.R.All ####
# confirm model by adding dropped terms back in 
CS.HIA.R.All<-glmmadmb(HIA.R.All~Rank*Preg.stage+Rank*Location+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data=CS_All )
Anova(CS.HIA.R.All)
coef(summary(CS.HIA.I.All))
write.csv((coef(summary(CS.HIA.R.All))), "CS.HIA.R.All.csv")
summary(CS.HIA.R.All)
ph <- lsmeans(CS.HIA.R.All, ~ Preg.stage|Rank)
summary(as.glht(pairs(ph), by = NULL))
ph <- lsmeans(CS.HIA.R.All, ~ Rank| Location)
summary(as.glht(pairs(ph), by = NULL))
plot(ph)
# looks like there is no difference between ranks during foraging but is at the burrow
# also looks like subordinate behaviour doesn't change across preg stage but control does


# Start with Olfactory.All locations ####

#check data is zero-inflated
z=summary(CS_All$Rank[CS_All$Olfactory.All==0])
nz=summary(CS_All$Rank[CS_All$Olfactory.All>0])
nz/z
dotchart(CS_All$Olfactory.All/CS_All$Focal.Length.Hour)
boxplot(CS_All$Olfactory.All/CS_All$Focal.Length.Hour~CS_All$Rank)
boxplot(CS_All$Olfactory.All/CS_All$Focal.Length.Hour~CS_All$Preg.stage)
# One big outlier in control

# start models - first check which distribution is best
p1<-glmmadmb(Olfactory.All~Rank*Preg.stage*Location+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson" , data= CS_All)
n1<-glmmadmb(Olfactory.All~Rank*Preg.stage*Location+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data= CS_All )
b1<-glmmadmb(Olfactory.All~Rank*Preg.stage*Location+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data= CS_All )
AICtab(p1,n1,b1)
# try without 3 way interaction
p2<-glmmadmb(Olfactory.All~Rank*Preg.stage+Rank*Location+Location*Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson" , data= CS_All)
n2<-glmmadmb(Olfactory.All~Rank*Preg.stage+Rank*Location+Location*Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data= CS_All )
b2<-glmmadmb(Olfactory.All~Rank*Preg.stage+Rank*Location+Location*Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data= CS_All )
#  try without interaction
p2<-glmmadmb(Olfactory.All~Rank+Location+Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson" , data= CS_All)
n2<-glmmadmb(Olfactory.All~Rank+Preg.stage+Location+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data= CS_All )
b2<-glmmadmb(Olfactory.All~Rank+Preg.stage+Location+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data= CS_All )
# only b2 works....

b3<-glmmadmb(Olfactory.All~Rank+Preg.stage+AM.PM+Location+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data= CS_All )
Anova(b3)
b4<-glmmadmb(Olfactory.All~Rank+Preg.stage+AM.PM+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data= CS_All )
Anova(b4)
b5<-glmmadmb(Olfactory.All~Rank+Preg.stage+AM.PM+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data= CS_All )
Anova(b5)
b6<-glmmadmb(Olfactory.All~Rank+Preg.stage+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data= CS_All )
Anova(b6)
b7<-glmmadmb(Olfactory.All~Rank+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data= CS_All )
Anova(b7)
AICtab(b1,b2,b3,b4,b5,b6,b7)

# Best model is b7 - CS.Olfactory.All ####
# confirm model by adding dropped terms back in 
CS.Olfactory.All<-glmmadmb(Olfactory.All~Rank+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data= CS_All )
Anova(CS.Olfactory.All)
summary(CS.Olfactory.All)
coef(summary(CS.Olfactory.All))
plot(cld(summary(glht(CS.Olfactory.All, linfct=mcp(Rank="Tukey")))))




################################################## Burrow focals only ####

# seperate out focals taken at the burrow only

# Create a datset for C vs S - burrow only ####
CS_Burrow <- CS_All %>% filter(Location != "Forage")
View(CS_Burrow)
str(CS_Burrow)

# check for overall colinearity
scatterplotMatrix(~Age+Total.Monthly.Rainfall+Rank+Preg.stage+AM.PM+Group.Size, data=CS_Burrow)
#age and rank might be correlated
age_rank<-lm(Age~Rank, data=CS_Burrow)
Anova(age_rank)
# age and rank are correlated

# Start with Prosocial.I.All - from the burrow only ####

#check data is zero-inflated
z=summary(CS_Burrow$Preg.stage[CS_Burrow$Prosocial.I.All==0])
nz=summary(CS_Burrow$Preg.stage[CS_Burrow$Prosocial.I.All>0])
nz/z
dotchart(CS_Burrow$Prosocial.I.All/CS_Burrow$Focal.Length.Hour)
boxplot(CS_Burrow$Prosocial.I.All/CS_Burrow$Focal.Length.Hour~CS_Burrow$Rank)
boxplot(CS_Burrow$Prosocial.I.All/CS_Burrow$Focal.Length.Hour~CS_Burrow$Preg.stage)

# start models - first check which distribution is best
p1<-glmmadmb(Prosocial.I.All~Rank*Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data=CS_Burrow )
n1<-glmmadmb(Prosocial.I.All~Rank*Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data=CS_Burrow )
b1<-glmmadmb(Prosocial.I.All~Rank*Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data=CS_Burrow )
AICtab(p1,n1,b1)
Anova(n1)
# n1 is the best

n2<-glmmadmb(Prosocial.I.All~Rank+Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data=CS_Burrow )
Anova(n2)
n3<-glmmadmb(Prosocial.I.All~Rank+Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data=CS_Burrow  )
Anova(n3)
n4<-glmmadmb(Prosocial.I.All~Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data=CS_Burrow  )
Anova(n4)
n5<-glmmadmb(Prosocial.I.All~Preg.stage+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data=CS_Burrow  )
Anova(n5)
n6<-glmmadmb(Prosocial.I.All~Preg.stage+Group.Size+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data=CS_Burrow  )
Anova(n6)
n7<-glmmadmb(Prosocial.I.All~Group.Size+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data=CS_Burrow )
Anova(n7)
AICtab(n1,n2,n3,n4,n5,n6,n7)

# Best model is n7 = CS.Prosocial.I.Burrow ####
# confirm model by adding dropped terms back in 
CS.Prosocial.I.Burrow<-glmmadmb(Prosocial.I.All~Group.Size+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data=CS_Burrow  )
Anova(CS.Prosocial.I.Burrow)
summary(CS.Prosocial.I.Burrow)
write.csv((coef(summary(CS.Prosocial.I.Burrow))), "CS.Prosocial.I.Burrow.csv")

# Start with Prosocial.R.All - from the burrow only ####

z=summary(CS_Burrow$Preg.stage[CS_Burrow$Prosocial.R.All==0])
nz=summary(CS_Burrow$Preg.stage[CS_Burrow$Prosocial.R.All>0])
nz/z
dotchart(CS_Burrow$Prosocial.R.All/CS_Burrow$Focal.Length.Hour)
boxplot(CS_Burrow$Prosocial.R.All/CS_Burrow$Focal.Length.Hour~CS_Burrow$Treatment)
boxplot(CS_Burrow$Prosocial.R.All/CS_Burrow$Focal.Length.Hour~CS_Burrow$Rank)
# one big outlier in control

# start models - first check which distribution is best
p1<-glmmadmb(Prosocial.R.All~Rank*Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = CS_Burrow)
n1<-glmmadmb(Prosocial.R.All~Rank*Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data = CS_Burrow)
b1<-glmmadmb(Prosocial.R.All~Rank*Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data = CS_Burrow )
AICtab(p1,n1,b1)
Anova(p1)
# p1 only model which will run

#try without interaction - as does not seem to be significant
p2<-glmmadmb(Prosocial.R.All~Rank+Preg.stage+AM.PMf+Group.Size+Total.Monthly.Rainfall+(1|FKPM.ID/FLitter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson" , data = CS_Burrow )
n2<-glmmadmb(Prosocial.R.All~Rank+Preg.stage+AM.PMf+Group.Size+Total.Monthly.Rainfall+(1|FKPM.ID/FLitter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = CS_Burrow  )
b2<-glmmadmb(Prosocial.R.All~Rank+Preg.stage+AM.PMf+Group.Size+Total.Monthly.Rainfall+(1|FKPM.ID/FLitter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1" , data = CS_Burrow )
AICtab(p2,n2,b2, p1)
# looks like n2 is best - b still doesnt run

Anova(n2)
n3<-glmmadmb(Prosocial.R.All~Rank+Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = CS_Burrow  )
Anova(n3)
n4<-glmmadmb(Prosocial.R.All~Rank+Preg.stage+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = CS_Burrow  )
Anova(n4)
n5<-glmmadmb(Prosocial.R.All~Rank+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = CS_Burrow  )
Anova(n5)
n6<-glmmadmb(Prosocial.R.All~Rank+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = CS_Burrow  )
Anova(n6)
n7<-glmmadmb(Prosocial.R.All~Rank+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data = CS_Burrow )
Anova(n7)
AICtab(n2,n3,n4,n5,n6,n7)

# Best model is n6 = CS.Prosocial.R.Burrow ####
# confirm model by adding dropped terms back in 
CS.Prosocial.R.Burrow<-glmmadmb(Prosocial.R.All~Rank+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = CS_Burrow  )
Anova(CS.Prosocial.R.Burrow)
summary(CS.Prosocial.R.Burrow)
write.csv((coef(summary(CS.Prosocial.R.Burrow))), "CS.Prosocial.R.Burrow.csv")




################################################## Forage focals only ####

# Create a datset for C vs S - forage only ####
CS_Forage <- CS_All %>% filter(Location != "Burrow")
View(CS_Forage)
str(CS_Forage)

# check for overall colinearity
scatterplotMatrix(~Age+Total.Monthly.Rainfall+Rank+Preg.stage+AM.PM+Group.Size, data=CS_Forage)
#age and rank might be correlated
age_rank<-lm(Age~Rank, data=CS_Forage)
Anova(age_rank)
# age and rank are correlated

# Start with Fcomp.R.All - from foraging only ####

#check data is zero-inflated
z=summary(CS_Forage$Preg.stage[CS_Forage$Fcomp.R==0])
nz=summary(CS_Forage$Preg.stage[CS_Forage$Fcomp.R>0])
nz/z
dotchart(CS_Forage$Fcomp.R/CS_Forage$Focal.Length.Hour)
boxplot(CS_Forage$Fcomp.R/CS_Forage$Focal.Length.Hour~CS_Forage$Treatment)
boxplot(CS_Forage$Fcomp.R/CS_Forage$Focal.Length.Hour~CS_Forage$Preg.stage)

# start models - first check which distribution is best
p1<-glmmadmb(Fcomp.R~Rank*Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = CS_Forage)
n1<-glmmadmb(Fcomp.R~rank*Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = CS_Forage )
b1<-glmmadmb(Fcomp.R~Rank*Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data = CS_Forage )
Anova(p1)
# only p1 works - try without interaction 
p2<-glmmadmb(Fcomp.R~Rank+Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = CS_Forage )
n1<-glmmadmb(Fcomp.R~Rank+Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = CS_Forage )
b1<-glmmadmb(Fcomp.R~Rank+Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data = CS_Forage )
AICtab(p1,n1,b1,p2)
# looks like p2 is best of ones that run

AICtab(p1,p2)
Anova(p2)
summary(p2)
p3<-glmmadmb(Fcomp.R~Rank+Preg.stage+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson" , data = CS_Forage)
Anova(p3)
# doesnt work well try removing next factor as well 
p4<-glmmadmb(Fcomp.R~Rank+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = CS_Forage )
Anova(p4)
p5<-glmmadmb(Fcomp.R~Rank+Group.Size+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = CS_Forage )
Anova(p5)
p6<-glmmadmb(Fcomp.R~Group.Size+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = CS_Forage )
# does not work with just group size - looking at summary looks like need to remove KMPID as random as random factor explains very little variance
p8<-glmmadmb(Fcomp.R~Group.Size+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson" , data = CS_Forage)
Anova(p8)
AICtab(p1,p2,p3,p4,p5,p6,p7,p8)

# Best model is p8 = CS.Fcomp.R.Forage ####
# confirm model by adding dropped terms back in 
CS.Fcomp.R.Forage<-glmmadmb(Fcomp.R~Group.Size+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = CS_Forage  )
Anova(CS.Fcomp.R.Forage)
summary(CS.Fcomp.R.Forage)
write.csv((coef(summary(CS.Fcomp.R.Forage))), "CS.Fcomp.R.Forage.csv")

# Start with Fcomp.I.All - from foraging only ####

#check data is zero-inflated
z=summary(CS_Forage$Rank[CS_Forage$Fcomp.I==0])
nz=summary(CS_Forage$Rank[CS_Forage$Fcomp.I>0])
nz/z
dotchart(CS_Forage$Fcomp.I/CS_Forage$Focal.Length.Hour)
boxplot(CS_Forage$Fcomp.I/CS_Forage$Focal.Length.Hour~CS_Forage$Rank)
boxplot(CS_Forage$Fcomp.I/CS_Forage$Focal.Length.Hour~CS_Forage$Preg.stage)

# start models - first check which distribution is best
p1<-glmmadmb(Fcomp.I~Rank*Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = CS_Forage)
n1<-glmmadmb(Fcomp.I~Rank*Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = CS_Forage )
b1<-glmmadmb(Fcomp.I~Rank*Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1" , data = CS_Forage)
AICtab(p1,n1,b1)
# none work- lets try without interaction
p1<-glmmadmb(Fcomp.I~Rank+Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = CS_Forage )
n1<-glmmadmb(Fcomp.I~Rank+Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data = CS_Forage)
b1<-glmmadmb(Fcomp.I~Rank+Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1" , data = CS_Forage)
AICtab(p1,n1,b1)
# still none work - try removing KMPID as a random factor
p1<-glmmadmb(Fcomp.I~Rank+Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = CS_Forage )
n1<-glmmadmb(Fcomp.I~Rank+Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data = CS_Forage)
b1<-glmmadmb(Fcomp.I~Rank+Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data = CS_Forage )
AICtab(n1,p1)
# hmm so NOTHING works - TRY WITHOUT KMP ID
p1<-glmmadmb(Fcomp.I~Rank+Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall++offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = CS_Forage )
n1<-glmmadmb(Fcomp.I~Rank+Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = CS_Forage )
b1<-glmmadmb(Fcomp.I~Rank+Preg.stage+AM.PM+Group.Size+Total.Monthly.Rainfall++offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1" , data = CS_Forage)
AICtab(n1,p1)
#n1 looks best 
Anova(n1)
n2<-glmmadmb(Fcomp.I~Rank+Preg.stage+AM.PM+Group.Size+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data = CS_Forage)
Anova(n2)
n3b<-glmmadmb(Fcomp.I~Rank+Preg.stage+Group.Size+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data = CS_Forage)
n3<-glmmadmb(Fcomp.I~Rank+Preg.stage+Group.Size+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data = CS_Forage)
Anova(n3)
n4<-glmmadmb(Fcomp.I~Rank+Group.Size+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data = CS_Forage)
Anova(n4)
AICtab(n1,p1, n2, n3, n4)
# have been able to add litter ID back in as a random factor

# Best model is n4 = CS.Fcomp.I.Forage ####
# confirm model by adding dropped terms back in 
CS.Fcomp.I.Forage<-glmmadmb(Fcomp.I~Rank+Group.Size+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data = CS_Forage)
Anova(CS.Fcomp.I.Forage)
summary(CS.Fcomp.I.Forage)
write.csv((coef(summary(CS.Fcomp.I.Forage))), "CS.Fcomp.I.Forage.csv")




################################################## Only dominant females ####
# Create a dataset only including dominant dams while foraging ####
C_Forage <- CS_Forage %>% filter(Rank != "S")
View(C_Forage)
str(C_Forage)

# check for overall colinearity
scatterplotMatrix(~Age+Total.Monthly.Rainfall+Preg.stage+AM.PM+Group.Size, data=C_Forage)

# Start with Fcomp.I - dominant dam only ####

#check data is zero-inflated
z=summary(C_Forage$Preg.stage[C_Forage$Fcomp.I==0])
nz=summary(C_Forage$Preg.stage[C_Forage$Fcomp.I>0])
nz/z
dotchart(C_Forage$Fcomp.I/C_Forage$Focal.length.Hour)
boxplot(C_Forage$Fcomp.I/C_Forage$Focal.Length.Hour~C_Forage$Preg.stage)
# quite alot of zeros

# start models - first check which distribution is best
p1<-glmmadmb(Fcomp.I~Preg.stage+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson" , data = C_Forage)
n1<-glmmadmb(Fcomp.I~Preg.stage+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = C_Forage )
b1<-glmmadmb(Fcomp.I~Preg.stage+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1" , data = C_Forage)
AICtab(p1,n1,b1)
# none work - try removing KMPID
p1<-glmmadmb(Fcomp.I~Preg.stage+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = C_Forage )
n1<-glmmadmb(Fcomp.I~Preg.stage+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data = C_Forage)
b1<-glmmadmb(Fcomp.I~Preg.stage+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1" , data = C_Forage)
AICtab(p1,n1,b1)
Anova(n1)
# n1 looks best

n2<-glmmadmb(Fcomp.I~Preg.stage+AM.PM+Age+Total.Monthly.Rainfall+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = C_Forage )
Anova(n2)
n3<-glmmadmb(Fcomp.I~Preg.stage+Age+Total.Monthly.Rainfall+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = C_Forage )
Anova(n3)
n4<-glmmadmb(Fcomp.I~Preg.stage+Age+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data = C_Forage)
# model will not work again..... maybe try with KMP.ID instead of Litter code
n5<-glmmadmb(Fcomp.I~Preg.stage+Age+(1|KPM.ID)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data = C_Forage)
# model now works including both random factors....
Anova(n5)
n6<-glmmadmb(Fcomp.I~Preg.stage+Age+(1|FKPM.ID/FLitter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" )
AICtab(n1,n2,n3,n4, n5, n6)
Anova(n6)

# Best model is n6 = C.Fcomp.I.Forage #### 
# confirm model by adding dropped terms back in 
C.Fcomp.I.Forage<-glmmadmb(Fcomp.I~Preg.stage+Age+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom" , data = C_Forage)
Anova(C.Fcomp.I.Forage)
summary(C.Fcomp.I.Forage)
ph <- lsmeans(C.Fcomp.I.Forage, ~ Preg.stage)
summary(as.glht(pairs(ph), by = NULL))
plot(cld(summary(glht(C.Fcomp.I.Forage, linfct=mcp(Preg.stage="Tukey")))))
summary(glht(C.Fcomp.I.Forage, linfct=mcp(Preg.stage="Tukey")))
write.csv((coef(summary(C.Fcomp.I.Forage))), "C.Fcomp.I.Forage.csv")




# Create a dataset only including dominant dams at the burrow ####
C_Burrow <- CS_Burrow %>% filter(Rank != "S")
View(C_Burrow)
str(C_Burrow)

# check for overall colinearity
scatterplotMatrix(~Age+Total.Monthly.Rainfall+Preg.stage+AM.PM+Group.Size, data=C_Burrow)

# Start with Prosocial.R.All ####

#check data is zero-inflated
z=summary(C_Burrow$Preg.stage[C_Burrow$Prosocial.I.All==0])
nz=summary(C_Burrow$Preg.stage[C_Burrow$Prosocial.I.All>0])
nz/z
dotchart(C_Burrow$Prosocial.I.All/C_Burrow$Focal.Length.Hour)
boxplot(C_Burrow$Prosocial.I.All/C_Burrow$Focal.Length.Hour~C_Burrow$Preg.stage)

# start models - first check which distribution is best
p1<-glmmadmb(Prosocial.R.All~Preg.stage+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="poisson", data = C_Burrow )
n1<-glmmadmb(Prosocial.R.All~Preg.stage+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = C_Burrow  )
b1<-glmmadmb(Prosocial.R.All~Preg.stage+AM.PM+Group.Size+Age+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom1", data = C_Burrow  )
AICtab(p1,n1,b1)
Anova(n1)
# n1 is best

n2<-glmmadmb(Prosocial.R.All~Preg.stage+Group.Size+Age+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = C_Burrow  )
Anova(n2)
n3<-glmmadmb(Prosocial.R.All~Preg.stage+Group.Size+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = C_Burrow  )
Anova(n3)
n4<-glmmadmb(Prosocial.R.All~Preg.stage+Total.Monthly.Rainfall+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = C_Burrow  )
Anova(n4)
n5<-glmmadmb(Prosocial.R.All~Preg.stage+(1|KPM.ID/Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = C_Burrow  )
AICtab(n1,n2,n3,n4,n5,p1)

# Best model is n5 - C.Prosocial.R.burrow ####
# confirm model by adding dropped terms back in 
C.Prosocial.R.burrow<-glmmadmb(Prosocial.R.All~Preg.stage+(1|Litter.Code)+offset(log(Focal.Length.Hour)),zeroInflation=TRUE, family="nbinom", data = C_Burrow  )
Anova(C.Prosocial.R.burrow)
summary(C.Prosocial.R.burrow)
summary(glht(C.Prosocial.R.burrow, linfct=mcp(Preg.stage="Tukey")))
coef(summary(C.Prosocial.R.burrow))
write.csv((coef(summary(C.Prosocial.R.burrow))), "C.Prosocial.R.burrow.csv")
