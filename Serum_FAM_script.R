
# Meerkat preg serum data + FAM - Control dominant vs subordinate dams
# 10/12/2019

# clear R of all objects
rm(list=ls())

# load libraries
library(tidyverse)
library(MASS)
library(nlme)
library(lsmeans)
library(car)
library(MuMIn)
library(multcomp)

# First set directory as downloaded zipped folder ####

# load serum dataset ####

# import and rename data file
library(readr)
Pregserum <- read_csv("Pregserum_dataset.csv", col_types = cols(A4.ng.ml = col_number(),
Age.months = col_number(), Preg.stage = col_factor(levels = c("EP","MP",
"LP", "PP")), Rank = col_factor(levels = c("D",
"S")), T.ng.ml = col_number(), Total.monthly.Rainfall = col_number()))

View(Pregserum)
str(Pregserum)

# Analysis for Androstenedione####

# create a datset for A4
df_ANDRO <- Pregserum %>% filter(A4.ng.ml != "na")%>% droplevels()
View(df_ANDRO)

# to convert variables to factors AND EXPLORE DATASET
df_ANDRO$Rank <- factor(df_ANDRO$Rank) 
df_ANDRO$ID <- factor(df_ANDRO$ID)
df_ANDRO$Preg.stage <- factor(df_ANDRO$Preg.stage)
df_ANDRO$LA4 <- log(df_ANDRO$A4.ng.ml+1)
df_ANDRO$Preg.stage= factor(df_ANDRO$Preg.stage, levels(df_ANDRO$Preg.stage) [c(1,2,3)])

#check for colinarity
vif(glmmPQL(A4.ng.ml~Rank*Preg.stage+Age.months+Total.monthly.Rainfall, random=  ~1|ID, family=Gamma(link='log'), data=df_ANDRO))

# RUN GLMM

A1<-glmmPQL(A4.ng.ml~Rank*Preg.stage+Age.months+Total.monthly.Rainfall, random=  ~1|ID, family=Gamma(link='log'), data=df_ANDRO)
Anova(A1)
qqnorm(resid(A1))
hist(resid(A1))
plot(resid(A1))
A2<-glmmPQL(A4.ng.ml~Rank+Preg.stage+Age.months+Total.monthly.Rainfall, random=  ~1|ID, family=Gamma(link='log'), data=df_ANDRO)
Anova(A2)
A3<-glmmPQL(A4.ng.ml~Rank+Preg.stage+Total.monthly.Rainfall, random=  ~1|ID, family=Gamma(link='log'), data=df_ANDRO)
Anova(gam3)
A4<-glmmPQL(A4.ng.ml~Rank+Preg.stage, random=  ~1|ID, family=Gamma(link='log'), data=df_ANDRO)
Anova(A4)
summary(A4)
qqnorm(resid(A4))
hist(resid(A4))
plot(resid(A4))
# check model by adding dropped terms back in and confirming they are none-significant

#POSTHOC 
APOSTHOC <- lsmeans(A4, ~ Preg.stage)
summary(as.glht(pairs(APOSTHOC), by = NULL))

# Analysis for Testosterone####

# create a datset for T
df_TESTO <- Pregserum %>% filter(T.ng.ml != "na")%>% droplevels()

# to convert variables to factors 
df_TESTO$Rank <- factor(df_TESTO$Rank) 
df_TESTO$ID <- factor(df_TESTO$ID)
df_TESTO$Preg.stage <- factor(df_TESTO$Preg.stage)
df_TESTO$Preg.stage= factor(df_TESTO$Preg.stage, levels(df_TESTO$Preg.stage) [c(1,2,3)])

#check for colinarity
vif(glmmPQL(T.ng.ml~Rank*Preg.stage+Age.months+Total.monthly.Rainfall, random=  ~1|ID, family=Gamma(link='log'), data=df_TESTO))

#run glmm

T1<-glmmPQL(T.ng.ml~Rank*Preg.stage+Age.months+Total.monthly.Rainfall, random=  ~1|ID, family=Gamma(link='log'), data=df_TESTO)
Anova(T1)
qqnorm(resid(T1))
hist(resid(T1))
plot(resid(T1))
T2<-glmmPQL(T.ng.ml~Rank+Preg.stage+Age.months+Total.monthly.Rainfall, random=  ~1|ID, family=Gamma(link='log'), data=df_TESTO)
Anova(T2)
T3<-glmmPQL(T.ng.ml~Rank+Preg.stage+Total.monthly.Rainfall, random=  ~1|ID, family=Gamma(link='log'), data=df_TESTO)
Anova(T3)
T4<-glmmPQL(T.ng.ml~Rank+Preg.stage, random=  ~1|ID, family=Gamma(link='log'), data=df_TESTO)
Anova(T4)
summary(T4)
qqnorm(resid(T4))
hist(resid(T4))
plot(resid(T4))

# check model by adding dropped terms back in and confirming they are none-significant

#POSTHOC
TPOSTHOC <- lsmeans(T4, ~ Preg.stage)
summary(as.glht(pairs(TPOSTHOC), by = NULL))


# Analysis for FAM ####

# load faecal dataset####

Faecal <- read_csv("FAM_datset.csv", col_types = cols(ng.g.feces = col_number(),Age.months = col_number(), Preg.Stage = col_factor(levels = c("MP","LP", "PP")), Preg.Treatment = col_factor(levels = c("Control","Subordinate")), AM.PM = col_factor(levels = c ("AM", "PM")), Total.monthly.Rainfall = col_number()))

View(Faecal)

#log FaM values
Faecal$lt=log(Faecal$ng.g.feces)

# explore dataset 
par(mfrow=c(2,2))
hist(Faecal$ng.g.feces)
hist(Faecal$lt)
hist(log10(Faecal$ng.g.feces))
hist(sqrt(Faecal$ng.g.feces))
par(mfrow=c(1,1)) 
shapiro.test(Faecal$lt)
shapiro.test(Faecal$ng.g.feces)

#check for colinarity
vif(lmer(lt~Preg.Stage*Preg.Treatment+Age.months+Total.monthly.Rainfall+AM.PM+(1|KMP.ID/Litter), data=Faecal))

#starting model 
m1<-lmer(lt~Preg.Stage*Preg.Treatment+Age.months+Total.monthly.Rainfall+AM.PM+(1|KMP.ID/Litter), data=Faecal)
Anova(m1)
res<-residuals(m1, type="response")
plot(exp(predict(m1)), res)
qqnorm(res)
m2<-lmer(lt~Preg.stage+Preg.Treatment+Age.months+Total.monthly.Rainfall+AM.PM+(1|KMP.ID/Litter), data=Faecal)
Anova(m2)
m3<-lmer(lt~Preg.stage+Preg.Treatment+Age.months+AM.PM+(1|KMP.ID/Litter), data=Faecal)
Anova(m3)
m4<-lmer(lt~Preg.stage+Preg.Treatment+Age.months+(1|KMP.ID/Litter), data=Faecal)
Anova(m4)
m5<-lmer(lt~Preg.stage+Preg.Treatment+(1|KMP.ID/Litter), data=Faecal)
Anova(m5)

#AICtab
AICc(m1,m2,m3,m4,m5)
# AIC confirms m5 is best

#m5 is best model - check residuals etc
m5<-lmer(lt~Preg.stage+Preg.Treatment+(1|KMP.ID/Litter), data=Faecal)
Anova(m5)
summary(m5)
res<-residuals(m5, type="response")
plot(exp(predict(m5)), res)
qqnorm(res)

# check model by adding dropped terms back in and confirming they are none-significant

#POSTHOC
Fposthoc <- lsmeans(m5, ~ Preg.stage)
summary(as.glht(pairs(Fposthoc), by = NULL))
