# Davies,C./Shearer,C.
# Meerkat manners - differences between control vs flut females

# 5/10/2021

# last modified 5/10/2021

#using R version 3.6.3

# clear R of all objects
rm(list=ls())

# load libraries
library(tidyverse)
#Load data and run exploratory analysis#############################################
#Re-name and organise####
# import and rename data file
######## FGc Dom Conrol vs. Dom Treated,  LP, vs. PP (2week PP)####
fcortdata<-read.csv("allpregcort_a.csv", header=T, na.strings=c(""," ","NA","na"))
head(fcortdata)
#View(fcortdata)
str(fcortdata)
#df_gc <- fcortdata %>% filter(NEW.DOSE..FGCM...ng.g. != "na")%>% droplevels()

fcortdata$fgcm.ng.g<-fcortdata$NEW.DOSE..FGCM...ng.g.
fcort_DCvsDT<- fcortdata[ which(fcortdata$rank2%in%c('DC','DF') ), ]
#rename
fcort_DCvsDT$Preg.Treatment<-fcort_DCvsDT$rank2
fcort_DCvsDT$Age.months<-fcort_DCvsDT$age.mo
fcort_DCvsDT$Total.monthly.Rainfall<-fcort_DCvsDT$rainfall..mm.
fcort_DCvsDT$KMP.ID<-fcort_DCvsDT$id
fcort_DCvsDT$Litter<-fcort_DCvsDT$litter

#OBS COUNTS by trt and stage
fcort_DCvsDTdf<- fcort_DCvsDT[ which(fcort_DCvsDT$Preg.Treatment%in%c('DF') ), ]
fcort_DCvsDTdc<- fcort_DCvsDT[ which(fcort_DCvsDT$Preg.Treatment%in%c('DC') ), ]
#EP: 8
count(fcort_DCvsDTdc[ which(fcort_DCvsDTdc$stage%in%c('(1) EP') ), ])
#MP: 32
count(fcort_DCvsDTdc[ which(fcort_DCvsDTdc$stage%in%c('(2) MP') ), ])
#LP: 45
count(fcort_DCvsDTdc[ which(fcort_DCvsDTdc$stage%in%c('(3) LP') ), ])
#PP: 49
count(fcort_DCvsDTdc[ which(fcort_DCvsDTdc$stage%in%c('(4) PP') ), ])

#EP: 1
count(fcort_DCvsDTdf[ which(fcort_DCvsDTdf$stage%in%c('(1) EP') ), ])
#MP: 0
count(fcort_DCvsDTdf[ which(fcort_DCvsDTdf$stage%in%c('(2) MP') ), ])
#LP: 18
count(fcort_DCvsDTdf[ which(fcort_DCvsDTdf$stage%in%c('(3) LP') ), ])
#PP: 22
count(fcort_DCvsDTdf[ which(fcort_DCvsDTdf$stage%in%c('(4) PP') ), ])

count(fcort_DCvsDTdf[ which(fcort_DCvsDTdf$AM.PM%in%c('AM') ), ])
count(fcort_DCvsDTdf[ which(fcort_DCvsDTdf$AM.PM%in%c('PM') ), ])
count(fcort_DCvsDTdc[ which(fcort_DCvsDTdc$AM.PM%in%c('PM') ), ])
count(fcort_DCvsDTdc[ which(fcort_DCvsDTdc$AM.PM%in%c('AM') ), ])

attach(fcort_DCvsDT)
lt=log(fgcm.ng.g)
hist(lt)
shapiro.test(lt)
#drop 1st two stages, not enough obs in DT
Preg.stage=factor(stage, levels(stage)[c(3,4)])
boxplot(lt~Preg.stage)
boxplot(lt~Preg.Treatment)
xtabs(~Preg.Treatment)
dotchart(fgcm.ng.g)

#check for colinarity want under 5 or under 10 for values
library(lme4)
vif(lmer(lt~Age.months+Total.monthly.Rainfall+AM.PM+Preg.stage*Preg.Treatment+(1|KMP.ID/Litter), fcort_DCvsDT))

library(lme4)
#starting model
m1<-lmer(lt~Age.months+Total.monthly.Rainfall+AM.PM+Preg.stage*Preg.Treatment+(1|KMP.ID/Litter), fcort_DCvsDT )
Anova(m1)
summary(m1)
res<-residuals(m1, type="response")
plot(exp(predict(m1)), res)
qqnorm(res)

m2<-lmer(lt~Age.months+Total.monthly.Rainfall+AM.PM+Preg.stage+Preg.Treatment+(1|KMP.ID/Litter), fcort_DCvsDT )
Anova(m2)
summary(m2)
# NOT TRUE FOR FGC THOUGH TRUE FOR FAM need to remove KMP.ID as a random factor- explains zero variance
m3<-lmer(lt~Total.monthly.Rainfall+Age.months+AM.PM+Preg.stage+(1|KMP.ID/Litter), fcort_DCvsDT )
Anova(m3)
summary(m3)

m4<-lmer(lt~Age.months+AM.PM+Preg.stage+(1|KMP.ID/Litter), fcort_DCvsDT )
Anova(m4)
summary(m4)

#MAM
m5<-lmer(lt~AM.PM+Preg.stage+(1|KMP.ID/Litter), fcort_DCvsDT )
Anova(m5)
summary(m5)

m6<-lmer(lt~Preg.stage+(1|KMP.ID/Litter), fcort_DCvsDT )
Anova(m6)
summary(m6)

m7<-lmer(lt~1+(1|KMP.ID/Litter), fcort_DCvsDT )
summary(m7)

AICtab(m1,m2,m3,m4,m5,m6,m7) # m5 is best

# check model by adding dropped terms back in and confirming they are none-significant

#pairs not needed since not more than two level comparison
#Fposthoc <- lsmeans(m5, ~ Preg.stage)
#summary(as.glht(pairs(Fposthoc), by = NULL))
#plot(Fposthoc)
