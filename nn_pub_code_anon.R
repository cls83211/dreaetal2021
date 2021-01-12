library(dplyr)
library(lubridate)
library(lmerTest)
library(car)

# read in analysis table:
analysis = readRDS("./data/nn_analysis_anon.rds")

head(analysis)

# check collinearity:
fm_temp = lm(nn_per_hr ~ 
               psex + 
               condition + 
               dom_m + 
               sub_f_preg + 
               grpsize, 
             data = analysis)

vif(fm_temp)

# run analysis:
fm_1 = lmer(nn_per_hr ~ 
              psex + 
              condition + 
              dom_m + 
              sub_f_preg + 
              grpsize + 
              (1|dominant_fem) + 
              (1|partner), 
            data = analysis2)

summary(fm_1)
