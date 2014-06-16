#IMPORTANT! THE DIFFERENT ANALYSES WILL REQUIRE DIFFERENT VERSIONS OF PACKAGES AND OF R - CHECK EACH FILE IN THE RE-RUN FO;DER FOR THE REQUIRED DEPENDENCIES

#this file combines all indivivual outputs together (see folder re-run/).
#it provides all material related to individual-level analyses presented in the paper.
#analyses-10-nlmer.R and analyses-20-bayesian.R will probably require many hours to run on a regular computer
#this file will only load RData files that already contain the outputs, so you will need R 2.15.3 running on your machine with the following packages in order to be able to fully replicate tables and figures:
#1. lme4 version 0.999999-0;
#2. R2jags version 0.04-01
#3. JAGS 3.3.0 (this is not a package, you have to install JAGS independently)

###########################################
# MODEL PERFORMANCE BASED ON LOG-LIKELIHOOD
###########################################
rm(list=ls())
load("re-run/analyses-10-nlmer.RData")
library(lme4)
anova(R1, R2)
anova(F1, F2) #F2 is the final parsimonious model
anova(F1, F3)
anova(F2, F4)

###################
# BAYESIAN ANALYSES
###################
rm(list=ls())
library(R2jags)
load("re-run/analyses-20-bayesian.RData")
#Topt range
round(djagsout['Topt',c('mean','2.5%','97.5%')] - 273.15)
#standard deviation of boTc according to JAGS
sqrt(djagsout['sigma2R[3,3]','mean'])
#how many families did not have their alphas and Ers within the estimated mean?
al  <-  djagsout[paste0("r[", 1:43, ",1]"),c("2.5%", "97.5%")] + djagsout["A", "mean"]
es  <-  djagsout[paste0("r[", 1:43, ",2]"),c("2.5%", "97.5%")] + djagsout['logit.Er','mean']
es  <-  djagsout['Ei','mean']/(1 + exp((-1)*es))
dim(al[al[,1] > 0.75 | al[,2] < 0.75,,drop=FALSE])[1] #8 families = 19%
dim(es[es[,1] > 0.7 | es[,2] < 0.6,,drop=FALSE])[1] #1 families = 2.3%
rm(al, es)


########################
# STANDARD RATES FROM SI
########################
rm(list=ls())
library(lme4)
load('re-run/analyses-a20-standardRates.RData')
pointEstimates  <-  coef(summary(St1))[,'Estimate']
stErrors        <-  coef(summary(St1))[,'Std. Error']
finalVals  <-  cbind(Estimates=pointEstimates, '2.5%'=(pointEstimates - stErrors*1.96), '97.5%'=(pointEstimates + stErrors*1.96))
finalVals  <-  round(finalVals, 2)
