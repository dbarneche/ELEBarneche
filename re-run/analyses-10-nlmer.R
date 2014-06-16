#you will need R 2.15.3 running on your machine with lme4 version 0.999999-0

#this file applies the non-linear functions defined in R/functions-analyses.R with all combinations of possible parameters.
#models are named following Table S1 in the online SI on Ecology Letters website
#when fitting all models, with the inactivation I(T) term (Eq. 2b) and associated parameters (Ei and Topt),
#rather than estimate Er directly, we instead estimated the transformed quantity logit.Er, where Er=Ei⁄(1+exp(-logit.Er))
#to ensure convergence and enforce the biologically realistic assumption that Topt exists, which requires that Ei > Er.
#random-effects model selection uses restricted maximum likelihood (default), whereas fixed-effects model selection
#uses maximum likelihood, i.e. REML=FALSE

library(lme4)
######################################################
# LOAD, PROCESS & ANALYSE METABOLIC RATES DATA IN lme4
######################################################
source("re-run/processData-10-metabolicRates.R") #the main dataset used in the paper is the object called sub.rate
source('R/functions-analyses.R')

#################################
## RANDOM-EFFECTS MODEL SELECTION
#################################

m.pars <- c(ln.boTc  =  -5.73,
			logit.Er =  -0.1,
			Ea       =  0,
			Ei       =  2,
			Topt     =  306,
			a        =  0.76)

R1  <-  nlmer(rate~FullModel(ln.boTc, logit.Er, Ea, Ei, Topt, a, exp(weight), avkT, temp)~(a + ln.boTc + logit.Er + Topt|family), sub.rate, start=m.pars)

R2  <-  nlmer(rate ~ FullModel(ln.boTc, logit.Er, Ea, Ei, Topt, a, exp(weight), avkT, temp) ~ (a + ln.boTc + logit.Er | family), sub.rate, start=m.pars)

################################
## FIXED-EFFECTS MODEL SELECTION
################################

F1  <-  nlmer(rate~FullModel(ln.boTc, logit.Er, Ea, Ei, Topt, a, exp(weight), avkT, temp)~(a + ln.boTc + logit.Er + Topt|family), sub.rate, REML=FALSE, start=m.pars)

F2 <- nlmer(rate ~ FullMinusEa(ln.boTc, logit.Er, Ei, Topt, a, exp(weight), temp) ~ (a + ln.boTc + logit.Er + Topt | family), sub.rate, REML=FALSE, start=m.pars[-3])

F3 <- nlmer(rate ~ FullMinusIT(ln.boTc, logit.Er, Ea, a, Ei=2.35, exp(weight), avkT, temp) ~ (a + ln.boTc + logit.Er | family), sub.rate, REML=FALSE, start=m.pars[-c(4,5)])

F4 <- nlmer(rate ~ FullMinusEaIT(ln.boTc, logit.Er, a, Ei=2.35, exp(weight), temp) ~ (a + ln.boTc + logit.Er | family), sub.rate, REML=FALSE, start=m.pars[-c(3,4,5)])

rm(list=ls()[!(ls() %in% c('R1', 'R2', 'F1', 'F2', 'F3', 'F4', 'sub.rate', 'rate'))])
save.image('re-run/analyses-10-nlmer.RData')
