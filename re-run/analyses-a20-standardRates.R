#you will need R 2.15.3 running on your machine with lme4 version 0.999999-0

library(lme4)
######################################################
# LOAD, PROCESS & ANALYSE METABOLIC RATES DATA IN lme4
######################################################
source("re-run/processData-20-standardRates.R")
source("R/functions-analyses.R")

St1 <- nlmer(rate~St(ln.boTc, Er, Ei, Topt, a, exp(weight), temp)~(a + ln.boTc + Er + Topt | family), standardRate, start=c(ln.boTc=-6.5, Er=0.6, Ei=2, Topt=306, a=0.76))

rm(list=ls()[!(ls() %in% c("standardRate", "St", "St1"))])
save.image("re-run/analyses-a20-standardRates.RData")
