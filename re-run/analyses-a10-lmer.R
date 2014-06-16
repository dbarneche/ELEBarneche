#you will need R 3.0.2 running on your machine with package R2jags version 0.04-01 and software JAGS 3.3.0

library(R2jags)
######################################################
# LOAD, PROCESS & ANALYSE METABOLIC RATES DATA IN lme4
######################################################
source("re-run/processData-10-metabolicRates.R")

#############################################
# BAYESIAN ANALYSIS FOR DAILY METABOLIC RATES
#############################################
set.seed(1)
mulJags       <-  list("ln.B"=sub.rate$rate, "ln.W"=sub.rate$weight, "FAM"=sub.rate$fam, "Te"=sub.rate$temp)
mulinits      <-  list(A=0.75, Er=0.65, boTc=-5, tauB=1)
mulfit        <-  jags(data=mulJags, inits=list(mulinits,mulinits,mulinits), parameters.to.save=c("A", "Er", "boTc", "sigma2B", "sigma2R", as.character(sapply(1:3, function(x)paste0("r[", 1:43, ",", x,"]")))), model.file="bugs/multinorm_boltzmann.bug", n.chains=3, n.iter=5e5, DIC=TRUE, n.thin=250)
mulfit        <-  autojags(mulfit, n.iter=5e5, n.thin=250, n.update=100)
muljagsout    <-  mulfit$BUGSoutput$summary #jags output
rm(list=ls()[!(ls() %in% c('muljagsout'))])

save.image('re-run/analyses-a10-lmer.RData')
