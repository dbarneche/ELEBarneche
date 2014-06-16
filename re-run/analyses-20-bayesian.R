#you will need R 3.0.2 running on your machine with package R2jags version 0.04-01 and software JAGS 3.3.0

#this script runs Bayesian analyses in JAGS 3.3.0 after finding the best model via nlmer function from package lme4
#it uses a trick to fit the temperature dependence term Er as follows:
#when fitting these models, with the inactivation I(T) term (Eq. 2b) and associated parameters (Ei and Topt),
#rather than estimate Er directly, we instead estimated the transformed quantity logit.Er, where Er=Ei⁄(1+exp(-logit.Er))
#to ensure convergence and enforce the biologically realistic assumption that Topt exists, which requires that Ei > Er.

library(R2jags)
source("re-run/processData-10-metabolicRates.R") #the main dataset used in the paper is the object called sub.rate

########################################################
# BAYESIAN ANALYSIS FOR DAILY METABOLIC RATES - logit.Er
########################################################
set.seed(1)
dJags       <-  list("ln.B"=sub.rate$rate, "ln.W"=sub.rate$weight, "FAM"=sub.rate$fam, "Te"=sub.rate$temp)
dinits1     <-  list(A=0.65, logit.Er=-1.5, Ei=1, Topt=295, tauB=3, boTc=-5.8)
dinits2     <-  list(A=0.75, logit.Er=-0.5, Ei=2, Topt=305, tauB=3, boTc=-5.8)
dinits3     <-  list(A=0.85, logit.Er=-2.5, Ei=3, Topt=315, tauB=3, boTc=-5.8)
dfit        <-  jags(data=dJags, inits=list(dinits1, dinits2, dinits3), parameters.to.save=c("A", "Er", "logit.Er", "Ei", "Topt", "boTc", "sigma2B", "sigma2R", as.character(sapply(1:4, function(x)paste0("r[", 1:43, ",", x,"]")))), model.file="bugs/multinorm_schoolfield.bug", n.chains=3, n.iter=5e5, DIC=TRUE, n.thin=250)
dfit        <-  autojags(dfit, n.iter=5e5, n.thin=250, n.update=100)
djagsout    <-  dfit$BUGSoutput$summary #jags output

rm(list=ls()[!(ls() %in% c('djagsout', 'dfit', 'rate', 'sub.rate'))])
save.image('re-run/analyses-20-bayesian.RData')
