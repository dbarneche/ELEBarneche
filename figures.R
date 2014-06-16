#this file will only load RData files that already contain the outputs, so you will need R 2.15.3 running on your machine with the following packages in order to be able to fully replicate the figures:
#1. lme4 version 0.999999-0;
#2. R2jags version 0.04-01
#3. JAGS 3.3.0 (this is not a package, you have to install JAGS independently)

library(lme4)
library(R2jags)
rm(list=ls())
source("R/import.R")
load("re-run/analyses-10-nlmer.RData")
load("re-run/analyses-20-bayesian.RData")
load("re-run/analyses-a10-lmer.RData")
load("re-run/analyses-a20-standardRates.RData")

#fig2
to.pdf(fig2(data=sub.rate, nlModel=djagsout),  fig.path(name="fig2.pdf"), width=11, height=5.5)

#figS1
to.pdf(figS1(data=sub.rate, JAGSmodel=djagsout, nlmerModel=F2), fig.path(name="figS1.pdf"), width=8, height=10)

#figS2
to.pdf(figS2(JAGSmatrix=dfit$BUGSoutput$sims.matrix), fig.path(name="figS2.pdf"), width=9, height=9)

#figS3
to.pdf(figS3(JAGSmatrix=dfit$BUGSoutput$sims.matrix), fig.path(name="figS3.pdf"), width=8.5, height=8)

#figS4
to.pdf(figS4(data=sub.rate, nlModel=djagsout, JAGSModel=muljagsout), fig.path(name="figS4.pdf"), width=11, height=10) 

#figS5
to.pdf(figS5(data=standardRate, nlmerModel=St1), fig.path(name="figS5.pdf"), width=11, height=5.5)
