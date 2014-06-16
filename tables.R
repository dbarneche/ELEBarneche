#this file will only load RData files that already contain the outputs, so you will need R 2.15.3 running on your machine with the following packages in order to be able to fully replicate the tables:
#1. lme4 version 0.999999-0;
#2. R2jags version 0.04-01
#3. JAGS 3.3.0 (this is not a package, you have to install JAGS independently)

###################
# TABLE 1 FROM JAGS
###################
library(R2jags)
rm(list=ls())
load("re-run/analyses-20-bayesian.RData")
table1  <-  djagsout[c('A', 
                       'Er',
                       'boTc',
                       'Topt',
                       'Ei'), 
c('mean', '2.5%', '97.5%')]

table1  <-  round(table1, 3)
print(table1)

##########################
# TABLE S1 MODEL SELECTION
##########################
library(lme4)
rm(list=ls())
load("re-run/analyses-10-nlmer.RData")
tableS1   <-  data.frame(Model=c('Stage 1', 'R1 Full', 'R2 Full - deltaTopt', 'Stage 2', 'F1 Full', 'F2 Full -Ea', 'F3 Full - Ei - Topt - deltaTopt', 'F4 F2 - Ei - Topt - deltaTopt'), 
                       DF=rep("",8),
                       Chisq=rep("",8),
                       P=rep("",8),
                       stringsAsFactors=FALSE)
tableS1[3,2:4]  <-  c(anova(R1,R2)[['Chi Df']][2], round(anova(R1,R2)[['Chisq']][2], 2), round(anova(R1,R2)[['Pr(>Chisq)']][2], 5))
tableS1[6,2:4]  <-  c(anova(F1,F2)[['Chi Df']][2], round(anova(F1,F2)[['Chisq']][2], 2), round(anova(F1,F2)[['Pr(>Chisq)']][2], 5))
tableS1[7,2:4]  <-  c(anova(F1,F3)[['Chi Df']][2], round(anova(F1,F3)[['Chisq']][2], 2), round(anova(F1,F3)[['Pr(>Chisq)']][2], 5))
tableS1[8,2:4]  <-  c(anova(F2,F4)[['Chi Df']][2], round(anova(F2,F4)[['Chisq']][2], 2), round(anova(F2,F4)[['Pr(>Chisq)']][2], 5))
print(tableS1)

###############################
# TABLE S2 ESTIMATES FROM NLMER
###############################
library(lme4)
rm(list=ls())
load("re-run/analyses-10-nlmer.RData")
tableS2   <-  coef(summary(F2))[c('a','logit.Er','ln.boTc','Topt','Ei'),]
stdevs    <-  as.numeric(attr(VarCorr(F2)[['family']], 'stddev'))[c(1,3,2,4)]
corrs     <-  attr(VarCorr(F2)[['family']], 'correlation')
corrs     <-  corrs[c(1,3,2,4),c(1,3,2,4)]
devs      <-  c(stdevs, corrs[col(corrs) < row(corrs)])
tableS2   <-  round(rbind(tableS2, matrix(c(devs, rep(0, length(devs)*2)), length(devs), 3)), 3)
rownames(tableS2)  <-  NULL
tableS2[tableS2==0]  <-  ""
tableS2  <-  as.data.frame(tableS2)
print(tableS2)

####################
# TABLE S3 FROM JAGS
####################
library(R2jags)
rm(list=ls())
load("re-run/analyses-20-bayesian.RData")
tableS3   <-  djagsout[c('A', 
                       'logit.Er',
                       'boTc',
                       'Topt',
                       'Ei', 
                       paste0('sigma2R[',1:4,',',1:4,']'), 
                       'sigma2R[1,2]',
                       'sigma2R[1,3]',
                       'sigma2R[1,4]',
                       'sigma2R[2,3]',
                       'sigma2R[2,4]',
                       'sigma2R[3,4]'), 
c('mean', '2.5%', '97.5%')]

tableS3  <-  round(tableS3, 3)
print(tableS3)