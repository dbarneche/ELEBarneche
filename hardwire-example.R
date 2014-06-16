#this file will only load RData files that already contain the outputs, and you will need R (any version compatible with R2jags version 0.04-01) and JAGS 3.3.0 (this is not a package, you have to install JAGS independently).

#the function hardWireBayesian propagates the individual-level estimated parameters to the community level
library(plyr)
#simulate a dataset
source('re-run/database-20-communityData.R')
source('R/import.R')
source('re-run/processData-30-SST_NPP.R')
load('re-run/analyses-20-bayesian.RData')

yjagsIter           <-  dfit$BUGSoutput$sims.matrix
yjagsIter[,'boTc']  <-  yjagsIter[,'boTc'] + log(365) #transform units to year rates on the natural log scale

communityTable  <-  hardWireBayesian(iterationMatrix=yjagsIter, community_dataset=data, temperature_dataset=sst, lnNPPvector=log(npp))
cjagsout        <-  t(apply(communityTable, 2, function(x)c("mean"=mean(x, na.rm=TRUE), "sd"=sd(x, na.rm=TRUE), quantile(x, probs=c(0.025, 0.25, 0.5, 0.75, 0.975), type=2, na.rm=TRUE))))

#cjagsout is a matrix, check it's row.names
rownames(cjagsout)

#take the ln.sizeCorrectedBiomass for example. You will notice that there are 49 values for total ln.sizeCorrectedBiomass, and also for each trophic group: ln.sizeCorrectedBiomass1, ln.sizeCorrectedBiomass2, ln.sizeCorrectedBiomass3, ln.sizeCorrectedBiomass4, ln.sizeCorrectedBiomass5, where 1:5 represent the different trophic groups from the community dataset. 

#For example, extract mean +- 95% CI ln.sizeCorrectedBiomass for trophic group 5
cjagsout[paste0('ln.sizeCorrectedBiomass5[', 1:49, ']'),c('mean', '2.5%', '97.5%')]