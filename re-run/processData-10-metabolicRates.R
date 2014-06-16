#brings in metabolic data with diet
load("re-run/database-10-metabolicRates.RData") 
rate      <-  metabolicRates[metabolicRates$activity == "routine", ]
rate      <-  rate[rate$stress == "none specified" | rate$stress == "", ]
rate      <-  rate[!is.na(rate$rate), c("species","family","rate","weight","temp","habitat","reef")]
rate      <-  rate[rate$weight > 0, ]
rate$temp <-  rate$temp + 273.15

AvOneOverkT  <-  1/(8.62e-5*rate$temp)
AvOneOverkT  <-  tapply(AvOneOverkT, rate$family, mean)
AvOneOverkT  <-  (1/(8.62e-5*293.15)) - AvOneOverkT
rate$avkT    <-  AvOneOverkT[match(rate$family, names(AvOneOverkT))]
rm(AvOneOverkT, metabolicRates)

#filter for families with at least 5 observations and 5 degree temperature range
fam5      <-  tapply(rate$temp, rate$family, function(x)c(max(x)-min(x), length(x)))
fam5      <-  names(fam5)[sapply(fam5,function(x)x[1]>=5 & x[2]>=5)]
sub.rate  <-  rate[rate$family %in% fam5,]
rm(fam5)

#transform rates for g C/m2/day
sub.rate$rate   <-  log(sub.rate$rate*(sub.rate$weight/1000)*0.009) #24/1000/32*12=0.0228125
sub.rate$weight <-  log(sub.rate$weight)
sub.rate        <-  sub.rate[!(sub.rate$family %in% c("Carangidae","Coryphaenidae")),]
sub.rate$fam    <-  as.numeric(as.factor(sub.rate$family))
sub.rate$rate.yr<-  log(exp(sub.rate$rate)*365)
