#########################
## bring fishbase data ##
#########################
library(XML)
library(Hmisc)

theurl   <- "http://www.fishbase.de/Topic/List.php?group=29"
pagetree <- htmlTreeParse(theurl, error=function(...){})
urltable <- pagetree$children$html$children$body$children$table

urls     <- vector(mode='character')
families <- vector(mode='character')
species  <- vector(mode='character')
for(i in 1:length(urltable$children[[2]])) {
  urls[i]     <- unlist(urltable$children[[2]][[i]])[7]	
  families[i] <- unlist(urltable$children[[2]][[i]][[5]])[5]
  species[i]  <- unlist(urltable$children[[2]][[i]])[9]	
}
urls    <-  sub("..",'http://www.fishbase.de',urls)


URLREF  <- vector(mode='character')
for(i in 1:length(urls)) {
  
  newurl  <-  urls[i]
  sptree  <-  htmlTreeParse(newurl, error=function(...){})
  sptable <-  sptree$children$html$children$body$children$table
  URL  <- vector(mode='character')
  for(a in 1:length(sptable$children[[2]])) {
    if(length(unlist(sptable$children[[2]][[a]])[6]) > 0) {
      URL[a]  <- unlist(sptable$children[[2]][[a]])[6]	
    }
  }
  URLREF  <-  c(URLREF, URL)
}

URLREF  <- unlist(lapply(URLREF, function(x){paste("http://www.fishbase.de", x, sep="")}))


fishbase.rate <- data.frame(family='',species='',rate=NA,rate20=NA,weight=NA,temp=NA,salinity=NA,activity='',stress='',stringsAsFactors = FALSE)[-1,]

for(i in 1:length(urls)) {
  
  theurl   <- urls[i]
  pagetree <- htmlTreeParse(theurl, error=function(...){})
  urltable <- pagetree$children$html$children$body$children$table
  
  dog <- data.frame(family='',species='',rate=NA,rate20=NA,weight=NA,temp=NA,salinity=NA,activity='',stress='',
                    stringsAsFactors = FALSE)[rep(1,length(urltable$children[[2]])),]
  
  
  for(j in 1:length(urltable$children[[2]])) {
    
    dog$family[j]  <- families[i]
    
    dog$species[j] <- species[i]
    
    if(length(urltable$children[[2]][[j]][[1]])>0) {
      dog[j,3] <- as.numeric(sub(',','',unlist(urltable$children[[2]][[j]][[1]][[1]])[4]))}
    
    for(k in 2:5) {
      if(length(urltable$children[[2]][[j]][[k]])>0) {
        dog[j,k+2] <- as.numeric(sub(',','',unlist(urltable$children[[2]][[j]][[k]][[1]])[2]))}}
    
    if(length(urltable$children[[2]][[j]][[6]])>0) {
      dog[j,8] <- unlist(urltable$children[[2]][[j]][[6]][[1]])[2]}
    
    if(length(urltable$children[[2]][[j]][[7]])>0) {
      dog[j,9] <- unlist(urltable$children[[2]][[j]][[7]][[1]])[2]}
    
  }
  fishbase.rate  <-  rbind(fishbase.rate, dog)
}

#################################
## get traits for each species ##
#################################
spp.met       <-  sort(unique(fishbase.rate$species))
fishbase.spp  <-  vector()
for(i in 1:length(spp.met)){
  vec  <-  unlist(strsplit(spp.met[i], " "))
  if(length(vec) == 2){
    fishbase.spp[i]  <-  paste(vec[1], "-", vec[2], sep="")
  } else {
    fishbase.spp[i]  <-  paste(vec[1], "-", vec[2], "+", vec[3], sep="")
  }
}

fishbase.spp  <-  paste("http://fishbase.de/summary/", fishbase.spp, ".html", sep="")
rm(i, vec)

diet_num  <-  vector(mode='character', length=length(fishbase.spp))
habitat   <-  vector(mode='character', length=length(fishbase.spp))
reef      <-  vector(mode='character', length=length(fishbase.spp))

for(j in seq_along(fishbase.spp)) { 
  theurl   <-  htmlTreeParse(fishbase.spp[j], error=function(...){})
  step1    <-  unlist(theurl$children$html$children$body)
  step2    <-  step1[grep("Based", step1)]
  step3    <-  step1[grep("Freshwater;", step1)]
  step4    <-  step1[grep("Marine;", step1)]
  step5    <-  step1[grep("reef-associated", step1, ignore.case=TRUE)]
  
  if(length(step2) != 0)
    diet_num[j]    <-  step2
  else
    diet_num[j]    <-  "missing_diet"
  
  if(length(step5) != 0)
    reef[j]    <-  "yes"
  else
    reef[j]    <-  "no"
  
  if(length(step3) != 0 & length(step4) != 0)
    habitat[j]    <-  "Both"
  if(length(step3) != 0 & length(step4) == 0)
    habitat[j]    <-  "Freshwater"	  	
  if(length(step3) == 0 & length(step4) != 0)
    habitat[j]    <-  "Marine"	  	
  if(length(step3) == 0 & length(step4) == 0)
    habitat[j]    <-  "None"	  	
}


re          <-  ".+[[:space:]]+([0-9.]+)[[:space:]].*"
diet_num    <-  as.numeric(unname(sub(re, "\\1", diet_num)))
diet.table  <-  data.frame(species=spp.met, diet_num=diet_num, habitat=habitat, reef=reef, stringsAsFactors=FALSE)

diet.table$diet[diet.table$diet_num >= 2   & diet.table$diet_num <  2.20]  <-  "H"
diet.table$diet[diet.table$diet_num >= 2.2 & diet.table$diet_num <  2.80]  <-  "O"
diet.table$diet[diet.table$diet_num >= 2.8 & diet.table$diet_num <  3.70]  <-  "C"
diet.table$diet[diet.table$diet_num >= 3.7]                                <-  "P"

unique(diet.table$diet_num) #good to go
fishbase.rate  <-  merge(fishbase.rate, diet.table, by=c("species","species"))
fishbase.rate  <-  fishbase.rate[,!names(fishbase.rate) %in% c("rate20","diet_num","diet")]

#################################################
## bind new metabolic rates data for reef fish ##
#################################################
#first, download the reef-fish metabolic rates .csv file provided in the online supporting information. Then follow the script as below.

reef.rates  <-  read.csv("data/ELEbarnecheST1.csv", header=TRUE, stringsAsFactors=FALSE, na.strings=c("",NA))

#convert ml or uL of 02/h to mg02/kg/h
reef.rates$rate[reef.rates$rate_unit=="mlO2_per_hour"]  <-  reef.rates$rate[reef.rates$rate_unit=="mlO2_per_hour"]*1.429 / (reef.rates$weight_mg[reef.rates$rate_unit=="mlO2_per_hour"]/1000000)
reef.rates$rate[reef.rates$rate_unit=="uLO2_per_h"]     <-  reef.rates$rate[reef.rates$rate_unit=="uLO2_per_h"]*1.429*10e-4 / (reef.rates$weight_mg[reef.rates$rate_unit=="uLO2_per_h"]/1000000)
reef.rates$rate_unit  <-  "mgO2_per_kg_per_h"
#convert mass from mg to grams
reef.rates$weight_mg  <-  reef.rates$weight_mg/1000
#standardize taxonomic names
reef.rates$family   <-  capitalize(reef.rates$family)
reef.rates$species  <-  capitalize(gsub("_"," ",reef.rates$species))
#assuming psu, ppt and ppm are the equivalente under tropical sealevel conditions
reef.rates$salinity[reef.rates$salinity=="field_seawater"]  <-  35
reef.rates$salinity  <-  as.numeric(reef.rates$salinity)
#specify types of stress
reef.rates$stress  <-  "none specified"
reef.rates$stress[reef.rates$dissolved_oxygen %in% 3:5]  <-  "hypoxia"
reef.rates$stress[!reef.rates$salinity %in% 30:35]       <-  "salinity"

reef.rates     <-  data.frame(species=reef.rates$species, family=reef.rates$family, rate=reef.rates$rate, 
                              weight=reef.rates$weight_mg, temp=reef.rates$temperature, salinity=reef.rates$salinity,
                              activity=reef.rates$rate_type, stress=reef.rates$stress, habitat="Marine",
                              reef="yes", stringsAsFactors=FALSE)
metabolicRates  <-  rbind(fishbase.rate, reef.rates)
rm(list=ls()[!(ls() %in% c("metabolicRates"))])
save.image("re-run/database-10-metabolicRates.RData")