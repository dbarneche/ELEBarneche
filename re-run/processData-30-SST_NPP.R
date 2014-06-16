#values in this netcdf have been converted...
#...to byte values to reduce the size, thus one...
#would need to rescale the pixel values (the...
#...values you extract) as follows:
#SST= pixel value *0.075 - 3

#######################################################################
# SST FILTERED FROM 1997 - 2007 TO CORRESPOND SEAWIFS NPP TEMPORAL SPAN
#######################################################################
sst   <-  read.csv('data/cortad_sst.csv', header=TRUE, stringsAsFactors=FALSE, row.names=1)
sst   <-  sst*0.075-3+273.15 #transforms to Kelvin
year  <-  unlist(lapply(names(sst), function(x){substr(x, 2, 5)}))
sst   <-  sst[,which(year > 1996 & year < 2008)] #filter for specific time window
rm(year)

#main Reference: http://www.jstor.org/stable/2838857
#behrenfeld MJ, Falkowski PG (1997) Photosynthetic Rates Derived from Satellite-Based Chlorophyll Concentration. Limnology and Oceanography, 42 (1), 1-20
#seaWIFS satellite data (1/6 degree of resolution) obtained from http://data.guillaumemaze.org/ocean_productivity
##########################################
# INCORPORATING NPP FROM SEAWIFS 1997-2007
##########################################
npp  <-  read.csv('data/seawifs_npp.csv', header=TRUE, stringsAsFactors=FALSE, row.names=1)
npp  <-  rowSums(npp)
