source('R/import.R')
#create a fake community-structure fish dataset to mimic the one presented at the paper
#150000 rows
#8 regions
#49 sites
#6000 transects
#1000 species
#5 trophic groups for each species
set.seed(1)
species  <-  sample(paste0('species_',1:1000), 150000, replace=TRUE)
regions  <-  sample(paste0('region_',1:8), 150000, replace=TRUE)
sites    <-  vector(mode='character', length=length(regions))
sites[regions=='region_1']  <-  'site_1'
sites[regions=='region_2']  <-  'site_2'
sites[regions=='region_3']  <-  sample(paste0('site_',3:6),   length(sites[regions=='region_3']), replace=TRUE)
sites[regions=='region_4']  <-  sample(paste0('site_',7:8),   length(sites[regions=='region_4']), replace=TRUE)
sites[regions=='region_5']  <-  sample(paste0('site_',9:33),  length(sites[regions=='region_5']), replace=TRUE)
sites[regions=='region_6']  <-  sample(paste0('site_',34:43), length(sites[regions=='region_6']), replace=TRUE)
sites[regions=='region_7']  <-  sample(paste0('site_',44:47), length(sites[regions=='region_7']), replace=TRUE)
sites[regions=='region_8']  <-  sample(paste0('site_',48:49), length(sites[regions=='region_8']), replace=TRUE)
trophnum        <-  sample(1:5, 1000, replace=TRUE)
names(trophnum) <-  unique(species)
trophnum        <-  trophnum[species]
transects       <-  sapply(unique(sites), function(x, data)paste0(x, '_', sample(1:100, length(data[data==x]), replace=TRUE)), data=sites)
transect_id     <-  vector(mode='character', length=length(regions))
for(j in unique(sites)) 
	transect_id[sites==j]  <-  transects[[j]] #tapply(transect_id, sites, function(x)length(unique(x)))
rm(transects)
biomass_g       <-  rgamma(150000,2,0.1)
abun            <-  rpois(150000, 10)
abun[abun==0]   <-  1
site.area       <-  rpois(49, 500000)
names(site.area)<-  unique(sites)
site.area       <-  site.area[sites]

data  <-  data.frame(region=regions,
					 sites=sites,
					 transect_id=transect_id,
					 species=species,
					 trophnum=trophnum,
					 abun=abun,
					 biomass_g=biomass_g,
					 site.area=site.area,
					 stringsAsFactors=FALSE)

rm(list=ls()[!(ls() %in% c('data'))])