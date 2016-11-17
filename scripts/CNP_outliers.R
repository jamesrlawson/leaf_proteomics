# looking for outliers in CNP data

require(readr)
require(dplyr)

source('scripts/transformations.R')

# leaf CNP

leaf_P <- read_csv('data/leaf_CNP/leaf_P.csv')
leaf_N <- read_csv('data/leaf_CNP/leaf_CN.csv')

hist(leaf_P$P) # some outliers apparent, above 2000 ppm P
hist(leaf_CN$N)
hist(leaf_CN$C)

  # plot leaf N values against total protein (only have values for eucs)
  
  plot(N ~ total_protein, merge(leaf_CN, protein_D14, by = 'sample'))
  
# soil CNP  

soil_P <- read_csv('data/leaf_CNP/soil_P.csv')
soil_N <- read_csv('data/leaf_CNP/soil_N.csv')

plot(soil_P$soil_P)
plot(soil_N$soil_N)
plot(soil_N$soil_C)

# plot values against modelled values

# get climate/site df

sample_locations <- read_csv('data/sample_locations.csv')
climate <- read_csv('data/discovery_site_climate.csv')

climate$Latitude <- round(climate$Latitude, 2)
climate$Longitude <- round(climate$Longitude, 2)
sample_locations$Latitude <- round(sample_locations$Latitude, 2)
sample_locations$Longitude <- round(sample_locations$Longitude, 2)

climate_locs <- merge(sample_locations, climate, by = c('Longitude', 'Latitude')) # we're going to add all the env and trait data into this df 'climate_locs'
climate_locs <- climate_locs[!duplicated(climate_locs$sample),]

# aggregate modelled soil CNP values by ID

replicates <- read_csv('output/replicates.csv')

soil_check <- merge(replicates, climate_locs, by = c('sample', 'Latitude', 'Longitude'), all=TRUE)
soil_check <- soil_check %>% group_by(ID) %>% summarise(soilN_mod = mean(soilN), soilP_mod = mean(soilP), soilC_mod = mean(soilC))
soil_check <- merge(soil_N, soil_check, by = 'ID')
soil_check <- merge(soil_P, soil_check, by = 'ID')

# plot measured vs modelled values

plot(soil_N ~ soilN_mod, soil_check)
plot(soil_P ~ soilP_mod, soil_check)
plot(soil_C ~ soilC_mod, soil_check)

# soil_P has outlier

unique(subset(merge(soil_check, replicates), soil_P >3000)$site_revised)
unique(subset(merge(soil_check, replicates), soil_P >3000)$ID)

# soil_C has outliers

unique(subset(merge(soil_check, replicates), soil_C >17)$site_revised)
unique(subset(merge(soil_check, replicates), soil_C >17)$ID)
unique(subset(merge(soil_check, replicates), site_revised == 'TAS_Mt_Wellington')$ID)

# ground litter        

litter_N <- read_csv('data/leaf_CNP/litter_N.csv')
litter_P <- read_csv('data/leaf_CNP/litter_P.csv')

hist(litter_N$litter_N)
hist(litter_N$litter_C)
hist(litter_P$litter_P)


# go back to leaf_P - check for outliers by plotting against soil_P

leaf_P_check <- merge(replicates, soil_check)
leaf_P_check <- merge(leaf_P_check, leaf_P)
plot(leaf_P_check$soil_P, leaf_P_check$P)

unique(subset(leaf_P_check, P > 2000)$site_revised)
unique(subset(leaf_P_check, P > 2000)$ID)

leaf_P_check2 <- merge(replicates, litter_P)
leaf_P_check2 <- merge(leaf_P_check2, leaf_P)
plot(leaf_P_check2$litter_P, leaf_P_check2$P)
