# the output of this script is several long format dataframes which combine environmental and trait data with protein amount data
# these dataframes can then be used in knitr .Rmd files found in the root directory to produce regression reports with scatterplots
# showing relationships between a target env variable or trait and all protein functional category amounts

# the following switches determine if protein data used is in mg_per_m2 or moles (used in binProteins.R)
if(!exists('mg_per_m2')){
  mg_per_m2 = TRUE
}

if(!exists('moles')) {
  moles = FALSE
}

# the following switches can be set to include or exclude traits/environmental variables which have missing data. 
# code will look for preexisting definitions set in knitr files first. Default is to skip these variables.

if(!exists('include_photosynthesis')) {
  include_photosynthesis = FALSE
}

if(!exists('include_chlorophyll')) {
  include_chlorophyll = FALSE
}

if(!exists('include_d13C')) {
  include_d13C = FALSE
}

if(!exists('include_leaf_N')) {
  include_leaf_N = FALSE
}

if(!exists('include_leaf_P')) {
  include_leaf_P = FALSE
}

if(!exists('include_soil_N')) {
  include_soil_N = FALSE
}

if(!exists('include_soil_P')) {
  include_soil_P = FALSE
}

source('scripts/binProteins.R')

sample_locations <- read_csv('data/sample_locations.csv')
climate <- read_csv('data/discovery_site_climate.csv')

# prentice alpha

alpha <- read.csv('data/alph.csv', header=T, stringsAsFactors = F)
x <- data.frame()

for(i in 70:112) {
  a <- data.frame(rep(i:i, 12))
  x <- rbind(a, x)
}

names(x) <- 'year'
alpha$year <- x[order(x$year),]
alpha$month <- NULL

x <- melt(alpha, id = 'year')
x <- ddply(x, .(variable), summarise, mean = mean(value))
climate$alpha <- x$mean
rm(a,x,alpha)

# climate

# round lats and longs by to capture spp. which are at close but non-identical locations
# merge climate data with sample locations data

climate$Latitude <- round(climate$Latitude, 2)
climate$Longitude <- round(climate$Longitude, 2)
sample_locations$Latitude <- round(sample_locations$Latitude, 2)
sample_locations$Longitude <- round(sample_locations$Longitude, 2)

climate_locs <- merge(sample_locations, climate, by = c('Longitude', 'Latitude')) # we're going to add all the env and trait data into this df 'climate_locs'
climate_locs <- climate_locs[!duplicated(climate_locs$sample),]
climate_locs <- climate_locs[climate_locs$sample %in% unique(protein_stand_D14$sample),]

climate_locs_test <- vector('numeric', length = 20) # create a vector to contain nrow(climate_locs) values for testing
n <- 1 # we'll iterate n every time a merge is conducted
climate_locs_test[n] <- nrow(climate_locs)

# climatic data calculated for one year preceding sample collection
recent_clim <- read_csv('data/recent_year_clim.csv')
recent_clim$date <- NULL
recent_clim$Latitude <- round(recent_clim$Latitude, 2)
recent_clim$Longitude <- round(recent_clim$Longitude, 2)

recent_clim_locs <- merge(sample_locations, recent_clim, by = c('Longitude', 'Latitude'))
recent_clim_locs <- recent_clim_locs[!duplicated(recent_clim_locs$sample),]
recent_clim_locs <- recent_clim_locs[recent_clim_locs$sample %in% unique(protein_stand_D14$sample),]

climate_locs <- merge(climate_locs, recent_clim_locs, by = c('sample', 'species', 'Longitude', 'Latitude'))

n <- n + 1 # we'll iterate n every time a merge is conducted
climate_locs_test[n] <- nrow(climate_locs)

# climatic data calculated for 30 days preceding sample collection

recent_clim <- read_csv('data/recent_month_clim.csv')
recent_clim$date <- NULL
recent_clim$Latitude <- round(recent_clim$Latitude, 2)
recent_clim$Longitude <- round(recent_clim$Longitude, 2)

recent_clim_locs <- merge(sample_locations, recent_clim, by = c('Longitude', 'Latitude'))
recent_clim_locs <- recent_clim_locs[!duplicated(recent_clim_locs$sample),]
recent_clim_locs <- recent_clim_locs[recent_clim_locs$sample %in% unique(protein_stand_D14$sample),]

rm(sample_locations,climate)

climate_locs <- merge(climate_locs, recent_clim_locs, by = c('sample', 'species', 'Longitude', 'Latitude', 'site_revised'))

n <- n + 1
climate_locs_test[n] <- nrow(climate_locs)

#climate_locs$tavg_recent <- round(climate_locs$tavg_recent , 2)

#climate_locs$prec <- log10(climate_locs$prec)
#climate_locs$soilN <- log10(climate_locs$soilN)
#climate_locs$prec_recent <- sqrt(climate_locs$prec_recent)

# LMA & LWC

LMA_LWC <- read_csv('data/LMA_LWC.csv')
LMA_LWC <- LMA_LWC[LMA_LWC$sample %in% climate_locs$sample,]
LMA_LWC$LMA_g_per_m2  <- as.numeric(LMA_LWC$LMA_g_per_m2)
LMA_LWC$LWC_percent  <- as.numeric(LMA_LWC$LWC_percent)

climate_locs <- merge(LMA_LWC, climate_locs)
rm(LMA_LWC)

n <- n + 1
climate_locs_test[n] <- nrow(climate_locs)

# photosynthesis

licor <- read.csv('data/licor/licor_data_canon2.csv', header=TRUE, stringsAsFactors = FALSE)
licor <- na.omit(licor)

licor <- licor[licor$sample %in% climate_locs$sample,]

licor$CO2_level <- 1
licor[licor$CO2S < 500.000,]$CO2_level <- 'photo_amb'
licor[licor$CO2S > 500.000,]$CO2_level <- 'photo_max'
licor <- tidyr::spread(licor, CO2_level, Photo)

licor$photo_max <- as.numeric(licor$photo_max)
licor$photo_amb <- as.numeric(licor$photo_amb)
licor$Cond <- as.numeric(licor$Cond)

#licor$photo_max <- NULL
licor$photo_amb <- NULL

licor <- na.omit(licor)

licor <- licor[licor$Cond > 0.05,]

licor <- dplyr::select(licor, photo_max, Cond, sample)

if(include_photosynthesis) {
  climate_locs <- merge(licor, climate_locs, all.y=TRUE, by = 'sample') # this is causing points to be deleted due to the na.omit(protein_climate_D14_stand) in the .Rmd's
  
  n <- n + 1
  climate_locs_test[n] <- nrow(climate_locs)
  
}

# leaf_CNP

leaf_CN <- read_csv('data/leaf_CNP/leaf_CN.csv')
leaf_CN <- leaf_CN[leaf_CN$sample %in% climate_locs$sample,]

if(include_leaf_N) {
  climate_locs <- merge(leaf_CN, climate_locs, all.y=TRUE, by = 'sample')
  climate_locs$N_per_area <- climate_locs$N * 10 * climate_locs$LMA_g_per_m2
  
  n <- n + 1
  climate_locs_test[n] <- nrow(climate_locs)
  
}

leaf_P <- read_csv('data/leaf_CNP/leaf_P.csv')
leaf_P <- leaf_P[leaf_P$sample %in% climate_locs$sample,]

if(include_leaf_P) {
  climate_locs <- merge(leaf_P, climate_locs, all.y=TRUE, by = 'sample')
  climate_locs$P_per_area <- climate_locs$P * 10 * climate_locs$LMA_g_per_m2
  
  n <- n + 1
  climate_locs_test[n] <- nrow(climate_locs)
  
}

# leaf age 

leaf_age <- read_csv('data/leaf_age.csv')
leaf_age <- leaf_age[!duplicated(leaf_age[,c('sample', 'leaf_age')]),]
leaf_age <- leaf_age[leaf_age$sample %in% climate_locs$sample,]
climate_locs <- merge(leaf_age, climate_locs, all.y=TRUE, by = 'sample')

climate_locs <- climate_locs[!duplicated(climate_locs$sample),]

n <- n + 1
climate_locs_test[n] <- nrow(climate_locs)

# canopy openness

gaps <- read_csv('data/sky_pics.csv')
gaps <- gaps[gaps$sample %in% climate_locs$sample,]
#gaps <- gaps[!duplicated(gaps[,c('sample', 'gap')]),]
gaps$gap <- as.numeric(gaps$gap)

climate_locs <- merge(gaps, climate_locs)
rm(gaps)

n <- n + 1
climate_locs_test[n] <- nrow(climate_locs)

# adjust gap fraction for leaf age (c.f. Reich et al. 2009) - averaged values for corgum and E. haemostoma (0.225)
climate_locs[climate_locs$leaf_age == 'mid',]$gap <- climate_locs[climate_locs$leaf_age == 'mid',]$gap * 0.8875
climate_locs[climate_locs$leaf_age == 'old',]$gap <- climate_locs[climate_locs$leaf_age == 'old',]$gap * 0.775

# irradiance

irradiance <- read_csv('data/irradiance.csv')

irradiance <- irradiance[irradiance$Longitude %in% climate_locs$Longitude & irradiance$Latitude %in% climate_locs$Latitude,]
irradiance <- irradiance[!duplicated(irradiance),]

climate_locs <- merge(irradiance, climate_locs, by = c('Latitude', 'Longitude'))

n <- n + 1
climate_locs_test[n] <- nrow(climate_locs)

# irradiance at leaf

climate_locs$leaf_rad <- climate_locs$irradiance * climate_locs$gap / 100  


# chlorophyll

chl <- na.omit(read_csv('data/chlorophyll.csv'))
chl <- chl[chl$sample %in% protein_D14$sample,]

if(include_chlorophyll) {
  climate_locs <- merge(chl, climate_locs, by = 'sample', all.y=TRUE)
  
  n <- n + 1
  climate_locs_test[n] <- nrow(climate_locs)
  
}


# add in replicate ID's

replicates <- read_csv('output/replicates.csv')
replicates$site_revised <- NULL

climate_locs <- merge(climate_locs, replicates, by = c('sample', 'Latitude', 'Longitude', 'leaf_age'))

n <- n + 1
climate_locs_test[n] <- nrow(climate_locs)

# calculate gap means and gap SE

gap_mean <- climate_locs %>% group_by(ID) %>% dplyr::summarise(gap_mean = mean(gap, na.rm=TRUE), gap_SE = SE(gap))
climate_locs <- merge(climate_locs, gap_mean)

n <- n + 1
climate_locs_test[n] <- nrow(climate_locs)

# calculate leafrad means and leafrad SE

leafrad_mean <- climate_locs %>% group_by(ID) %>% dplyr::summarise(leafrad_mean = mean(leaf_rad, na.rm=TRUE), leafrad_SE = SE(leaf_rad))
climate_locs <- merge(climate_locs, leafrad_mean)

n <- n + 1
climate_locs_test[n] <- nrow(climate_locs)

# calculate LMA mean and LMA SE

LMA_mean <- climate_locs %>% group_by(ID) %>% dplyr::summarise(LMA_mean = mean(LMA_g_per_m2, na.rm=TRUE), LMA_SE = SE(LMA_g_per_m2))
climate_locs <- merge(climate_locs, LMA_mean)

n <- n + 1
climate_locs_test[n] <- nrow(climate_locs)

# calculate N_per_area mean and SE

if(include_leaf_N) {
  Narea_mean <- climate_locs %>% group_by(ID) %>% dplyr::summarise(Narea_mean = mean(N_per_area, na.rm=TRUE), Narea_SE = SE(N_per_area), Narea_CV = CV(N_per_area))
  climate_locs <- merge(Narea_mean, climate_locs, all.y = TRUE)
  
  n <- n + 1
  climate_locs_test[n] <- nrow(climate_locs)
}

# soil and litter data

soil_N <- read_csv('data/leaf_CNP/soil_N.csv')
soil_P <- read_csv('data/leaf_CNP/soil_P.csv') %>% 
  filter(soil_P < 1500)

soil_N <- soil_N[soil_N$ID %in% climate_locs$ID,]

if(include_soil_N) {
  climate_locs <- merge(soil_N, climate_locs, by = 'ID', all.y = TRUE)
  
  n <- n + 1
  climate_locs_test[n] <- nrow(climate_locs)
  
}

if(include_soil_P) {
  climate_locs <- merge(soil_P, climate_locs, by = 'ID', all.y = TRUE)
  
  n <- n + 1
  climate_locs_test[n] <- nrow(climate_locs)
  
}

# delta C 13

if(include_d13C) {
  
  d13C <- read_csv('data/d13C.csv')
  d13C <- d13C[d13C$sample %in% protein_D14$sample,]
  climate_locs <- merge(d13C, climate_locs, by = 'sample', all.y=TRUE)
  
  id <- unique(climate_locs$ID)
  
  # fill in values so d13C for all leaf ages of a branch = mid value
  
  for(i in 1:length(id)) {
    
    for(j in 1:3) {
      
      temp <- climate_locs[climate_locs$ID %in% id[i] & climate_locs$biological_rep %in% j,]  
      
      if(any(temp$leaf_age %in% 'mid')) {
        
        temp[temp$leaf_age %in% 'mid',]$d13C
        
        climate_locs[climate_locs$ID %in% id[i] & climate_locs$biological_rep %in% j,]$d13C <- temp[temp$leaf_age %in% 'mid',]$d13C
        
      }
      
    }
    
  }
  
  n <- n + 1
  climate_locs_test[n] <- nrow(climate_locs)
  
}

climate_locs$ID <- NULL
# climate_locs <- na.omit(climate_locs)

climate_locs <- climate_locs[!climate_locs$leaf_age %in% 'sen',]

n <- n + 1
climate_locs_test[n] <- nrow(climate_locs)

# these are the df's that are used in the knitr report

protein_D14 <- na.omit(protein_D14)
protein_stand_D14 <- na.omit(protein_stand_D14)

# cleanup

rm(chl, gap_mean, irradiance, leaf_age, leaf_CN, leaf_P, leafrad_mean, licor, LMA_mean, Narea_mean, 
   recent_clim, recent_clim_locs, replicates, soil_N, soil_P, d13C)

gc(verbose=FALSE)