# the output of this script is several long format dataframes which combine environmental and trait data with protein amount data
# these dataframes can then be used in knitr .Rmd files found in the root directory to produce regression reports with scatterplots
# showing relationships between a target env variable or trait and all protein functional category amounts

require(ggplot2)
require(vegan)
library(readr)
library(knitr)
require(plyr)
require(reshape2)
require(dplyr)
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

# climatic data calculated for one year preceding sample collection
recent_clim <- read_csv('data/recent_clim.csv')
recent_clim$date <- NULL
recent_clim$Latitude <- round(recent_clim$Latitude, 2)
recent_clim$Longitude <- round(recent_clim$Longitude, 2)

recent_clim_locs <- merge(sample_locations, recent_clim, by = c('Longitude', 'Latitude'))
recent_clim_locs <- recent_clim_locs[!duplicated(recent_clim_locs$sample),]
recent_clim_locs <- recent_clim_locs[recent_clim_locs$sample %in% unique(protein_stand_D14$sample),]

rm(sample_locations,climate)

climate_locs <- merge(climate_locs, recent_clim_locs, by = c('sample', 'species', 'Longitude', 'Latitude'))

#climate_locs$tavg_recent <- round(climate_locs$tavg_recent , 2)

climate_locs$prec <- log10(climate_locs$prec)
#climate_locs$soilN <- log10(climate_locs$soilN)
climate_locs$prec_recent <- sqrt(climate_locs$prec_recent)

# soil and litter data
  
  # add in replicate ID's
  
#  replicates <- read_csv('output/replicates.csv')

#soil_N <- read_csv('data/leaf_CNP/soil_N.csv')
#soil_P <- read_csv('data/leaf_CNP/soil_P.csv') %>% 
#          filter(soil_P < 1500)

#  climate_locs <- merge(climate_locs, replicates, by = c('sample', 'Latitude', 'Longitude'))
  
#  climate_locs <- merge(climate_locs, soil_N, by = 'ID')
#  climate_locs <- merge(climate_locs, soil_P, by = 'ID')

#  plot(soilN ~ soil_N, climate_locs)
#  plot(soilP ~ soil_P, subset(climate_locs, soil_P < 1500))
  
# canopy openness

gaps <- read_csv('data/sky_pics.csv')
gaps <- gaps[gaps$sample %in% climate_locs$sample,]
gaps$gap <- as.numeric(gaps$gap)

climate_locs <- merge(gaps, climate_locs)
rm(gaps)

# LMA & LWC

LMA_LWC <- read_csv('data/LMA_LWC.csv')
LMA_LWC <- LMA_LWC[LMA_LWC$sample %in% climate_locs$sample,]
LMA_LWC$LMA_g_per_m2  <- as.numeric(LMA_LWC$LMA_g_per_m2)
LMA_LWC$LWC_percent  <- as.numeric(LMA_LWC$LWC_percent)

climate_locs <- merge(LMA_LWC, climate_locs)
rm(LMA_LWC)


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

#licor$photo_max <- NULL
licor$photo_amb <- NULL

licor <- na.omit(licor)

licor <- licor[licor$Cond > 0.05,]

#climate_locs <- merge(licor, climate_locs, all.y=TRUE, by = 'sample') # this is causing points to be deleted due to the na.omit(protein_climate_D14_stand) in the .Rmd's

# leaf_CN

leaf_CN <- read_csv('data/leaf_CNP/leaf_CN.csv')
leaf_CN <- leaf_CN[leaf_CN$sample %in% climate_locs$sample,]
climate_locs <- merge(leaf_CN, climate_locs, all.y=TRUE, by = 'sample')


climate_locs$N_per_area <- climate_locs$N * 10 * climate_locs$LMA_g_per_m2

# leaf age

leaf_age <- read_csv('data/leaf_age.csv')
leaf_age <- leaf_age[!duplicated(leaf_age[,c('sample', 'leaf_age')]),]
leaf_age <- leaf_age[leaf_age$sample %in% climate_locs$sample,]
climate_locs <- merge(leaf_age, climate_locs, all.y=TRUE, by = 'sample')


## these are the df's used in most of the knitr reports

# relative [protein]
protein_stand_D14_long <- melt(protein_stand_D14)
names(protein_stand_D14_long) <- c('sample', 'bin_arch_name', 'sum')

protein_climate_D14_stand <- merge(climate_locs, protein_stand_D14_long, by = 'sample')
protein_climate_D14_stand$bin_arch_name <- as.character(protein_climate_D14_stand$bin_arch_name)

# absolute [protein]
protein_D14_long <- melt(protein_D14)
names(protein_D14_long) <- c('sample', 'bin_arch_name', 'sum')

protein_climate_D14 <- merge(climate_locs, protein_D14_long, by = 'sample')
protein_climate_D14$bin_arch_name <- as.character(protein_climate_D14$bin_arch_name)


# relative [protein] separated by leaf age

leaf_age <- read_csv('data/leaf_age.csv')
leaf_age <- leaf_age[!duplicated(leaf_age[,c('sample', 'leaf_age')]),]
protein_stand_D14_age <- merge(leaf_age, protein_stand_D14, by = 'sample')

protein_stand_D14_new <- subset(protein_stand_D14_age, leaf_age == "new")
protein_stand_D14_mid <- subset(protein_stand_D14_age, leaf_age == "mid")
protein_stand_D14_old <- subset(protein_stand_D14_age, leaf_age == "old")

protein_climate_D14_new_stand <- merge(climate_locs, melt(protein_stand_D14_new), by = 'sample')
names(protein_climate_D14_new_stand)[names(protein_climate_D14_new_stand) == 'variable'] <- 'bin_arch_name'
names(protein_climate_D14_new_stand)[names(protein_climate_D14_new_stand) == 'value'] <- 'sum'
protein_climate_D14_new_stand$bin_arch_name <- as.character(protein_climate_D14_new_stand$bin_arch_name)

protein_climate_D14_mid_stand <- merge(climate_locs, melt(protein_stand_D14_mid), by = 'sample')
names(protein_climate_D14_mid_stand)[names(protein_climate_D14_mid_stand) == 'variable'] <- 'bin_arch_name'
names(protein_climate_D14_mid_stand)[names(protein_climate_D14_mid_stand) == 'value'] <- 'sum'
protein_climate_D14_mid_stand$bin_arch_name <- as.character(protein_climate_D14_mid_stand$bin_arch_name)

protein_climate_D14_old_stand <- merge(climate_locs, melt(protein_stand_D14_old), by = 'sample')
names(protein_climate_D14_old_stand)[names(protein_climate_D14_old_stand) == 'variable'] <- 'bin_arch_name'
names(protein_climate_D14_old_stand)[names(protein_climate_D14_old_stand) == 'value'] <- 'sum'
protein_climate_D14_old_stand$bin_arch_name <- as.character(protein_climate_D14_old_stand$bin_arch_name)


# absolute [protein] separated by leaf age

leaf_age <- read_csv('data/leaf_age.csv')
leaf_age <- leaf_age[unique(leaf_age$sample) %in% protein_D14$sample,]
protein_D14_age <- merge(leaf_age, protein_D14, by = 'sample')

protein_D14_new <- subset(protein_D14_age, leaf_age == "new")
protein_D14_mid <- subset(protein_D14_age, leaf_age == "mid")
protein_D14_old <- subset(protein_D14_age, leaf_age == "old")

protein_climate_D14_new <- merge(climate_locs, melt(protein_D14_new), by = 'sample')
names(protein_climate_D14_new)[names(protein_climate_D14_new) == 'variable'] <- 'bin_arch_name'
names(protein_climate_D14_new)[names(protein_climate_D14_new) == 'value'] <- 'sum'
protein_climate_D14_new$bin_arch_name <- as.character(protein_climate_D14_new$bin_arch_name)

protein_climate_D14_mid <- merge(climate_locs, melt(protein_D14_mid), by = 'sample')
names(protein_climate_D14_mid)[names(protein_climate_D14_mid) == 'variable'] <- 'bin_arch_name'
names(protein_climate_D14_mid)[names(protein_climate_D14_mid) == 'value'] <- 'sum'
protein_climate_D14_mid$bin_arch_name <- as.character(protein_climate_D14_mid$bin_arch_name)

protein_climate_D14_old <- merge(climate_locs, melt(protein_D14_old), by = 'sample')
names(protein_climate_D14_old)[names(protein_climate_D14_old) == 'variable'] <- 'bin_arch_name'
names(protein_climate_D14_old)[names(protein_climate_D14_old) == 'value'] <- 'sum'
protein_climate_D14_old$bin_arch_name <- as.character(protein_climate_D14_old$bin_arch_name)




