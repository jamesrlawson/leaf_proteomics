require(readr)
require(dplyr)
require(tidyr)

replicates <- read_csv('data/misc/replicates_raw.csv')
replicates$Latitude <- round(replicates$Latitude, 2)
replicates$Longitude <- round(replicates$Longitude, 2)

replicates_short <- select(replicates, species_confirmed, Latitude, Longitude)
replicates_short <- replicates_short[!duplicated(replicates_short),]
replicates_short$ID <- paste(replicates_short$species_confirmed, with(replicates_short, ave(species_confirmed, FUN = seq_along)), sep = "_") #adds a sequence number for multiple values within datasets so I can spread long to wide without aggregating

replicates <- merge(replicates, replicates_short, by = c('species_confirmed', 'Latitude', 'Longitude'))

names(replicates)[4] <- 'sample'

write_csv(replicates, 'data/misc/replicates.csv')