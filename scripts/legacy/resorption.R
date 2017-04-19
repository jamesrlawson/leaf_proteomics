require(readr)
require(dplyr)
require(tidyr)

replicates <- read_csv('data/replicates.csv')
replicates$Latitude <- round(replicates$Latitude, 2)
replicates$Longitude <- round(replicates$Longitude, 2)

replicates_short <- select(replicates, species_confirmed, Latitude, Longitude)
replicates_short <- replicates_short[!duplicated(replicates_short),]
replicates_short$ID <- paste(replicates_short$species_confirmed, with(replicates_short, ave(species_confirmed, FUN = seq_along)), sep = "_") #adds a sequence number for multiple values within datasets so I can spread long to wide without aggregating

replicates <- merge(replicates, replicates_short, by = c('species_confirmed', 'Latitude', 'Longitude'))

names(replicates)[4] <- 'sample'

#write_csv(replicates, 'output/replicates.csv')


leaf_CN <- read_csv('data/leaf_CNP/leaf_CN.csv')


blah <- merge(replicates, leaf_CN, by = 'sample') 
#blah <- merge(blah, leaf_P, by = 'sample')

blax_N <- blah %>% group_by(ID, biological_rep, leaf_age) %>% 
                   summarise(mean_N = mean(N, na.rm=TRUE)) %>%
                   spread(key = leaf_age, value = mean_N) %>%
                   mutate(resorption = sen/old)



leaf_P <- read_csv('data/leaf_CNP/leaf_P.csv')

blah_P <- merge(replicates, leaf_P, by = 'sample')

blax_P <- blah_P %>% group_by(ID, biological_rep, leaf_age) %>% 
  summarise(mean_P = mean(P, na.rm=TRUE)) %>%
  spread(key = leaf_age, value = mean_P) %>%
  mutate(resorption = sen/old)

names(blax_P)[7] <- 'resorption_P'
names(blax_N)[7] <- 'resorption_N'

d <- merge(blax_P, blax_N, by = 'ID', all=TRUE)

plot(resorption_N ~ resorption_P, d)

hist(log10(d$resorption_P))
hist(log10(d$resorption_N))

plot(log10(resorption_N) ~ log10(resorption_P), d, ylab = 'resorption_N', xlab = 'resorption_P')

plot(resorption_N ~ resorption_P, d, ylab = 'resorption_N', xlab = 'resorption_P')


#plot(resorption_N ~ soil_N, merge(d, climate_locs), ylab = 'N_resorption', xlab = 'soil_N')
#plot(resorption_P ~ soil_P, merge(d, climate_locs), ylab = 'P_resorption', xlab = 'soil_P')

  
#plot(N ~ soil_N, climate_locs)


bl <- merge(protein_D14, replicates, by = 'sample')

ba <- merge(blax_N, merge(protein_stand_D14, replicates, by = 'sample'), by = 'ID')

plot(ba$stress ~ ba$resorption_N, ylab = 'stress protein (mg/m2)', xlab = 'N resorption (senecent leaf N / old leaf N)')
abline(lm(ba$stress ~ ba$resorption_N))
summary(lm(ba$stress ~ ba$resorption_N))

                                                        