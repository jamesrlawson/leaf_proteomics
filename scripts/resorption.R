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

leaf_CN <- read_csv('data/leaf_CN.csv')


blah <- merge(replicates, leaf_CN, by = 'sample') 
#blah <- merge(blah, leaf_P, by = 'sample')

plot(blah$S.x, blah$S.y, ylab = 'leaf S (ppm)', xlab= 'leaf S (per mille)')

blax_N <- blah %>% group_by(ID, biological_rep, leaf_age) %>% 
                   summarise(mean_N = mean(N, na.rm=TRUE)) %>%
                   spread(key = leaf_age, value = mean_N) %>%
                   mutate(resorption = sen/old)



leaf_P <- read_csv('data/leaf_P.csv')

blah_P <- merge(replicates, leaf_P, by = 'sample')

blax_P <- blah_P %>% group_by(ID, biological_rep, leaf_age) %>% 
  summarise(mean_P = mean(P, na.rm=TRUE)) %>%
  spread(key = leaf_age, value = mean_P) %>%
  mutate(resorption = sen/old)

names(blax_P)[7] <- 'resorption_N'
names(blax_N)[7] <- 'resorption_P'

d <- merge(blax_P, blax_N, by = 'ID', all=TRUE)

plot(resorption_N ~ resorption_P, d)

hist(log10(d$resorption_P))
hist(log10(d$resorption_N))

plot(log10(resorption_N) ~ log10(resorption_P), d)



blah <- blah[-337,]

b <- select(blah, sample, leaf_age, N) %>%
  spread(key = leaf_age, value = N) %>%
  mutate(resorption = sen/old)
                 
c <- spread(blah_P, key = leaf_age, value = P) %>%
  mutate(resorption = sen/old)

d <- merge(b, c, by = 'sample')



                                                        
                                                        