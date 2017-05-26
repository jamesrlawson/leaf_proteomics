# protein summary stats

l <- select(protein_D14_age, total_protein, Photosystems, Calvin_cycle, Photorespiration, Light_reactions)
l$all_photosynth <- l$Light_reactions + l$Calvin_cycle + l$Photorespiration
l <- melt(l)
l <- ddply(l, .(variable), summarise, mean = mean(value, na.rm=TRUE), SD = sd(value, na.rm=TRUE))


# within species / site combination
# - leaf age
# - biological rep

require(dplyr)

# total CV
CV(protein_D14_age$total_protein)

blah <- merge(protein_D14_age, replicates, by = 'sample')

# interspecific variation (CV of ID means)

blax <- blah %>% group_by(ID) %>%
  summarise(mean = CV(total_protein)) 
CV(blax$mean) 

# proportion of total variation accounted for by leaf age


blaz <- blah %>% group_by(ID, leaf_age.x) %>%
  summarise(mean = mean(total_protein)) %>%
  group_by(ID) %>%
  summarise(CV = CV(mean))
mean(blaz$CV, na.rm=TRUE) 


# proportion of total variation accounted for by ID (i.e species)
# this is the average error bar size

blax <- blah %>% group_by(ID) %>%
  summarise(CV = CV(total_protein))
mean(blax$CV)  # should this be total CV of everything, or CV of the ID means?



# proportion of total variation accounted for by intraspecific variation

blaz <- blah %>% group_by(ID, biological_rep) %>%
  summarise(mean = mean(total_protein)) %>%
  group_by(ID) %>%
  summarise(CV = CV(mean))
mean(blaz$CV, na.rm=TRUE) 
