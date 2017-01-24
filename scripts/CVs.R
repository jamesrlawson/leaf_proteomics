# within species / site combination
# - leaf age
# - biological rep

require(dplyr)

# total CV
CV(protein_D14_age$Photosystems)

blah <- merge(protein_D14_age, replicates, by = 'sample')

# proportion of total variation accounted for by ID (i.e species)

blax <- blah %>% group_by(ID) %>%
                 summarise(CV = CV(Photosystems))
mean(blax$CV)  # should this be total CV of everything, or CV of the ID means?


# proportion of total variation accounted for by leaf age


blaz <- blah %>% group_by(ID, leaf_age.x) %>%
  summarise(mean = mean(Photosystems)) %>%
  group_by(ID) %>%
  summarise(CV = CV(mean))
mean(blaz$CV, na.rm=TRUE) 



# proportion of total variation accounted for by intraspecific variation

blaz <- blah %>% group_by(ID, biological_rep) %>%
  summarise(mean = mean(Photosystems)) %>%
  group_by(ID) %>%
  summarise(CV = CV(mean))
mean(blaz$CV, na.rm=TRUE) 
