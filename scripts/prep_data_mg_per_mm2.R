replicates <- read_csv('output/replicates.csv')
replicates <- replicates[replicates$sample %in% protein_D14_age$sample,]

replicates <- merge(replicates, read_csv('data/lineage.csv'))

climate_locs$biological_rep <- NULL

data <- merge(protein_D14_age, climate_locs)
#data$ID <- NULL
data <- merge(data, replicates, by = c('sample', 'Latitude', 'Longitude', 'leaf_age', 'site_revised', 'species_confirmed', 'date'))
data <- data[!duplicated(data$sample),]

data <- filter(data, ID != 'melpal_106')

#data <- na.omit(data)

# calculate total_protein means and SE

total_protein_means <- data %>% group_by(ID) %>% summarise(total_protein_mean = mean(total_protein, na.rm=TRUE),
                                                           total_protein_SE = SE(total_protein))
data <- merge(total_protein_means, data)