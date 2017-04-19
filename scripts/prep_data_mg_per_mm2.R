replicates <- read_csv('output/replicates.csv')
replicates <- replicates[replicates$sample %in% protein_D14$sample,]

replicates <- merge(replicates, read_csv('data/lineage.csv'))

climate_locs$biological_rep <- NULL

data <- merge(protein_D14, climate_locs)
#data$ID <- NULL
data <- merge(data, replicates, by = c('sample', 'Latitude', 'Longitude', 'leaf_age', 'site_revised', 'species_confirmed', 'date'))
data <- data[!duplicated(data$sample),]

data <- filter(data, ID != 'melpal_106')

#data <- na.omit(data)

# calculate total_protein means and SE

total_protein_means <- data %>% group_by(ID) %>% summarise(total_protein_mean = mean(total_protein, na.rm=TRUE),
                                                           total_protein_SE = SE(total_protein))
data <- merge(total_protein_means, data)

data$electron_transport <- data$electron_transport_minATPsynth + data$ATP_synthase_chloroplastic

data$LHC <- data$LHC_I + data$LHC_II

data$LHCI_per_PSI <- data$LHC_I / data$PSI_min_LHCI
data$LHCI_per_PSII <- data$LHC_II / data$PSII_min_LHCII