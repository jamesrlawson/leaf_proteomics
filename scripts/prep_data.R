# called by scripts/functions which use ID's to aggregate data 

replicates <- read_csv('data/misc_data/replicates.csv')
replicates <- replicates[replicates$sample %in% protein_D14$sample,]

replicates <- merge(replicates, read_csv('data/leaf_data/misc/lineage.csv'))

climate_locs$biological_rep <- NULL

data <- merge(protein_stand_D14, climate_locs)
#data$ID <- NULL
data <- merge(data, replicates, by = c('sample', 'Latitude', 'Longitude', 'leaf_age', 'site_revised', 'species_confirmed', 'date'))
data <- data[!duplicated(data$sample),]

data <- filter(data, ID != 'melpal_106')

#data <- na.omit(data)

# calculate total_protein means and SE

total_protein_means <- data %>% group_by(ID) %>% dplyr::summarise(total_protein_mean = mean(total_protein, na.rm=TRUE),
                                                           total_protein_SE = SE(total_protein),
                                                           total_protein_CV = CV(total_protein))
data <- merge(total_protein_means, data)


if(include_chlorophyll) {
  data$Cl_per_LHC <- data$mg_Cl_total_per_m2 / data$LHC
}

data$prec <- data$prec/1000



data$jmax <- data$cytochrome_b6f 

data$jmax_over_vcmax <- data$jmax / data$calvin_cycle