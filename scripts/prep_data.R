replicates <- read_csv('output/replicates.csv')
replicates <- replicates[replicates$sample %in% protein_stand_D14_age$sample,]

replicates <- merge(replicates, read_csv('data/lineage.csv'))

climate_locs$biological_rep <- NULL

data <- merge(protein_stand_D14_age, climate_locs)
#data$ID <- NULL
data <- merge(data, replicates, by = c('sample', 'Latitude', 'Longitude', 'leaf_age', 'site_revised', 'species_confirmed', 'date'))
data <- data[!duplicated(data$sample),]

data <- filter(data, ID != 'melpal_106')

data <- na.omit(data)

if(any(data$ID %in% 'corgum_47')) {
  data[data$ID == 'corgum_47',]$ID <- 'corgum_46'
}

corcit <- data.frame(data = NA, leaf_age = 'mid', biological_rep = 2, ID = 'corcit_43')
cortes <- data.frame(data = NA, leaf_age = 'old', biological_rep = 3, ID = 'cortes_51')
eucmed <- data.frame(data = NA, leaf_age = 'new', biological_rep = 3, ID = 'eucmed_65')
eucdel <- data.frame(data = NA, leaf_age = 'old', biological_rep = 3, ID = 'eucdel_58')
eucglo <- data.frame(data = NA, leaf_age = 'new', biological_rep = 2, ID = 'eucglo_63')
eucpun <- data.frame(data = NA, leaf_age = 'mid', biological_rep = 3, ID = 'eucpun_71')
eucspa <- data.frame(data = NA, leaf_age = 'old', biological_rep = 2, ID = 'eucspa_75')
eucten_new <- data.frame(data = NA, leaf_age = 'new', biological_rep = 2, ID = 'eucten_77')
eucten_old <- data.frame(data = NA, leaf_age = 'old', biological_rep = 3, ID = 'eucten_77')
corexi <- data.frame(data = NA, leaf_age = 'new', biological_rep = 3, ID = 'corexi_45')
corpoc <- data.frame(data = NA, leaf_age = 'old', biological_rep = 2, ID = 'corpoc_49')
eucdum <- data.frame(data = NA, leaf_age = 'new', biological_rep = 1, ID = 'eucdum_60')
euchae <- data.frame(data = NA, leaf_age = 'mid', biological_rep = 1, ID = 'euchae_64')
eucrub <- data.frame(data = NA, leaf_age = 'new', biological_rep = 3, ID = 'eucrub_73')

#eucspa doesnt have biological rep 3 (fixed) - and YG118 is missing from protein_D14_age (!) this needs further investigation
#biological rep needs to be checked for added lines

data <- data %>% 
  bind_rows(corcit) %>% 
  bind_rows(cortes) %>% 
  bind_rows(eucmed) %>%
  bind_rows(eucdel) %>%
  bind_rows(eucglo) %>%
  bind_rows(eucpun) %>%
  bind_rows(eucspa) %>%
  bind_rows(eucten_new) %>%
  bind_rows(eucten_old) %>%
  bind_rows(corexi) %>%
  bind_rows(corpoc) %>%
  bind_rows(eucdum) %>%
  bind_rows(euchae) %>%
  bind_rows(eucrub)

rm(corcit,cortes,eucmed,eucdel,eucglo,eucpun,eucspa,eucten_new,eucten_old,corexi,corpoc,eucdum,euchae,eucrub)

data$data <- NULL

#data <- na.omit(data)
