source('scripts/transformations.R')

replicates <- read_csv('output/replicates.csv')
replicates <- replicates[replicates$sample %in% protein_stand_D14_age$sample,]

data <- merge(protein_stand_D14_age, climate_locs)
data <- merge(data, replicates, by = c('sample', 'Latitude', 'Longitude', 'leaf_age'), all=TRUE)

data$leaf_age <- factor(data$leaf_age, levels = c('new', 'mid', 'old'))

#boxplot(total_protein ~ leaf_age, data)


total_protein <- data %>% select(total_protein, leaf_age, biological_rep, ID) %>%
  mutate(total_protein_stand = NA) %>%
  filter(!ID %in% 'melpal_106')

total_protein[total_protein$ID == 'corgum_47',]$ID <- 'corgum_46'

corcit <- data.frame(total_protein = NA, leaf_age = 'mid', biological_rep = 2, ID = 'corcit_43')
cortes <- data.frame(total_protein = NA, leaf_age = 'old', biological_rep = 3, ID = 'cortes_51')
eucmed <- data.frame(total_protein = NA, leaf_age = 'new', biological_rep = 3, ID = 'eucmed_65')
eucdel <- data.frame(total_protein = NA, leaf_age = 'old', biological_rep = 3, ID = 'eucdel_58')
eucglo <- data.frame(total_protein = NA, leaf_age = 'new', biological_rep = 3, ID = 'eucglo_63')
eucpun <- data.frame(total_protein = NA, leaf_age = 'mid', biological_rep = 3, ID = 'eucpun_71')
eucspa <- data.frame(total_protein = NA, leaf_age = 'mid', biological_rep = 3, ID = 'eucspa_75')
eucten_new <- data.frame(total_protein = NA, leaf_age = 'new', biological_rep = 2, ID = 'eucten77')
eucten_old <- data.frame(total_protein = NA, leaf_age = 'old', biological_rep = 3, ID = 'eucten77')

#eucspa doesnt have biological rep 3 - and YG118 is missing from protein_stand_D14_age (!) this needs further investigation
#biological rep needs to be checked for added lines

total_protein <- total_protein %>% 
  bind_rows(corcit) %>% 
  bind_rows(cortes) %>% 
  bind_rows(eucmed) %>%
  bind_rows(eucdel) %>%
  bind_rows(eucglo) %>%
  bind_rows(eucpun) %>%
  bind_rows(eucspa)

total_protein <- total_protein[order(total_protein[,4], total_protein[,2], total_protein[,3]),]

for(j in 0:2) {
  
  for(i in unique(total_protein$ID)) {
    
    x <- seq((j*3+1), (j+1)*3)
    
    total_protein[total_protein$ID %in% i,]$total_protein_stand[x] <- total_protein[total_protein$ID %in% i,]$total_protein[x] / subset(total_protein[total_protein$ID %in% i,], leaf_age == 'new')$total_protein
    
    #   browser()
  }
  
  browser()
  
}