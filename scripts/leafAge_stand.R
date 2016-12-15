#source('scripts/transformations.R')

# Leaf age plots

replicates <- read_csv('output/replicates.csv')
replicates <- replicates[replicates$sample %in% protein_stand_D14_age$sample,]

data <- merge(protein_stand_D14_age, climate_locs)
data <- merge(data, replicates, by = c('sample', 'Latitude', 'Longitude', 'leaf_age', 'site_revised'))
data <- data[!duplicated(data$sample),]

data[data$ID == 'corgum_47',]$ID <- 'corgum_46'

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

# TOTAL PROTEIN

total_protein <- data %>% select(total_protein, leaf_age, biological_rep, ID) %>%
  mutate(total_protein_stand = NA) %>%
  filter(!ID %in% 'melpal_106')

total_protein <- total_protein[order(total_protein[,4], total_protein[,2], total_protein[,3]),]

for(j in 0:2) {
  
  for(i in unique(total_protein$ID)) {
    
    x <- seq((j*3+1), (j+1)*3)
    
    total_protein[total_protein$ID %in% i,]$total_protein_stand[x] <- total_protein[total_protein$ID %in% i,]$total_protein[x] / subset(total_protein[total_protein$ID %in% i,], leaf_age == 'new')$total_protein
    
  }
  
}

check <- total_protein %>% 
  filter(leaf_age == 'new') %>% 
  group_by(ID, leaf_age) %>%
  summarise(mean_new = mean(total_protein_stand, na.rm=TRUE))

  # PLOT
  
  data$leaf_age <- factor(data$leaf_age, levels = c('new', 'mid', 'old'))
  
  boxplot(total_protein ~ leaf_age, data)
  
  total_protein$leaf_age <- factor(total_protein$leaf_age, levels = c('new', 'mid', 'old'))
  
  boxplot(total_protein_stand ~ leaf_age, total_protein)
  
  summary(aov(total_protein_stand ~ leaf_age, subset(total_protein, leaf_age != 'new')))
  

# PHOTOSYSTEMS
  
  Photosystems <- data %>% select(Photosystems, leaf_age, biological_rep, ID) %>%
    mutate(Photosystems_stand = NA) %>%
    filter(!ID %in% 'melpal_106')
  
  Photosystems <- Photosystems[order(Photosystems[,4], Photosystems[,2], Photosystems[,3]),]
  
  for(j in 0:2) {
    
    for(i in unique(Photosystems$ID)) {
      
      x <- seq((j*3+1), (j+1)*3)
      
      Photosystems[Photosystems$ID %in% i,]$Photosystems_stand[x] <- Photosystems[Photosystems$ID %in% i,]$Photosystems[x] / subset(Photosystems[Photosystems$ID %in% i,], leaf_age == 'new')$Photosystems
      
    }
    
  }
  
  check <- Photosystems %>% 
    filter(leaf_age == 'new') %>% 
    group_by(ID, leaf_age) %>%
    summarise(mean_new = mean(Photosystems_stand, na.rm=TRUE))

  # PLOT
  
  data$leaf_age <- factor(data$leaf_age, levels = c('new', 'mid', 'old'))
  
  boxplot(Photosystems ~ leaf_age, data)
  
  Photosystems$leaf_age <- factor(Photosystems$leaf_age, levels = c('new', 'mid', 'old'))
  
  boxplot(Photosystems_stand ~ leaf_age, Photosystems)
  
  summary(aov(Photosystems_stand ~ leaf_age, subset(Photosystems, leaf_age != 'new')))
  
  # PSI
  
  PSI <- data %>% select(PSI, leaf_age, biological_rep, ID) %>%
    mutate(PSI_stand = NA) %>%
    filter(!ID %in% 'melpal_106')
  
  PSI <- PSI[order(PSI[,4], PSI[,2], PSI[,3]),]
  
  for(j in 0:2) {
    
    for(i in unique(PSI$ID)) {
      
      x <- seq((j*3+1), (j+1)*3)
      
      PSI[PSI$ID %in% i,]$PSI_stand[x] <- PSI[PSI$ID %in% i,]$PSI[x] / subset(PSI[PSI$ID %in% i,], leaf_age == 'new')$PSI
      
    }
    
  }
  
  check <- PSI %>% 
    filter(leaf_age == 'new') %>% 
    group_by(ID, leaf_age) %>%
    summarise(mean_new = mean(PSI_stand, na.rm=TRUE))
  
  # PLOT
  
  data$leaf_age <- factor(data$leaf_age, levels = c('new', 'mid', 'old'))
  
  boxplot(PSI ~ leaf_age, data)
  
  PSI$leaf_age <- factor(PSI$leaf_age, levels = c('new', 'mid', 'old'))
  
  boxplot(PSI_stand ~ leaf_age, PSI)
  
  summary(aov(PSI_stand ~ leaf_age, subset(PSI, leaf_age != 'new')))
  
  
  # electron_transport_minATPsynth
  
  electron_transport_minATPsynth <- data %>% select(electron_transport_minATPsynth, leaf_age, biological_rep, ID) %>%
    mutate(electron_transport_minATPsynth_stand = NA) %>%
    filter(!ID %in% 'melpal_106')
  
  electron_transport_minATPsynth <- electron_transport_minATPsynth[order(electron_transport_minATPsynth[,4], electron_transport_minATPsynth[,2], electron_transport_minATPsynth[,3]),]
  
  for(j in 0:2) {
    
    for(i in unique(electron_transport_minATPsynth$ID)) {
      
      x <- seq((j*3+1), (j+1)*3)
      
      electron_transport_minATPsynth[electron_transport_minATPsynth$ID %in% i,]$electron_transport_minATPsynth_stand[x] <- electron_transport_minATPsynth[electron_transport_minATPsynth$ID %in% i,]$electron_transport_minATPsynth[x] / subset(electron_transport_minATPsynth[electron_transport_minATPsynth$ID %in% i,], leaf_age == 'new')$electron_transport_minATPsynth
      
    }
    
  }
  
  check <- electron_transport_minATPsynth %>% 
    filter(leaf_age == 'new') %>% 
    group_by(ID, leaf_age) %>%
    summarise(mean_new = mean(electron_transport_minATPsynth_stand, na.rm=TRUE))
  
  # PLOT
  
  data$leaf_age <- factor(data$leaf_age, levels = c('new', 'mid', 'old'))
  
  boxplot(electron_transport_minATPsynth ~ leaf_age, data)
  
  electron_transport_minATPsynth$leaf_age <- factor(electron_transport_minATPsynth$leaf_age, levels = c('new', 'mid', 'old'))
  
  boxplot(electron_transport_minATPsynth_stand ~ leaf_age, electron_transport_minATPsynth)
  
  summary(aov(electron_transport_minATPsynth_stand ~ leaf_age, subset(electron_transport_minATPsynth, leaf_age != 'new')))
  
  
  # Cytochrome_b6f
  
  Cytochrome_b6f <- data %>% select(Cytochrome_b6f, leaf_age, biological_rep, ID) %>%
    mutate(Cytochrome_b6f_stand = NA) %>%
    filter(!ID %in% 'melpal_106')
  
  Cytochrome_b6f <- Cytochrome_b6f[order(Cytochrome_b6f[,4], Cytochrome_b6f[,2], Cytochrome_b6f[,3]),]
  
  for(j in 0:2) {
    
    for(i in unique(Cytochrome_b6f$ID)) {
      
      x <- seq((j*3+1), (j+1)*3)
      
      Cytochrome_b6f[Cytochrome_b6f$ID %in% i,]$Cytochrome_b6f_stand[x] <- Cytochrome_b6f[Cytochrome_b6f$ID %in% i,]$Cytochrome_b6f[x] / subset(Cytochrome_b6f[Cytochrome_b6f$ID %in% i,], leaf_age == 'new')$Cytochrome_b6f
      
    }
    
  }
  
  check <- Cytochrome_b6f %>% 
    filter(leaf_age == 'new') %>% 
    group_by(ID, leaf_age) %>%
    summarise(mean_new = mean(Cytochrome_b6f_stand, na.rm=TRUE))
  
  # PLOT
  
  data$leaf_age <- factor(data$leaf_age, levels = c('new', 'mid', 'old'))
  
  boxplot(Cytochrome_b6f ~ leaf_age, data)
  
  Cytochrome_b6f$leaf_age <- factor(Cytochrome_b6f$leaf_age, levels = c('new', 'mid', 'old'))
  
  boxplot(Cytochrome_b6f_stand ~ leaf_age, Cytochrome_b6f)
  
  summary(aov(Cytochrome_b6f_stand ~ leaf_age, subset(Cytochrome_b6f, leaf_age != 'new')))
  
  
  
  # Calvin_cycle
  
  Calvin_cycle <- data %>% select(Calvin_cycle, leaf_age, biological_rep, ID) %>%
    mutate(Calvin_cycle_stand = NA) %>%
    filter(!ID %in% 'melpal_106')
  
  Calvin_cycle <- Calvin_cycle[order(Calvin_cycle[,4], Calvin_cycle[,2], Calvin_cycle[,3]),]
  
  for(j in 0:2) {
    
    for(i in unique(Calvin_cycle$ID)) {
      
      x <- seq((j*3+1), (j+1)*3)
      
      Calvin_cycle[Calvin_cycle$ID %in% i,]$Calvin_cycle_stand[x] <- Calvin_cycle[Calvin_cycle$ID %in% i,]$Calvin_cycle[x] / subset(Calvin_cycle[Calvin_cycle$ID %in% i,], leaf_age == 'new')$Calvin_cycle
      
    }
    
  }
  
  check <- Calvin_cycle %>% 
    filter(leaf_age == 'new') %>% 
    group_by(ID, leaf_age) %>%
    summarise(mean_new = mean(Calvin_cycle_stand, na.rm=TRUE))
  
  # PLOT
  
  data$leaf_age <- factor(data$leaf_age, levels = c('new', 'mid', 'old'))
  
  boxplot(Calvin_cycle ~ leaf_age, data)
  
  Calvin_cycle$leaf_age <- factor(Calvin_cycle$leaf_age, levels = c('new', 'mid', 'old'))
  
  boxplot(Calvin_cycle_stand ~ leaf_age, Calvin_cycle)
  
  summary(aov(Calvin_cycle_stand ~ leaf_age, subset(Calvin_cycle, leaf_age != 'new')))