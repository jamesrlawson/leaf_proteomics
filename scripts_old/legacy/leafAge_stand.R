source('scripts/transformations.R')

# Leaf age plots

source('scripts/prep_data.R')

# TOTAL PROTEIN

total_protein <- data %>% dplyr::select(total_protein, leaf_age, biological_rep, ID) %>%
  mutate(total_protein_stand = NA) %>%
  filter(!ID %in% 'melpal_106')

total_protein <- total_protein[order(total_protein[,4], total_protein[,2], total_protein[,3]),]

for(j in 0:2) {
  
  for(i in unique(total_protein$ID)) {
    
    x <- seq((j*3+1), (j+1)*3)
    
    total_protein[total_protein$ID %in% i,]$total_protein_stand[x] <- total_protein[total_protein$ID %in% i,]$total_protein[x] / subset(total_protein[total_protein$ID %in% i,], leaf_age == 'new')$total_protein
    
   #browser()
    
  }
  
}

check <- total_protein %>% 
  filter(leaf_age == 'new') %>% 
  group_by(ID, leaf_age) %>%
  summarise(mean_new = mean(total_protein_stand, na.rm=TRUE))

  # PLOT
  
  data$leaf_age <- factor(data$leaf_age, levels = c('new', 'mid', 'old'))
  
  boxplot(total_protein ~ leaf_age, data, main = "Total protein vs leaf age", ylab = "total protein (mg/mm2)")
  
  total_protein$leaf_age <- factor(total_protein$leaf_age, levels = c('new', 'mid', 'old'))
  
  boxplot(total_protein_stand ~ leaf_age, total_protein, main = "Total protein vs leaf age, rel. to newest leaf", ylab = "total protein (rel. to newest leaf)")
  
  summary(aov(total_protein ~ leaf_age, subset(total_protein, leaf_age != 'new')))
  
  summary(aov(total_protein_stand ~ leaf_age, subset(total_protein, leaf_age != 'new')))
  

# PHOTOSYSTEMS
  
  Photosystems <- data %>% dplyr::select(Photosystems, leaf_age, biological_rep, ID) %>%
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
  
  boxplot(Photosystems ~ leaf_age, data, main = "Photosystems protein vs leaf age", ylab = "Photosystems protein (mg/mm2)")
  
  Photosystems$leaf_age <- factor(Photosystems$leaf_age, levels = c('new', 'mid', 'old'))
  
  boxplot(Photosystems_stand ~ leaf_age, Photosystems, main = "Photosystems protein vs leaf age, rel. to newest leaf", ylab = "Photosystems protein (rel. to newest leaf)")
  
  summary(aov(Photosystems ~ leaf_age, subset(Photosystems, leaf_age != 'new')))
  
  summary(aov(Photosystems_stand ~ leaf_age, subset(Photosystems, leaf_age != 'new')))
  
  # PSI
  
  PSI <- data %>% dplyr::select(PSI, leaf_age, biological_rep, ID) %>%
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
  
  boxplot(PSI ~ leaf_age, data, main = "PSI protein vs leaf age", ylab = "PSI protein (mg/mm2)")
  
  PSI$leaf_age <- factor(PSI$leaf_age, levels = c('new', 'mid', 'old'))
  
  boxplot(PSI_stand ~ leaf_age, PSI, main = "PSI protein vs leaf age, rel. to newest leaf", ylab = "PSI protein (rel. to newest leaf)")
  
  summary(aov(PSI_stand ~ leaf_age, subset(PSI, leaf_age != 'new')))
  
  
  # electron_transport_minATPsynth
  
  electron_transport_minATPsynth <- data %>% dplyr::select(electron_transport_minATPsynth, leaf_age, biological_rep, ID) %>%
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
  
  boxplot(electron_transport_minATPsynth ~ leaf_age, data, main = "Electron transport protein vs leaf age", ylab = "Electron transport protein (mg/mm2)")
  
  electron_transport_minATPsynth$leaf_age <- factor(electron_transport_minATPsynth$leaf_age, levels = c('new', 'mid', 'old'))
  
  boxplot(electron_transport_minATPsynth_stand ~ leaf_age, electron_transport_minATPsynth, main = "Electron transport protein vs leaf age, rel. to newest leaf", ylab = "Electron transport protein (rel. to newest leaf)")
  
  summary(aov(electron_transport_minATPsynth ~ leaf_age, subset(electron_transport_minATPsynth, leaf_age != 'new')))
  
  summary(aov(electron_transport_minATPsynth_stand ~ leaf_age, subset(electron_transport_minATPsynth, leaf_age != 'new')))
  
  
  # Cytochrome_b6f
  
  Cytochrome_b6f <- data %>% dplyr::select(Cytochrome_b6f, leaf_age, biological_rep, ID) %>%
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
  
  boxplot(Cytochrome_b6f ~ leaf_age, data, main = "Cytochrome b6f protein vs leaf age", ylab = "Cytochrome b6f transport protein (mg/mm2)")
  
  Cytochrome_b6f$leaf_age <- factor(Cytochrome_b6f$leaf_age, levels = c('new', 'mid', 'old'))
  
  boxplot(Cytochrome_b6f_stand ~ leaf_age, Cytochrome_b6f, main = "Cytochrome b6f protein vs leaf age, rel. to newest leaf", ylab = "Cytochrome b6f protein (rel. to newest leaf)")
  
  summary(aov(Cytochrome_b6f_stand ~ leaf_age, subset(Cytochrome_b6f, leaf_age != 'new')))
  
  
  
  # Calvin_cycle
  
  Calvin_cycle <- data %>% dplyr::select(Calvin_cycle, leaf_age, biological_rep, ID) %>%
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
  
  boxplot(Calvin_cycle ~ leaf_age, data, main = "Calvin cycle protein vs leaf age", ylab = "Calvin cycle transport protein (mg/mm2)")
  
  Calvin_cycle$leaf_age <- factor(Calvin_cycle$leaf_age, levels = c('new', 'mid', 'old'))
  
  boxplot(Calvin_cycle_stand ~ leaf_age, Calvin_cycle, main = "Calvin cycle protein vs leaf age, rel. to newest leaf", ylab = "Calvin cycle protein (rel. to newest leaf)")
  
  summary(aov(Calvin_cycle ~ leaf_age, subset(Calvin_cycle, leaf_age != 'new')))
  
  summary(aov(Calvin_cycle_stand ~ leaf_age, subset(Calvin_cycle, leaf_age != 'new')))
  