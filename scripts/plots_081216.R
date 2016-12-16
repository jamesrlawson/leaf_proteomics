# protein functional categories vs env gradients

# aggregated by species, site

require(plyr)
require(ggplot2)

source('scripts/transformations.R')

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

# total protein

agg_plot(data, 'total_protein', 'prec', logx = TRUE, labs = c('Total protein vs MAP', 'Mean annual precip (log10)', 'Species mean of [Total protein]'))
agg_plot(data, 'total_protein', 'gap', logx = FALSE, labs = c('Total protein vs canopy gap fraction', 'Canopy gap fraction (%) ', 'Species mean of [Total protein] '))
agg_plot(data, 'total_protein', 'tavg', logx = FALSE, labs = c('Total protein vs MAT', 'Mean annual temp (degC)', 'Species mean of [Total protein]'))

# Calvin cycle

agg_plot(data, 'Calvin_cycle', 'prec', logx = TRUE, labs = c('Calvin Cycle vs MAP', 'Mean annual precip (log10)', 'Species mean of [Calvin cycle proteins] (rel)'))
agg_plot(data, 'Calvin_cycle', 'gap', logx = FALSE, labs = c('Calvin Cycle vs canopy gap fraction', 'Canopy gap fraction (%) ', 'Species mean of [Calvin cycle proteins] (rel)'))
agg_plot(data, 'Calvin_cycle', 'tavg', logx = FALSE, labs = c('Calvin Cycle vs MAT', 'Mean annual temp (degC)', 'Species mean of [Calvin cycle proteins] (rel)'))

# Photosystems

agg_plot(data, 'electron_transport_minATPsynth', 'prec', logx = TRUE, labs = c('Electron transport vs MAP', 'Mean annual precip (log10)', 'Species mean of [Electron transport proteins] (rel)'))
agg_plot(data, 'electron_transport_minATPsynth', 'gap', logx = FALSE, labs = c('Electron transport vs canopy gap fraction', 'Canopy gap fraction (%) ', 'Species mean of [Electron transport proteins] (rel)'))
agg_plot(data, 'electron_transport_minATPsynth', 'tavg', logx = FALSE, labs = c('Electron transport vs MAT', 'Mean annual temp (degC)', 'Species mean of [Electron transport proteins] (rel)'))

# Photosystems

agg_plot(data, 'Photosystems', 'prec', logx = TRUE, labs = c('Photosystems vs MAP', 'Mean annual precip (log10)', 'Species mean of [Photosystems proteins] (rel)'))
agg_plot(data, 'Photosystems', 'gap', logx = FALSE, labs = c('Photosystems vs canopy gap fraction', 'Canopy gap fraction (%) ', 'Species mean of [Photosystems proteins] (rel)'))
agg_plot(data, 'Photosystems', 'tavg', logx = FALSE, labs = c('Photosystems vs MAT', 'Mean annual temp (degC)', 'Species mean of [Photosystems proteins] (rel)'))









# total protein

total_protein_means <- ddply(data, .(ID), summarise, total_protein_mean = mean(total_protein, na.rm=TRUE), total_protein_se = SE(total_protein))
total_protein_means <- merge(data, total_protein_means, by = c('ID'))
total_protein_means <- total_protein_means[!duplicated(total_protein_means[,c('total_protein_mean')]),]

errorbar_width <- (max(data$prec, na.rm=TRUE) - min(data$tavg, na.rm=TRUE)) / 50
p <- ggplot(total_protein_means, aes(y = total_protein_mean, x = prec)) + geom_point(size = 3)
p <- p + geom_errorbar(aes(ymin = total_protein_mean - total_protein_se, ymax = total_protein_mean + total_protein_se), width = errorbar_width)
p <- p + geom_smooth(method = 'lm', se = F) + xlab('Mean annual precip (log)') + ylab('Species mean of [Calvin cycle proteins] (rel)')
p

summary(lm(total_protein ~ prec, total_protein_means))

errorbar_width <- (max(data$gap, na.rm=TRUE) - min(data$tavg, na.rm=TRUE)) / 50
p <- ggplot(total_protein_means, aes(y = total_protein_mean, x = gap)) + geom_point(size = 3)
p <- p + geom_errorbar(aes(ymin = total_protein_mean - total_protein_se, ymax = total_protein_mean + total_protein_se), width = errorbar_width)
p <- p + geom_smooth(method = 'lm', se = F) + xlab('Canopy gap fraction (%)') + ylab('Species mean of [Calvin cycle proteins] (rel)')
p

summary(lm(total_protein ~ gap, total_protein_means))

errorbar_width <- (max(data$tavg, na.rm=TRUE) - min(data$tavg, na.rm=TRUE)) / 50
p <- ggplot(total_protein_means, aes(y = total_protein_mean, x = tavg)) + geom_point(size = 3)
p <- p + geom_errorbar(aes(ymin = total_protein_mean - total_protein_se, ymax = total_protein_mean + total_protein_se), width = errorbar_width)
p <- p + geom_smooth(method = 'lm', se = F) + xlab('Mean annual temperature (deg C)') + ylab('Species mean of [Calvin cycle proteins] (rel)')
p

summary(lm(total_protein ~ tavg, total_protein_means))


# Calvin cycle

Calvin_cycle_means <- ddply(data, .(ID), summarise, Calvin_cycle_mean = mean(Calvin_cycle, na.rm=TRUE), Calvin_cycle_se = SE(Calvin_cycle))
Calvin_cycle_means <- merge(data, Calvin_cycle_means, by = c('ID'))
Calvin_cycle_means <- Calvin_cycle_means[!duplicated(Calvin_cycle_means[,c('Calvin_cycle_mean')]),]

errorbar_width <- (max(data$prec, na.rm=TRUE) - min(data$prec, na.rm=TRUE)) / 50
p <- ggplot(Calvin_cycle_means, aes(y = Calvin_cycle_mean, x = prec)) + geom_point(size = 3)
p <- p + geom_errorbar(aes(ymin = Calvin_cycle_mean - Calvin_cycle_se, ymax = Calvin_cycle_mean + Calvin_cycle_se), width = 0.02)
p <- p + geom_smooth(method = 'lm', se = F) + xlab('Mean annual precip (log)') + ylab('Species mean of [Calvin cycle proteins] (rel)')
p


summary(lm(Calvin_cycle_mean ~ prec, Calvin_cycle_means))

errorbar_width <- (max(data$gap, na.rm=TRUE) - min(data$prec, na.rm=TRUE)) / 50
p <- ggplot(Calvin_cycle_means, aes(y = Calvin_cycle_mean, x = gap)) + geom_point(size = 3)
p <- p + geom_errorbar(aes(ymin = Calvin_cycle_mean - Calvin_cycle_se, ymax = Calvin_cycle_mean + Calvin_cycle_se), width = 0.7)
p <- p + geom_smooth(method = 'lm', se = F) + xlab('Canopy gap fraction (%)') + ylab('Species mean of [Calvin cycle proteins] (rel)')
p

summary(lm(Calvin_cycle_mean ~ gap, Calvin_cycle_means))

errorbar_width <- (max(data$tavg, na.rm=TRUE) - min(data$prec, na.rm=TRUE)) / 50
p <- ggplot(Calvin_cycle_means, aes(y = Calvin_cycle_mean, x = tavg)) + geom_point(size = 3)
p <- p + geom_errorbar(aes(ymin = Calvin_cycle_mean - Calvin_cycle_se, ymax = Calvin_cycle_mean + Calvin_cycle_se), width = 0.7)
p <- p + geom_smooth(method = 'lm', se = F) + xlab('Mean annual temperature (deg C)') + ylab('Species mean of [Calvin cycle proteins] (rel)')
p

summary(lm(Calvin_cycle_mean ~ tavg, Calvin_cycle_means))



# electron transport

electron_transport_minATPsynth_means <- ddply(data, .(ID), summarise, electron_transport_minATPsynth_mean = mean(electron_transport_minATPsynth, na.rm=TRUE), electron_transport_minATPsynth_se = SE(electron_transport_minATPsynth))
electron_transport_minATPsynth_means <- merge(data, electron_transport_minATPsynth_means, by = c('ID'))
electron_transport_minATPsynth_means <- electron_transport_minATPsynth_means[!duplicated(electron_transport_minATPsynth_means[,c('electron_transport_minATPsynth_mean')]),]

errorbar_width <- (max(data$prec, na.rm=TRUE) - min(data$prec, na.rm=TRUE)) / 50
p <- ggplot(electron_transport_minATPsynth_means, aes(y = electron_transport_minATPsynth_mean, x = prec)) + geom_point(size = 3)
p <- p + geom_errorbar(aes(ymin = electron_transport_minATPsynth_mean - electron_transport_minATPsynth_se, ymax = electron_transport_minATPsynth_mean + electron_transport_minATPsynth_se), width = 0.02)
p <- p + geom_smooth(method = 'lm', se = F) + xlab('Mean annual precip (log)') + ylab('Species mean of [electron transport proteins] (rel)')
p

summary(lm(electron_transport_minATPsynth_mean ~ prec, electron_transport_minATPsynth_means))

errorbar_width <- (max(data$gap, na.rm=TRUE) - min(data$prec, na.rm=TRUE)) / 50
p <- ggplot(electron_transport_minATPsynth_means, aes(y = electron_transport_minATPsynth_mean, x = gap)) + geom_point(size = 3)
p <- p + geom_errorbar(aes(ymin = electron_transport_minATPsynth_mean - electron_transport_minATPsynth_se, ymax = electron_transport_minATPsynth_mean + electron_transport_minATPsynth_se), width = 0.7)
p <- p + geom_smooth(method = 'lm', se = F) + xlab('Canopy gap fraction (%)') + ylab('Species mean of [electron transport proteins] (rel)')
p

summary(lm(electron_transport_minATPsynth_mean ~ gap, electron_transport_minATPsynth_means))

errorbar_width <- (max(data$tavg, na.rm=TRUE) - min(data$prec, na.rm=TRUE)) / 50
p <- ggplot(electron_transport_minATPsynth_means, aes(y = electron_transport_minATPsynth_mean, x = tavg)) + geom_point(size = 3)
p <- p + geom_errorbar(aes(ymin = electron_transport_minATPsynth_mean - electron_transport_minATPsynth_se, ymax = electron_transport_minATPsynth_mean + electron_transport_minATPsynth_se), width = 0.7)
p <- p + geom_smooth(method = 'lm', se = F) + xlab('Mean annual temperature (degC)') + ylab('Species mean of [electron transport proteins] (rel)')
p

summary(lm(electron_transport_minATPsynth_mean ~ tavg, electron_transport_minATPsynth_means))


# photosystems

Photosystems_means <- ddply(data, .(ID), summarise, Photosystems_mean = mean(Photosystems, na.rm=TRUE), Photosystems_se = SE(Photosystems))
Photosystems_means <- merge(data, Photosystems_means, by = c('ID'))
Photosystems_means <- Photosystems_means[!duplicated(Photosystems_means[,c('Photosystems_mean')]),]

errorbar_width <- (max(data$tavg, na.rm=TRUE) - min(data$prec, na.rm=TRUE)) / 50
p <- ggplot(Photosystems_means, aes(y = Photosystems_mean, x = prec)) + geom_point(size = 3)
p <- p + geom_errorbar(aes(ymin = Photosystems_mean - Photosystems_se, ymax = Photosystems_mean + Photosystems_se), width = 0.02)
p <- p + geom_smooth(method = 'lm', se = F) + xlab('Mean annual precip (log)') + ylab('Species mean of [Photosystems proteins] (rel)')
p

summary(lm(Photosystems_mean ~ prec, Photosystems_means))

errorbar_width <- (max(data$gap, na.rm=TRUE) - min(data$gap, na.rm=TRUE)) / 50
p <- ggplot(Photosystems_means, aes(y = Photosystems_mean, x = gap)) + geom_point(size = 3)
p <- p + geom_errorbar(aes(ymin = Photosystems_mean - Photosystems_se, ymax = Photosystems_mean + Photosystems_se), width = 0.7)
p <- p + geom_smooth(method = 'lm', se = F) + xlab('Canopy gap fraction (%)') + ylab('Species mean of [Photosystems proteins] (rel)')
p

summary(lm(Photosystems_mean ~ gap, Photosystems_means))

# PSI

PSI_means <- ddply(data, .(ID), summarise, PSI_mean = mean(PSI, na.rm=TRUE), PSI_se = SE(PSI))
PSI_means <- merge(data, PSI_means, by = c('ID'))
PSI_means <- PSI_means[!duplicated(PSI_means[,c('PSI_mean')]),]

errorbar_width <- (max(data$tavg, na.rm=TRUE) - min(data$prec, na.rm=TRUE)) / 50
p <- ggplot(PSI_means, aes(y = PSI_mean, x = prec)) + geom_point(size = 3)
p <- p + geom_errorbar(aes(ymin = PSI_mean - PSI_se, ymax = PSI_mean + PSI_se), width = 0.02)
p <- p + geom_smooth(method = 'lm', se = F) + xlab('Mean annual precip (log)') + ylab('Species mean of [PSI proteins] (rel)')
p

summary(lm(PSI_mean ~ prec, PSI_means))

p <- ggplot(PSI_means, aes(y = PSI_mean, x = gap)) + geom_point(size = 3)
p <- p + geom_errorbar(aes(ymin = PSI_mean - PSI_se, ymax = PSI_mean + PSI_se), width = 0.7)
p <- p + geom_smooth(method = 'lm', se = F) + xlab('Canopy gap fraction (%)') + ylab('Species mean of [PSI proteins] (rel)')
p

summary(lm(PSI_mean ~ gap, PSI_means))


# leaf age plots

data$leaf_age <- factor(data$leaf_age, levels = c("new", "mid", "old"))

boxplot(total_protein ~ leaf_age, data, ylab = "[total protein] (mg/m2)")

summary(aov(total_protein ~ leaf_age, data))

boxplot(Photosystems ~ leaf_age, data, ylab = "[Photosystems] (rel)")

summary(aov(Photosystems ~ leaf_age, data))

boxplot(electron_transport_minATPsynth ~ leaf_age, data, ylab = "[electron transport] (rel)")

summary(aov(electron_transport_minATPsynth ~ leaf_age, data))

boxplot(Calvin_cycle ~ leaf_age, data, ylab = "[Calvin cycle] (rel)")

summary(aov(Calvin_cycle ~ leaf_age, data))


# influence of lineage

cor <- c('corcit', 'corcla', 'corexi', 'corgum', 'corint', 'corpoc', 'corter', 'cortes')
euc <- unique(data$species)[!unique(data$species) %in% cor & !unique(data$species) == 'angcos']

data$lineage <- NA
data[data$species %in% cor,]$lineage <- 'Corymbia'
data[data$species %in% euc,]$lineage <- 'Eucalyptus'


p <- ggplot(Calvin_cycle_means, aes(y = Calvin_cycle_mean, x = prec)) + geom_point(size = 3, aes(color = lineage))
p <- p + geom_errorbar(aes(ymin = Calvin_cycle_mean - Calvin_cycle_se, ymax = Calvin_cycle_mean + Calvin_cycle_se), width = 0.02)
p <- p + geom_smooth(data = Calvin_cycle_means[!is.na(Calvin_cycle_means$lineage),], method = 'lm', aes(color = lineage), se = F) + xlab('Mean annual precip') + ylab('Species mean of Calvin cycle (rel)')
p

p <- ggplot(Calvin_cycle_means, aes(y = Calvin_cycle_mean, x = 100 - gap)) + geom_point(size = 3, aes(color = lineage))
p <- p + geom_errorbar(aes(ymin = Calvin_cycle_mean - Calvin_cycle_se, ymax = Calvin_cycle_mean + Calvin_cycle_se), width = 0.7)
p <- p + geom_smooth(data = Calvin_cycle_means[!is.na(Calvin_cycle_means$lineage),], method = 'lm', aes(color = lineage), se = F) + xlab('Canopy cover') + ylab('Species mean of Calvin cycle (rel)')
p

p <- ggplot(electron_transport_minATPsynth_means, aes(y = electron_transport_minATPsynth_mean, x = prec)) + geom_point(size = 3, aes(color = lineage))
p <- p + geom_errorbar(aes(ymin = electron_transport_minATPsynth_mean - electron_transport_minATPsynth_se, ymax = electron_transport_minATPsynth_mean + electron_transport_minATPsynth_se), width = 0.02)
p <- p + geom_smooth(data = electron_transport_minATPsynth_means[!is.na(electron_transport_minATPsynth_means$lineage),], method = 'lm', aes(color = lineage), se = F) + xlab('Mean annual precip') + ylab('Species mean of electron transport (rel)')
p

p <- ggplot(electron_transport_minATPsynth_means, aes(y = electron_transport_minATPsynth_mean, x = 100 - gap)) + geom_point(size = 3, aes(color = lineage))
p <- p + geom_errorbar(aes(ymin = electron_transport_minATPsynth_mean - electron_transport_minATPsynth_se, ymax = electron_transport_minATPsynth_mean + electron_transport_minATPsynth_se), width = 0.7)
p <- p + geom_smooth(data = electron_transport_minATPsynth_means[!is.na(electron_transport_minATPsynth_means$lineage),], method = 'lm', aes(color = lineage), se = F) + xlab('Canopy cover') + ylab('Species mean of electron transport (rel)')
p

p <- ggplot(Photosystems_means, aes(y = Photosystems_mean, x = prec)) + geom_point(size = 3, aes(color = lineage))
p <- p + geom_errorbar(aes(ymin = Photosystems_mean - Photosystems_se, ymax = Photosystems_mean + Photosystems_se), width = 0.02)
p <- p + geom_smooth(data = Photosystems_means[!is.na(Photosystems_means$lineage),], method = 'lm', aes(color = lineage), se = F) + xlab('Mean annual precip') + ylab('Species mean of Photosystems (rel)')
p

p <- ggplot(Photosystems_means, aes(y = Photosystems_mean, x =  gap)) + geom_point(size = 3, aes(color = lineage))
p <- p + geom_errorbar(aes(ymin = Photosystems_mean - Photosystems_se, ymax = Photosystems_mean + Photosystems_se), width = 0.7)
p <- p + geom_smooth(data = Photosystems_means[!is.na(Photosystems_means$lineage),], method = 'lm', aes(color = lineage), se = F) + xlab('Canopy gap fraction') + ylab('Species mean of Photosystems (rel)')
p


# leaf age plots standardised to new leaf of biological rep

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


# variance partitioning

require(vegan)

Photosystems_means.varpart <- varpart(Photosystems_means$Photosystems_mean,
                                      ~soilN,
                                      ~gap,
                                      ~prec,
                                      ~tavg,
                                      data = Photosystems_means)
Photosystems_means.varpart
plot(Photosystems_means.varpart)



plot(total_protein ~ tavg, data)

