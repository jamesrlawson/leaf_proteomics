# protein functional categories vs env gradients

# aggregated by species, site

require(plyr)
require(ggplot2)

source('scripts/transformations.R')

data <- na.omit(merge(protein_stand_D14_age, climate_locs))

data$Latitude <- round(data$Latitude, 2)
data$Longitude <- round(data$Longitude, 2)

# Calvin cycle

Calvin_cycle_means <- ddply(data, .(species, Longitude, Latitude), summarise, Calvin_cycle_mean = mean(Calvin_cycle), Calvin_cycle_se = SE(PSII))
Calvin_cycle_means <- merge(data, Calvin_cycle_means, by = c('species','Longitude', 'Latitude'))
Calvin_cycle_means <- Calvin_cycle_means[!duplicated(Calvin_cycle_means[,c('Longitude', 'Latitude', 'species', 'leaf_age')]),]

p <- ggplot(Calvin_cycle_means, aes(y = Calvin_cycle_mean, x = prec)) + geom_point(size = 3)
p <- p + geom_errorbar(aes(ymin = Calvin_cycle_mean - Calvin_cycle_se, ymax = Calvin_cycle_mean + Calvin_cycle_se), width = 0.02)
p <- p + geom_smooth(method = 'lm', se = F) + xlab('Mean annual precip (log)') + ylab('Species mean of [Calvin cycle proteins] (rel)')
p

summary(lm(Calvin_cycle_mean ~ prec, Calvin_cycle_means))

p <- ggplot(Calvin_cycle_means, aes(y = Calvin_cycle_mean, x = gap)) + geom_point(size = 3)
p <- p + geom_errorbar(aes(ymin = Calvin_cycle_mean - Calvin_cycle_se, ymax = Calvin_cycle_mean + Calvin_cycle_se), width = 0.7)
p <- p + geom_smooth(method = 'lm', se = F) + xlab('Canopy gap fraction (%)') + ylab('Species mean of [Calvin cycle proteins] (rel)')
p

summary(lm(Calvin_cycle_mean ~ gap, Calvin_cycle_means))

# electron transport

electron_transport_minATPsynth_means <- ddply(data, .(species, Longitude, Latitude), summarise, electron_transport_minATPsynth_mean = mean(electron_transport_minATPsynth), electron_transport_minATPsynth_se = SE(PSII))
electron_transport_minATPsynth_means <- merge(data, electron_transport_minATPsynth_means, by = c('species','Longitude', 'Latitude'))
electron_transport_minATPsynth_means <- electron_transport_minATPsynth_means[!duplicated(electron_transport_minATPsynth_means[,c('Longitude', 'Latitude', 'species', 'leaf_age')]),]

p <- ggplot(electron_transport_minATPsynth_means, aes(y = electron_transport_minATPsynth_mean, x = prec)) + geom_point(size = 3)
p <- p + geom_errorbar(aes(ymin = electron_transport_minATPsynth_mean - electron_transport_minATPsynth_se, ymax = electron_transport_minATPsynth_mean + electron_transport_minATPsynth_se), width = 0.02)
p <- p + geom_smooth(method = 'lm', se = F) + xlab('Mean annual precip (log)') + ylab('Species mean of [electron transport proteins] (rel)')
p

summary(lm(electron_transport_minATPsynth_mean ~ prec, electron_transport_minATPsynth_means))

p <- ggplot(electron_transport_minATPsynth_means, aes(y = electron_transport_minATPsynth_mean, x = gap)) + geom_point(size = 3)
p <- p + geom_errorbar(aes(ymin = electron_transport_minATPsynth_mean - electron_transport_minATPsynth_se, ymax = electron_transport_minATPsynth_mean + electron_transport_minATPsynth_se), width = 0.7)
p <- p + geom_smooth(method = 'lm', se = F) + xlab('Canopy gap fraction (%)') + ylab('Species mean of [electron transport proteins] (rel)')
p

summary(lm(electron_transport_minATPsynth_mean ~ gap, electron_transport_minATPsynth_means))

# photosystems

Photosystems_means <- ddply(data, .(species, Longitude, Latitude), summarise, Photosystems_mean = mean(Photosystems), Photosystems_se = SE(PSII))
Photosystems_means <- merge(data, Photosystems_means, by = c('species','Longitude', 'Latitude'))
Photosystems_means <- Photosystems_means[!duplicated(Photosystems_means[,c('Longitude', 'Latitude', 'species', 'leaf_age')]),]

p <- ggplot(Photosystems_means, aes(y = Photosystems_mean, x = prec)) + geom_point(size = 3)
p <- p + geom_errorbar(aes(ymin = Photosystems_mean - Photosystems_se, ymax = Photosystems_mean + Photosystems_se), width = 0.02)
p <- p + geom_smooth(method = 'lm', se = F) + xlab('Mean annual precip (log)') + ylab('Species mean of [Photosystems proteins] (rel)')
p

summary(lm(Photosystems_mean ~ prec, Photosystems_means))

p <- ggplot(Photosystems_means, aes(y = Photosystems_mean, x = gap)) + geom_point(size = 3)
p <- p + geom_errorbar(aes(ymin = Photosystems_mean - Photosystems_se, ymax = Photosystems_mean + Photosystems_se), width = 0.7)
p <- p + geom_smooth(method = 'lm', se = F) + xlab('Canopy gap fraction (%)') + ylab('Species mean of [Photosystems proteins] (rel)')
p

summary(lm(Photosystems_mean ~ gap, Photosystems_means))



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

