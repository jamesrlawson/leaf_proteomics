## running some models!

require('MuMIn')
require(vegan)
require(lme4)
require(plot3D)
require(plyr)
source('scripts/transformations.R')

### plots and models, data aggregated by species (within site) ##

data <- na.omit(merge(protein_stand_D14_age, climate_locs))
data <- na.omit(merge(protein_stand_D14_age, climate_locs))


# adjust gap fraction for leaf age (c.f. Reich et al. 2009) 
data[data$leaf_age == 'mid',]$gap <- data[data$leaf_age == 'mid',]$gap * 0.83
data[data$leaf_age == 'old',]$gap <- data[data$leaf_age == 'old',]$gap * 0.66

# scatter plots

scatter3D(x = data$prec, y = 100-data$gap, z = data$Light_reactions, bty='g', phi=10, highlight.3d=TRUE, pch = 20, col = 'blue',
          xlab = 'prec', ylab = 'cover', zlab = 'Light reactions (rel)')

plot(x = 100-data$gap, y = data$Light_reactions, xlab = 'canopy cover (%)', ylab = 'light reactions (rel)')

plot(x = data$prec, y = data$Light_reactions, xlab = 'precip', ylab = 'light reactions (rel)')

data$leaf_age = factor(data$leaf_age,c("new","mid","old"))
plot(x = as.factor(data$leaf_age), y = data$Light_reactions, xlab = 'leaf age', ylab = 'light reactions (rel)')


####### regression plots with species as point + error bars for intraspecific var

SE <- function(x) {
  x <- sd(x) / sqrt(length(x))
}

# round lats and longs by to capture spp. which are at close but non-identical locations

#data$Latitude <- round(data$Latitude, 2)
#data$Longitude <- round(data$Longitude, 2)

# add in lineages
cor <- c('corcit', 'corcla', 'corexi', 'corgum', 'corint', 'corpoc', 'corter', 'cortes')
euc <- unique(data$species)[!unique(data$species) %in% cor & !unique(data$species) == 'angcos']

data$lineage <- NA
data[data$species %in% cor,]$lineage <- 'Corymbia'
data[data$species %in% euc,]$lineage <- 'Eucalyptus'
#data[data$species == 'angcos',]$lineage <- 'ang'

# find aggregated species means (one per 'site')

LR_means <- ddply(data, .(species, Longitude, Latitude), summarise, LR_mean = mean(PSII), LR_se = SE(PSII))
LR_means <- merge(data, LR_means, by = c('species','Longitude', 'Latitude'))
LR_means <- LR_means[!duplicated(LR_means[,c('Longitude', 'Latitude', 'species', 'leaf_age')]),]

# canopy cover and precip vs LR_means

p <- ggplot(LR_means, aes(y = LR_mean, x = 100-gap)) + geom_point(size = 3, aes(color = leaf_age))
p <- p + geom_errorbar(aes(ymin = LR_mean - LR_se, ymax = LR_mean + LR_se), width = 0.7)
p <- p + geom_smooth(method = 'lm', se = F, aes(color = leaf_age)) + xlab('Canopy Cover') + ylab('Species mean of light reactions (rel)')
p

summary(lm(LR_mean ~ gap, LR_means))

p <- ggplot(LR_means, aes(y = LR_mean, x = prec)) + geom_point(size = 3)
p <- p + geom_errorbar(aes(ymin = LR_mean - LR_se, ymax = LR_mean + LR_se), width = 0.02)
p <- p + geom_smooth(method = 'lm', se = F) + xlab('Mean annual precip (log)') + ylab('Species mean of [light reactions proteins] (rel)')
p

summary(lm(LR_mean ~ prec, LR_means))

# lineage plots

p <- ggplot(LR_means, aes(y = LR_mean, x = prec)) + geom_point(size = 3, aes(color = lineage))
p <- p + geom_errorbar(aes(ymin = LR_mean - LR_se, ymax = LR_mean + LR_se), width = 0.02)
p <- p + geom_smooth(data = LR_means[!is.na(LR_means$lineage),], method = 'lm', aes(color = lineage), se = F) + xlab('Mean annual precip') + ylab('Species mean of light reactions (rel)')
p

p <- ggplot(LR_means, aes(y = LR_mean, x = 100 - gap)) + geom_point(size = 3, aes(color = lineage))
p <- p + geom_errorbar(aes(ymin = LR_mean - LR_se, ymax = LR_mean + LR_se), width = 0.7)
p <- p + geom_smooth(data = LR_means[!is.na(LR_means$lineage),], method = 'lm', aes(color = lineage), se = F) + xlab('Canopy cover') + ylab('Species mean of light reactions (rel)')
p

### mixed models ###

lm1 <- lmer(Light_reactions ~ prec + gap + (1|species), data = data, REML=F)
lm2 <- lmer(Light_reactions ~ prec + (1|species), data = data, REML=F)
lm3 <- lmer(Light_reactions ~ gap + (1|species), data = data, REML=F)
lm4 <- lmer(Light_reactions ~ (1|species), data = data, REML=F)

lm5 <- lmer(Light_reactions ~ prec + gap + leaf_age + (1|species), data = data, REML=F)
lm6 <- lmer(Light_reactions ~ gap + leaf_age + (1|species), data = data, REML=F)
lm7 <- lmer(Light_reactions ~ leaf_age + (1|species), data = data, REML=F)


AIC(lm1,lm2,lm3,lm4,lm5,lm6,lm7)
anova(lm1,lm2,lm3,lm4,lm5,lm6,lm7)


summary(lm5)


# variance partitioning in linear models

data.varpart <- varpart(data$Light_reactions,
                        ~ prec,
                        ~ gap,
                        #  ~ tavg,
                        #  ~ soilN,
                        #  ~ alpha,
                        ~ leaf_age,
                        #  ~ species,
                        data = data)
data.varpart
plot(data.varpart)

# permanovas between lineages #

data_photosynth <- data[data$lineage %in% c('Eucalyptus', 'Corymbia'),]
data_photosynth.dist <- data_photosynth[,c('ATP_synthase_chloroplastic', 'Cytochrome_b6f', 'Calvin_cycle', 'Rubisco', 'PSI', 'PSII', 'electron_transport_minATPsynth')]
data_photosynth.dist <- dist(data_photosynth.dist)
adonis(data_photosynth.dist ~ data_photosynth$lineage)
rm(data_photosynth, data_photosynth.dist)
