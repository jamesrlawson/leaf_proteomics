#WUE

include_photosynthesis=TRUE
source('scripts/transformations_photomax.R')
source('scripts/prep_data_mg_per_mm2.R')
#source('scripts/prep_data.R')

means <- data %>% group_by(ID) %>% dplyr::summarise(calv_mean = mean(calvin_cycle, na.rm=TRUE),
                                                    calv_SE = SE(calvin_cycle),
                                                    rubisco_mean = mean(Rubisco, na.rm=TRUE),
                                                    rubisco_SE = SE(Rubisco),
                                                    phot_mean = mean(Photosystems, na.rm=TRUE),
                                                    Cond_mean = mean(Cond, na.rm=TRUE),
                                                    photo_max_mean = mean(photo_max, na.rm=TRUE),
                                                    WUE_SE = SE(photo_max/Cond))
data <-  full_join(data, means)
data <- distinct(data, calv_mean, phot_mean, leafrad_mean, Cond_mean, photo_max_mean, WUE_SE, .keep_all=TRUE)

data$log10prec <- log10(data$prec)
data$WUE_mean <- data$photo_max_mean/data$Cond_mean

data <- data[data$WUE_mean > 1,]

plot(data$WUE_mean ~ data$calv_mean)
summary(lm(data$WUE_mean ~ data$calv_mean))
abline(lm(data$WUE_mean ~ data$calv_mean))

plot(data$WUE_mean ~ data$total_protein_mean)
summary(lm(data$WUE_mean ~ data$total_protein_mean))
abline(lm(data$WUE_mean ~ data$total_protein_mean))

plot(data$WUE_mean ~ data$LMA_mean)
summary(lm(data$WUE_mean ~ data$LMA_mean))
abline(lm(data$WUE_mean ~ data$LMA_mean))


ggplot(data, aes(x = calv_mean, y = WUE_mean)) + geom_point(size=2) +
  geom_smooth(method='lm', se=FALSE, colour='black', size=0.7) +
  geom_errorbar(aes(ymin = WUE_mean - WUE_SE, ymax = WUE_mean + WUE_SE)) +
  geom_errorbarh(aes(xmin = calv_mean - calv_SE, xmax = calv_mean + calv_SE)) +
  xlab('Calvin cycle (mg/m2)') + ylab('Photosynthetic Water-use Efficiency (Amax/Gs)') +
 # xlim(0,max(data$calv_mean)) +
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.line = element_line(colour = "black"),
        text = element_text(size = 18))

ggplot(data, aes(x = rubisco_mean, y = WUE_mean)) + geom_point(size=2) +
  geom_smooth(method='lm', se=FALSE, colour='black', size=0.7) +
  geom_errorbar(aes(ymin = WUE_mean - WUE_SE, ymax = WUE_mean + WUE_SE)) +
  geom_errorbarh(aes(xmin = rubisco_mean - rubisco_SE, xmax = rubisco_mean + rubisco_SE)) +
  xlab('Rubisco (mg/m2)') + ylab('Photosynthetic Water-use Efficiency (Amax/Gs)') +
  # xlim(0,max(data$calv_mean)) +
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.line = element_line(colour = "black"),
        text = element_text(size = 18))
                    
bla <- lm(calv_mean ~ tavg + log10prec, data)
summary(bla)
visreg2d(bla, "tavg", "log10prec", plot.type="image",cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)


include_photosynthesis=TRUE
source('scripts/transformations.R')
source('scripts/prep_data_mg_per_mm2.R')

plot(data$photo_amb)

data <- data[data$photo_amb > 10,]

bla1 <- lm(photo_amb ~ Cond + calvin_cycle, data)
summary(bla1)
bla2 <- lm(photo_amb ~ Cond * calvin_cycle, data)
summary(bla2)
bla3 <- lm(photo_amb ~ Cond, data)
summary(bla3)
bla4 <- lm(photo_amb ~ calvin_cycle, data)
summary(bla4)

#require(MuMIn)
AICc(bla1,bla2,bla3,bla4)

require(visreg)
visreg2d(bla2, "Cond", "calvin_cycle", plot.type="image",cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
plot(photo_amb ~ calvin_cycle, data)
plot(photo_amb ~ Cond, data)


include_photosynthesis=TRUE
source('scripts/transformations_photomax.R')
source('scripts/prep_data_mg_per_mm2.R')

plot(data$photo_max)

data <- data[data$photo_max > 10,]

bla1 <- lm(photo_max_mean ~ Cond_mean + calv_mean, data)
summary(bla1)
bla2 <- lm(photo_max_mean ~ Cond_mean * calv_mean, data)
summary(bla2)
bla3 <- lm(photo_max_mean ~ Cond_mean, data)
summary(bla3)
bla4 <- lm(photo_max_mean ~ calv_mean, data)
summary(bla4)

#require(MuMIn)
AICc(bla1,bla2,bla3,bla4)

require(visreg)
visreg2d(bla1, "Cond_mean", "calv_mean", plot.type="image",cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
plot(photo_max ~ calvin_cycle, data)
plot(photo_max ~ Cond, data)
plot(calv_mean ~ Cond, data)
ggplot(data, aes(y = calvin_cycle, x = Cond, colour = photo_max)) + geom_point()  + scale_colour_gradientn('Mean annual irradiance @ leaf (MJ/m2/yr)', colours=rainbow(2))
