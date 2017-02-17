

require(lme4)
#source('scripts/transformations.R')
source('scripts/prep_data_.R')


data_means <- data %>% group_by(ID) %>% summarise(mean = mean(Calvin_cycle, na.rm=TRUE),
                                                          SE = SE(Calvin_cycle))
data_means <- merge(data, data_means)
data_means <- distinct(data_means, ID, .keep_all = TRUE)

#total_protein_means <- data %>% group_by(ID) %>% summarise(total_protein_mean = mean(total_protein, na.rm=TRUE),
#                                                           total_protein_SE = SE(total_protein))
#data_means <- merge(total_protein_means, data_means)


i <- varpart(data_means$mean,
             ~log10(leafrad_mean),
      #       ~log10(gap_mean),
      #       ~log10(prec),
            ~tavg,
      #       ~total_protein_mean,
       #      ~leaf_age,
             data = data_means)
i 
plot(i)

summary(lm(mean ~ log10(leafrad_mean), data_means))
summary(lm(mean ~ log10(gap_mean), data_means))
summary(lm(mean ~ log10(prec), data_means))
summary(lm(mean ~ leaf_age, data_means))
summary(lm(mean ~ tavg, data_means))

h <- lm(mean ~ gap_mean, data_means)
i <- lm(mean ~ total_protein_mean, data_means)
j <- lm(mean ~ total_protein_mean + gap_mean, data_means)
k <- lm(mean ~ total_protein_mean * gap_mean, data_means)
AICc(h,i,j,k)
dredge(k)
summary(k)

j <- lm(scale(mean) ~ scale(total_protein_mean) + scale(gap_mean), data_means)
summary(j)

k <-  lm(scale(mean) ~ scale(total_protein_mean) * scale(gap_mean), data_means)
#summary(k)

h <- lm(scale(mean) ~ scale(leafrad_mean), data_means)


errorbar_width <- (max(data_means$leafrad_mean, na.rm=TRUE) - min(data_means$leafrad_mean, na.rm=TRUE)) / 50

g <- ggplot(data_means, aes(y = mean, x = log10(gap_mean))) + geom_point() + geom_smooth(method = 'lm') +xlab('gap mean') +ylab('Calvin_cycle mg/m2')
g <- g + geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), alpha = 0.8)

g

plot(tavg ~ total_protein_mean, data_means)
abline(lm(tavg ~ total_protein_mean, data_means))
plot(mean ~ total_protein_mean, data_means, ylab = 'Calvin_cycle proteins (mg/m2)')
abline(lm(mean ~ total_protein_mean, data_means))
summary(lm(mean ~ total_protein_mean, data_means))


plot(total_protein_mean ~ log10(leafrad_mean), data_means)
abline(lm(total_protein_mean ~ log10(leafrad_mean), data_means))
summary(lm(total_protein_mean ~ log10(leafrad_mean), data_means))

plot(tavg ~ log10(leafrad_mean), data_means)
cor.test(data_means$tavg, log10(data_means$leafrad_mean))

plot(log10(prec) ~ log10(leafrad_mean), data_means)
cor.test(log10(data_means$prec), log10(data_means$leafrad_mean))

plot(log10(prec) ~ tavg, data_means)
cor.test(log10(data_means$prec), data_means$tavg)


plot(log10(prec) ~ log10(gap_mean), data_means)
cor.test(log10(data_means$prec), data_means$gap_mean)

plot(tavg ~ log10(gap_mean), data_means)
cor.test(data_means$tavg, log10(data_means$gap_mean))

plot(log10(gap_mean) ~ log10(leafrad_mean), data_means)
cor.test(log10(data_means$gap_mean), log10(data_means$leafrad_mean))



# modelled reduction in Calvin_cycle proteins across gradients

p <- lm(mean ~ log10(leafrad_mean), data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(log10(data_means$leafrad_mean))
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(log10(data_means$leafrad_mean))
  
(y - x)/y # modelled reduction in Calvin_cycle protein across gradient of gap_mean
  
p <- lm(mean ~ log10(gap_mean), data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(log10(data_means$gap_mean))
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(log10(data_means$gap_mean))

(y - x)/y # modelled reduction in Calvin_cycle protein across gradient of leafrad_mean





##### rubisco standardised by Calvin_cycle amounts #####
# all being equal, calvin cycle amounts should stay constant for a given rubisco activity
# so to respond to reduced CO2 availability at lower rainfall, rubisco proportion of calvin cycle should increase

source('scripts/prep_data_.R')

data$rubisco_stand <- data$Rubisco / data$Calvin_cycle

data_means <- data %>% group_by(ID) %>% summarise(mean = mean(rubisco_stand, na.rm=TRUE),
                                                  SE = SE(rubisco_stand))
data_means <- merge(data, data_means)
data_means <- distinct(data_means, ID, .keep_all = TRUE)

plot(rubisco_stand ~ gap, data_means)
abline(lm(rubisco_stand ~ gap, data_means))
summary(lm(rubisco_stand ~ gap, data_means))
