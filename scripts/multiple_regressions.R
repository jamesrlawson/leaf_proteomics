

require(lme4)
source('scripts/transformations.R')
source('scripts/prep_data_.R')


data_means <- data %>% group_by(ID) %>% summarise(mean = mean(Photosystems, na.rm=TRUE),
                                                          SE = SE(Photosystems))
data_means <- merge(data, data_means)
data_means <- distinct(data_means, ID, .keep_all = TRUE)

#total_protein_means <- data %>% group_by(ID) %>% summarise(total_protein_mean = mean(total_protein, na.rm=TRUE),
#                                                           total_protein_SE = SE(total_protein))
#data_means <- merge(total_protein_means, data_means)


i <- varpart(data_means$mean,
             ~log10(leafrad_mean),
             ~tavg,
       #      ~total_protein_mean,
             ~leaf_age,
             ~log10(gap_mean),
             data = data_means)
i 
plot(i)

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

g <- ggplot(data_means, aes(y = mean, x = log10(gap_mean))) + geom_point() + geom_smooth(method = 'lm') +xlab('gap mean') +ylab('Photosystems mg/m2')
g <- g + geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), alpha = 0.8)

g

plot(tavg ~ total_protein_mean, data_means)
abline(lm(tavg ~ total_protein_mean, data_means))
plot(mean ~ total_protein_mean, data_means, ylab = 'Photosystems proteins (mg/m2)')
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



# modelled reduction in photosystems proteins across gradients

p <- lm(mean ~ log10(leafrad_mean), data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(log10(data_means$leafrad_mean))
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(log10(data_means$leafrad_mean))
  
(y - x)/y # modelled reduction in Photosystems protein across gradient of gap_mean
  
p <- lm(mean ~ log10(gap_mean), data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(log10(data_means$gap_mean))
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(log10(data_means$gap_mean))

(y - x)/y # modelled reduction in Photosystems protein across gradient of leafrad_mean


