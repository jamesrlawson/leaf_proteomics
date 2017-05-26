

require(lme4)
source('scripts/transformations.R')
source('scripts/prep_data.R')


data_means <- data %>% dplyr::group_by(ID) %>% dplyr::summarise(mean = mean(Photosystems, na.rm=TRUE),
                                                          SE = SE(Photosystems),
                                                          N_per_area_mean = mean(N_per_area, na.rm=TRUE),
                                                          LMA_mean = mean(LMA_g_per_m2, na.rm=TRUE))
data_means <- merge(data, data_means)
data_means <- distinct(data_means, ID, .keep_all = TRUE)

#total_protein_means <- data %>% group_by(ID) %>% summarise(total_protein_mean = mean(total_protein, na.rm=TRUE),
#                                                           total_protein_SE = SE(total_protein))
#data_means <- merge(total_protein_means, data_means)


i <- varpart(data_means$mean,
             ~leafrad_mean,
             ~gap_mean,
            ~log10(prec),
       #     ~tavg,
     # ~pdmt,
      #       ~total_protein_mean,
             ~leaf_age,
             data = data_means)
i 
plot(i)

summary(lm(mean ~ gap_mean, data_means))
summary(lm(mean ~ leafrad_mean, data_means))
summary(lm(mean ~ log10(prec), data_means))
summary(lm(mean ~ pdmt, data_means))

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

g <- ggplot(data_means, aes(y = mean, x = gap_mean)) + geom_point() + geom_smooth(method = 'lm') +xlab('gap mean') +ylab('Photosystems mg/m2')
g <- g + geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), alpha = 0.8)

g

plot(tavg ~ total_protein_mean, data_means)
abline(lm(tavg ~ total_protein_mean, data_means))
plot(mean ~ total_protein_mean, data_means, ylab = 'Photosystems proteins (mg/m2)')
abline(lm(mean ~ total_protein_mean, data_means))
summary(lm(mean ~ total_protein_mean, data_means))


plot(total_protein_mean ~ leafrad_mean, data_means)
abline(lm(total_protein_mean ~ leafrad_mean, data_means))
summary(lm(total_protein_mean ~ leafrad_mean, data_means))

plot(tavg ~ leafrad_mean, data_means)
cor.test(data_means$tavg, data_means$leafrad_mean)

plot(log10(prec) ~ leafrad_mean, data_means)
cor.test(log10(data_means$prec), data_means$leafrad_mean)

plot(log10(prec) ~ tavg, data_means)
cor.test(log10(data_means$prec), data_means$tavg)


plot(log10(prec) ~ gap_mean, data_means)
cor.test(log10(data_means$prec), data_means$gap_mean)

plot(tavg ~ gap_mean, data_means)
cor.test(data_means$tavg, data_means$gap_mean)

plot(gap_mean ~ leafrad_mean, data_means)
cor.test(data_means$gap_mean, data_means$leafrad_mean)



# modelled change in total protein across gradients


p <- lm(total_protein_mean ~ log10(prec), data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(log10(data_means$prec))
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(log10(data_means$prec))

(x - y)/x

p <- lm(total_protein_mean ~ tavg, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$tavg)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$tavg)

(x - y)/x



# modelled change in Photosystems proteins across env gradients

# leafrad_mean

p <- lm(mean ~ leafrad_mean, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$leafrad_mean)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$leafrad_mean)
  
(x - y)/x

# modelled change in Photosystems protein across gradient of gap_mean
  
p <- lm(mean ~ gap_mean, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$gap_mean)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$gap_mean)

(x - y)/x


# modelled change in Photosystems protein across gradient of tavg

p <- lm(mean ~ tavg, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$tavg)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$tavg)

(x - y)/x

# modelled change in Photosystems protein across gradient of prec

p <- lm(mean ~ log10(prec), data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(log10(data_means$prec))
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(log10(data_means$prec))

(x - y)/x

# modelled change in Photosystems protein across gradient of prec in driest month

p <- lm(mean ~ pdmt, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$pdmt)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$pdmt)

(x - y)/x


# correlation table

cors <- dplyr::select(data_means, prec, pdmt, tavg, gap_mean, leafrad_mean)

cors$prec <- log10(cors$prec)
cors$gap_mean <- log10(cors$gap_mean)
cors$leafrad_mean <- log10(cors$leafrad_mean)

envcors <- data.frame(cor(cors))

write_csv(envcors, 'output/env_cors.csv')

envcors.pca <- prcomp(envcors, scale=TRUE, center=TRUE)

plot(envcors.pca)
summary(envcors.pca)
biplot(envcors.pca)


# Photosystems with PSI and PSII added

data_means <- data %>% dplyr::group_by(ID) %>% dplyr::summarise(mean = mean(Photosystems, na.rm=TRUE),
                                                                SE = SE(Photosystems),
                                                                meanPSI = mean(PSI, na.rm=TRUE),
                                                                meanPSII = mean(PSII, na.rm=TRUE))
data_means <- merge(data, data_means)
data_means <- distinct(data_means, ID, .keep_all = TRUE)


g <- ggplot(data_means, aes(y = mean, x = gap_mean)) + geom_point() + geom_smooth(method = 'lm', se=FALSE) +xlab('gap mean') +ylab('Photosystems mg/m2')
g <- g + geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), alpha = 0.8)
g <- g + geom_smooth(aes(y = meanPSI, x = gap_mean), method = 'lm', se = FALSE, colour = 'red')
g <- g + geom_smooth(aes(y = meanPSII, x = gap_mean), method = 'lm', se = FALSE, colour = 'green')

g


## Photosystems vs total protein, colour scaled

source('scripts/prep_data_mg_per_mm2.R')

x <- ggplot(data, aes(x = total_protein, y = Photosystems)) + geom_point(aes(colour = gap)) + scale_colour_gradientn('canopy openness', colours=rainbow(2))
x <- x + ylab('Photosystem protein abundance (mg/mm2)') + xlab('Total protein abundance (mm2/mg)')
x <- x + theme_bw() + theme(axis.line = element_line(colour = "black"),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.border = element_blank(),
                            panel.background = element_blank(),
                            legend.position="top"
)
x

x <- ggplot(data, aes(x = total_protein, y = Calvin_cycle)) + geom_point(aes(colour = log10(pdmt))) + scale_colour_gradientn('Mean precip in driest month (log10)', colours=rainbow(2))
x <- x + ylab('Calvin cycle protein abundance (mg/mm2)') + xlab('Total protein abundance (mm2/mg)')
x <- x + theme_bw() + theme(axis.line = element_line(colour = "black"),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.border = element_blank(),
                            panel.background = element_blank(),
                            legend.position="top"
)
x

# total protein varparts

total_protein.varpart <- varpart(data$total_protein,
                                 ~N_per_area,
                                 ~LMA_g_per_m2,
                                 ~log10(pwmt),
                                 ~tavg,
                                 data = data)
total_protein.varpart
plot(total_protein.varpart)

bla <- lm(total_protein ~ LMA_g_per_m2 + tavg, data)
bla1 <- lm(total_protein ~ LMA_g_per_m2 * tavg, data)
bla2 <- lm(total_protein ~ LMA_g_per_m2, data)
bla3 <- lm(total_protein ~ tavg, data)
bla4 <- lm(total_protein ~ tavg + log10(pwmt), data) 
bla5 <- lm(total_protein ~ LMA_g_per_m2 + tavg + log10(pwmt), data) 
bla6 <- lm(total_protein ~ LMA_g_per_m2 + tavg * log10(pwmt), data) 
bla7 <- lm(total_protein ~ LMA_g_per_m2 * tavg * log10(pwmt), data) 
  
AICc(bla,bla1,bla2,bla3,bla4,bla5,bla6,bla7)

summary(lm(scale(total_protein_mean) ~ scale(LMA_g_per_m2) + scale(tavg) + scale(log(pwmt)), data))
summary(lm(scale(total_protein_mean) ~ scale(LMA_g_per_m2) +  scale(log(pwmt)) * scale(tavg) , data))


total_protein.varpart <- varpart(data_means$total_protein_mean,
                                 ~N_per_area_mean,
                                 ~LMA_mean,
                                 ~log10(pwmt),
                                 ~tavg,
                                 data = data_means)
total_protein.varpart
plot(total_protein.varpart)

cor(data_means$total_protein_mean, data_means$N_per_area_mean)
cor(data_means$total_protein_mean, data_means$LMA_mean)

##### rubisco standardised by Photosystems amounts #####
# all being equal, calvin cycle amounts should stay constant for a given rubisco activity
# so to respond to reduced CO2 availability at lower rainfall, rubisco proportion of calvin cycle should increase

source('scripts/prep_data_mg_per_mm2.R')

data$rubisco_stand <- data$Rubisco / data$Photosystems

data_means <- data %>% group_by(ID) %>% summarise(mean = mean(rubisco_stand, na.rm=TRUE),
                                                  SE = SE(rubisco_stand))
data_means <- merge(data, data_means)
data_means <- distinct(data_means, ID, .keep_all = TRUE)

plot(rubisco_stand ~ gap, data_means)
abline(lm(rubisco_stand ~ gap, data_means))
summary(lm(rubisco_stand ~ gap, data_means))





