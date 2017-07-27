??plot3D

require(plot3D)


source('scripts/prep_data.R')

bla <- points3D(data$LMA_g_per_m2, data$gap, data$Photosystems, theta=50, phi =20)
bla <- points3D(data$LMA_g_per_m2, data$gap, data$calvin_cycle, theta=50, phi =20)

source('scripts/prep_data_mg_per_mm2.R')

bla <- points3D(data$LMA_g_per_m2, data$leaf_rad, data$Photosystems, theta=50, phi =20)
bla <- points3D(data$LMA_g_per_m2, data$leaf_rad, data$calvin_cycle, theta=50, phi =20)


bla



###########################

library(car)
source('scripts/prep_data.R')


blah <- gather(data, key = 'funccat', value = 'protein_amount', Photosystems, calvin_cycle)

scatter3d(x = data$LMA_g_per_m2, y = data$Photosystems, z = data$leaf_rad, surface=TRUE, fit='linear')

scatter3d(x = data$LMA_g_per_m2, y = data$calvin_cycle, z = data$leaf_rad, surface=TRUE, fit='linear')

source('scripts/prep_data_mg_per_mm2.R')


blah <- gather(data, key = 'funccat', value = 'protein_amount', Photosystems, calvin_cycle)

scatter3d(x = data$LMA_g_per_m2, y = data$Photosystems, z = data$leaf_rad, surface=TRUE, fit='linear')

scatter3d(x = data$LMA_g_per_m2, y = data$calvin_cycle, z = data$leaf_rad, surface=TRUE, fit='linear')


##########################









source('scripts/prep_data.R')
data_means <- data %>% dplyr::group_by(ID) %>% dplyr::summarise(phot_mean = mean(Photosystems, na.rm=TRUE),
                                                                phot_SE = SE(Photosystems),
                                                                calv_mean = mean(calvin_cycle, na.rm=TRUE),
                                                                calv_SE = mean(calvin_cycle))
data_means <- merge(data, data_means)
data_means <- distinct(data_means, ID, .keep_all = TRUE) %>%
              gather(key = 'funccat', value = 'protein_amount', Photosystems, calvin_cycle)


scatter3d(x = data_means$LMA_mean, y = data_means$phot_mean, z = data_means$leafrad_mean, surface=TRUE, fit='linear')
scatter3d(x = data_means$LMA_mean, y = data_means$calv_mean, z = data_means$leafrad_mean, surface=TRUE, fit='linear')


 data_means <- data %>% dplyr::group_by(ID) %>% dplyr::summarise(phot_mean = mean(Photosystems, na.rm=TRUE),
                                                                phot_SE = SE(Photosystems),
                                                                calv_mean = mean(calvin_cycle, na.rm=TRUE),
                                                                calv_SE = mean(calvin_cycle))
data_means <- merge(data, data_means)
data_means <- distinct(data_means, ID, .keep_all = TRUE) %>%
  gather(key = 'funccat', value = 'protein_amount', Photosystems, calvin_cycle)


scatter3d(z = data_means$LMA_mean, y = data_means$phot_mean, x = data_means$leafrad_mean, surface=TRUE, fit='linear', model.summary=TRUE)
scatter3d(z = data_means$LMA_mean, y = data_means$calv_mean, x = data_means$leafrad_mean, surface=TRUE, fit='linear')

scatter3d(z = data_means$LMA_mean, y = data_means$protein_amount, x = data_means$leafrad_mean, groups = as.factor(data_means$funccat), surface=TRUE, fogtype = 'linear', parallel=FALSE)

scatter3d(protein_amount ~ leafrad_mean + LMA_mean | as.factor(funccat), data = data_means, surface=TRUE, fit = 'additive', fogtype = 'linear')



source('scripts/prep_data_mg_per_mm2.R')
dep_means <- group_by(data, ID) %>%
  dplyr::summarise(other_calv = mean(calvin_cycle, na.rm=TRUE) - mean(Rubisco, na.rm=TRUE),
                   rub = mean(Rubisco, na.rm=TRUE),
                   phot = mean(Photosystems_min_LHC, na.rm=TRUE),
                   LHC_ = mean(LHC, na.rm=TRUE)) %>%
  full_join(data, by = 'ID') %>%
  distinct(ID, .keep_all=TRUE)

a <- lm(rub ~ leafrad_mean, dep_means)
plot(rub ~ leafrad_mean, dep_means)
b <- lm(rub ~ tavg, dep_means)
plot(rub ~ tavg, dep_means)
ab <- lm(rub ~ leafrad_mean + tavg, dep_means)
axb <- lm(rub ~ leafrad_mean * tavg, dep_means)
MuMIn::AICc(a,b,ab,axb)



scatter3d(z = dep_means$tavg, y = dep_means$rub, x = dep_means$leafrad_mean, fit='additive')






source('scripts/prep_data_mg_per_mm2.R')

scatter3d(x = data$LMA_g_per_m2, y = data$calvin_cycle, z = data$leaf_rad, fit='smooth')


require(scatterplot3d)
source('scripts/prep_data.R')

with(data, {
  scatterplot3d(LMA_g_per_m2, data$leaf_rad, Photosystems,
                color = 'blue', pch=19,
                type = 'h', lty.hplot=3,
                main = 'Photosystem protein-trait-environment relationships',
                xlab = 'LMA', 
                ylab = 'irradiance',
                zlab = 'Photosystems proportion')
})