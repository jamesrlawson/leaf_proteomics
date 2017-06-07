install.packages('visreg')
require(visreg)

source('scripts/prep_data.R')

means <- data %>% group_by(ID) %>% dplyr::summarise(calv_mean = mean(calvin_cycle, na.rm=TRUE),
                                                    phot_mean = mean(Photosystems, na.rm=TRUE)) 
  data <-  full_join(data, means)
  data <- distinct(data, calv_mean, phot_mean, gap_mean, .keep_all=TRUE)

fit_calv_rel <- lm(calv_mean ~ tavg * gap_mean, data=data)
summary(fit_calv_rel)
fit_phot_rel <- lm(phot_mean ~ tavg * gap_mean, data=data)
summary(fit_phot_rel)

visreg2d(fit_calv_rel, "tavg", "gap_mean", plot.type="image")
visreg2d(fit_phot_rel, "tavg", "gap_mean", plot.type="image")

source('scripts/prep_data_mg_per_mm2.R')

means <- data %>% group_by(ID) %>% dplyr::summarise(calv_mean = mean(calvin_cycle, na.rm=TRUE),
                                                    phot_mean = mean(Photosystems, na.rm=TRUE)) 
data <-  full_join(data, means)
data <- distinct(data, calv_mean, phot_mean, gap_mean, .keep_all=TRUE)

fit_calv_abs <- lm(calvin_cycle ~ tavg + prec, data=data)
summary(fit_calv_abs)
fit_phot_abs <- lm(Photosystems ~ tavg * prec, data=data)
summary(fit_phot_abs)

visreg2d(fit_calv_abs, "tavg", "prec", plot.type="image")
visreg2d(fit_phot_abs, "tavg", "prec", plot.type="image")

fit_total <- lm(total_protein_mean ~ scale(tavg) * scale(prec), data)
summary(fit_total)

visreg2d(fit_total, xvar = 'tavg', yvar = 'prec', plot.type='image')
