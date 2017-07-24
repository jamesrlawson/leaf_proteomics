include_leaf_N = TRUE

source('scripts/transformations.R')
source('scripts/prep_data_mg_per_mm2.R')

plot(data$calvin_cycle ~ data$N_per_area)
summary(lm(calvin_cycle ~ N_per_area, data))

data_out <- subset(data, N_per_area < 5000)
data_out <- subset(data_out, calvin_cycle < 30000)

plot(data_out$calvin_cycle ~ data_out$N_per_area)
summary(lm(calvin_cycle ~ N_per_area, data_out))

require(MuMIn)
options(na.action = "na.fail")

summary(lm(calvin_cycle ~ N_per_area, data_out))
bla <- lm(calvin_cycle ~ N_per_area + LMA_g_per_m2 + leaf_rad + tavg + log10(prec), data_out)
bla.dredge <- dredge(bla, extra = 'R^2')

plot(data_out$cytochrome_b6f ~ data_out$N_per_area)
summary(lm(cytochrome_b6f ~ N_per_area, data_out))
bla <- lm(cytochrome_b6f ~ N_per_area + LMA_g_per_m2 + leaf_rad + tavg + log10(prec), data_out)
bla.dredge <- dredge(bla, extra = 'R^2')

plot(data_out$Photosystems ~ data_out$N_per_area)
summary(lm(Photosystems ~ N_per_area, data_out))
bla <- lm(Photosystems ~ N_per_area + LMA_g_per_m2 + leaf_rad + tavg + log10(prec), data_out)
bla.dredge <- dredge(bla, extra = 'R^2')




require(visreg)

######### per leaf area protein amounts ########


source('scripts/prep_data_mg_per_mm2.R')

means <- data %>% group_by(ID) %>% dplyr::summarise(calv_mean = mean(calvin_cycle, na.rm=TRUE),
                                                    rub_mean = mean(Rubisco, na.rm=TRUE),
                                                    phot_mean = mean(Photosystems, na.rm=TRUE),
                                                    cyt_mean = mean(cytochrome_b6f, na.rm=TRUE)) 
data <-  full_join(data, means)
data <- distinct(data, cyt_mean, Narea_mean, rub_mean, calv_mean, phot_mean, leafrad_mean, .keep_all=TRUE)

data$log10prec <- log10(data$prec)

data$Jmax_per_Vcmax <- data$cyt_mean /data$rub_mean
Jmax_per_Vcmax.lm <- lm(Jmax_per_Vcmax ~ Narea_mean + LMA_mean + leafrad_mean + tavg + log10(prec), data)
Jmax_per_Vcmax.dredge <- dredge(Jmax_per_Vcmax.lm, extra = 'R^2', m.max=2)
Jmax_per_Vcmax.lm <- lm(Jmax_per_Vcmax ~ tavg + log10prec, data)

visreg::visreg2d(Jmax_per_Vcmax.lm, "tavg", "log10prec", plot.type="image" ,cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

summary(bla)

plot(calv_mean ~ Narea_mean, data)
summary(lm(calv_mean ~ Narea_mean, data))

plot(rub_mean ~ Narea_mean, data)
summary(lm(rub_mean ~ Narea_mean, data))

plot(cyt_mean ~ Narea_mean, data)
summary(lm(cyt_mean ~ Narea_mean, data))

plot(phot_mean ~ Narea_mean, data)
summary(lm(phot_mean ~ Narea_mean, data))

require(MuMIn)
options(na.action = "na.fail")

bla <- lm(phot_mean ~ Narea_mean + LMA_mean + leafrad_mean + tavg + log10(prec), data)
bla.dredge <- dredge(bla, extra = 'R^2')
bla.phot <- lm(phot_mean ~ Narea_mean + leafrad_mean, data)
summary(bla.phot)

bla <- lm(cyt_mean ~ Narea_mean + LMA_mean + leafrad_mean + tavg + log10(prec), data)
bla.dredge <- dredge(bla, extra = 'R^2')
bla.cyt <- lm(cyt_mean ~ Narea_mean + tavg, data)
summary(bla.cyt) 


bla <- lm(calv_mean ~ Narea_mean + LMA_mean + leafrad_mean + tavg + log10(prec), data)
bla.dredge <- dredge(bla, extra = 'R^2')
bla.calv <- lm(calv_mean ~ Narea_mean, data)
summary(bla.calv)
