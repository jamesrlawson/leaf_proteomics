#########

#MULTIPLE REGRESSIONS#

# depvars: per leaf area protein amounts
# predictors: tavg, log10prec, leafrad_mean, LMA_g_per_m2, total_protein_mean

source('scripts/transformations.R')
source('scripts/prep_data_mg_per_mm2.R')

require(visreg)

means <- data %>% dplyr::group_by(ID) %>% dplyr::summarise(calv_mean = mean(calvin_cycle, na.rm=TRUE),
                                                           phot_mean = mean(Photosystems, na.rm=TRUE)) 
data <-  dplyr::full_join(data, means)
data <- dplyr::distinct(data, calv_mean, phot_mean, leafrad_mean, .keep_all=TRUE)

data$log10prec <- log10(data$prec)

require(MuMIn)
options(na.action = "na.fail")


# for total protein

total_prot <- lm(total_protein_mean ~ tavg + log10prec, data)
summary(total_prot)
visreg2d(total_prot, 'tavg', 'log10prec', plot.type ='image', cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

data$leaf_protein_content <- data$total_protein_mean / data$LMA_mean / 1000
protein_per_LMA <- lm(leaf_protein_content ~ log10prec * tavg, data)
summary(protein_per_LMA)
visreg2d(protein_per_LMA, 'tavg', 'log10prec', plot.type='image', cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

visreg2d_(protein_per_LMA, 'tavg', 'log10prec', plot.type='image', cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
plot(data$tavg, data$log10prec)

LMA <- lm(LMA_mean ~ tavg * log10prec, data)
summary(LMA)
visreg2d(LMA, 'tavg', 'log10prec', plot.type='image', cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

# for Calvin cycle protein only

calv <- lm(calv_mean ~ tavg + log10prec, data)
summary(calv)
visreg2d(calv, 'tavg', 'log10prec', plot.type ='image', cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

data$leaf_protein_content_calv <- data$calv_mean / data$LMA_mean / 1000
calv_per_LMA <- lm(leaf_protein_content_calv ~ log10prec * tavg, data)
summary(calv_per_LMA)
visreg2d(calv_per_LMA, 'tavg', 'log10prec', plot.type='image', cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

LMA <- lm(LMA_mean ~ tavg * log10prec, data)
summary(LMA)
visreg2d(LMA, 'tavg', 'log10prec', plot.type='image', cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)



# for photosystems proteins only

PS <- lm(phot_mean ~ tavg * log10prec, data)
summary(PS)

visreg2d(PS, 'tavg', 'log10prec', plot.type ='image', cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

data$leaf_protein_content_PS <- data$phot_mean / data$LMA_mean / 1000
PS_per_LMA <- lm(leaf_protein_content_PS ~ log10prec * tavg, data)
summary(PS_per_LMA)
visreg2d(PS_per_LMA, 'tavg', 'log10prec', plot.type='image', cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

LMA <- lm(LMA_mean ~ tavg * log10prec, data)
summary(LMA)
visreg2d(LMA, 'tavg', 'log10prec', plot.type='image', cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)


