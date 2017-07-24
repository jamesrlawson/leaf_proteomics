source('scripts/transformations.R')

require(visreg)

######### per leaf area protein amounts ########


source('scripts/prep_data_mg_per_mm2.R')

means <- data %>% group_by(ID) %>% dplyr::summarise(calv_mean = mean(calvin_cycle, na.rm=TRUE),
                                                    phot_mean = mean(Photosystems, na.rm=TRUE)) 
data <-  full_join(data, means)
data <- distinct(data, calv_mean, phot_mean, leafrad_mean, .keep_all=TRUE)

data$log10prec <- log10(data$prec)

require(MuMIn)
options(na.action = "na.fail")

global_calv_abs <- lm(calv_mean ~ log10prec * leafrad_mean * tavg, data)
global_calv_abs.dredge <- dredge(global_calv_abs, extra = c('R^2'))                        
global_calv_abs.dredge[1,]

calv_abs <- lm(calv_mean ~ tavg + log10prec, data)
summary(calv_abs)
visreg2d(calv_abs, "tavg", "log10prec", plot.type="image" ,cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

data$calv_frac <- data$calv_mean / data$LMA_mean / 1000
data$calv_frac

calv_abs_frac <- lm(calv_frac ~ tavg * log10prec, data)
summary(calv_abs_frac)

visreg2d(calv_abs_frac, "tavg", "log10prec", plot.type="image" ,cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)




global_phot_abs <- lm(phot_mean ~ log10prec * leafrad_mean * tavg, data)
global_phot_abs.dredge <- dredge(global_phot_abs, extra = c('R^2'))
global_phot_abs.dredge[1,]

phot_abs <- lm(phot_mean ~ tavg * log10prec, data)
summary(phot_abs)
visreg2d(phot_abs, "tavg", "log10prec", plot.type="image",cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

######### relative protein amounts ########


source('scripts/prep_data.R')

means <- data %>% group_by(ID) %>% dplyr::summarise(calv_mean = mean(calvin_cycle, na.rm=TRUE),
                                                    phot_mean = mean(Photosystems, na.rm=TRUE)) 
data <-  full_join(data, means)
data <- distinct(data, calv_mean, phot_mean, leafrad_mean, .keep_all=TRUE)
data$log10prec <- log10(data$prec)


global_calv_rel <- lm(calv_mean ~ log10prec * gap_mean * leafrad_mean * tavg, data)
global_calv_rel.dredge <- dredge(global_calv_rel, extra = c('R^2'))                        
global_calv_rel.dredge[1,]

calv_rel <- lm(calv_mean ~ leafrad_mean * tavg, data)
summary(calv_rel)
visreg2d(calv_rel, "tavg", "leafrad_mean", plot.type="image",cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

global_phot_rel <- lm(phot_mean ~ log10prec * tavg * gap_mean * leafrad_mean, data)
global_phot_rel.dredge <- dredge(global_phot_rel, extra = c('R^2'))
global_phot_rel.dredge[1,]     

phot_rel <- lm(phot_mean ~ tavg * leafrad_mean, data)
summary(phot_rel)
visreg2d(phot_rel, "tavg", "leafrad_mean", plot.type="image",cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
plot(leafrad_mean ~ tavg, data)

######### total protein ########

global_totalprot_rel <- lm(total_protein_mean ~ log10prec * leafrad_mean * tavg * LMA_mean, data)
global_totalprot_rel.dredge <- dredge(global_totalprot_rel, extra = c('R^2'))                        
global_totalprot_rel.dredge[1,]

total_prot <- lm(total_protein_mean ~ leafrad_mean * tavg, data)
summary(total_prot)
total_prot <- lm(total_protein_mean ~ leafrad_mean + tavg, data)
summary(total_prot)
visreg2d(total_prot, 'tavg', 'leafrad_mean', plot.type ='image')

total_prot <- lm(total_protein_mean ~ log10prec * tavg, data)
summary(total_prot)
total_prot <- lm(total_protein_mean ~ tavg + log10prec, data)
summary(total_prot)
visreg2d(total_prot, 'tavg', 'log10prec', plot.type ='image')

###


visreg2d(calv_rel, "tavg", "leafrad_mean", plot.type="image", xlim=c(0,30), ylim=c(30,100))
points(data$tavg, data$leafrad_mean)

visreg2d(total_prot, 'tavg', 'log10prec', plot.type ='image', cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

require(car)

phot_rel <- lm(phot_mean ~ leafrad_mean * tavg, data)
summary(phot_rel)

avPlots(phot_rel)
visreg2d(phot_rel, 'tavg', 'leafrad_mean', plot.type='persp',  col=c("red"), border="yellow",
         theta=120, phi=10, r=10)
visreg2d(phot_rel, 'tavg', 'leafrad_mean', plot.type='image', border="yellow")


#########

#MULTIPLE REGRESSIONS#

# depvars: per leaf area protein amounts
# predictors: tavg, log10prec, leafrad_mean, LMA_g_per_m2, total_protein_mean

include_leaf_N=TRUE

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

calv_global <- lm(calv_mean ~ tavg + log10prec + leafrad_mean + LMA_mean, data)
calv_global.dredge <- dredge(calv_global, extra = c('adjR^2'))
phot_global <- lm(phot_mean ~ tavg + log10prec + leafrad_mean + LMA_mean, data)
phot_global.dredge <- dredge(phot_global, extra = c('adjR^2'))

phot_3d <- lm(phot_mean ~ LMA_mean + total_protein_mean, data)
summary(phot_3d)
visreg2d(phot_3d, 'LMA_mean', 'total_protein_mean', plot.type='image', cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

calv_3d <- lm(calv_mean ~ LMA_mean + total_protein_mean, data)
summary(calv_3d)
visreg2d(calv_3d, 'LMA_mean', 'total_protein_mean', plot.type='image')


LMA_3d <- lm(LMA_mean ~ tavg * log10prec, data)
summary(LMA_3d)
visreg2d(LMA_3d, 'tavg', 'log10prec', plot.type='image', cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

LMA_3d <- lm(total_protein_mean ~ LMA_mean, data)
summary(LMA_3d)
visreg2d(LMA_3d, 'LMA_mean', 'log10prec', plot.type='image', cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)


data$leaf_protein_content <- data$total_protein_mean / data$LMA_mean / 1000

protein_per_LMA <- lm(leaf_protein_content ~ log10prec * tavg * leafrad_mean, data)
protein_per_LMA.dredge <- dredge(protein_per_LMA, extra = c('adjR^2'))

protein_per_LMA <- lm(leaf_protein_content ~ log10prec * tavg, data)
summary(protein_per_LMA)
visreg2d(protein_per_LMA, 'tavg', 'log10prec', plot.type='image', cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

blah <- lm(calv_mean ~ leaf_protein_content + LMA_mean, data)
summary(blah)
visreg2d(blah, 'LMA_mean', 'leaf_protein_content', plot.type='image', cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)


plot(total_protein_mean ~ leafrad_mean, data)
summary(lm(total_protein_mean ~ leafrad_mean, data))

agg_plot_save(proportion = FALSE, depvar = 'total_protein', indepvar = 'leafrad_mean', indepvarType = 'leafrad', logx = FALSE, 
              labs=c('Total leaf protein', 'Irradiance (MJ/m2/year)'), outDir = 'output/figures/20170704/tiff', goldenRatio = FALSE, fileType = 'tiff')


plot(total_protein_mean ~ log10prec, data)
abline(lm(total_protein_mean ~ log10prec, data))
summary(lm(total_protein_mean ~ log10prec, data))

agg_plot_save(proportion = FALSE, depvar = 'total_protein', indepvar = 'prec', indepvarType = 'standard', logx = TRUE, 
              labs=c('Total leaf protein', 'Mean annual precipitation (mm/year)'), outDir = 'output/figures/20170704/tiff', goldenRatio = FALSE, fileType = 'tiff')


plot(total_protein_mean ~ tavg, data)
abline(lm(total_protein_mean ~ tavg, data))
summary(lm(total_protein_mean ~ tavg, data))

agg_plot_save(proportion = FALSE, depvar = 'total_protein', indepvar = 'tavg', indepvarType = 'standard', logx = FALSE, 
              labs=c('Total leaf protein', 'Mean annual temperature (oC)'), outDir = 'output/figures/20170704/tiff', goldenRatio = FALSE, fileType = 'tiff')

agg_plot_save(proportion = FALSE, depvar = 'total_protein', indepvar = 'LMA_mean', indepvarType = 'LMA', logx = FALSE, 
              labs=c('Total leaf protein', 'Leaf mass per area (g/m2)'), outDir = 'output/figures/20170704/tiff', goldenRatio = FALSE, fileType = 'tiff')


source('scripts/prep_data.R')

means <- data %>% dplyr::group_by(ID) %>% dplyr::summarise(calv_mean = mean(calvin_cycle, na.rm=TRUE),
                                                           phot_mean = mean(Photosystems, na.rm=TRUE),
                                                           phot_SE = SE(Photosystems),
                                                           calv_SE = SE(calvin_cycle)) 
data <-  dplyr::full_join(data, means)
data <- dplyr::distinct(data, calv_mean, phot_mean, leafrad_mean, .keep_all=TRUE)

data$log10prec <- log10(data$prec)

data$leaf_protein_content <- data$total_protein_mean / data$LMA_mean / 1000

plot(phot_mean ~ leaf_protein_content, data)
summary(lm(phot_mean ~ leaf_protein_content, data))
abline(lm(phot_mean ~ leaf_protein_content, data))

plot(calv_mean ~ leaf_protein_content, data)
summary(lm(calv_mean ~ leaf_protein_content, data))
abline(lm(calv_mean ~ leaf_protein_content, data))


blah <- lm(Photosystems ~ N_per_area, data)
summary(blah)

