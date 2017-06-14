require(visreg)

######### per leaf area protein amounts ########


source('scripts/prep_data_mg_per_mm2.R')

means <- data %>% group_by(ID) %>% dplyr::summarise(calv_mean = mean(calvin_cycle, na.rm=TRUE),
                                                    phot_mean = mean(Photosystems, na.rm=TRUE)) 
data <-  full_join(data, means)
data <- distinct(data, calv_mean, phot_mean, gap_mean, .keep_all=TRUE)

data$log10prec <- log10(data$prec)

require(MuMIn)

global_calv_abs <- lm(calv_mean ~ log10prec * gap_mean * tavg, data)
global_calv_abs.dredge <- dredge(global_calv_abs, extra = c('R^2'))                        
global_calv_abs.dredge[1,]

calv_abs <- lm(calv_mean ~ tavg + log10prec, data)
summary(calv_abs)
visreg2d(calv_abs, "tavg", "log10prec", plot.type="image" ,cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)


global_phot_abs <- lm(phot_mean ~ log10prec * gap_mean * tavg, data)
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
data <- distinct(data, calv_mean, phot_mean, gap_mean, .keep_all=TRUE)
data$log10prec <- log10(data$prec)


global_calv_rel <- lm(calv_mean ~ tavg * gap_mean * tavg, data)
global_calv_rel.dredge <- dredge(global_calv_rel, extra = c('R^2'))                        
global_calv_rel.dredge[1,]

calv_rel <- lm(calv_mean ~ log10prec * tavg, data)
summary(calv_rel)
visreg2d(calv_rel, "tavg", "log10prec", plot.type="image",cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

global_phot_rel <- lm(phot_mean ~ log10prec * tavg * gap_mean, data)
global_phot_rel.dredge <- dredge(global_phot_rel, extra = c('R^2'))
global_phot_rel.dredge[1,]     

phot_rel <- lm(phot_mean ~ log10prec * tavg, data)
summary(phot_rel)
visreg2d(phot_rel, "tavg", "log10prec", plot.type="image",cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
plot(gap_mean ~ tavg, data)

######### total protein ########

global_totalprot_rel <- lm(total_protein_mean ~ log10prec * gap_mean * tavg, data)
global_totalprot_rel.dredge <- dredge(global_calv_rel, extra = c('R^2'))                        
global_totalprot_rel.dredge[1,]

total_prot <- lm(total_protein_mean ~ gap_mean * tavg, data)
summary(total_prot)
total_prot <- lm(total_protein_mean ~ gap_mean + tavg, data)
summary(total_prot)
visreg2d(total_prot, 'tavg', 'gap_mean', plot.type ='image')

total_prot <- lm(total_protein_mean ~ log10prec * tavg, data)
summary(total_prot)
total_prot <- lm(total_protein_mean ~ log10prec + tavg, data)
summary(total_prot)
visreg2d(total_prot, 'tavg', 'log10prec', plot.type ='image')

###


visreg2d(calv_rel, "tavg", "gap_mean", plot.type="image", xlim=c(0,30), ylim=c(30,100))
points(data$tavg, data$gap_mean)

visreg2d(total_prot, 'tavg', 'log10prec', plot.type ='image', cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)

require(car)

phot_rel <- lm(phot_mean ~ leafrad_mean * tavg, data)
summary(phot_rel)

avPlots(phot_rel)
visreg2d(phot_rel, 'tavg', 'leafrad_mean', plot.type='persp',  col=c("red"), border="yellow",
         theta=120, phi=10, r=10)
visreg2d(phot_rel, 'tavg', 'leafrad_mean', plot.type='image', border="yellow")
