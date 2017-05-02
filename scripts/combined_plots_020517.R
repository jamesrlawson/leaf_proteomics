
source('scripts/functions.R')


agg_plot_save_combined(indepvar = 'prec', logx = TRUE, proportion = FALSE, indepvarType = 'standard', 
              labs =  c('Protein amount', 'Mean annual precip. (m/yr)'), outDir = 'output/figures/20170501/svg', goldenRatio = FALSE, fileType = 'svg')

agg_plot_save_combined(indepvar = 'prec', logx = TRUE, proportion = TRUE, indepvarType = 'standard', 
              labs =  c('Protein amount', 'Mean annual precip. (m/yr)'), outDir = 'output/figures/20170501/svg', goldenRatio = FALSE, fileType = 'svg')


agg_plot_save_combined(indepvar = 'gap_mean', logx = FALSE, proportion = FALSE, indepvarType = 'gap', 
              labs =  c('Protein amount', 'Canopy openness (%)'), outDir = 'output/figures/20170501/svg', goldenRatio = FALSE, fileType = 'svg')

agg_plot_save_combined(indepvar = 'gap_mean', logx = FALSE, proportion = TRUE, indepvarType = 'gap', 
              labs =  c('Protein amount', 'Canopy openness (%)'), outDir = 'output/figures/20170501/svg', goldenRatio = FALSE, fileType = 'svg')


agg_plot_save_combined(indepvar = 'leafrad_mean', logx = FALSE, proportion = FALSE, indepvarType = 'leafrad', 
                       labs =  c('Protein amount', 'Irradiance (MJ/m2/year)'), outDir = 'output/figures/20170501/svg', goldenRatio = FALSE, fileType = 'svg')

agg_plot_save_combined(indepvar = 'leafrad_mean', logx = FALSE, proportion = TRUE, indepvarType = 'leafrad', 
                       labs =  c('Protein amount', 'Irradiance (MJ/m2/year)'), outDir = 'output/figures/20170501/svg', goldenRatio = FALSE, fileType = 'svg')

agg_plot_save_combined(indepvar = 'tavg', logx = FALSE, proportion = FALSE, indepvarType = 'standard', 
                       labs =  c('Protein amount', 'Mean annual temperature (oC)'), outDir = 'output/figures/20170501/svg', goldenRatio = FALSE, fileType = 'svg')

agg_plot_save_combined(indepvar = 'tavg', logx = FALSE, proportion = TRUE, indepvarType = 'standard', 
                       labs =  c('Protein amount', 'Mean annual temperature (oC)'), outDir = 'output/figures/20170501/svg', goldenRatio = FALSE, fileType = 'svg')

agg_plot_save_combined(indepvar = 'LMA_mean', logx = FALSE, proportion = FALSE, indepvarType = 'LMA', 
                       labs =  c('Protein amount', 'LMA (g/m2)'), outDir = 'output/figures/20170501/svg', goldenRatio = FALSE, fileType = 'svg')

agg_plot_save_combined(indepvar = 'LMA_mean', logx = FALSE, proportion = TRUE, indepvarType = 'LMA', 
                       labs =  c('Protein amount', 'LMA (g/m2)'), outDir = 'output/figures/20170501/svg', goldenRatio = FALSE, fileType = 'svg')


agg_plot_save_combined(indepvar = 'total_protein_mean', logx = FALSE, proportion = FALSE, indepvarType = 'total_protein', 
                       labs =  c('Protein amount', 'total protein (mg/m2)'), outDir = 'output/figures/20170501/svg', goldenRatio = FALSE, fileType = 'svg', total_prot = TRUE)

agg_plot_save_combined(indepvar = 'total_protein_mean', logx = FALSE, proportion = TRUE, indepvarType = 'total_protein', 
                       labs =  c('Protein amount', 'total protein (mg/m2)'), outDir = 'output/figures/20170501/svg', goldenRatio = FALSE, fileType = 'svg', total_prot = TRUE)


# modelled changes in protein amounts across env gradients

# proportional data first
#source('scripts/transformations.R')
source('scripts/prep_data.R')

data_means <- data %>% dplyr::group_by(ID) %>% dplyr::summarise(phot_mean = mean(Photosystems, na.rm=TRUE),
                                                                phot_SE = SE(Photosystems),
                                                                calv_mean = mean(calvin_cycle, na.rm=TRUE),
                                                                calv_SE = mean(calvin_cycle))
data_means <- merge(data, data_means)
data_means <- distinct(data_means, ID, .keep_all = TRUE)

# total protein

p <- lm(total_protein_mean ~ log10(prec), data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(log10(data_means$prec))
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(log10(data_means$prec))

total_protein_propChange_prec <- (x - y)/x

p <- lm(total_protein_mean ~ leafrad_mean, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$leafrad_mean)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$leafrad_mean)

total_protein_propChange_leafrad <- (x - y)/x

p <- lm(total_protein_mean ~ gap_mean, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$gap_mean)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$gap_mean)

total_protein_propChange_gap <- (x - y)/x

p <- lm(total_protein_mean ~ tavg, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$tavg)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$tavg)

total_protein_propChange_tavg <- (x - y)/x

p <- lm(total_protein_mean ~ LMA_mean, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$LMA_mean)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$LMA_mean)

total_protein_propChange_LMA <- (x - y)/x

# photosystems

p <- lm(phot_mean ~ log10(prec), data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(log10(data_means$prec))
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(log10(data_means$prec))

phot_propChange_prec <- (x - y)/x

p <- lm(phot_mean ~ leafrad_mean, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$leafrad_mean)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$leafrad_mean)

phot_propChange_leafrad <- (x - y)/x

p <- lm(phot_mean ~ gap_mean, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$gap_mean)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$gap_mean)

phot_propChange_gap <- (x - y)/x

p <- lm(phot_mean ~ tavg, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$tavg)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$tavg)

phot_propChange_tavg <- (x - y)/x

p <- lm(phot_mean ~ LMA_mean, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$LMA_mean)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$LMA_mean)

phot_propChange_LMA <- (x - y)/x

p <- lm(phot_mean ~ total_protein_mean, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$total_protein_mean)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$total_protein_mean)

phot_propChange_totalprot <- (x - y)/x

# calvin cycle

p <- lm(calv_mean ~ log10(prec), data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(log10(data_means$prec))
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(log10(data_means$prec))

calv_propChange_prec <- (x - y)/x

p <- lm(calv_mean ~ leafrad_mean, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$leafrad_mean)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$leafrad_mean)

calv_propChange_leafrad <- (x - y)/x

p <- lm(calv_mean ~ gap_mean, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$gap_mean)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$gap_mean)

calv_propChange_gap <- (x - y)/x

p <- lm(calv_mean ~ tavg, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$tavg)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$tavg)

calv_propChange_tavg <- (x - y)/x

p <- lm(calv_mean ~ LMA_mean, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$LMA_mean)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$LMA_mean)

calv_propChange_LMA <- (x - y)/x

p <- lm(calv_mean ~ total_protein_mean, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$total_protein_mean)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$total_protein_mean)

calv_propChange_totalprot <- (x - y)/x
  

  # mg/m2 data next
  
source('scripts/prep_data_mg_per_mm2.R')

data_means <- data %>% dplyr::group_by(ID) %>% dplyr::summarise(phot_mean = mean(Photosystems, na.rm=TRUE),
                                                                phot_SE = SE(Photosystems),
                                                                calv_mean = mean(calvin_cycle, na.rm=TRUE),
                                                                calv_SE = mean(calvin_cycle))
data_means <- merge(data, data_means)
data_means <- distinct(data_means, ID, .keep_all = TRUE)

# total protein

p <- lm(total_protein_mean ~ log10(prec), data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(log10(data_means$prec))
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(log10(data_means$prec))

total_protein_absChange_prec <- (x - y)/x

p <- lm(total_protein_mean ~ leafrad_mean, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$leafrad_mean)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$leafrad_mean)

total_protein_absChange_leafrad <- (x - y)/x

p <- lm(total_protein_mean ~ gap_mean, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$gap_mean)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$gap_mean)

total_protein_absChange_gap <- (x - y)/x

p <- lm(total_protein_mean ~ tavg, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$tavg)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$tavg)

total_protein_absChange_tavg <- (x - y)/x

p <- lm(total_protein_mean ~ LMA_mean, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$LMA_mean)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$LMA_mean)

total_protein_absChange_LMA <- (x - y)/x

# photosystems

p <- lm(phot_mean ~ log10(prec), data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(log10(data_means$prec))
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(log10(data_means$prec))

phot_absChange_prec <- (x - y)/x

p <- lm(phot_mean ~ leafrad_mean, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$leafrad_mean)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$leafrad_mean)

phot_absChange_leafrad <- (x - y)/x

p <- lm(phot_mean ~ gap_mean, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$gap_mean)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$gap_mean)

phot_absChange_gap <- (x - y)/x

p <- lm(phot_mean ~ tavg, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$tavg)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$tavg)

phot_absChange_tavg <- (x - y)/x

p <- lm(phot_mean ~ LMA_mean, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$LMA_mean)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$LMA_mean)

phot_absChange_LMA <- (x - y)/x

p <- lm(phot_mean ~ total_protein_mean, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$total_protein_mean)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$total_protein_mean)

phot_absChange_totalprot <- (x - y)/x

# calvin cycle

p <- lm(calv_mean ~ log10(prec), data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(log10(data_means$prec))
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(log10(data_means$prec))

calv_absChange_prec <- (x - y)/x

p <- lm(calv_mean ~ leafrad_mean, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$leafrad_mean)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$leafrad_mean)

calv_absChange_leafrad <- (x - y)/x

p <- lm(calv_mean ~ gap_mean, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$gap_mean)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$gap_mean)

calv_absChange_gap <- (x - y)/x

p <- lm(calv_mean ~ tavg, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$tavg)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$tavg)

calv_absChange_tavg <- (x - y)/x

p <- lm(calv_mean ~ LMA_mean, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$LMA_mean)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$LMA_mean)

calv_absChange_LMA <- (x - y)/x

p <- lm(calv_mean ~ total_protein_mean, data_means)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data_means$total_protein_mean)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data_means$total_protein_mean)

calv_absChange_totalprot <- (x - y)/x





source('scripts/prep_data_mg_per_mm2.R')

datax <- select(data, total_protein, leaf_rad, Photosystems, calvin_cycle) %>%
  gather(key = 'funccat', value = 'protein_amount', -total_protein, -leaf_rad)

j <- ggplot(datax, aes(y = protein_amount, x = total_protein, shape = funccat, colour = leaf_rad)) + geom_point() +
   scale_colour_gradientn('Mean annual irradiance @ leaf (MJ/m2/yr)', colours=rainbow(2))
j

cor(data$calvin_cycle, data$total_protein)
cor(data$Photosystems, data$total_protein)











## Photosystems & calvin_cycle vs total protein, colour scaled

source('scripts/transformations.R')
source('scripts/prep_data_mg_per_mm2.R')

tiff('C:/Users/James/Desktop/stuff/PEPMOB/D14/output/figures/20170501/tiff/photosys_vs_totalprotein.tiff', height = 6, width = 6.5, units = 'in', res =300)
x <- ggplot(data, aes(x = total_protein, y = Photosystems)) + geom_point(aes(colour = leaf_rad)) + scale_colour_gradientn('Mean annual irradiance @ leaf (MJ/m2/yr)', colours=rainbow(2))
x <- x + ylab('Photosystem protein abundance (mg/m2)') + xlab('Total protein abundance (mg/m2)') #+ xlim(0,150000)
x <- x + theme_classic() + theme(text = element_text(size = 17),
                                 legend.title=element_text(size=14),
                                 legend.position="top"
)
x
dev.off()

cor(data$Photosystems, data$total_protein)

tiff('C:/Users/James/Desktop/stuff/PEPMOB/D14/output/figures/20170501/tiff/calvin_vs_totalprotein.tiff', height = 6, width = 6.5, units = 'in', res =300)
x <- ggplot(data, aes(x = total_protein, y = calvin_cycle)) + geom_point(aes(colour = leaf_rad)) + scale_colour_gradientn('Mean annual irradiance @ leaf (MJ/m2/yr)', colours=rainbow(2))
x <- x + ylab('Calvin cycle protein abundance (mg/m2)') + xlab('Total protein abundance (mg/m2)') #+ xlim(0,150000)
x <- x + theme_classic() + theme(text = element_text(size = 17),
                                 legend.title=element_text(size=15),
                                 legend.position="top"
)
x
dev.off()

cor(data$calvin_cycle, data$total_protein)