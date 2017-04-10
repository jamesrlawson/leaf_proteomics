source('scripts/transformations.R')

agg_plot_save(depvar = 'LMA_g_per_m2', indepvar = 'leafrad_mean', logx = FALSE, proportion = FALSE, indepvarType = 'leafrad', 
              labs =  c('LMA (g per m2)', 'Mean annual irradiance @ leaf (MJ/m2/yr)'), outDir = 'output/figures/20170329/tiff', goldenRatio = TRUE)

agg_plot_save(depvar = 'LMA_g_per_m2', indepvar = 'gap_mean', logx = FALSE, proportion = FALSE, indepvarType = 'gap', 
              labs =  c('LMA (g per m2)', 'Canopy openness (%)'), outDir = 'output/figures/20170329/tiff', goldenRatio = TRUE)

agg_plot_save(depvar = 'LMA_g_per_m2', indepvar = 'N_per_area', logx = FALSE, proportion = FALSE, indepvarType = 'Narea', 
              labs =  c('LMA (g per m2)', 'N per area'), outDir = 'output/figures/20170329/tiff', goldenRatio = TRUE)

agg_plot_save(depvar = 'LMA_g_per_m2', indepvar = 'tavg', logx = FALSE, proportion = FALSE, indepvarType = 'standard', 
              labs =  c('LMA (g per m2)', 'Mean annual temperature (oC)'), outDir = 'output/figures/20170329/tiff', goldenRatio = TRUE)

agg_plot_save(depvar = 'LMA_g_per_m2', indepvar = 'prec', logx = TRUE, proportion = FALSE, indepvarType = 'standard', 
              labs =  c('LMA (g per m2)', 'Mean annual precipitation (mm)'), outDir = 'output/figures/20170329/tiff', goldenRatio = TRUE)

# relative and proportional funccats vs LMA

agg_plot_save(depvar = 'Photosystems', indepvar = 'LMA_mean', logx = FALSE, proportion = FALSE, indepvarType = 'LMA', 
              labs =  c('Photosystems', 'Leaf mass per area (g / m2)'), outDir = 'output/figures/20170329/tiff', goldenRatio = FALSE)

agg_plot_save(depvar = 'Photosystems', indepvar = 'LMA_mean', logx = FALSE, proportion = TRUE, indepvarType = 'LMA', 
              labs =  c('Photosystems', 'Leaf mass per area (g / m2)'), outDir = 'output/figures/20170329/tiff', goldenRatio = FALSE)

agg_plot_save(depvar = 'Calvin_cycle', indepvar = 'LMA_mean', logx = FALSE, proportion = FALSE, indepvarType = 'LMA', 
              labs =  c('Calvin cycle', 'Leaf mass per area (g / m2)'), outDir = 'output/figures/20170329/tiff', goldenRatio = FALSE)

agg_plot_save(depvar = 'Calvin_cycle', indepvar = 'LMA_mean', logx = FALSE, proportion = TRUE, indepvarType = 'LMA', 
              labs =  c('Calvin cycle', 'Leaf mass per area (g / m2)'), outDir = 'output/figures/20170329/tiff', goldenRatio = FALSE)


# total protein trends (LMA, Narea, tavg, prec, gap,leafrad)

agg_plot_save(depvar = 'total_protein', indepvar = 'LMA_mean', logx = FALSE, proportion = FALSE, indepvarType = 'LMA', 
              labs =  c('Protein / leaf area', 'Leaf mass per area (g per m2)'), outDir = 'output/figures/20170329/tiff', goldenRatio = FALSE)

agg_plot_save(depvar = 'total_protein', indepvar = 'Narea_mean', logx = FALSE, proportion = FALSE, indepvarType = 'Narea', 
              labs =  c('Protein / leaf area', 'N per area (mg per m2)'), outDir = 'output/figures/20170329/tiff', goldenRatio = FALSE)

agg_plot_save(depvar = 'total_protein', indepvar = 'tavg', logx = FALSE, proportion = FALSE, indepvarType = 'standard', 
              labs =  c('Protein / leaf area', 'Mean annual temperature (oC)'), outDir = 'output/figures/20170329/tiff', goldenRatio = TRUE)

agg_plot_save(depvar = 'total_protein', indepvar = 'prec', logx = TRUE, proportion = FALSE, indepvarType = 'standard', 
              labs =  c('Protein / leaf area', 'Mean annual precipitation (mm)'), outDir = 'output/figures/20170329/tiff', goldenRatio = TRUE)

agg_plot_save(depvar = 'total_protein', indepvar = 'gap_mean', logx = FALSE, proportion = FALSE, indepvarType = 'gap', 
              labs =  c('Protein / leaf area', 'Canopy openness (%)'), outDir = 'output/figures/20170329/tiff', goldenRatio = TRUE)

agg_plot_save(depvar = 'total_protein', indepvar = 'leafrad_mean', logx = FALSE, proportion = FALSE, indepvarType = 'leafraf', 
              labs =  c('Protein / leaf area', 'Mean annual irradiance @ leaf (MJ/m2/yr)'), outDir = 'output/figures/20170329/tiff', goldenRatio = TRUE)


# photosystems, calv and photoresp vs env gradients plots

agg_plot_save(depvar = 'Photosystems', indepvar = 'tavg', logx = FALSE, proportion = TRUE, indepvarType = 'standard', 
              labs =  c('Photosystems protein', 'Mean annual temperature (oC)'), outDir = 'output/figures/20170329/tiff', goldenRatio = TRUE)

agg_plot_save(depvar = 'Photosystems', indepvar = 'gap_mean', logx = FALSE, proportion = TRUE, indepvarType = 'gap', 
              labs =  c('Photosystems protein', 'Canopy openness (%)'), outDir = 'output/figures/20170329/tiff', goldenRatio = TRUE)

agg_plot_save(depvar = 'Photosystems', indepvar = 'leafrad_mean', logx = FALSE, proportion = TRUE, indepvarType = 'leafrad', 
              labs =  c('Photosystems protein', 'Mean annual irradiance @ leaf (MJ/m2/yr)'), outDir = 'output/figures/20170329/tiff', goldenRatio = TRUE)

agg_plot_save(depvar = 'Photosystems', indepvar = 'prec', logx = TRUE, proportion = TRUE, indepvarType = 'standard', 
              labs =  c('Photosystems protein', 'Mean annual precipitation (mm)'), outDir = 'output/figures/20170329/tiff', goldenRatio = TRUE)



agg_plot_save(depvar = 'Calvin_cycle', indepvar = 'tavg', logx = FALSE, proportion = TRUE, indepvarType = 'standard', 
              labs =  c('Calvin cycle protein', 'Mean annual temperature (oC)'), outDir = 'output/figures/20170329/tiff', goldenRatio = TRUE)

agg_plot_save(depvar = 'Calvin_cycle', indepvar = 'gap_mean', logx = FALSE, proportion = TRUE, indepvarType = 'gap', 
              labs =  c('Calvin cycle protein', 'Canopy openness (%)'), outDir = 'output/figures/20170329/tiff', goldenRatio = TRUE)

agg_plot_save(depvar = 'Calvin_cycle', indepvar = 'leafrad_mean', logx = FALSE, proportion = TRUE, indepvarType = 'leafrad', 
              labs =  c('Calvin cycle protein', 'Mean annual irradiance @ leaf (MJ/m2/yr)'), outDir = 'output/figures/20170329/tiff', goldenRatio = TRUE)

agg_plot_save(depvar = 'Calvin_cycle', indepvar = 'prec', logx = TRUE, proportion = TRUE, indepvarType = 'standard', 
              labs =  c('Calvin cycle protein', 'Mean annual precipitation (mm)'), outDir = 'output/figures/20170329/tiff', goldenRatio = TRUE)


agg_plot_save(depvar = 'Photorespiration', indepvar = 'tavg', logx = FALSE, proportion = TRUE, indepvarType = 'standard', 
              labs =  c('Photorespiration protein', 'Mean annual temperature (oC)'), outDir = 'output/figures/20170329/tiff', goldenRatio = TRUE)

agg_plot_save(depvar = 'Photorespiration', indepvar = 'gap_mean', logx = FALSE, proportion = TRUE, indepvarType = 'gap', 
              labs =  c('Photorespiration protein', 'Canopy openness (%)'), outDir = 'output/figures/20170329/tiff', goldenRatio = TRUE)

agg_plot_save(depvar = 'Photorespiration', indepvar = 'leafrad_mean', logx = FALSE, proportion = TRUE, indepvarType = 'leafrad', 
              labs =  c('Photorespiration protein', 'Mean annual irradiance @ leaf (MJ/m2/yr)'), outDir = 'output/figures/20170329/tiff', goldenRatio = TRUE)

agg_plot_save(depvar = 'Photorespiration', indepvar = 'prec', logx = TRUE, proportion = TRUE, indepvarType = 'standard', 
              labs =  c('Photorespiration protein', 'Mean annual precipitation'), outDir = 'output/figures/20170329/tiff', goldenRatio = TRUE)



agg_plot_save(depvar = 'electron_transport_minATPsynth', indepvar = 'tavg', logx = FALSE, proportion = TRUE, indepvarType = 'standard', 
              labs =  c('Electron transport protein', 'Mean annual temperature (oC)'), outDir = 'output/figures/20170329/tiff', goldenRatio = TRUE)

agg_plot_save(depvar = 'electron_transport_minATPsynth', indepvar = 'tmax_recent_month', logx = FALSE, proportion = TRUE, indepvarType = 'standard', 
              labs =  c('Electron transport protein', 'Mean max daily temperature in last month (oC)'), outDir = 'output/figures/20170329/tiff', goldenRatio = TRUE)


agg_plot_save(depvar = 'Cytochrome_b6f', indepvar = 'tavg', logx = FALSE, proportion = TRUE, indepvarType = 'standard', 
              labs =  c('Cytochrome b6f / leaf area', 'Mean annual temperature (oC)'), outDir = 'output/figures/20170329/tiff', goldenRatio = TRUE)



# total protein CV vs Narea CV
source('scripts/prep_data.R')

TPCVvsNareaCV <- ggplot(data, aes(x = Narea_CV, y = total_protein_CV)) + geom_point(size = 2) + geom_smooth(method='lm', se=FALSE, colour = 'black', size=0.5) +
  theme_classic() +
  theme(text = element_text(size = 18)) +
  xlab('Variation in N per area (CV)') +
  ylab('Variation in total protein (CV)')
TPCVvsNareaCV


## Photosystems & Calvin_cycle vs total protein, colour scaled

source('scripts/transformations.R')
source('scripts/prep_data_mg_per_mm2.R')

tiff('C:/Users/James/Desktop/stuff/PEPMOB/D14/output/figures/20170329/tiff/photosys_vs_totalprotein.tiff', height = 6, width = 6.5, units = 'in', res =300)
x <- ggplot(data, aes(x = total_protein, y = Photosystems)) + geom_point(aes(colour = gap)) + scale_colour_gradientn('Canopy openness (%)', colours=rainbow(2))
x <- x + ylab('Photosystem protein abundance (mg/m2)') + xlab('Total protein abundance (mg/m2)') #+ xlim(0,150000)
x <- x + theme_classic() + theme(text = element_text(size = 17),
                                 legend.title=element_text(size=14),
                                 legend.position="top"
)
x
dev.off()

cor(data$Photosys, data$total_protein)

tiff('C:/Users/James/Desktop/stuff/PEPMOB/D14/output/figures/20170329/tiff/calvin_vs_totalprotein.tiff', height = 6, width = 6.5, units = 'in', res =300)
x <- ggplot(data, aes(x = total_protein, y = Calvin_cycle)) + geom_point(aes(colour = leaf_rad)) + scale_colour_gradientn('Mean annual irradiance @ leaf (MJ/m2/yr)', colours=rainbow(2))
x <- x + ylab('Calvin cycle protein abundance (mg/m2)') + xlab('Total protein abundance (mg/m2)') #+ xlim(0,150000)
x <- x + theme_classic() + theme(text = element_text(size = 17),
                                 legend.title=element_text(size=15),
                                 legend.position="top"
)
x
dev.off()

cor(data$Calvin_cycle, data$total_protein)

# means and ranges

source('scripts/prep_data_mg_per_mm2.R')

calvmean <- data %>% group_by(ID) %>% summarise(calvmean = mean(Calvin_cycle, na.rm=TRUE), 
                                                photmean = mean(Photosystems, na.rm=TRUE),
                                                photorespmean = mean(Photorespiration, na.rm=TRUE))

max(calvmean$calvmean) /min(calvmean$calvmean)
mean(calvmean$calvmean, na.rm=TRUE)
sd(calvmean$calvmean, na.rm=TRUE)


max(calvmean$photmean) / min(calvmean$photmean)
mean(calvmean$photmean, na.rm=TRUE)
sd(calvmean$photmean, na.rm=TRUE)

max(calvmean$photorespmean) / min(calvmean$photorespmean)
min(calvmean$photorespmean)
max(calvmean$photorespmean)
mean(calvmean$photorespmean, na.rm=TRUE)
sd(calvmean$photorespmean, na.rm=TRUE)

max(data$total_protein_mean) / min(data$total_protein_mean)
mean(data$total_protein_mean, na.rm=TRUE)
sd(data$total_protein_mean, na.rm=TRUE)



min((data$total_protein_mean / data$LMA_g_per_m2) / 1000)
max((data$total_protein_mean / data$LMA_g_per_m2) / 1000)
mean((data$total_protein_mean / data$LMA_g_per_m2) / 1000)
sd((data$total_protein_mean / data$LMA_g_per_m2) / 1000)






agg_plot_save(depvar = 'mitochondrial_electron_transport_ATP_synthesis', indepvar = 'tavg', logx = FALSE, proportion = TRUE, indepvarType = 'standard', 
              labs =  c('mitochondrial_electron_transport_ATP_synthesis', 'Mean annual temperature (oC)'), outDir = 'output/figures/20170329/tiff', goldenRatio = TRUE)

agg_plot_save(depvar = 'mitochondrial_electron_transport_ATP_synthesis', indepvar = 'tavg', logx = FALSE, proportion = FALSE, indepvarType = 'standard', 
              labs =  c('mitochondrial_electron_transport_ATP_synthesis', 'Mean annual temperature (oC)'), outDir = 'output/figures/20170329/tiff', goldenRatio = TRUE)



