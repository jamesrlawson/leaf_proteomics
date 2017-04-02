source('scripts/transformations.R')

agg_plot_save(depvar = 'LMA_g_per_m2', indepvar = 'leafrad_mean', logx = FALSE, proportion = FALSE, indepvarType = 'leafrad', 
              labs =  c('LMA (g per m2)', 'leaf level irradiance'), outDir = 'output/figures/20170329/tiff', goldenRatio = FALSE)

agg_plot_save(depvar = 'LMA_g_per_m2', indepvar = 'gap_mean', logx = FALSE, proportion = FALSE, indepvarType = 'gap', 
              labs =  c('LMA (g per m2)', 'canopy gap fraction (%)'), outDir = 'output/figures/20170329/tiff', goldenRatio = FALSE)

agg_plot_save(depvar = 'LMA_g_per_m2', indepvar = 'N_per_area', logx = FALSE, proportion = FALSE, indepvarType = 'LMA', 
              labs =  c('LMA (g per m2)', 'N per area'), outDir = 'output/figures/20170329/tiff', goldenRatio = FALSE)

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


