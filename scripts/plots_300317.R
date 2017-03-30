source('scripts/transformations.R')

agg_plot_save(depvar = 'LMA_g_per_m2', indepvar = 'leafrad_mean', logx = FALSE, proportion = FALSE, indepvarType = 'leafrad', 
              labs =  c('LMA (g per m2)', 'leaf level irradiance'), outDir = 'output/figures/20170329/tiff', goldenRatio = FALSE)

agg_plot_save(depvar = 'LMA_g_per_m2', indepvar = 'gap_mean', logx = FALSE, proportion = FALSE, indepvarType = 'gap', 
              labs =  c('LMA (g per m2)', 'canopy gap fraction (%)'), outDir = 'output/figures/20170329/tiff', goldenRatio = FALSE)


agg_plot_save(depvar = 'LMA_g_per_m2', indepvar = 'N_per_area', logx = FALSE, proportion = FALSE, indepvarType = 'LMA', 
              labs =  c('LMA (g per m2)', 'N per area'), outDir = 'output/figures/20170329/tiff', goldenRatio = FALSE)

agg_plot_save(depvar = 'total_protein', indepvar = 'LMA_mean', logx = FALSE, proportion = FALSE, indepvarType = 'LMA', 
              labs =  c('total protein', 'Leaf mass per area (g / m2)'), outDir = 'output/figures/20170329/tiff', goldenRatio = FALSE)

# relative and proportional funccats vs LMA

agg_plot_save(depvar = 'Photosystems', indepvar = 'LMA_mean', logx = FALSE, proportion = FALSE, indepvarType = 'LMA', 
              labs =  c('Photosystems', 'Leaf mass per area (g / m2)'), outDir = 'output/figures/20170329/tiff', goldenRatio = FALSE)

agg_plot_save(depvar = 'Photosystems', indepvar = 'LMA_mean', logx = FALSE, proportion = TRUE, indepvarType = 'LMA', 
              labs =  c('Photosystems', 'Leaf mass per area (g / m2)'), outDir = 'output/figures/20170329/tiff', goldenRatio = FALSE)

agg_plot_save(depvar = 'Calvin_cycle', indepvar = 'LMA_mean', logx = FALSE, proportion = FALSE, indepvarType = 'LMA', 
              labs =  c('Calvin cycle', 'Leaf mass per area (g / m2)'), outDir = 'output/figures/20170329/tiff', goldenRatio = FALSE)

agg_plot_save(depvar = 'Calvin_cycle', indepvar = 'LMA_mean', logx = FALSE, proportion = TRUE, indepvarType = 'LMA', 
              labs =  c('Calvin cycle', 'Leaf mass per area (g / m2)'), outDir = 'output/figures/20170329/tiff', goldenRatio = FALSE)

agg_plot_save(depvar = 'Photorespiration', indepvar = 'LMA_mean', logx = FALSE, proportion = FALSE, indepvarType = 'LMA', 
              labs =  c('Photorespiration', 'Leaf mass per area (g / m2)'), outDir = 'output/figures/20170329/tiff', goldenRatio = FALSE)

agg_plot_save(depvar = 'Photorespiration', indepvar = 'LMA_mean', logx = FALSE, proportion = TRUE, indepvarType = 'LMA', 
              labs =  c('Photorespiration', 'Leaf mass per area (g / m2)'), outDir = 'output/figures/20170329/tiff', goldenRatio = FALSE)

agg_plot_save(depvar = 'protein', indepvar = 'LMA_mean', logx = FALSE, proportion = FALSE, indepvarType = 'LMA', 
              labs =  c('Protein synth', 'Leaf mass per area (g / m2)'), outDir = 'output/figures/20170329/tiff', goldenRatio = FALSE)

agg_plot_save(depvar = 'protein', indepvar = 'LMA_mean', logx = FALSE, proportion = TRUE, indepvarType = 'LMA', 
              labs =  c('Protein synth', 'Leaf mass per area (g / m2)'), outDir = 'output/figures/20170329/tiff', goldenRatio = FALSE)

agg_plot_save(depvar = 'ATP_synthase_chloroplastic', indepvar = 'LMA_mean', logx = FALSE, proportion = FALSE, indepvarType = 'LMA', 
              labs =  c('ATP_synthase_chloroplastic (mg per m2)', 'Leaf mass per area (g / m2)'), outDir = 'output/figures/20170329/tiff', goldenRatio = FALSE)

agg_plot_save(depvar = 'ATP_synthase_chloroplastic', indepvar = 'LMA_mean', logx = FALSE, proportion = TRUE, indepvarType = 'LMA', 
              labs =  c('ATP_synthase_chloroplastic (mg per m2)', 'Leaf mass per area (g / m2)'), outDir = 'output/figures/20170329/tiff', goldenRatio = FALSE)


# total protein trends (LMA, Narea, tavg, prec, gap,leafrad)

agg_plot_save(depvar = 'total_protein', indepvar = 'LMA_mean', logx = FALSE, proportion = FALSE, indepvarType = 'LMA', 
              labs =  c('Protein / leaf area', 'Leaf mass per area (g per m2)'), outDir = 'output/figures/20170329/tiff', goldenRatio = FALSE)

agg_plot_save(depvar = 'total_protein', indepvar = 'Narea_mean', logx = FALSE, proportion = FALSE, indepvarType = 'Narea', 
              labs =  c('Protein / leaf area', 'N per area (mg per m2)'), outDir = 'output/figures/20170329/tiff', goldenRatio = FALSE)

agg_plot_save(depvar = 'total_protein', indepvar = 'tavg', logx = FALSE, proportion = FALSE, indepvarType = 'standard', 
              labs =  c('Protein / leaf area', 'Mean annual temperature (oC)'), outDir = 'output/figures/20170329/tiff', goldenRatio = TRUE)

agg_plot_save(depvar = 'total_protein', indepvar = 'prec', logx = TRUE, proportion = FALSE, indepvarType = 'standard', 
              labs =  c('Protein / leaf area', 'Mean annual precipitation'), outDir = 'output/figures/20170329/tiff', goldenRatio = TRUE)

agg_plot_save(depvar = 'total_protein', indepvar = 'gap_mean', logx = FALSE, proportion = FALSE, indepvarType = 'gap', 
              labs =  c('Protein / leaf area', 'Canopy openness (%)'), outDir = 'output/figures/20170329/tiff', goldenRatio = TRUE)

agg_plot_save(depvar = 'total_protein', indepvar = 'leafrad_mean', logx = FALSE, proportion = FALSE, indepvarType = 'leafraf', 
              labs =  c('Protein / leaf area', 'Mean annual irradiance @ leaf (MJ/m2/yr)'), outDir = 'output/figures/20170329/tiff', goldenRatio = TRUE)


