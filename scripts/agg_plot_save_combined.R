
source('scripts/functions.R')


agg_plot_save_combined(indepvar = 'prec', logx = TRUE, proportion = FALSE, indepvarType = 'standard', 
              labs =  c('Protein amount', 'Mean annual precip. (mm)'), outDir = 'output/figures/20170501/tiff', goldenRatio = FALSE)

agg_plot_save_combined(indepvar = 'prec', logx = TRUE, proportion = TRUE, indepvarType = 'standard', 
              labs =  c('Protein amount', 'Mean annual precip. (mm)'), outDir = 'output/figures/20170501/tiff', goldenRatio = FALSE)


agg_plot_save_combined(indepvar = 'gap_mean', logx = TRUE, proportion = FALSE, indepvarType = 'gap', 
              labs =  c('Protein amount', 'Canopy openness (%)'), outDir = 'output/figures/20170501/tiff', goldenRatio = FALSE)

agg_plot_save_combined(indepvar = 'gap_mean', logx = TRUE, proportion = TRUE, indepvarType = 'gap', 
              labs =  c('Protein amount', 'Canopy openness (%)'), outDir = 'output/figures/20170501/tiff', goldenRatio = FALSE)


agg_plot_save_combined(indepvar = 'leafrad_mean', logx = FALSE, proportion = FALSE, indepvarType = 'leafrad', 
                       labs =  c('Protein amount', 'Irradiance (MJ/m2/year)'), outDir = 'output/figures/20170501/tiff', goldenRatio = FALSE)

agg_plot_save_combined(indepvar = 'leafrad_mean', logx = FALSE, proportion = TRUE, indepvarType = 'leafrad', 
                       labs =  c('Protein amount', 'Irradiance (MJ/m2/year)'), outDir = 'output/figures/20170501/tiff', goldenRatio = FALSE)

agg_plot_save_combined(indepvar = 'tavg', logx = FALSE, proportion = FALSE, indepvarType = 'standard', 
                       labs =  c('Protein amount', 'Mean annual temperature (oC)'), outDir = 'output/figures/20170501/tiff', goldenRatio = FALSE)

agg_plot_save_combined(indepvar = 'tavg', logx = FALSE, proportion = TRUE, indepvarType = 'standard', 
                       labs =  c('Protein amount', 'Mean annual temperature (oC)'), outDir = 'output/figures/20170501/tiff', goldenRatio = FALSE)



agg_plot_save_combined(indepvar = 'LMA_mean', logx = FALSE, proportion = FALSE, indepvarType = 'LMA', 
                       labs =  c('Protein amount', 'LMA (g/m2)'), outDir = 'output/figures/20170501/tiff', goldenRatio = FALSE)

agg_plot_save_combined(indepvar = 'LMA_mean', logx = FALSE, proportion = TRUE, indepvarType = 'LMA', 
                       labs =  c('Protein amount', 'LMA (g/m2)'), outDir = 'output/figures/20170501/tiff', goldenRatio = FALSE)




