source('scripts/transformations.R')

agg_plot_save(depvar = 'LHC1', indepvar = 'gap_mean', logx = FALSE, proportion = TRUE, indepvarType = 'gap', 
              labs =  c('LHC1', 'Canopy openness (%)'), outDir = 'output/figures/20170410/tiff', goldenRatio = TRUE)

agg_plot_save(depvar = 'LHC1', indepvar = 'gap_mean', logx = FALSE, proportion = FALSE, indepvarType = 'gap', 
              labs =  c('LHC1', 'Canopy openness (%)'), outDir = 'output/figures/20170410/tiff', goldenRatio = TRUE)


agg_plot_save(depvar = 'LHC2', indepvar = 'gap_mean', logx = FALSE, proportion = TRUE, indepvarType = 'gap', 
              labs =  c('LHC2', 'Canopy openness (%)'), outDir = 'output/figures/20170410/tiff', goldenRatio = TRUE)

agg_plot_save(depvar = 'LHC2', indepvar = 'gap_mean', logx = FALSE, proportion = FALSE, indepvarType = 'gap', 
              labs =  c('LHC2', 'Canopy openness (%)'), outDir = 'output/figures/20170410/tiff', goldenRatio = TRUE)


agg_plot_save(depvar = 'PSI_min_LHC1', indepvar = 'gap_mean', logx = FALSE, proportion = TRUE, indepvarType = 'gap', 
              labs =  c('PSI_min_LHC1', 'Canopy openness (%)'), outDir = 'output/figures/20170410/tiff', goldenRatio = TRUE)

agg_plot_save(depvar = 'PSI_min_LHC1', indepvar = 'gap_mean', logx = FALSE, proportion = FALSE, indepvarType = 'gap', 
              labs =  c('PSI_min_LHC1', 'Canopy openness (%)'), outDir = 'output/figures/20170410/tiff', goldenRatio = TRUE)


agg_plot_save(depvar = 'PSII_min_LHC2', indepvar = 'gap_mean', logx = FALSE, proportion = TRUE, indepvarType = 'gap', 
              labs =  c('PSII_min_LHC2', 'Canopy openness (%)'), outDir = 'output/figures/20170410/tiff', goldenRatio = TRUE)

agg_plot_save(depvar = 'PSII_min_LHC2', indepvar = 'gap_mean', logx = FALSE, proportion = FALSE, indepvarType = 'gap', 
              labs =  c('PSII_min_LHC2', 'Canopy openness (%)'), outDir = 'output/figures/20170410/tiff', goldenRatio = TRUE)


agg_plot_save(depvar = 'LHC1_per_PSI', indepvar = 'gap_mean', logx = FALSE, proportion = TRUE, indepvarType = 'gap', 
              labs =  c('LHC1_per_PSI', 'Canopy openness (%)'), outDir = 'output/figures/20170410/tiff', goldenRatio = TRUE)

agg_plot_save(depvar = 'LHC2_per_PSII', indepvar = 'gap_mean', logx = FALSE, proportion = TRUE, indepvarType = 'gap', 
              labs =  c('LHC2_per_PSII', 'Canopy openness (%)'), outDir = 'output/figures/20170410/tiff', goldenRatio = TRUE)



agg_plot_save(depvar = 'LHC1', indepvar = 'PSI_min_LHC1', logx = FALSE, proportion = FALSE, indepvarType = 'standard', 
              labs =  c('LHC1', 'PSI'), outDir = 'output/figures/20170410/tiff', goldenRatio = TRUE)

agg_plot_save(depvar = 'LHC2', indepvar = 'PSII_min_LHC2', logx = FALSE, proportion = FALSE, indepvarType = 'standard', 
              labs =  c('LHC2', 'PSII'), outDir = 'output/figures/20170410/tiff', goldenRatio = TRUE)



