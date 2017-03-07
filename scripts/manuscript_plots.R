# plots for total protein

source('scripts/transformations.R')

# total protein

agg_plot_save('total_protein', 'pwmt', logx=TRUE, proportion = FALSE, indepvarType = 'standard', labs = c('total protein abundance', 'Mean precipitation in the wettest month (mm)'))
agg_plot_save('total_protein', 'tavg', logx=FALSE, proportion = FALSE, indepvarType = 'standard', labs = c('total protein abundance', expression(paste("Mean annual temperature ",degree,"C", sep = ""))))
agg_plot_save('total_protein', 'leafrad_mean', logx=FALSE, proportion = FALSE, indepvarType = 'leafrad', labs =  c('total protein abundance', expression(paste('Mean annual irradiance'~(MJ~m^{-2}~yr^{-1})), sep="")))
agg_plot_save('total_protein', 'gap_mean', logx=FALSE, proportion = FALSE, indepvarType = 'gap', labs =  c('total protein abundance', 'Canopy openness (%)'))

# photosystems

agg_plot_save('Photosystems', 'prec', logx=TRUE, proportion = TRUE, indepvarType = 'standard', labs = c('Calvin cycle / photosystems abundance', 'Mean annual precipitation (mm)'))
agg_plot_save('Photosystems', 'tavg', logx=FALSE, proportion = TRUE, indepvarType = 'standard', labs = c('photosystems protein abundance', expression(paste("Mean annual temperature ",degree,"C", sep = ""))))
agg_plot_save('Photosystems', 'leafrad_mean', logx=FALSE, proportion = TRUE, indepvarType = 'leafrad', labs =  c('photosystems protein abundance', expression(paste('Mean annual irradiance'~(MJ~m^{-2}~yr^{-1})), sep="")))
agg_plot_save('Photosystems', 'gap_mean', logx=FALSE, proportion = TRUE, indepvarType = 'gap', labs =  c('photosystems protein abundance', 'Canopy openness (%)'))


# calvin cycle

agg_plot_save('Calvin_cycle', 'pdmt', logx=FALSE, proportion = TRUE, indepvarType = 'standard', labs = c('Calvin cycle protein abundance', 'Mean precipitation in the driest month (mm)'))
agg_plot_save('Calvin_cycle', 'tavg', logx=FALSE, proportion = TRUE, indepvarType = 'standard', labs = c('Calvin cycle protein abundance', expression(paste("Mean annual temperature ",degree,"C", sep = ""))))
agg_plot_save('Calvin_cycle', 'leafrad_mean', logx=FALSE, proportion = TRUE, indepvarType = 'leafrad', labs =  c('Calvin cycle protein abundance', expression(paste('Mean annual irradiance'~(MJ~m^{-2}~yr^{-1})), sep="")))
agg_plot_save('Calvin_cycle', 'gap_mean', logx=FALSE, proportion = TRUE, indepvarType = 'gap', labs =  c('Calvin cycle protein abundance', 'Canopy openness (%)'))

# photosystems per calvin cycle

agg_plot_save('calv_per_photo', 'pdmt', logx=FALSE, proportion = TRUE, indepvarType = 'standard', labs = c('Calvin cycle / photosystems abundance', 'Mean precipitation in the driest month (mm)'))
agg_plot_save('calv_per_photo', 'prec', logx=TRUE, proportion = TRUE, indepvarType = 'standard', labs = c('Calvin cycle / photosystems abundance', 'Mean annual precipitation (mm)'))
agg_plot_save('calv_per_photo', 'tavg', logx=FALSE, proportion = TRUE, indepvarType = 'standard', labs = c('Calvin cycle / photosystems abundance', expression(paste("Mean annual temperature ",degree,"C", sep = ""))))
agg_plot_save('calv_per_photo', 'leafrad_mean', logx=FALSE, proportion = TRUE, indepvarType = 'leafrad', labs =  c('Calvin cycle / photosystems abundance', expression(paste('Mean annual irradiance'~(MJ~m^{-2}~yr^{-1})), sep="")))
agg_plot_save('calv_per_photo', 'gap_mean', logx=FALSE, proportion = TRUE, indepvarType = 'gap', labs =  c('Calvin cycle / photosystems abundance', 'Canopy openness (%)'))

