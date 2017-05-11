
source('scripts/functions.R')


agg_plot_save_combined(indepvar = 'prec', logx = TRUE, proportion = FALSE, indepvarType = 'standard', 
              labs =  c('Protein amount', 'Mean annual precip. (m/yr)'), outDir = 'output/figures/20170501/tiff', goldenRatio = FALSE, fileType = 'tiff')

agg_plot_save_combined(indepvar = 'prec', logx = TRUE, proportion = TRUE, indepvarType = 'standard', 
              labs =  c('Protein amount', 'Mean annual precip. (m/yr)'), outDir = 'output/figures/20170501/tiff', goldenRatio = FALSE, fileType = 'tiff')


agg_plot_save_combined(indepvar = 'gap_mean', logx = FALSE, proportion = FALSE, indepvarType = 'gap', 
              labs =  c('Protein amount', 'Canopy openness (%)'), outDir = 'output/figures/20170501/tiff', goldenRatio = FALSE, fileType = 'tiff')

agg_plot_save_combined(indepvar = 'gap_mean', logx = FALSE, proportion = TRUE, indepvarType = 'gap', 
              labs =  c('Protein amount', 'Canopy openness (%)'), outDir = 'output/figures/20170501/tiff', goldenRatio = FALSE, fileType = 'tiff')


agg_plot_save_combined(indepvar = 'leafrad_mean', logx = FALSE, proportion = FALSE, indepvarType = 'leafrad', 
                       labs =  c('Protein amount', 'Irradiance (MJ/m2/year)'), outDir = 'output/figures/20170501/tiff', goldenRatio = FALSE, fileType = 'tiff')

agg_plot_save_combined(indepvar = 'leafrad_mean', logx = FALSE, proportion = TRUE, indepvarType = 'leafrad', 
                       labs =  c('Protein amount', 'Irradiance (MJ/m2/year)'), outDir = 'output/figures/20170501/tiff', goldenRatio = FALSE, fileType = 'tiff')

agg_plot_save_combined(indepvar = 'tavg', logx = FALSE, proportion = FALSE, indepvarType = 'standard', 
                       labs =  c('Protein amount', 'Mean annual temperature (oC)'), outDir = 'output/figures/20170501/tiff', goldenRatio = FALSE, fileType = 'tiff')

agg_plot_save_combined(indepvar = 'tavg', logx = FALSE, proportion = TRUE, indepvarType = 'standard', 
                       labs =  c('Protein amount', 'Mean annual temperature (oC)'), outDir = 'output/figures/20170501/tiff', goldenRatio = FALSE, fileType = 'tiff')

agg_plot_save_combined(indepvar = 'LMA_mean', logx = FALSE, proportion = FALSE, indepvarType = 'LMA', 
                       labs =  c('Protein amount', 'LMA (g/m2)'), outDir = 'output/figures/20170501/tiff', goldenRatio = FALSE, fileType = 'tiff')

agg_plot_save_combined(indepvar = 'LMA_mean', logx = FALSE, proportion = TRUE, indepvarType = 'LMA', 
                       labs =  c('Protein amount', 'LMA (g/m2)'), outDir = 'output/figures/20170501/tiff', goldenRatio = FALSE, fileType = 'tiff')


agg_plot_save_combined(indepvar = 'total_protein_mean', logx = FALSE, proportion = FALSE, indepvarType = 'total_protein', 
                       labs =  c('Protein amount', 'total protein (mg/m2)'), outDir = 'output/figures/20170501/tiff', goldenRatio = FALSE, fileType = 'tiff', total_prot = TRUE)

agg_plot_save_combined(indepvar = 'total_protein_mean', logx = FALSE, proportion = TRUE, indepvarType = 'total_protein', 
                       labs =  c('Protein amount', 'total protein (mg/m2)'), outDir = 'output/figures/20170501/tiff', goldenRatio = FALSE, fileType = 'tiff', total_prot = TRUE)


# modelled changes in protein amounts across env gradients

effect_size_all   <-  rbind(effect_size(proportion=TRUE, depvar = 'phot_mean'),
                                 effect_size(proportion=FALSE, depvar = 'phot_mean'),
                                 effect_size(proportion=TRUE, depvar = 'calv_mean'),
                                 effect_size(proportion=FALSE, depvar = 'calv_mean'),
                                 effect_size(proportion=TRUE, depvar = 'total_protein_mean'),
                                 effect_size(proportion=FALSE, depvar = 'total_protein_mean'))
write_csv(effect_size_all, 'output/effect_size_all.csv')




## Photosystems & calvin_cycle vs total protein, colour scaled

source('scripts/prep_data_mg_per_mm2.R')

bla <- gather(data, key = 'funccat', value = 'protein_amount', Photosystems, calvin_cycle, -total_protein)

ggplot(bla, aes(x = total_protein, y = protein_amount)) + geom_point(aes(shape = funccat, colour = leaf_rad)) + scale_colour_gradientn('Mean annual irradiance @ leaf (MJ/m2/yr)', colours=rainbow(2))



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