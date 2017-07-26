


aggplot_save_combined_4 <- function(proportion, indepvar, logx = FALSE, indepvarType, labs, outDir, fileType = 'tiff', goldenRatio = FALSE, yaxis = FALSE) {
  
  # data - the output of scripts/prep_data.R or sources/prep_data_mg_per_mm2.R (contains all climate/env variables and protein amounts) 
  # depvar - the dependent variable (e.g. 'Photosystems')  
  # indepvar - the independent variable (e.g. 'prec')  # logx - should x be logged? (x axis to be printed as log scaled)
  # proportion - does the y axis label describe proportions or mg/mm2 ?
  # indepvarType - standard, gap, leafrad or LMA (to use appropriate horizontal error bars - standard uses none)  
  # labs - should contain:   
  #   a pretty form of the dependent variable (labs[1])  
  #   the x axis label (labs[2]) 
  # outdir - directory to output to 
  # golden ratio - plot dimensions golden ratio? 
  
  require(lazyeval) 
  if(proportion) { # source data conditional on proportion 
    source('scripts/prep_data.R') 
  } else { 
    source('scripts/prep_data_mg_per_mm2.R') 
  } 
  
  if(proportion) { # set labels conditional on proportion
    type = '_proportion_' 
  } else { 
    type = '_mg-per-m2_' 
  } 
  
  
  dep_means <- group_by(data, ID) %>% 
    dplyr::summarise(#calv = mean(calvin_cycle, na.rm=TRUE),     
                     phot = mean(Photosystems, na.rm=TRUE),
                     etrans = mean(cytochrome_b6f, na.rm=TRUE)) %>%   
    gather(key = 'funccat', value = 'protein_mean', -ID)
  
  dep_means <- group_by(data, ID) %>% 
    dplyr::summarise(#calv = SE(calvin_cycle),     
                     phot = SE(Photosystems),
                     etrans = SE(cytochrome_b6f)) %>%
    gather(key = 'funccat', value = 'protein_SE', -ID) %>%  
    full_join(dep_means, by = c('ID', 'funccat')) %>%  
    full_join(data, by = 'ID') %>% 
    distinct(ID, funccat, .keep_all=TRUE) 
  
  if(logx) { # get significance TRUE/FALSE for model, depending on if logx or not
  #  sig_calv <- summary(lm(dep_means[dep_means$funccat %in% 'calv',]$protein_mean ~ log10(dep_means[dep_means$funccat %in% 'calv',][[indepvar]])))$coefficients[,4][2] < 0.05  
    sig_phot <- summary(lm(dep_means[dep_means$funccat %in% 'phot',]$protein_mean ~ log10(dep_means[dep_means$funccat %in% 'phot',][[indepvar]])))$coefficients[,4][2] < 0.05 
    sig_etrans <- summary(lm(dep_means[dep_means$funccat %in% 'etrans',]$protein_mean ~ log10(dep_means[dep_means$funccat %in% 'etrans',][[indepvar]])))$coefficients[,4][2] < 0.05 
    
  } else {
  #  sig_calv <- summary(lm(dep_means[dep_means$funccat %in% 'calv',]$protein_mean ~ dep_means[dep_means$funccat %in% 'calv',][[indepvar]]))$coefficients[,4][2] < 0.05 
    sig_phot <- summary(lm(dep_means[dep_means$funccat %in% 'phot',]$protein_mean ~ dep_means[dep_means$funccat %in% 'phot',][[indepvar]]))$coefficients[,4][2] < 0.05
    sig_etrans <- summary(lm(dep_means[dep_means$funccat %in% 'etrans',]$protein_mean ~ dep_means[dep_means$funccat %in% 'etrans',][[indepvar]]))$coefficients[,4][2] < 0.05
    
  }
  
  dep_means$sig <- NA 
 # dep_means[dep_means$funccat == 'calv',]$sig <- sig_calv                                
  dep_means[dep_means$funccat == 'phot',]$sig <- sig_phot 
  dep_means[dep_means$funccat == 'etrans',]$sig <- sig_etrans 
  
  
  # directory to save to
  mainDir <- getwd() 
  
  if (!file.exists(outDir)){ 
    dir.create(file.path(mainDir, outDir), showWarnings = FALSE, recursive = TRUE)
  }
  
  if(goldenRatio) { 
    Width = 6*1.618 
  } else {
    Width = 7.5 
  } 
  
  
  
  if(fileType == 'svg'){
    svg(paste(outDir, '/',indepvar, type, '.svg', sep = ""), height = 6, width = Width)
  } else {
    if(fileType == 'tiff') {
      tiff(paste(outDir, '/', indepvar, type, '.tiff', sep = ""), height = 6, width = Width, units = 'in', res = 300)
    }
  }
  
  # plot object
  
  p <- ggplot(data = dep_means, aes(y = protein_mean, x = dep_means[[indepvar]])) + geom_point(aes(colour = funccat, shape = funccat),size = 1.5)
  
  p <- p + geom_smooth(data = dep_means[dep_means$sig == TRUE,], aes(y = protein_mean, x = get(indepvar), colour = funccat), method = 'lm', se = F, size = 0.5)
  
  p <- p + scale_colour_manual(values = c('blue','red', 'black'))
  
  if(proportion) {
    p <- p + xlab(labs[2]) +  ylab(paste(labs[1], ' (proportion)', sep = ""))
  } else {
    p <- p + xlab(labs[2]) +  ylab(bquote(.(labs[1]) ~ '(mg / ' ~ m^{2} ~ ')'))
    
  }
  
  p <- p + expand_limits(y=0)
  
  p <- p + theme_classic()
  p <- p + theme(legend.title=element_blank(),
                 legend.position="none",
                 axis.line = element_line(colour = "black"),
                 text = element_text(size = 18))
  
  if(!yaxis) {
    p <- p + theme(axis.title.y=element_blank())
  }
  
  number_ticks <- function(n) {function(limits) pretty(limits, n)}
  
  if(logx) {
    
    errorbar_width <- (log10(max(data[[indepvar]], na.rm=TRUE)) - log10(min(data[[indepvar]], na.rm=TRUE))) / 50
    p <- p + geom_errorbar(aes(colour = funccat, ymin = protein_mean - protein_SE, ymax = protein_mean + protein_SE), width = errorbar_width, alpha = 0.3, size = 0.5)
    
    if(indepvarType == 'gap') {
      
      p <- p + geom_errorbarh(aes(xmin = gap_mean - gap_SE, xmax = gap_mean + gap_SE, colour = funccat), alpha = 0.3, size = 0.5)
      
    } else { if(indepvarType == 'leafrad') {
      
      p <- p + geom_errorbarh(aes(xmin = leafrad_mean - leafrad_SE, xmax = leafrad_mean + leafrad_SE, colour = funccat), alpha = 0.3, size = 0.5)
      
    }
    }
    
    p <- p + scale_x_continuous(breaks=scales::pretty_breaks(n=10),trans='log10')
    
  } else {
    
    errorbar_width <- (max(data[[indepvar]], na.rm=TRUE) - min(data[[indepvar]], na.rm=TRUE)) / 50
    p <- p + geom_errorbar(aes(colour = funccat, ymin = protein_mean - protein_SE, ymax = protein_mean + protein_SE), width = errorbar_width, alpha = 0.3, size = 0.5)
    
    if(indepvarType == 'gap') {
      
      p <- p + geom_errorbarh(aes(xmin = gap_mean - gap_SE, xmax = gap_mean + gap_SE, colour = funccat), alpha = 0.3, size = 0.5)
      
    } else { if(indepvarType == 'leafrad') {
      
      p <- p + geom_errorbarh(aes(xmin = leafrad_mean - leafrad_SE, xmax = leafrad_mean + leafrad_SE, colour = funccat), alpha = 0.3, size = 0.5)
      
      
    } else { if(indepvarType == 'LMA') {
      
      p <- p + geom_errorbarh(aes(xmin = LMA_mean - LMA_SE, xmax = LMA_mean + LMA_SE, colour = funccat), alpha = 0.3, size = 0.5)
      
    } else { if(indepvarType == 'Narea') {
      
      p <- p + geom_errorbarh(aes(xmin = Narea_mean - Narea_SE, xmax = Narea_mean + Narea_SE, colour = funccat), alpha = 0.3, size = 0.5)
      
    } else { if(indepvarType == 'total_protein') {
      
      p <- p + geom_errorbarh(aes(xmin = total_protein_mean - total_protein_SE, xmax = total_protein_mean + total_protein_SE, colour = funccat), alpha = 0.3, size = 0.5)
      
    }
      
    }
      
    }
      
    }
      
    }
    
    p <- p + scale_x_continuous(breaks=scales::pretty_breaks(n=10))
    
  }
  
  print(p)
  
  dev.off()
  
}



aggplot_save_combined_4(indepvar = 'prec', logx = TRUE, proportion = FALSE, indepvarType = 'standard', 
                        labs =  c('Protein amount', 'Mean annual precip. (m/yr)'), outDir = 'output/figures/20170724/tiff', goldenRatio = FALSE, fileType = 'tiff')

aggplot_save_combined_4(indepvar = 'prec', logx = TRUE, proportion = TRUE, indepvarType = 'standard', 
                        labs =  c('Protein amount', 'Mean annual precip. (m/yr)'), outDir = 'output/figures/20170724/tiff', goldenRatio = FALSE, fileType = 'tiff')


aggplot_save_combined_4(indepvar = 'gap_mean', logx = FALSE, proportion = FALSE, indepvarType = 'gap', 
                        labs =  c('Protein amount', 'Canopy openness (%)'), outDir = 'output/figures/20170724/tiff', goldenRatio = FALSE, fileType = 'tiff')

aggplot_save_combined_4(indepvar = 'gap_mean', logx = FALSE, proportion = TRUE, indepvarType = 'gap', 
                        labs =  c('Protein amount', 'Canopy openness (%)'), outDir = 'output/figures/20170724/tiff', goldenRatio = FALSE, fileType = 'tiff')


aggplot_save_combined_4(indepvar = 'leafrad_mean', logx = FALSE, proportion = FALSE, indepvarType = 'leafrad', 
                        labs =  c('Protein amount', 'Irradiance (MJ/m2/year)'), outDir = 'output/figures/20170724/tiff', goldenRatio = FALSE, fileType = 'tiff')

aggplot_save_combined_4(indepvar = 'leafrad_mean', logx = FALSE, proportion = TRUE, indepvarType = 'leafrad', 
                        labs =  c('Protein amount', 'Irradiance (MJ/m2/year)'), outDir = 'output/figures/20170724/tiff', goldenRatio = FALSE, fileType = 'tiff')

aggplot_save_combined_4(indepvar = 'tavg', logx = FALSE, proportion = FALSE, indepvarType = 'standard', 
                        labs =  c('Protein amount', 'Mean annual temperature (oC)'), outDir = 'output/figures/20170724/tiff', goldenRatio = FALSE, fileType = 'tiff')

aggplot_save_combined_4(indepvar = 'tavg', logx = FALSE, proportion = TRUE, indepvarType = 'standard', 
                        labs =  c('Protein amount', 'Mean annual temperature (oC)'), outDir = 'output/figures/20170724/tiff', goldenRatio = FALSE, fileType = 'tiff')

aggplot_save_combined_4(indepvar = 'LMA_mean', logx = FALSE, proportion = FALSE, indepvarType = 'LMA', 
                        labs =  c('Protein amount', 'LMA (g/m2)'), outDir = 'output/figures/20170724/tiff', goldenRatio = FALSE, fileType = 'tiff')

aggplot_save_combined_4(indepvar = 'LMA_mean', logx = FALSE, proportion = TRUE, indepvarType = 'LMA', 
                        labs =  c('Protein amount', 'LMA (g/m2)'), outDir = 'output/figures/20170724/tiff', goldenRatio = FALSE, fileType = 'tiff')


aggplot_save_combined_4(indepvar = 'total_protein_mean', logx = FALSE, proportion = FALSE, indepvarType = 'total_protein', 
                        labs =  c('Protein amount', 'total protein (mg/m2)'), outDir = 'output/figures/20170724/tiff', goldenRatio = FALSE, fileType = 'tiff')

aggplot_save_combined_4(indepvar = 'total_protein_mean', logx = FALSE, proportion = TRUE, indepvarType = 'total_protein', 
                        labs =  c('Protein amount', 'total protein (mg/m2)'), outDir = 'output/figures/20170724/tiff', goldenRatio = FALSE, fileType = 'tiff')


include_leaf_N=TRUE
source('scripts/transformations.R')

aggplot_save_combined_4(indepvar = 'Narea_mean', logx = FALSE, proportion = FALSE, indepvarType = 'Narea', 
                        labs =  c('Protein amount', 'Leaf N (mg/m2)'), outDir = 'output/figures/20170706/tiff', goldenRatio = FALSE, fileType = 'tiff')

aggplot_save_combined_4(indepvar = 'Narea_mean', logx = FALSE, proportion = TRUE, indepvarType = 'Narea', 
                        labs =  c('Protein amount', 'Leaf N (mg/m2)'), outDir = 'output/figures/20170706/tiff', goldenRatio = FALSE, fileType = 'tiff')
