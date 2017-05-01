
# this function outputs a pretty aggregated plot to your directory of choice, with stacked points for different funccats

agg_plot_save_combined <- function(proportion, indepvar, logx = FALSE, indepvarType, labs, outDir, fileType = 'tiff', goldenRatio = FALSE) {
  
  # data - the output of scripts/prep_data.R or sources/prep_data_mg_per_mm2.R (contains all climate/env variables and protein amounts)
  # depvar - the dependent variable (e.g. 'Photosystems')
  # indepvar - the independent variable (e.g. 'prec')
  # logx - should x be logged? (x axis to be printed as log scaled)
  # proportion - does the y axis label describe proportions or mg/mm2 ?
  # indepvarType - standard, gap or leafrad (to use appropriate horizontal error bars - standard uses none)
  # labs - should contain: 
  #   a pretty form of the dependent variable (labs[1])
  #   the x axis label (labs[2])
  # outdir - directory to output to
  # golden ratio - plot dimensions golden ratio?
  
  require(lazyeval)
  
  if(proportion) {
    source('scripts/prep_data.R')
  } else {
    source('scripts/prep_data_mg_per_mm2.R')
  }
  
  if(proportion) {
    type = '_proportion_'
  } else {
    type = '_mg-per-m2_'
  }
  
  dep_means <- group_by(data, ID) %>%
    dplyr::summarise(calv = mean(calvin_cycle, na.rm=TRUE),
                     phot = mean(Photosystems, na.rm=TRUE)) %>%
    gather(key = 'funccat', value = 'protein_mean', -ID)
  
  dep_means <- group_by(data, ID) %>%
    dplyr::summarise(calv = SE(calvin_cycle),
                     phot = SE(Photosystems)) %>%
    gather(key = 'funccat', value = 'protein_SE', -ID) %>%
    full_join(dep_means, by = c('ID', 'funccat')) %>%
    full_join(data, by = 'ID') %>%
    distinct(ID, funccat, .keep_all=TRUE)
  
  if(logx) {
    model_calv <- summary(lm(dep_means[dep_means$funccat %in% 'calv',]$protein_mean ~ log10(dep_means[dep_means$funccat %in% 'calv',][[indepvar]])))
    sig_calv <- summary(lm(dep_means[dep_means$funccat %in% 'calv',]$protein_mean ~ log10(dep_means[dep_means$funccat %in% 'calv',][[indepvar]])))$coefficients[,4][2] < 0.05
    model_phot <- summary(lm(dep_means[dep_means$funccat %in% 'phot',]$protein_mean ~ log10(dep_means[dep_means$funccat %in% 'phot',][[indepvar]])))
    sig_phot <- summary(lm(dep_means[dep_means$funccat %in% 'phot',]$protein_mean ~ log10(dep_means[dep_means$funccat %in% 'phot',][[indepvar]])))$coefficients[,4][2] < 0.05
    
  } else {
    model_calv <- summary(lm(dep_means[dep_means$funccat %in% 'calv',]$protein_mean ~ dep_means[dep_means$funccat %in% 'calv',][[indepvar]]))
    sig_calv <- summary(lm(dep_means[dep_means$funccat %in% 'calv',]$protein_mean ~ dep_means[dep_means$funccat %in% 'calv',][[indepvar]]))$coefficients[,4][2] < 0.05
    model_phot <- summary(lm(dep_means[dep_means$funccat %in% 'phot',]$protein_mean ~ dep_means[dep_means$funccat %in% 'phot',][[indepvar]]))
    sig_phot <- summary(lm(dep_means[dep_means$funccat %in% 'phot',]$protein_mean ~ dep_means[dep_means$funccat %in% 'phot',][[indepvar]]))$coefficients[,4][2] < 0.05
  }
  
  
  # directory to save to
  
  
  mainDir <- getwd()

 if (!file.exists(outDir)){
    dir.create(file.path(mainDir, outDir), showWarnings = FALSE, recursive = TRUE)
  }
  
  if(goldenRatio) {
    Width = 6*1.618
  } else {
    Width = 6
  }
  
  
  if(fileType == 'svg'){
    svg(paste(outDir, '/',indepvar, type, '.svg', sep = ""), height = 6, width = Width)
  } else {
    if(fileType == 'tiff') {
      tiff(paste(outDir, '/', indepvar, type, '.tiff', sep = ""), height = 6, width = Width, units = 'in', res = 300)
    }
  }
  
  # plot object
  
  dep_means$sig <- NA
  dep_means[dep_means$funccat == 'calv',]$sig <- sig_calv
  dep_means[dep_means$funccat == 'phot',]$sig <- sig_phot
  
  p <- ggplot(data = dep_means, aes(y = protein_mean, x = dep_means[[indepvar]])) + geom_point(aes(colour = funccat),size = 2)
  
  p <- p + geom_smooth(data = dep_means[dep_means$sig == TRUE,], aes(y = protein_mean, x = get(indepvar), colour = funccat), method = 'lm', se = F, size = 0.5)
  
  if(proportion) {
    p <- p + xlab(labs[2]) +  ylab(paste(labs[1], ' (proportion)', sep = ""))
  } else {
    p <- p + xlab(labs[2]) +  ylab(bquote(.(labs[1]) ~ '(mg / ' ~ m^{2} ~ ')'))
    
  }
  
  #p <- p + expand_limits(y=0,x=0)
  
  p <- p + theme_classic()
  p <- p + theme(legend.title=element_blank(),
                 legend.position="bottom",
                 axis.line = element_line(colour = "black"),
                 text = element_text(size = 18))
  
  number_ticks <- function(n) {function(limits) pretty(limits, n)}
  
  if(logx) {
    
    errorbar_width <- (log10(max(data[[indepvar]], na.rm=TRUE)) - log10(min(data[[indepvar]], na.rm=TRUE))) / 50
    p <- p + geom_errorbar(aes(colour = funccat, ymin = protein_mean - protein_SE, ymax = protein_mean + protein_SE), width = errorbar_width, alpha = 0.8)
    
    if(indepvarType == 'gap') {
      
      p <- p + geom_errorbarh(aes(xmin = gap_mean - gap_SE, xmax = gap_mean + gap_SE), alpha = 0.6)
      
    } else { if(indepvarType == 'leafrad') {
      
      p <- p + geom_errorbarh(aes(xmin = leafrad_mean - leafrad_SE, xmax = leafrad_mean + leafrad_SE), alpha = 0.6)
      
    }
    }
    
    p <- p + scale_x_continuous(breaks=scales::pretty_breaks(n=10),trans='log10')
    
  } else {
    
    errorbar_width <- (max(data[[indepvar]], na.rm=TRUE) - min(data[[indepvar]], na.rm=TRUE)) / 50
    p <- p + geom_errorbar(aes(colour = funccat, ymin = protein_mean - protein_SE, ymax = protein_mean + protein_SE), width = errorbar_width, alpha = 0.8)
    
    if(indepvarType == 'gap') {
      
      p <- p + geom_errorbarh(aes(xmin = gap_mean - gap_SE, xmax = gap_mean + gap_SE), alpha = 0.6)
      
    } else { if(indepvarType == 'leafrad') {
      
      p <- p + geom_errorbarh(aes(xmin = leafrad_mean - leafrad_SE, xmax = leafrad_mean + leafrad_SE), alpha = 0.6)
      
      
    } else { if(indepvarType == 'LMA') {
      
      p <- p + geom_errorbarh(aes(xmin = LMA_mean - LMA_SE, xmax = LMA_mean + LMA_SE), alpha = 0.6)
      
    } else { if(indepvarType == 'Narea') {
      
      p <- p + geom_errorbarh(aes(xmin = Narea_mean - Narea_SE, xmax = Narea_mean + Narea_SE), alpha = 0.6)
      
    }
      
    }
      
    }
      
    }
    
    p <- p + scale_x_continuous(breaks=scales::pretty_breaks(n=10))
    
  }
  
  print(p)
  
  dev.off()
  
}


agg_plot_save(indepvar = 'prec', logx = TRUE, proportion = FALSE, indepvarType = 'standard', 
              labs =  c('Protein amount', 'Mean annual precip. (mm)'), outDir = 'output/figures/20170501/tiff', goldenRatio = TRUE)

agg_plot_save(indepvar = 'prec', logx = TRUE, proportion = TRUE, indepvarType = 'standard', 
              labs =  c('Protein amount', 'Mean annual precip. (mm)'), outDir = 'output/figures/20170501/tiff', goldenRatio = TRUE)


agg_plot_save(indepvar = 'gap_mean', logx = TRUE, proportion = FALSE, indepvarType = 'gap', 
              labs =  c('Protein amount', 'Canopy openness (%)'), outDir = 'output/figures/20170501/tiff', goldenRatio = TRUE)

agg_plot_save(indepvar = 'gap_mean', logx = TRUE, proportion = TRUE, indepvarType = 'gap', 
              labs =  c('Protein amount', 'Canopy openness (%)'), outDir = 'output/figures/20170501/tiff', goldenRatio = TRUE)

