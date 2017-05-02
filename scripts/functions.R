# sourced by binProteins.R

require(lazyeval)
require(readr)
require(ggplot2)
require(vegan)
library(readr)
require(plyr)
require(reshape2)
require(dplyr)

SE <- function(x) {
  x <- sd(x, na.rm=TRUE) / sqrt(length(x))
}

CV <- function(x){
  sqrt(var(x))/mean(x)
}

getTotalProtein <- function(protein_samples) {
  
  protein_samples_long <- melt(protein_samples)
  total_protein <- ddply(protein_samples_long, .(variable), summarise, total_protein = sum(value, na.rm=TRUE))
  names(total_protein) <- c('sample', 'total_protein')                       
  
  return(total_protein)
  
}

getProteinBins <- function(protein_samples, mercator) {
  
  protein_samples$Protein <- tolower(protein_samples$Protein)
  mercator$IDENTIFIER <- tolower(mercator$IDENTIFIER)
  
  protein_samples$BINCODE <- NA
  protein_samples$NAME <- NA
  
  for(i in 1:length(protein_samples$Protein)) {
    merc_row <- grep(protein_samples$Protein[i], mercator$IDENTIFIER, fixed = TRUE)
    protein_samples$BINCODE[i] <- paste0(mercator[merc_row,]$BINCODE, collapse = ", ")
    protein_samples$NAME[i] <- paste0(mercator[merc_row,]$NAME, collapse = ", ")
  }
  
  return(protein_samples)
  
}


top2 <- function(x) {
  sum(sort(x, decreasing = TRUE)[1:2], na.rm=TRUE)
}

top3 <- function(x) {
  mean(sort(x, decreasing = TRUE)[1:3], na.rm=TRUE)
}

top2avg <- function(x) {
  mean(sort(x, decreasing = TRUE)[1:2], na.rm=TRUE)
}



agg_plot <- function(data, depvar, indepvar, logx = FALSE, labs) {
  
  dep_means <- data %>%
    group_by(ID) %>%
    summarise_(mean = interp(~mean(var, na.rm=TRUE), var = as.name(depvar)),
               SE = lazyeval::interp(~SE(var), var = as.name(depvar))) %>%
    full_join(data, by = 'ID')  %>%
    distinct(mean, .keep_all=TRUE) 
  
  if(logx) {
    
    model <- summary(lm(dep_means$mean ~ log10(dep_means[[indepvar]])))
    
  } else {
    
    model <- summary(lm(dep_means$mean ~ dep_means[[indepvar]]))
    
  }
  
  
              
  p <- ggplot(dep_means, aes(y = mean, x = dep_means[[indepvar]])) + geom_point(size = 2)
  p <- p + geom_smooth(method = 'lm', se = F)
  p <- p + xlab(labs[1]) + ggtitle(depvar) + ylab(paste('Species mean of [', depvar, ' proteins] (proportion)', sep = ""))
  p <- p + expand_limits(y=0)
  
    p <- p + annotate("text" , x = max(dep_means[[indepvar]])*0.3, y = (max(dep_means$mean)+max(dep_means$SE)), label = paste('R2 =', round(model$r.squared,3), sep = " "))
    p <- p + annotate("text" , x = max(dep_means[[indepvar]])*0.3, y = (max(dep_means$mean)+max(dep_means$SE))*0.9, label = paste('pval =', round(model$coefficients[,4][2],2), sep = " "))
  
  
  p <- p + theme_bw()
  p <- p + theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), 
                 legend.title=element_blank(),
                 legend.position="bottom",
                 axis.line = element_line(colour = "black"))
  
  number_ticks <- function(n) {function(limits) pretty(limits, n)}
  
  if(logx) {
    
    errorbar_width <- (log10(max(data[[indepvar]], na.rm=TRUE)) - log10(min(data[[indepvar]], na.rm=TRUE))) / 50
    p <- p + geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = errorbar_width, alpha = 0.8)
    
    p <- p + scale_x_continuous(breaks=scales::pretty_breaks(n=10),trans='log10')
    
    # print(summary(lm(dep_means$mean ~ log10(dep_means[[indepvar]]))))
    
  } else {
    
    errorbar_width <- (max(data[[indepvar]], na.rm=TRUE) - min(data[[indepvar]], na.rm=TRUE)) / 50
    p <- p + geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = errorbar_width, alpha = 0.8)
    
    p <- p + scale_x_continuous(breaks=scales::pretty_breaks(n=10))
    
    # print(summary(lm(dep_means$mean ~ dep_means[[indepvar]])))
    
  }
  
  print(p)
  
  
}


regression_output_agg <- function(df, indepvar, logx) {
  
  # called by regression_agg, creates table of model stats
  
  df <- df[names(df) %in% names(protein_D14[,2:ncol(protein_D14)]) | names(df) %in% indepvar]
  
  a <- data.frame()
  my.list <- vector("list", length(df)-2)
  
  for(i in 1:ncol(df)) {
    
    protname <- names(df[i])
    
    if(logx) {
      
      model <- summary(lm(df[[protname]] ~ log10(df[[indepvar]])))
      
    } else {
      
      model <- summary(lm(df[[protname]] ~ df[[indepvar]]))
      
    }
    
    model.stats <-  data.frame(cbind(protname,round(model$r.squared,3),
                                     round(model$coefficients[,4][2],2)))
    
    names(model.stats) <- c('submodel','R2','pval') 
    rownames(model.stats) <- NULL
    
 #   browser()
    
    my.list[[i]] <- model.stats
    
  }
  
  a <- rbind(a, do.call(rbind, my.list))
  
  a <-  filter(a, submodel != indepvar)
  
  a$p.adj <- p.adjust(as.numeric(as.character(a$pval)), method = "BH")
  
  a$R2 <- signif(as.numeric(as.character(a$R2)), 3)
  a$pval <- signif(as.numeric(as.character(a$pval)), 5)
  a$p.adj <- signif(as.numeric(as.character(a$p.adj)), 5)
  a$submodel <- as.character(a$submodel)
  
  a <- a[order(a$R2, decreasing=T),]
  
  return(a)
  
}


regression_agg <- function(data, indepvar, logx = FALSE) {
  
  # recombobulates data to aggregate by ID
  # then calls regression_output_agg to output table
  
  a <- data.frame()
  my.list <- vector("list", length(names(protein_stand_D14[,2:ncol(protein_D14)])))
  
  for(i in names(protein_stand_D14[,2:ncol(protein_stand_D14)])) {
    
    depvar = i 
    
    dep_means <- data %>%
      group_by(ID) %>%
      summarise_(mean = interp(~mean(var, na.rm=TRUE), var = as.name(depvar)))
    
    dep_means$bin_arch_name <- depvar
    
    dep_means <- na.omit(dep_means)
  
    my.list[[i]] <- dep_means
    
  }
  
  a <- rbind(a, do.call(rbind, my.list))
  
  a <- tidyr::spread(a, bin_arch_name, mean)
  
  climate_locs$ID <- NULL
  
  clim  <- merge(climate_locs, replicates, by = c('sample', 'Latitude', 'Longitude'))
  

  bla <- regression_output_agg(merge(a, clim, by = 'ID'), indepvar, logx=logx)
  
  return(bla)
  
}


agg_plot_gap <- function(data, depvar, indepvar, logx=FALSE, labs) {
  
  dep_means <- data %>%
    group_by(ID) %>%
    summarise_(mean = interp(~mean(var, na.rm=TRUE), var = as.name(depvar)),
               SE = lazyeval::interp(~SE(var), var = as.name(depvar))) %>%
    full_join(data, by = 'ID')  %>%
    distinct(mean, .keep_all=TRUE) 
  
  if(logx) {
    
    model <- summary(lm(dep_means$mean ~ log10(dep_means[[indepvar]])))
    
  } else {
    
    model <- summary(lm(dep_means$mean ~ dep_means[[indepvar]]))
    
  }
  
  p <- ggplot(dep_means, aes(y = mean, x = dep_means[[indepvar]])) + geom_point(size = 2)
  p <- p + geom_smooth(method = 'lm', se = F)
  p <- p + xlab(labs[1]) + ggtitle(depvar) + ylab(paste('Species mean of [', depvar, ' proteins] (proportion)', sep = ""))
  p <- p + expand_limits(y=0)
  p <- p + annotate("text" , x = max(dep_means[[indepvar]])*0.4, y = (max(dep_means$mean)+max(dep_means$SE)), label = paste('R2 =', round(model$r.squared,3), sep = " "))
  p <- p + annotate("text" , x = max(dep_means[[indepvar]])*0.4, y = (max(dep_means$mean)+max(dep_means$SE))*0.9, label = paste('pval =', round(model$coefficients[,4][2],2), sep = " "))
  p <- p + theme_bw()
  p <- p + theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), 
                 legend.title=element_blank(),
                 legend.position="bottom",
                 axis.line = element_line(colour = "black"))
  
  number_ticks <- function(n) {function(limits) pretty(limits, n)}
  
  if(logx) {
    
    errorbar_width <- (log10(max(data[[indepvar]], na.rm=TRUE)) - log10(min(data[[indepvar]], na.rm=TRUE))) / 50
    p <- p + geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = errorbar_width, alpha = 0.8)
    p <- p + geom_errorbarh(aes(xmin = gap_mean - gap_SE, xmax = gap_mean + gap_SE), alpha = 0.6)
    p <- p + scale_x_continuous(breaks=scales::pretty_breaks(n=10),trans='log10')
    
    # print(summary(lm(dep_means$mean ~ log10(dep_means[[indepvar]]))))
    
  } else {
    
    errorbar_width <- (max(data[[indepvar]], na.rm=TRUE) - min(data[[indepvar]], na.rm=TRUE)) / 50
    p <- p + geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = errorbar_width, alpha = 0.8)
    p <- p + geom_errorbarh(aes(xmin = gap_mean - gap_SE, xmax = gap_mean + gap_SE), alpha = 0.6)
    
    p <- p + scale_x_continuous(breaks=scales::pretty_breaks(n=10))
    
    # print(summary(lm(dep_means$mean ~ dep_means[[indepvar]])))
    
  }
  

  print(p)
  
  
}


agg_plot_leafrad <- function(data, depvar, indepvar, logx=FALSE, labs) {
  
  dep_means <- data %>%
    group_by(ID) %>%
    summarise_(mean = interp(~mean(var, na.rm=TRUE), var = as.name(depvar)),
               SE = lazyeval::interp(~SE(var), var = as.name(depvar))) %>%
    full_join(data, by = 'ID')  %>%
    distinct(mean, .keep_all=TRUE) 
  dep_means <- na.omit(dep_means)
  
  if(logx) {
    
    model <- summary(lm(dep_means$mean ~ log10(dep_means[[indepvar]])))
    
  } else {
    
    model <- summary(lm(dep_means$mean ~ dep_means[[indepvar]]))
    
  }
  
  p <- ggplot(dep_means, aes(y = mean, x = dep_means[[indepvar]])) + geom_point(size = 2)
  p <- p + geom_smooth(method = 'lm', se = F)
  p <- p + xlab(labs[1]) + ggtitle(depvar) + ylab(paste('Species mean of [', depvar, ' proteins] (proportion)', sep = ""))
  p <- p + expand_limits(y=0)
  p <- p + annotate("text" , x = max(dep_means[[indepvar]])*0.33, y = (max(dep_means$mean)+max(dep_means$SE)), label = paste('R2 =', round(model$r.squared,3), sep = " "))
  p <- p + annotate("text" , x = max(dep_means[[indepvar]])*0.33, y = (max(dep_means$mean)+max(dep_means$SE))*0.9, label = paste('pval =', round(model$coefficients[,4][2],2), sep = " "))
  p <- p + theme_bw()
  p <- p + theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), 
                 legend.title=element_blank(),
                 legend.position="bottom",
                 axis.line = element_line(colour = "black"))
  
  number_ticks <- function(n) {function(limits) pretty(limits, n)}
  
  if(logx) {
    
    errorbar_width <- (log10(max(data[[indepvar]], na.rm=TRUE)) - log10(min(data[[indepvar]], na.rm=TRUE))) / 50
    p <- p + geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = errorbar_width, alpha = 0.8)
    p <- p + geom_errorbarh(aes(xmin = leafrad_mean - leafrad_SE, xmax = leafrad_mean + leafrad_SE), alpha = 0.6)
    p <- p + scale_x_continuous(breaks=scales::pretty_breaks(n=10),trans='log10')
    
    # print(summary(lm(dep_means$mean ~ log10(dep_means[[indepvar]]))))
    
  } else {
    
    errorbar_width <- (max(data[[indepvar]], na.rm=TRUE) - min(data[[indepvar]], na.rm=TRUE)) / 50
    p <- p + geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = errorbar_width, alpha = 0.8)
    p <- p + geom_errorbarh(aes(xmin = leafrad_mean - leafrad_SE, xmax = leafrad_mean + leafrad_SE), alpha = 0.6)
    
    p <- p + scale_x_continuous(breaks=scales::pretty_breaks(n=10))
    
    # print(summary(lm(dep_means$mean ~ dep_means[[indepvar]])))
    
  }
  
  print(p)
  
  
}





# this function outputs a pretty aggregated plot to your directory of choice

agg_plot_save <- function(proportion, depvar, indepvar, logx = FALSE, indepvarType, labs, outDir, fileType = 'tiff', goldenRatio = FALSE) {
  
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

  dep_means <- data %>%
    group_by(ID) %>%
    summarise_(mean = interp(~mean(var, na.rm=TRUE), var = as.name(depvar)),
               SE = lazyeval::interp(~SE(var), var = as.name(depvar))) %>%
    full_join(data, by = 'ID')  %>%
    distinct(mean, .keep_all=TRUE) 
  
  if(logx) {
    
    model <- summary(lm(dep_means$mean ~ log10(dep_means[[indepvar]])))
    
  } else {
    
    model <- summary(lm(dep_means$mean ~ dep_means[[indepvar]]))
    
  }
  
  
  # directory to save to
  
  
  mainDir <- getwd()
  outDir
  
  if (!file.exists(outDir)){
    dir.create(file.path(mainDir, outDir), showWarnings = FALSE, recursive = TRUE)
  }
  
  if(goldenRatio) {
    Width = 6*1.618
  } else {
    Width = 6
  }
    
  
  if(fileType == 'svg'){
    svg(paste(outDir, '/', depvar, '_vs_', indepvar, type, '_R2-', round(model$r.squared,3), '_pval-', round(model$coefficients[,4][2],2),'.svg', sep = ""), height = 6, width = Width)
  } else {
    if(fileType == 'tiff') {
      tiff(paste(outDir, '/', depvar, '_vs_', indepvar, type, '_R2-', round(model$r.squared,3), '_pval-', round(model$coefficients[,4][2],2),'.tiff', sep = ""), height = 6, width = Width, units = 'in', res = 300)
    }
  }
    
  # plot object
  
  p <- ggplot(dep_means, aes(y = mean, x = dep_means[[indepvar]])) + geom_point(size = 2)
  
  if(model$coefficients[,4][2] < 0.05) {
    p <- p + geom_smooth(method = 'lm', se = F, colour = 'black', size = 0.5)
  }
  
  if(proportion) {
    p <- p + xlab(labs[2]) +  ylab(paste(labs[1], ' (proportion)', sep = ""))
  } else {
    p <- p + xlab(labs[2]) +  ylab(bquote(.(labs[1]) ~ '(mg / ' ~ m^{2} ~ ')'))
   
  }
  
  #bquote('Assimilation ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')  
  
  p <- p + expand_limits(y=0,x=0)
  
  p <- p + theme_classic()
  p <- p + theme(legend.title=element_blank(),
                 legend.position="bottom",
                 axis.line = element_line(colour = "black"),
                 text = element_text(size = 18))
  
  number_ticks <- function(n) {function(limits) pretty(limits, n)}
  
  if(logx) {
    
    errorbar_width <- (log10(max(data[[indepvar]], na.rm=TRUE)) - log10(min(data[[indepvar]], na.rm=TRUE))) / 50
    p <- p + geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = errorbar_width, alpha = 0.8)
    
    if(indepvarType == 'gap') {
      
      p <- p + geom_errorbarh(aes(xmin = gap_mean - gap_SE, xmax = gap_mean + gap_SE), alpha = 0.6)
      
    } else { if(indepvarType == 'leafrad') {
      
      p <- p + geom_errorbarh(aes(xmin = leafrad_mean - leafrad_SE, xmax = leafrad_mean + leafrad_SE), alpha = 0.6)
      
    }
    }
    
    p <- p + scale_x_continuous(breaks=scales::pretty_breaks(n=10),trans='log10')
    
  } else {
    
    errorbar_width <- (max(data[[indepvar]], na.rm=TRUE) - min(data[[indepvar]], na.rm=TRUE)) / 50
    p <- p + geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = errorbar_width, alpha = 0.8)
    
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



# this function outputs a pretty aggregated plot to your directory of choice, with stacked points for different funccats

agg_plot_save_combined <- function(proportion, indepvar, logx = FALSE, indepvarType, labs, outDir, fileType = 'tiff', goldenRatio = FALSE, total_prot = FALSE) {
  
  # data - the output of scripts/prep_data.R or sources/prep_data_mg_per_mm2.R (contains all climate/env variables and protein amounts)
  # depvar - the dependent variable (e.g. 'Photosystems')
  # indepvar - the independent variable (e.g. 'prec')
  # logx - should x be logged? (x axis to be printed as log scaled)
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
  
  if(proportion | total_prot) { # if proportional or total_protein to be indepvar, only aggregate calvin cycle and photosystems
    
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
    
  } else { # else add in total protein
    
    dep_means <- group_by(data, ID) %>%
      dplyr::summarise(calv = mean(calvin_cycle, na.rm=TRUE),
                       phot = mean(Photosystems, na.rm=TRUE),
                       total_prot = mean(total_protein, na.rm=TRUE)) %>%
      gather(key = 'funccat', value = 'protein_mean', -ID)
    
    dep_means <- group_by(data, ID) %>%
      dplyr::summarise(calv = SE(calvin_cycle),
                       phot = SE(Photosystems),
                       total_prot = SE(total_protein)) %>%
      gather(key = 'funccat', value = 'protein_SE', -ID) %>%
      full_join(dep_means, by = c('ID', 'funccat')) %>%
      full_join(data, by = 'ID') %>%
      distinct(ID, funccat, .keep_all=TRUE)
    
  }
  
  
  if(logx) { # get significance TRUE/FALSE for model, depending on if logx or not
    sig_calv <- summary(lm(dep_means[dep_means$funccat %in% 'calv',]$protein_mean ~ log10(dep_means[dep_means$funccat %in% 'calv',][[indepvar]])))$coefficients[,4][2] < 0.05
    sig_phot <- summary(lm(dep_means[dep_means$funccat %in% 'phot',]$protein_mean ~ log10(dep_means[dep_means$funccat %in% 'phot',][[indepvar]])))$coefficients[,4][2] < 0.05
    if(!proportion && !total_prot) {
      sig_total_prot <- summary(lm(dep_means[dep_means$funccat %in% 'total_prot',]$protein_mean ~ log10(dep_means[dep_means$funccat %in% 'total_prot',][[indepvar]])))$coefficients[,4][2] < 0.05
    }  
  } else {
    sig_calv <- summary(lm(dep_means[dep_means$funccat %in% 'calv',]$protein_mean ~ dep_means[dep_means$funccat %in% 'calv',][[indepvar]]))$coefficients[,4][2] < 0.05
    sig_phot <- summary(lm(dep_means[dep_means$funccat %in% 'phot',]$protein_mean ~ dep_means[dep_means$funccat %in% 'phot',][[indepvar]]))$coefficients[,4][2] < 0.05
    if(!proportion && !total_prot) {
      sig_total_prot <- summary(lm(dep_means[dep_means$funccat %in% 'total_prot',]$protein_mean ~ dep_means[dep_means$funccat %in% 'total_prot',][[indepvar]]))$coefficients[,4][2] < 0.05
    }  
  }
  
  dep_means$sig <- NA
  dep_means[dep_means$funccat == 'calv',]$sig <- sig_calv
  dep_means[dep_means$funccat == 'phot',]$sig <- sig_phot
  if(!proportion && !total_prot) {
    dep_means[dep_means$funccat == 'total_prot',]$sig <- sig_total_prot
  }
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
  
  p <- p + scale_colour_manual(values = c('red','forestgreen','blue'))
  
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