# sourced by binProteins.R

require(readr)

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

multipleBins <- function(protein_samples) {
  
  protein_samples <- cbind(protein_samples, str_split_fixed(protein_samples$BINCODE, ", ", 5))
  
  bins <- protein_samples[,317:322]
  bins[] <- lapply(bins,as.character)
  bins[,2:6][!(bins[,2:6]) == ""] <- 1
  bins[,2:6] <- lapply(bins[,2:6], as.numeric)
  bins$bins <- rowSums(bins[,2:6], na.rm=TRUE)
  
  protein_samples$num_bins <- bins$bins
  
  for(i in unique(protein_samples$num_bins)) {
    protein_samples[,2:315][which(protein_samples1$num_bins == i),] <- protein_samples1[,2:315][which(protein_samples1$num_bins == i),] / i
  }
  
  return(protein_samples[,1:317])
  
}

populateProteinBins <- function(protein_samples, bin_arch.list) {
  
  protein_samples$bin_arch <- NA
  
  for(i in 1:length(bin_arch.list)) {
    
    x <-  bin_arch.list[[i]][1]
    rows1 <- nrow(protein_samples[grep(bin_arch.list[[i]][1], protein_samples$BINCODE, fixed=TRUE),])
    rows2 <- nrow(protein_samples[grep(bin_arch.list[[i]][2], protein_samples$BINCODE, fixed=TRUE),])
    
    if(nrow(protein_samples[grep(bin_arch.list[[i]][1], protein_samples$BINCODE, fixed=TRUE),]) > 0) {
      protein_samples[grep(bin_arch.list[[i]][1], protein_samples$BINCODE, fixed=TRUE),]$bin_arch <- bin_arch.list[[i]][1]
    }
    
    if(nrow(protein_samples[grep(bin_arch.list[[i]][2], protein_samples$BINCODE, fixed=TRUE),]) > 0) {
      protein_samples[grep(bin_arch.list[[i]][2], protein_samples$BINCODE, fixed=TRUE),]$bin_arch <- bin_arch.list[[i]][1]
      
    } 
    
  }
  
  
  protein_bins <- melt(protein_samples, id = c('Protein', 'bin_arch'))
  names(protein_bins) <- c('Protein', 'bin_arch', 'sample', 'value')
  
  # protein_bins1 <- ddply(protein_bins, .(bin_arch, sample), summarise, sum = sum(as.numeric(value), na.rm=TRUE), nprot=length(value))
  
  protein_bins <- protein_bins %>% dplyr::group_by(bin_arch, sample) %>% dplyr::summarise(sum = sum(as.numeric(value), na.rm=TRUE), nprot=length(value)) %>% ungroup(.)
  
  protein_bins <- na.omit(protein_bins)
  
  # recombobulate arch categories (where lower levels need to be incorporated into an higher level)
  protein_bins$sample <- as.character(protein_bins$sample)
  
  protein_bins[protein_bins$bin_arch == '_1.1.1_',]$sum <-     protein_bins[protein_bins$bin_arch == '_1.1.1_',]$sum +
    protein_bins[protein_bins$bin_arch == '_1.1.1.1_',]$sum 
  
  protein_bins[protein_bins$bin_arch == '_1.1.2_',]$sum <-     protein_bins[protein_bins$bin_arch == '_1.1.2_',]$sum +
    protein_bins[protein_bins$bin_arch == '_1.1.2.1_',]$sum 
  
  protein_bins[protein_bins$bin_arch == '_1.1_',]$sum <-     protein_bins[protein_bins$bin_arch == '_1.1_',]$sum +
    protein_bins[protein_bins$bin_arch == '_1.1.1_',]$sum +
    protein_bins[protein_bins$bin_arch == '_1.1.2_',]$sum +
    protein_bins[protein_bins$bin_arch == '_1.1.3_',]$sum +
    protein_bins[protein_bins$bin_arch == '_1.1.4_',]$sum +
    protein_bins[protein_bins$bin_arch == '_1.1.5_',]$sum +
    protein_bins[protein_bins$bin_arch == '_1.1.6_',]$sum 
  
  protein_bins[protein_bins$bin_arch == '_1.3_',]$sum <-     protein_bins[protein_bins$bin_arch == '_1.3_',]$sum +
    protein_bins[protein_bins$bin_arch == '_1.3.1_',]$sum +
    protein_bins[protein_bins$bin_arch == '_1.3.2_',]$sum
  
  protein_bins[protein_bins$bin_arch == '_29_',]$sum <- protein_bins[protein_bins$bin_arch == '_29_',]$sum +
    protein_bins[protein_bins$bin_arch == '_29.6_',]$sum
  
  protein_bins[protein_bins$bin_arch == '_20_',]$sum <- protein_bins[protein_bins$bin_arch == '_20_',]$sum +
    protein_bins[protein_bins$bin_arch == '_20.1_',]$sum +
    protein_bins[protein_bins$bin_arch == '_20.2_',]$sum +
    protein_bins[protein_bins$bin_arch == '_20.2.1_',]$sum
  
  #depends if we want to include heat stress with abiotic stress
  protein_bins[protein_bins$bin_arch == '_20.2_',]$sum <- protein_bins[protein_bins$bin_arch == '_20.2_',]$sum +
    protein_bins[protein_bins$bin_arch == '_20.2.1_',]$sum 
  
  # add some summed bins for special cases
  
  rubisco <- protein_bins[protein_bins$bin_arch %in% c('_1.3.1_', '_1.3.2_'),]
  rubisco <- ddply(na.omit(rubisco), .(sample), summarise, sum = sum(sum))
  rubisco$bin_arch <- '_1.3.1_+_1.3.2_'
  
  photosystems <- protein_bins[protein_bins$bin_arch %in% c('_1.1.1_', '_1.1.2_'),]
  photosystems <- ddply(na.omit(photosystems), .(sample), summarise, sum = sum(sum))
  photosystems$bin_arch <- '_1.1.1_+_1.1.2_'
  
  TCA_orgtrans <- protein_bins[protein_bins$bin_arch %in% c('_8.1_', '_8.2_'),]
  TCA_orgtrans <- ddply(na.omit(TCA_orgtrans), .(sample), summarise, sum = sum(sum))
  TCA_orgtrans$bin_arch <- '_8.1_+_8.2_'
  
  redox <- protein_bins[protein_bins$bin_arch %in% c('_21_', '_26.9_'),]
  redox <- ddply(na.omit(redox), .(sample), summarise, sum = sum(sum))
  redox$bin_arch <- '_21_+_26.9_'
  
  electron_transport_minATPsynth <- protein_bins[protein_bins$bin_arch %in% c('_1.1.3_', '_1.1.5_', '_1.1.6_'),]
  electron_transport_minATPsynth <- ddply(na.omit(electron_transport_minATPsynth), .(sample), summarise, sum = sum(sum))
  electron_transport_minATPsynth$bin_arch <- '_1.1.3_+_1.1.5_+_1.1.6_'
  
  CHO_metabolism <- protein_bins[protein_bins$bin_arch %in% c('_2_', '_3_'),]
  CHO_metabolism <- ddply(na.omit(CHO_metabolism), .(sample), summarise, sum = sum(sum))
  CHO_metabolism$bin_arch <- '_2_+_3_'
  
  LHC2 <- protein_bins[protein_bins$bin_arch %in% c('_1.1.1.1_'),]
  LHC2 <- ddply(na.omit(LHC2), .(sample), summarise, sum = sum(sum))
  LHC2 <- dplyr::arrange(LHC2, sample)
  PSII_min_LHC2 <- protein_bins[protein_bins$bin_arch %in% c('_1.1.1_'),]
  PSII_min_LHC2 <- ddply(na.omit(PSII_min_LHC2), .(sample), summarise, sum = sum(sum)) 
  PSII_min_LHC2 <- dplyr::arrange(PSII_min_LHC2, sample)
  PSII_min_LHC2$sum <- PSII_min_LHC2$sum - LHC2$sum # subtract sum of LHC2 protein amounts
  PSII_min_LHC2$bin_arch <- '_1.1.1_-_1.1.1.1_'
  
  LHC1 <- protein_bins[protein_bins$bin_arch %in% c('_1.1.2.1_'),]
  LHC1 <- ddply(na.omit(LHC1), .(sample), summarise, sum = sum(sum))
  LHC1 <- dplyr::arrange(LHC1, sample)
  PSI_min_LHC1 <- protein_bins[protein_bins$bin_arch %in% c('_1.1.2_'),]
  PSI_min_LHC1 <- ddply(na.omit(PSI_min_LHC1), .(sample), summarise, sum = sum(sum))
  PSI_min_LHC1 <- dplyr::arrange(PSI_min_LHC1, sample)
  PSI_min_LHC1$sum <- PSI_min_LHC1$sum - LHC1$sum # subtract sum of LHC1 protein amounts
  PSI_min_LHC1$bin_arch <- '_1.1.2_-_1.1.2.1_'
  
  
  protein_bins <- rbind(protein_bins[,c('bin_arch','sample','sum')], rubisco, photosystems, TCA_orgtrans, redox, electron_transport_minATPsynth, CHO_metabolism, PSII_min_LHC2, PSI_min_LHC1)
  
  rm(rubisco,photosystems,TCA_orgtrans,redox, electron_transport_minATPsynth, CHO_metabolism, PSII_min_LHC2, PSI_min_LHC1)
  protein_bins$sample <- as.character(protein_bins$sample)
  
  protein_bins[!protein_bins$sample %in% c('BINCODE','NAME'),]
  
  # merge in bin names from mercator_bins
  
  protein_bins$bin_arch <- gsub(" ", "", protein_bins$bin_arch) # remove _'s
  
  unique(protein_bins$bin_arch)[!unique(protein_bins$bin_arch) %in% mercator_bins$bin_arch]
  
  protein_bins <- protein_bins[!protein_bins$bin_arch %in% c('_1.3.1_', '1.3.2_'),]
  
  protein_bins$bin_arch_name <- NA
  
  
  for(i in 1:length(unique(protein_bins$bin_arch))) {
    #x <- paste(unique(protein_bins$bin_arch)[i], "_", sep = "")
    x <- unique(protein_bins$bin_arch)[i]
    y <- mercator_bins[grep(x, mercator_bins$bin_arch, fixed =TRUE),]$bin_name
    
    if (x %in% mercator_bins$bin_arch) {
      protein_bins[grep(x, paste(protein_bins$bin_arch, "_", sep = ""), fixed=TRUE),]$bin_arch_name <- y[1]
    } else {
      protein_bins[grep(x, paste(protein_bins$bin_arch, "_", sep = ""), fixed=TRUE),]$bin_arch_name <- NA
    }
    # browser()
  }
  
  return(protein_bins)
  
}



regression_output <- function(df, indepvar, logx) {
  
  #  df <- na.omit(merge(protein_stand_D14_age, climate_locs, all=TRUE))
  
  a <- data.frame()
  my.list <- vector("list", length(protein_stand_D14_age)-2)
  
  for(i in 3:length(protein_stand_D14_age)) {
    
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
    
    my.list[[i]] <- model.stats
    
  }
  
  a <- rbind(a, do.call(rbind, my.list))
  a$p.adj <- p.adjust(as.numeric(as.character(a$pval)), method = "BH")
  
  a$R2 <- signif(as.numeric(as.character(a$R2)), 3)
  a$pval <- signif(as.numeric(as.character(a$pval)), 5)
  a$p.adj <- signif(as.numeric(as.character(a$p.adj)), 5)
  a$submodel <- as.character(a$submodel)
  
  a <- a[order(a$R2, decreasing=T),]
  
  return(a)
  
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


populateProteinBins_mean <- function(protein_samples, bin_arch.list) {
  
  protein_samples$bin_arch <- NA
  
  for(i in 1:length(bin_arch.list)) {
    
    x <-  bin_arch.list[[i]][1]
    rows1 <- nrow(protein_samples[grep(bin_arch.list[[i]][1], protein_samples$BINCODE, fixed=TRUE),])
    rows2 <- nrow(protein_samples[grep(bin_arch.list[[i]][2], protein_samples$BINCODE, fixed=TRUE),])
    
    if(nrow(protein_samples[grep(bin_arch.list[[i]][1], protein_samples$BINCODE, fixed=TRUE),]) > 0) {
      protein_samples[grep(bin_arch.list[[i]][1], protein_samples$BINCODE, fixed=TRUE),]$bin_arch <- bin_arch.list[[i]][1]
    }
    
    if(nrow(protein_samples[grep(bin_arch.list[[i]][2], protein_samples$BINCODE, fixed=TRUE),]) > 0) {
      protein_samples[grep(bin_arch.list[[i]][2], protein_samples$BINCODE, fixed=TRUE),]$bin_arch <- bin_arch.list[[i]][1]
      
    } 
    
  }
  
  
  protein_bins <- melt(protein_samples, id = c('Protein', 'bin_arch'))
  names(protein_bins) <- c('Protein', 'bin_arch', 'sample', 'value')
  
  protein_bins <- ddply(protein_bins, .(bin_arch, sample), summarise, sum = mean(as.numeric(value), na.rm=TRUE), nprot=length(value)) # changed to mean, still called sum for hacking ease
  protein_bins <- na.omit(protein_bins)
  
  # recombobulate arch categories (where lower levels need to be incorporated into an higher level)
  protein_bins$sample <- as.character(protein_bins$sample)
  
  protein_bins[protein_bins$bin_arch == '_1.1_',]$sum <-     protein_bins[protein_bins$bin_arch == '_1.1_',]$sum +
    protein_bins[protein_bins$bin_arch == '_1.1.1_',]$sum +
    protein_bins[protein_bins$bin_arch == '_1.1.2_',]$sum +
    protein_bins[protein_bins$bin_arch == '_1.1.3_',]$sum +
    protein_bins[protein_bins$bin_arch == '_1.1.4_',]$sum +
    protein_bins[protein_bins$bin_arch == '_1.1.5_',]$sum +
    protein_bins[protein_bins$bin_arch == '_1.1.6_',]$sum 
  
  protein_bins[protein_bins$bin_arch == '_1.3_',]$sum <-     protein_bins[protein_bins$bin_arch == '_1.3_',]$sum +
    protein_bins[protein_bins$bin_arch == '_1.3.1_',]$sum +
    protein_bins[protein_bins$bin_arch == '_1.3.2_',]$sum
  
  protein_bins[protein_bins$bin_arch == '_29_',]$sum <- protein_bins[protein_bins$bin_arch == '_29_',]$sum +
    protein_bins[protein_bins$bin_arch == '_29.6_',]$sum
  
  protein_bins[protein_bins$bin_arch == '_20_',]$sum <- protein_bins[protein_bins$bin_arch == '_20_',]$sum +
    protein_bins[protein_bins$bin_arch == '_20.1_',]$sum +
    protein_bins[protein_bins$bin_arch == '_20.2_',]$sum +
    protein_bins[protein_bins$bin_arch == '_20.2.1_',]$sum
  
  #depends if we want to include heat stress with abiotic stress
  protein_bins[protein_bins$bin_arch == '_20.2_',]$sum <- protein_bins[protein_bins$bin_arch == '_20.2_',]$sum +
    protein_bins[protein_bins$bin_arch == '_20.2.1_',]$sum 
  
  # add some summed bins for special cases
  
  rubisco <- protein_bins[protein_bins$bin_arch %in% c('_1.3.1_', '_1.3.2_'),]
  rubisco <- ddply(na.omit(rubisco), .(sample), summarise, sum = sum(sum))
  rubisco$bin_arch <- '_1.3.1_+_1.3.2_'
  
  photosystems <- protein_bins[protein_bins$bin_arch %in% c('_1.1.1_', '_1.1.2_'),]
  photosystems <- ddply(na.omit(photosystems), .(sample), summarise, sum = sum(sum))
  photosystems$bin_arch <- '_1.1.1_+_1.1.2_'
  
  TCA_orgtrans <- protein_bins[protein_bins$bin_arch %in% c('_8.1_', '_8.2_'),]
  TCA_orgtrans <- ddply(na.omit(TCA_orgtrans), .(sample), summarise, sum = sum(sum))
  TCA_orgtrans$bin_arch <- '_8.1_+_8.2_'
  
  redox <- protein_bins[protein_bins$bin_arch %in% c('_21_', '_26.9_'),]
  redox <- ddply(na.omit(redox), .(sample), summarise, sum = sum(sum))
  redox$bin_arch <- '_21_+_26.9_'
  
  electron_transport_minATPsynth <- protein_bins[protein_bins$bin_arch %in% c('_1.1.3_', '_1.1.5_', '_1.1.6_'),]
  electron_transport_minATPsynth <- ddply(na.omit(electron_transport_minATPsynth), .(sample), summarise, sum = sum(sum))
  electron_transport_minATPsynth$bin_arch <- '_1.1.3_+_1.1.5_+_1.1.6_'
  
  protein_bins <- rbind(protein_bins[,c('bin_arch','sample','sum')], rubisco, photosystems, TCA_orgtrans, redox, electron_transport_minATPsynth)
  rm(rubisco,photosystems,TCA_orgtrans,redox, electron_transport_minATPsynth)
  protein_bins$sample <- as.character(protein_bins$sample)
  
  # merge in bin names from mercator_bins
  
  protein_bins$bin_arch <- gsub(" ", "", protein_bins$bin_arch) # remove _'s
  
  unique(protein_bins$bin_arch)[!unique(protein_bins$bin_arch) %in% mercator_bins$bin_arch]
  
  protein_bins <- protein_bins[!protein_bins$bin_arch %in% c('_1.3.1_', '1.3.2_'),]
  
  protein_bins$bin_arch_name <- NA
  
  
  for(i in 1:length(unique(protein_bins$bin_arch))) {
    #x <- paste(unique(protein_bins$bin_arch)[i], "_", sep = "")
    x <- unique(protein_bins$bin_arch)[i]
    y <- mercator_bins[grep(x, mercator_bins$bin_arch, fixed =TRUE),]$bin_name
    
    if (x %in% mercator_bins$bin_arch) {
      protein_bins[grep(x, paste(protein_bins$bin_arch, "_", sep = ""), fixed=TRUE),]$bin_arch_name <- y[1]
    } else {
      protein_bins[grep(x, paste(protein_bins$bin_arch, "_", sep = ""), fixed=TRUE),]$bin_arch_name <- NA
    }
    # browser()
  }
  
  return(protein_bins)
  
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








agg_plot_ <- function(data, depvar, indepvar, logx = FALSE, labs) {
  
  dep_means <- data %>%
    group_by(ID) %>%
    summarise_(mean = interp(~mean(var, na.rm=TRUE), var = as.name(depvar)),
               SE = lazyeval::interp(~SE(var), var = as.name(depvar))) %>%
    full_join(data, by = 'ID')  %>%
    distinct(mean, .keep_all=TRUE) 
  
  p <- ggplot(dep_means, aes(y = SE, x = dep_means[[indepvar]])) + geom_point(size = 2)
  p <- p + geom_smooth(method = 'lm', se = F)
  p <- p + xlab(labs[1]) + ggtitle(depvar) + ylab(paste('Species mean of [', depvar, ' proteins] (proportion)', sep = ""))
  p <- p + expand_limits(y=0)
  p <- p + theme_bw()
  p <- p + theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), 
                 legend.title=element_blank(),
                 legend.position="bottom",
                 axis.line = element_line(colour = "black"))
  
  number_ticks <- function(n) {function(limits) pretty(limits, n)}
  
  if(logx) {
    
   # errorbar_width <- (log10(max(data[[indepvar]], na.rm=TRUE)) - log10(min(data[[indepvar]], na.rm=TRUE))) / 50
  #  p <- p + geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = errorbar_width, alpha = 0.8)
    
    p <- p + scale_x_continuous(breaks=scales::pretty_breaks(n=10),trans='log10')
    
    # print(summary(lm(dep_means$mean ~ log10(dep_means[[indepvar]]))))
    
  } else {
    
   # errorbar_width <- (max(data[[indepvar]], na.rm=TRUE) - min(data[[indepvar]], na.rm=TRUE)) / 50
    #p <- p + geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = errorbar_width, alpha = 0.8)
    
    p <- p + scale_x_continuous(breaks=scales::pretty_breaks(n=10))
    
    # print(summary(lm(dep_means$mean ~ dep_means[[indepvar]])))
    
  }
  
  print(p)
  
  
}


regression_agg_ <- function(data, indepvar, logx = FALSE) {
  
  # recombobulates data to aggregate by ID
  # then calls regression_output_agg to output table
  
  a <- data.frame()
  my.list <- vector("list", length(names(protein_D14[,2:ncol(protein_D14)])))
  
  for(i in names(protein_D14[,2:ncol(protein_D14)])) {
    
    depvar = i 
    
    dep_means <- data %>%
      group_by(ID) %>%
      summarise_(mean = interp(~SE(var), var = as.name(depvar)))
    
    dep_means$bin_arch_name <- depvar
    
    my.list[[i]] <- dep_means
    
  }
  
  a <- rbind(a, do.call(rbind, my.list))
  
  a <- tidyr::spread(a, bin_arch_name, mean)
  
  climate_locs$ID <- NULL
  
  clim  <- merge(climate_locs, replicates, by = c('sample', 'Latitude', 'Longitude'))
  
  # browser()
  
  #bla <- regression_output_agg(merge(a, clim, by = 'ID'), indepvar, logx=TRUE)
  bla <- regression_output_agg(merge(a, clim, by = 'ID'), indepvar, logx=logx)
  
  
  return(bla)
  
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

