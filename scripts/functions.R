# sourced by binProteins.R

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
  
  protein_bins <- ddply(protein_bins, .(bin_arch, sample), summarise, sum = sum(as.numeric(value), na.rm=TRUE), nprot=length(value))
  protein_bins <- na.omit(protein_bins)
  
  # recombobulate arch categories
  protein_bins$sample <- as.character(protein_bins$sample)
  
  protein_bins[protein_bins$bin_arch == '_1.1_',]$sum <-     protein_bins[protein_bins$bin_arch == '_1.1_',]$sum +
    protein_bins[protein_bins$bin_arch == '_1.1.1_',]$sum +
    protein_bins[protein_bins$bin_arch == '_1.1.2_',]$sum +
    protein_bins[protein_bins$bin_arch == '_1.1.3_',]$sum +
    protein_bins[protein_bins$bin_arch == '_1.1.4_',]$sum
  
  protein_bins[protein_bins$bin_arch == '_1.3_',]$sum <-     protein_bins[protein_bins$bin_arch == '_1.3_',]$sum +
    protein_bins[protein_bins$bin_arch == '_1.3.1_',]$sum +
    protein_bins[protein_bins$bin_arch == '_1.3.2_',]$sum
  
  protein_bins[protein_bins$bin_arch == '_29_',]$sum <- protein_bins[protein_bins$bin_arch == '_29_',]$sum +
    protein_bins[protein_bins$bin_arch == '_29.6_',]$sum
  
  protein_bins[protein_bins$bin_arch == '_20_',]$sum <- protein_bins[protein_bins$bin_arch == '_20_',]$sum +
    protein_bins[protein_bins$bin_arch == '_20.1_',]$sum +
    protein_bins[protein_bins$bin_arch == '_20.2_',]$sum +
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

regression_output <- function(df, x, y) { 
  
  # x is the dependent variable (whatever climate variable), y is 'sum'
  x <- deparse(substitute(x))
  y <- deparse(substitute(y))
  
  a <- data.frame()
  my.list <- vector("list", length(unique(df$bin_arch_name)))
  
  for(i in 1:length(unique(df$bin_arch_name))) {
    sub <- unique(df$bin_arch_name)[i]
    df_sub <- subset(df, bin_arch_name == sub)
    
    model <- lm(paste(y, '~', x, sep = ""), data = df_sub)
    
  #  model.stats <-  data.frame(cbind(sub,summary(model)$r.squared, summary(model)$coefficients[,4][2]))
    
    model.stats <-  data.frame(cbind(sub,round(summary(model)$r.squared,3),
                                     round(summary(model)$coefficients[,4][2],2)))
    
    
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
  
  browser()
}


