# sensitivity analysis 
  
  sequence <- seq(0, 1, by = 0.01)
  
  def <- data.frame()
  my.list <- vector("list", length(sequence))
  
  for (j in 1:length(sequence)) {
    
    blz <- sequence[j]
    
    source('scripts/transformations_quick.R')
    
    source('scripts/prep_data.R')
    
    # aggregate data
    
    # alter gap fraction for mid and old leaves by i or 1-((1-i)/2)
    
    data[data$leaf_age == 'mid',]$gap <- data[data$leaf_age == 'mid',]$gap * (1-(1-sequence[j])/2)
    data[data$leaf_age == 'old',]$gap <- data[data$leaf_age == 'old',]$gap * sequence[j]
    
    blah.lm <- lm(Photosystems ~ gap, data)
    
    model.stats <- cbind(sequence[j], blah.lm$coefficients[2], summary(blah.lm)$r.squared)
    
    names(model.stats) <- c('fraction_light_reaching_oldest_leaf','slope','R2') 
    rownames(model.stats) <- NULL
    
    my.list[[j]] <- model.stats
    
   # browser()
    
  }
  
  def <- rbind(def, do.call(rbind, my.list))
  names(def) <- c('fraction_light_reaching_oldest_leaf','slope','R2') 
  #return(def)

  plot(1-def$R2 ~ def$fraction_light_reaching_oldest_leaf, xlab = 'fraction of light reaching oldest leaf', ylab = '1-R2')
  plot(def$slope ~ def$fraction_light_reaching_oldest_leaf, xlab = 'fraction of light reaching oldest leaf', ylab = 'slope')

  
  staz <- blah  %>% summarise(mean = mean(1-def$fraction_light_reaching_oldest_leaf, na.rm=TRUE), 
                                                   SE = SE(1-def$fraction_light_reaching_oldest_leaf),
                                                   n = length(def$fraction_light_reaching_oldest_leaf)) %>%  
    mutate(E = qt(.975, df=n−1)∗SE)

  
  
# compare slopes of Photosystems ~ Reich adjusted and raw gap fraction
  
  source('scripts/transformations_quick.R')
  
  source('scripts/prep_data.R')
  
  data_raw <- data
  data_raw$model <- 'raw'
  
  raw.lm <- lm(Photosystems ~ gap, data_raw)
  
  data_reich <- data
  data_reich$model <- 'reich'
  
  data_reich[data_reich$leaf_age == 'mid',]$gap <- data_reich[data_reich$leaf_age == 'mid',]$gap * (1-0.115)
  data_reich[data_reich$leaf_age == 'old',]$gap <- data_reich[data_reich$leaf_age == 'old',]$gap * (1-0.23)
  
  reich.lm <- lm(Photosystems ~ gap, data_reich)
  
  data_reich_raw <- rbind(data_raw, data_reich)
  
  reich_raw.lm <- lm(Photosystems ~ gap*model, data_reich_raw)
  anova(reich_raw.lm)
  

  
# ggplots
  
  
  R2 <- ggplot(def, aes(x = fraction_light_reaching_oldest_leaf, y = 1-R2))
  R2 <- R2 + geom_smooth(method = 'loess', se=FALSE)
  R2 <- R2 + geom_vline(xintercept = 1-0.225, colour = 'red') 
  R2 <- R2 + geom_vline(xintercept = def[1-def$R2 == min(1-def$R2),]$fraction_light_reaching_oldest_leaf + staz$E)
  R2 <- R2 + geom_vline(xintercept = def[1-def$R2 == min(1-def$R2),]$fraction_light_reaching_oldest_leaf - staz$E)
  R2 <- R2 + theme_bw()
  R2
  
  slope <- ggplot(def, aes(x = fraction_light_reaching_oldest_leaf, y = slope))
  slope <- slope + geom_smooth(method = 'loess', se=FALSE)
  slope <- slope + geom_vline(xintercept = 1-0.225, colour = 'red') 
  slope <- slope + geom_hline(yintercept = confint(zero.lm)[2])
  slope <- slope + geom_hline(yintercept = confint(zero.lm)[4])

  slope <- slope + geom_vline(xintercept = def[def$slope == min(def$slope),]$fraction_light_reaching_oldest_leaf + staz$E) 
  slope <- slope + geom_vline(xintercept = def[def$slope == min(def$slope),]$fraction_light_reaching_oldest_leaf - staz$E)
  slope <- slope + theme_bw()
  slope 
  
  
