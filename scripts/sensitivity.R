# sensitivity analysis 
  
  sequence <- seq(1, 0, by = 0.02)
  
  def <- data.frame()
  my.list <- vector("list", length(sequence))
  
  for (j in 1:51) {
    
    source('scripts/transformations_quick.R')
    
    source('scripts/prep_data.R')
    
    # aggregate data
    
    # alter gap fraction for mid and old leaves by i or 1-((1-i)/2)
    
    data[data$leaf_age == 'mid',]$gap <- data[data$leaf_age == 'mid',]$gap * (1-(1-sequence[j])/2)
    data[data$leaf_age == 'old',]$gap <- data[data$leaf_age == 'old',]$gap * sequence[j]
    
    blah.lm <- lm(Photosystems ~ gap, data)
    
    model.stats <- cbind(sequence[j], blah.lm$coefficients[2], round(summary(blah.lm)$r.squared,3))
    
    names(model.stats) <- c('pc_light_reduction_at_oldest_leaf','slope','R2') 
    rownames(model.stats) <- NULL
    
    my.list[[j]] <- model.stats
    
  }
  
  def <- rbind(def, do.call(rbind, my.list))
  names(def) <- c('pc_light_reduction_at_oldest_leaf','slope','R2') 
  #return(def)

  plot(1-def$R2 ~ def$pc_light_reduction_at_oldest_leaf, xlab = 'fraction of light reaching oldest leaf', ylab = '1-R2')
  plot(def$slope ~ def$pc_light_reduction_at_oldest_leaf, xlab = 'fraction of light reaching oldest leaf', ylab = 'slope')



