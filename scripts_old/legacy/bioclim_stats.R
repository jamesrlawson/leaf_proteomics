source('scripts/prep_data.R')

data_means <- data %>% group_by(ID) %>% summarise(mean = mean(Photosystems, na.rm=TRUE),
                                                  SE = SE(Photosystems))
data_means <- merge(data, data_means)
data_means <- distinct(data_means, ID, .keep_all = TRUE)

my.list <- vector("list", length(48:66))
a <- data.frame()

for(i in 48:66) {
  
  protname <- names(data_means[i])
 
  model <- summary(lm(data_means$mean ~ data_means[[i]]))
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

outs <- a
