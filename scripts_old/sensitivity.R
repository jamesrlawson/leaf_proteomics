# sensitivity analysis 

sequence <- seq(0, 1, by = 0.005)

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

plot(def$R2 ~ def$fraction_light_reaching_oldest_leaf, xlab = 'fraction of light reaching oldest leaf', ylab = 'R2')
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

data_reich[data_reich$leaf_age == 'mid',]$gap <- data_reich[data_reich$leaf_age == 'mid',]$gap * (1-0.1125)
data_reich[data_reich$leaf_age == 'old',]$gap <- data_reich[data_reich$leaf_age == 'old',]$gap * (1-0.225)

reich.lm <- lm(Photosystems ~ gap, data_reich)

data_bestfit <- data
data_bestfit$model <- 'bestfit'

data_bestfit[data_bestfit$leaf_age == 'mid',]$gap <- data_bestfit[data_bestfit$leaf_age == 'mid',]$gap * (1-0.165)
data_bestfit[data_bestfit$leaf_age == 'old',]$gap <- data_bestfit[data_bestfit$leaf_age == 'old',]$gap * (1-0.33)

data_reich_best_raw <- rbind(data_raw, data_reich)
data_reich_best_raw <- rbind(data_reich_best_raw, data_bestfit)


reich_best_raw.lm <- lm(Photosystems ~ gap*model, data_reich_best_raw)
anova(reich_best_raw.lm)



# ggplots


R2 <- ggplot(def, aes(x = fraction_light_reaching_oldest_leaf, y = R2))
R2 <- R2 + geom_point()
R2 <- R2 + geom_vline(xintercept = 1-0.225, colour = 'red') 
R2 <- R2 + geom_vline(xintercept = def[def$R2 == max(def$R2),]$fraction_light_reaching_oldest_leaf + staz$E)
R2 <- R2 + geom_vline(xintercept = def[def$R2 == max(def$R2),]$fraction_light_reaching_oldest_leaf - staz$E)
R2 <- R2 + geom_vline(xintercept = def[def$R2 == max(def$R2),]$fraction_light_reaching_oldest_leaf, color = 'orange')

R2 <- R2 + theme_bw()
R2

slope <- ggplot(def, aes(x = fraction_light_reaching_oldest_leaf, y = slope))
slope <- slope + geom_smooth(method = 'loess', se=FALSE)
slope <- slope + geom_vline(xintercept = 1-0.225, colour = 'red') 
slope <- slope + geom_hline(yintercept = confint(zero.lm)[2], color = 'darkgreen')
slope <- slope + geom_hline(yintercept = confint(zero.lm)[4], color = 'darkgreen')

slope <- slope + geom_vline(xintercept = def[def$slope == min(def$slope),]$fraction_light_reaching_oldest_leaf + staz$E) 
slope <- slope + geom_vline(xintercept = def[def$slope == min(def$slope),]$fraction_light_reaching_oldest_leaf - staz$E)
slope <- slope + geom_vline(xintercept = def[def$R2 == max(def$R2),]$fraction_light_reaching_oldest_leaf, color = 'orange')

slope <- slope + theme_bw()
slope 












# compare slopes of Photosystems ~ Reich adjusted and raw gap fraction

source('scripts/transformations_quick.R')

source('scripts/prep_data.R')


data_nolight <- data
data_nolight[data_nolight$leaf_age == 'mid',]$gap <- data_nolight[data_nolight$leaf_age == 'mid',]$gap * (0)
data_nolight[data_nolight$leaf_age == 'old',]$gap <- data_nolight[data_nolight$leaf_age == 'old',]$gap * (0.5)

data_75 <- data
data_75[data_75$leaf_age == 'mid',]$gap <- data_75[data_75$leaf_age == 'mid',]$gap * (0.25)
data_75[data_75$leaf_age == 'old',]$gap <- data_75[data_75$leaf_age == 'old',]$gap * (1-(1-0.25)/2)

data_50 <- data
data_50[data_50$leaf_age == 'mid',]$gap <- data_50[data_50$leaf_age == 'mid',]$gap * (0.5)
data_50[data_50$leaf_age == 'old',]$gap <- data_50[data_50$leaf_age == 'old',]$gap * (1-(1-0.5)/2)

data_25 <- data
data_25[data_25$leaf_age == 'mid',]$gap <- data_25[data_25$leaf_age == 'mid',]$gap * (0.75)
data_25[data_25$leaf_age == 'old',]$gap <- data_25[data_25$leaf_age == 'old',]$gap * (1-(1-0.75)/2)


plot(Photosystems ~ gap, data_nolight)
plot(Photosystems ~ gap, data_75)
plot(Photosystems ~ gap, data_50)
plot(Photosystems ~ gap, data_25)
plot(Photosystems ~ gap, data)

