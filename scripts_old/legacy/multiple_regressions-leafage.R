source('scripts/transformations.R')
source('scripts/prep_data.R')

require(vegan)
require(MuMIn)
require(ggplot2)

data$leaf_age <- factor(data$leaf_age, levels = c('new','mid','old'))

# Photosystems vs leaf age, multiple regression models with gap fraction as covariate

  bla1 <- lm(Photosystems ~ gap, data)
  
  bla2 <- lm(Photosystems ~ leaf_age + gap , data) # these results are quoted in euc MS
  summary(bla2)
  anova(bla2) 
  
  bla3 <- aov(Photosystems ~ gap, data_means)
  
  bla4 <- lm(Photosystems ~ leaf_age * gap, data)
  
  bla5 <- lm(Photosystems ~ leaf_age, data = )
  
  AICc(bla1,bla2,bla3,bla4,bla5)
  
  boxplot(data$Photosystems ~ data$leaf_age)
  
  p <- ggplot(data, aes(y = Photosystems, x = leaf_age)) + geom_boxplot()
  p <- p + expand_limits(y=0)
  p <- p + theme_bw()
  p <- p + xlab('leaf age') + ylab('Photosystems protein abundance (proportion)')
  p <- p + theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), 
                 legend.title=element_blank(),
                 legend.position="bottom",
                 axis.line = element_line(colour = "black"))
  p
  
  jou <- varpart(data$Photosystems, 
                 ~leaf_age,
            #     ~log10(leaf_rad),
            #     ~log10(prec),
                 ~log10(gap),
                 data = data)
  jou                         # this is quoted in results
  plot(jou)


# Calvin cycle  
  

data$leaf_age <- factor(data$leaf_age, levels = c('new','mid','old'))

bla1 <- lm(Calvin_cycle ~ gap, data)

#bla2 <- lm(Calvin_cycle ~ gap + leaf_age, data[!data$leaf_age %in% 'new',])
bla2 <- lm(Calvin_cycle ~ gap + leaf_age, data)

summary(bla2)
anova(bla2) 

bla3 <- lm(Calvin_cycle ~ gap, data)

bla4 <- lm(Calvin_cycle ~ gap * leaf_age, data)

bla5 <- lm(Calvin_cycle ~ leaf_age, data)

AICc(bla1,bla2,bla3,bla4, bla5)

boxplot(data$Calvin_cycle ~ data$leaf_age)

p <- ggplot(data, aes(y = Calvin_cycle, x = leaf_age)) + geom_boxplot()
p <- p + expand_limits(y=0)
p <- p + xlab('leaf age') + ylab('Calvin cycle protein abundance (proportion)')
p <- p + theme_bw()
p <- p + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), 
               legend.title=element_blank(),
               legend.position="bottom",
               axis.line = element_line(colour = "black"))
p

jou <- varpart(data$Calvin_cycle, 
               ~leaf_age,
               ~gap,
           #    ~log10(prec),
          #     ~leaf_rad,
               data = data)
jou
plot(jou)

plot(Calvin_cycle ~ pdmt, data)


# percent changes with leaf age

# total_protein

bla1 <- lm(total_protein ~ gap, data)

#bla2 <- lm(total_protein ~ gap + leaf_age, data[!data$leaf_age %in% 'new',])
bla2 <- lm(total_protein ~ gap + leaf_age, data)



bla3 <- lm(total_protein ~ gap, data)

bla4 <- lm(total_protein ~ gap * leaf_age, data)

bla5 <- lm(total_protein ~ leaf_age, data)

summary(bla5)
anova(bla5) 

AICc(bla1,bla2,bla3,bla4, bla5)

boxplot(data$total_protein ~ data$leaf_age)

p <- ggplot(data, aes(y = total_protein, x = leaf_age)) + geom_boxplot()
p <- p + expand_limits(y=0)
p <- p + xlab('leaf age') + ylab('Total protein abundance (mg / m2)')
p <- p + theme_bw()
p <- p + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), 
               legend.title=element_blank(),
               legend.position="bottom",
               axis.line = element_line(colour = "black"))
p

p <- lm(total_protein ~ leaf_age, data) 
boxplot(total_protein ~ leaf_age, data)

tp_new <- mean(data[data$leaf_age == 'new',]$total_protein, na.rm=TRUE)
tp_mid <- mean(data[data$leaf_age == 'mid',]$total_protein, na.rm=TRUE)
tp_old <- mean(data[data$leaf_age == 'old',]$total_protein, na.rm=TRUE)

tp_new_mid <- (tp_mid - tp_new) / tp_new
tp_mid_old <- (tp_old - tp_mid) / tp_mid
tp_new_old <- (tp_old - tp_new) / tp_new
tp_new_mid
tp_mid_old
tp_new_old

tp_leafage <- aov(total_protein ~ leaf_age, data)
TukeyHSD(tp_leafage)



# Photosystems

p <- lm(Photosystems ~ leaf_age, data) 
boxplot(Photosystems ~ leaf_age, data)

ps_new <- mean(data[data$leaf_age == 'new',]$Photosystems, na.rm=TRUE)
ps_mid <- mean(data[data$leaf_age == 'mid',]$Photosystems, na.rm=TRUE)
ps_old <- mean(data[data$leaf_age == 'old',]$Photosystems, na.rm=TRUE)

ps_new_mid <- (ps_mid - ps_new) / ps_new
ps_mid_old <- (ps_old - ps_mid) / ps_mid
ps_new_old <- (ps_old - ps_new) / ps_new
ps_new_mid
ps_mid_old
ps_new_old


ps_leafage <- lm(Photosystems ~ leaf_age, data)
TukeyHSD(ps_leafage)


# Calvin_cycle

p <- lm(Calvin_cycle ~ leaf_age, data) 
boxplot(Calvin_cycle ~ leaf_age, data)

calv_new <- mean(data[data$leaf_age == 'new',]$Calvin_cycle, na.rm=TRUE)
calv_mid <- mean(data[data$leaf_age == 'mid',]$Calvin_cycle, na.rm=TRUE)
calv_old <- mean(data[data$leaf_age == 'old',]$Calvin_cycle, na.rm=TRUE)

calv_new_mid <- (calv_mid - calv_new) / calv_new
calv_mid_old <- (calv_old - calv_mid) / calv_mid
calv_new_old <- (calv_old - calv_new) / calv_new
calv_new_mid
calv_mid_old
calv_new_old


calv_leafage <- aov(Calvin_cycle ~ leaf_age, data)
TukeyHSD(calv_leafage)



source('scripts/prep_data_mg_per_mm2.R')
plot(Calvin_cycle ~ total_protein, data)
plot(Photosystems ~ total_protein, data)


cor(data$Calvin_cycle,data$total_protein)
cor(data$Photosystems,data$total_protein)









