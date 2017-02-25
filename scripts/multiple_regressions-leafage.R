source('scripts/prep_data.R')


leafrad_mean <- data %>% group_by(ID, leaf_age) %>% summarise(leafrad_mean = mean(leaf_rad, na.rm=TRUE), leafrad_SE = SE(leaf_rad))
leafrad_mean <- merge(data, leafrad_mean, by = c('leaf_age', 'ID'))


data$leafrad_mean <- leafrad_mean$leafrad_mean
data$leafrad_SE <- leafrad_mean$leafrad_SE


data_means <- data %>% group_by(ID, leaf_age) %>% summarise(mean = mean(Photosystems, na.rm=TRUE),
                                                  SE = SE(Photosystems))
data_means <- merge(data, data_means)
data_means <- distinct(data_means, ID, leaf_age, .keep_all = TRUE)

data_means$leaf_age <- factor(data_means$leaf_age, levels = c('new','mid','old'))

bla1 <- lm(mean ~ as.factor(leaf_age), data_means)
anova(bla1)

bla2 <- lm(mean ~ leafrad_mean + as.factor(leaf_age), data_means)
summary(bla2)
anova(bla2) 

bla3 <- lm(mean ~ leafrad_mean, data_means)

AICc(bla1,bla2,bla3)

boxplot(data_means$mean ~ as.factor(data_means$leaf_age))

jou <- varpart(data_means$mean, 
               ~leaf_age,
               ~leafrad_mean,
            #   ~gap_mean,
               data = data_means)
jou
plot(jou)











data$leaf_age <- factor(data$leaf_age, levels = c('new','mid','old'))

bla1 <- lm(Photosystems ~ gap, data)
anova(bla1)

bla2 <- lm(Photosystems ~ gap + leaf_age, data)
summary(bla2)
anova(bla2) 

bla3 <- lm(Photosystems ~ gap, data_means)

bla4 <- lm(Photosystems ~ gap * leaf_age, data)


AICc(bla1,bla2,bla3,bla4)

boxplot(data$Photosystems ~ data$leaf_age)

jou <- varpart(data$Photosystems, 
               ~leaf_age,
          #     ~log10(leaf_rad),
          #     ~log10(prec),
               ~log10(gap),
               data = data)
jou
plot(jou)



plot(Photosystems ~ gap, data)
p <- lm(Photosystems ~ gap, data)
abline(p)
summary(p)

y = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * min(data$gap)
x = as.numeric(coef(p)[1]) + as.numeric(coef(p)[2]) * max(data$gap)

(x - y)/x




data$leaf_age <- factor(data$leaf_age, levels = c('new','mid','old'))

bla1 <- lm(Calvin_cycle ~ gap, data)
anova(bla1)
summary(bla1)

bla2 <- lm(Calvin_cycle ~ gap + leaf_age, data)
summary(bla2)
anova(bla2) 

bla3 <- lm(Calvin_cycle ~ gap, data_means)

bla4 <- lm(Calvin_cycle ~ gap * leaf_age, data)


AICc(bla1,bla2,bla3,bla4)

boxplot(data$Calvin_cycle ~ data$leaf_age)

jou <- varpart(data$Calvin_cycle, 
               ~leaf_age,
               ~gap,
               ~log10(prec),
               ~log10(leaf_rad),
               data = data)
jou
plot(jou)

plot(Calvin_cycle ~ pdmt, data)


# percent changes with leaf age

# total_protein

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

tp_leafage <- lm(total_protein ~ leaf_age, data)
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


calv_leafage <- lm(Calvin_cycle ~ leaf_age, data)
TukeyHSD(calv_leafage)