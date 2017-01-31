replicates <- read_csv('output/replicates.csv')

replicates <- merge(replicates, read_csv('data/lineage.csv'), by = 'species_confirmed')

blah <- merge(replicates, protein_stand_D14_age)

blah <- blah[!blah$lineage %in% 'Angophora',]

blah$lineage <- as.factor(blah$lineage)

#boxplot(total_protein ~ lineage, blah)

bla <- aov(total_protein ~ lineage, blah)
  
anova(bla)
require('multcomp')

tuk <- glht(bla, linfct = mcp(lineage = "Tukey"))
summary(tuk)          # standard display
tuk.cld <- cld(tuk)   # letter-based display
opar <- par(mai=c(1,1,1.5,1))
plot(tuk.cld)
par(opar)



source('scripts/prep_data.R')
require(ggplot2)
require(vegan)
library(readr)
library(knitr)
require(plyr)
require(reshape2)
require(dplyr)
require(lazyeval)

#data[data$lineage %in% 'Angophora',]$lineage <- NA

agg_plot_lineage <- function(data, depvar, indepvar, logx = FALSE, labs) {
  
  dep_means <- data %>%
    group_by(ID) %>%
    summarise_(mean = interp(~mean(var, na.rm=TRUE), var = as.name(depvar)),
               SE = lazyeval::interp(~SE(var), var = as.name(depvar))) %>%
    full_join(data, by = 'ID')  %>%
    distinct(mean, .keep_all=TRUE) 
  
  p <- ggplot(dep_means, aes(y = mean, x = dep_means[[indepvar]])) + geom_point(size = 2, aes(colour = lineage))
  p <- p + geom_smooth(aes(colour = lineage), method = 'lm', se = F)
  
  p <- p + xlab(labs[1]) + ggtitle(depvar) + ylab(paste('Species mean of [', depvar, ' proteins] (proportion)', sep = ""))
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







agg_plot_lineage(data, depvar = 'Photosystems', indepvar = 'leaf_rad', logx=FALSE, labs = 'irradiance')


