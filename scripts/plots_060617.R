require(readr)

leafN_7_16 <- read_csv('data/leaf_data/CNP/raw/CN_Leaf_Weights all trays (7-16)_combined.csv')

leafN <- read_csv('data/leaf_data/CNP/leaf_CN.csv')

LMA <- read_csv('data/leaf_data/misc/LMA_LWC.csv')

LMA_leafN <- merge(LMA, leafN, by = 'sample')

LMA_leafN$N_area <- LMA_leafN$N * 10 * LMA_leafN$LMA_g_per_m2

hist(LMA_leafN$N_area)

plot(LMA_leafN$N_area ~ LMA_leafN$LMA_g_per_m2)

include_leaf_N=TRUE

source('scripts/transformations.R')
source('scripts/prep_data.R')

datax <- data

ggplot(data, aes(y = total_protein_mean, x = gap_mean)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = total_protein_mean - total_protein_SE, ymax = total_protein_mean + total_protein_SE)) +
  geom_smooth(method = 'lm')
 


ggplot(data, aes(y = total_protein_mean, x = Narea_mean)) + 
  geom_point(size = 2) + 
  geom_errorbar(aes(ymin = total_protein_mean - total_protein_SE, ymax = total_protein_mean + total_protein_SE)) +
  geom_errorbarh(aes(xmin = Narea_mean - Narea_SE, xmax = Narea_mean + Narea_SE)) +
  geom_smooth(method = 'lm',se=FALSE) +
  ylab('Leaf protein (mg/m2)') + xlab('Leaf nitrogen (mg/m2)') +
  theme_classic() +
  theme(legend.title=element_blank(),
               legend.position="bottom",
               axis.line = element_line(colour = "black"),
               text = element_text(size = 18))

ggplot(data, aes(y = total_protein, x = N_per_area)) + 
  geom_point(size = 2) + 
  geom_smooth(method = 'lm',se=FALSE) +
  ylab('Leaf protein (mg/m2)') + xlab('Leaf nitrogen (mg/m2)') +
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.position="bottom",
        axis.line = element_line(colour = "black"),
        text = element_text(size = 18))


bla <- filter(datax, N_per_area > 4000 & total_protein < 30000) %>%
       select(sample, ID, total_protein, N_per_area, site_revised)


ggplot(data, aes(y = total_protein, x = tavg)) + 
  geom_point(size = 2) + 
#  geom_errorbar(aes(ymin = total_protein_mean - total_protein_SE, ymax = total_protein_mean + total_protein_SE)) +
  geom_smooth(method = 'lm',se=FALSE) +
  ylab('Leaf protein (mg/m2)') + xlab('Mean annual temperature (oC)') +
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.position="bottom",
        axis.line = element_line(colour = "black"),
        text = element_text(size = 18))

ggplot(data, aes(y = total_protein, x = prec)) + 
  geom_point(size = 2) + 
  #  geom_errorbar(aes(ymin = total_protein_mean - total_protein_SE, ymax = total_protein_mean + total_protein_SE)) +
  geom_smooth(method = 'lm',se=FALSE) +
  ylab('Leaf protein (mg/m2)') + xlab('Mean annual preciptation (mm)') +
  scale_x_continuous(breaks=scales::pretty_breaks(n=10),trans='log10') +
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.position="bottom",
        axis.line = element_line(colour = "black"),
        text = element_text(size = 18))


source('scripts/prep_data_mg_per_mm2.R')

bla <- select(data, total_protein, Rubisco, Photosystems, calvin_cycle)

bla$calv <- bla$calvin_cycle - bla$Rubisco


blax <- gather(bla, key = 'funcat', value = 'amount', - total_protein)

ggplot(blax, aes(y = amount, x = total_protein)) + geom_point(aes(colour = funcat, shape = funcat))
