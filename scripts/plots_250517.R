source('scripts/prep_data_mg_per_mm2.R')

bla <- select(data, LMA_g_per_m2, total_protein, Photosystems, calvin_cycle, leaf_rad) %>%
       gather(key = 'funccat', value = 'protein_amount', -LMA_g_per_m2, -total_protein, -leaf_rad)

ggplot(bla, aes(y = protein_amount, x = LMA_g_per_m2)) + geom_point(aes(colour = funccat, shape = funccat)) + geom_smooth(aes(colour = funccat), method = lm, se=FALSE) +
scale_colour_manual(values = c('red','green'))

ggplot(bla, aes(y = protein_amount, x = total_protein)) + geom_point(aes(shape = funccat, colour=funccat)) + geom_smooth(aes(colour = funccat), method = lm, se=FALSE) +
   scale_colour_manual(values = c('red','green'))


summary(lm(Photosystems ~ LMA_g_per_m2, data = data))
summary(lm(calvin_cycle ~ LMA_g_per_m2, data = data))


ggplot(data, aes(y = tavg, x = prec)) + geom_point(size = 3, shape = 17) + scale_x_continuous(breaks=scales::pretty_breaks(n=10),trans='log10') +
  xlab('Mean annual precipiation (m/yr)') + ylab('Mean annual temperature (oC)') +
  theme_classic() +
  theme(legend.title=element_blank(),
               legend.position="bottom",
               axis.line = element_line(colour = "black"),
               text = element_text(size = 18))

ggplot(data, aes(y = total_protein_mean, x = prec)) + geom_point(size = 2) + scale_x_continuous(breaks=scales::pretty_breaks(n=10),trans='log10') +
  xlab('Mean annual precipiation (m/yr)') + ylab('Leaf protein (mg/m2)') +
  geom_errorbar(aes(ymin = total_protein_mean - total_protein_SE, ymax = total_protein_mean + total_protein_SE), width = 0, alpha = 0.8) +
  geom_smooth(method = 'lm', se=FALSE, colour = 'black') +
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.position="bottom",
        axis.line = element_line(colour = "black"),
        text = element_text(size = 18))

ggplot(data, aes(y = total_protein_mean, x = tavg)) + geom_point(size = 2) +
  geom_errorbar(aes(ymin = total_protein_mean - total_protein_SE, ymax = total_protein_mean + total_protein_SE), width = 0, alpha = 0.8) +
  geom_smooth(method = 'lm', se=FALSE, colour = 'black') +
  xlab('Mean annual temperature (oC)') + ylab('Leaf protein (mg/m2)') +
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.position="bottom",
        axis.line = element_line(colour = "black"),
        text = element_text(size = 18))





