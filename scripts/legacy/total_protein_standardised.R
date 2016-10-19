# not really sure what this stuff is or if it's redundant

# TOTAL PROTEIN #

protein_samples <- read_csv('../data/D14_protein_sites.csv')
protein_samples_long <- melt(protein_samples)

total_protein <- ddply(protein_samples_long, .(variable), summarise, total_protein = sum(value, na.rm=TRUE))
names(total_protein) <- c('sample', 'total_protein')                       

protein_bins_spread <- na.omit(protein_bins[,c('sample','bin_arch_name','sum')])
protein_bins_spread <- spread(protein_bins_spread, key = bin_arch_name, value=sum)
protein_bins_spread <- protein_bins_spread[-1,]

protein_stand <- merge(protein_bins_spread, total_protein, by = 'sample')
protein_stand[,2:25] <- protein_stand[,2:25]/protein_stand$total_protein

protein_stand <- melt(protein_stand)
names(protein_stand) <- c('sample','bin_arch_name', 'sum')
protein_stand <- merge(protein_stand, climate_locs, by = 'sample')
protein_stand$bin_arch_name <- as.character(protein_stand$bin_arch_name)

tavg_lms <- regression_output(na.omit(protein_stand), tavg, sum)
tavg_lms <- tavg_lms[order(tavg_lms$R2, decreasing =TRUE),]
write.csv(tavg_lms, '../output/MAT_models.csv')

prec_lms <- regression_output(na.omit(protein_stand), prec, sum)
prec_lms <- prec_lms[order(prec_lms$R2, decreasing =TRUE),]
write.csv(prec_lms, '../output/MAT_models.csv')


protein_stand_spread <- spread(protein_stand, bin_arch_name, sum)

plot(PSII ~ tavg, data = protein_stand_spread)
abline(lm(PSII~ tavg, protein_stand_spread))
summary(lm(PSII ~ tavg + I(tavg^2) , protein_stand_spread))

plot(total_protein ~ prec, data = protein_stand_spread)
abline(lm(total_protein~ prec, protein_stand_spread))
summary(lm(total_protein ~ prec, protein_stand_spread))


protein_stand_spread$arid <- protein_stand_spread$prec/protein_stand_spread$tavg

b <- ggplot(protein_stand_spread, aes(x = prec, y = PSII))
b <- b + geom_point(alpha = 0.5, col = 'red')
b <- b + geom_smooth(method = lm,formula = y ~ x + I(x^2))
b

protein_stand_spread$PSII_max <- protein_stand_spread$PSII/max(protein_stand_spread$PSII)

#  svg(sprintf("output/figures/linear/%s.svg", proteinname), width = 2.795, height = 2.329, pointsize=8)

pdf('../output/figures/PSII.pdf', width = 4, height = 2.88, pointsize=10)
b <- ggplot(protein_stand_spread, aes(x = prec, y = PSII*100, col = tavg))
#b <- b + geom_point(alpha = 0.5)
#b <- b + scale_colour_gradient(low = 'green', high='red', name = "Mean annual temperature (deg C)")
b <- b + geom_smooth(method = lm,formula = y ~ x + I(x^2), size = 0.3, alpha = 0.3)
b <- b + xlab('Average annual precipitation (mm)') 
b <- b + ylab('% Photosystem II')
b <- b + theme_bw()
b <- b + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), 
               legend.position="top",
               axis.line = element_line(colour = "black"))
b
dev.off()













