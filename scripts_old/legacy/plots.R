# produces density distributions for important protein functional categories
# compares field vs. aspinwall data

source('../scripts/transformations.R') # ?

require(gridExtra)

# [protein] unstandardised by [total protein]

protein_climate_spread <-protein_climate[,c('sample', 'bin_arch_name','sum')]
protein_climate_spread <- na.omit(protein_climate_spread)
protein_climate_spread <- spread(protein_climate_spread, key = bin_arch_name, value=sum)
#protein_climate_spread[,2:23] <- protein_climate_spread[,2:23]/protein_climate_spread$Rubisco

protein_climate_melt <- melt(protein_climate_spread)
names(protein_climate_melt) <- c('sample', 'bin_arch_name', 'sum')
protein_climate_melt <- merge(protein_climate_melt, climate_locs, by = 'sample')
protein_climate_melt$bin_arch_name <- as.character(protein_climate_melt$bin_arch_name)

tavg_lms <- regression_output(na.omit(protein_climate_melt), tavg, sum)
tavg_lms <- tavg_lms[order(tavg_lms$R2, decreasing =TRUE),]
write.csv(tavg_lms, 'output/MAT_models.csv')

prec_lms <- regression_output(na.omit(protein_climate_melt), prec, sum)
prec_lms <- prec_lms[order(prec_lms$R2, decreasing =TRUE),]
write.csv(prec_lms, 'output/MAP_models.csv')

protein_climate_spread <- spread(protein_climate_stand, bin_arch_name, sum)

plot(log(Rubisco) ~ prec, data = protein_climate_spread)
abline(lm(log(Rubisco)~ prec, protein_climate_spread))
summary(lm(log(Rubisco) ~ prec, protein_climate_spread))

plot(lipid_metabolism ~ prec, data = protein_climate_spread)
abline(lm(lipid_metabolism ~ prec, protein_climate_spread))
summary(lm(lipid_metabolism ~ prec, protein_climate_spread))

######### histograms #########
 # standardised by Rubisco #

require(ggplot2)

#protein_climate_stand <- melt(protein_climate_spread)
#names(protein_climate_stand) <- c('sample','bin_arch_name','sum')



protein_climate_RbcStand <-protein_climate[,c('sample', 'bin_arch_name','sum')]
protein_climate_RbcStand <- na.omit(protein_climate_RbcStand)
protein_climate_RbcStand <- spread(protein_climate_RbcStand, key = bin_arch_name, value=sum)
protein_climate_RbcStand[,2:26] <- protein_climate_RbcStand[,2:26]/protein_climate_RbcStand$Rubisco


big4 <- protein_climate_RbcStand[,c('PSI','PSII','ATP_synthase_chloroplastic', 'electron_transport_minATPsynth')]

big4 <- melt(big4)
names(big4) <- c('bin_arch_name', 'sum')

x <- ggplot(big4, aes(x = sum, fill = bin_arch_name, color = bin_arch_name))
x <- x + geom_density(alpha = 0.2, stat = 'density')
#x <- x + geom_histogram(alpha = 0.2, binwidth = 0.05)

x <- x + xlab('[protein] (mg / m2)')
x <- x + theme_bw()
x <- x + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), 
               legend.title=element_blank(),
               legend.position="bottom",
               axis.line = element_line(colour = "black"))
plot(x)


blah <- ddply(melt(protein_climate_RbcStand[,c('Rubisco','Calvin_cycle','Cytochrome_b6f','PSI','PSII','ATP_synthase_chloroplastic', 'electron_transport_minATPsynth')]), 
              .(variable), 
              summarise,
              mean = mean(value),
              CV = CV(value))

write_csv(blah, '../output/D14_proteins_standRbc.csv')


### aspinwall vs field [protein]
aspinwall <- read_csv('data/aspinwall_preheatwave.csv')
aspinwall_long <- melt(aspinwall, id = c('id'))
names(aspinwall_long) <- c('sample', 'bin_arch_name', 'sum')
aspinwall_long$bin_arch_name <- as.character(aspinwall_long$bin_arch_name)
big4_asp <- aspinwall_long[aspinwall_long$bin_arch_name %in% c('Rubisco','PSI','PSII','ATP_synthase_chloroplastic'),]
big4_asp$study <- 'aspinwall'

big4$study <- 'D14'
big4_d14_asp <- rbind(big4[,c('sample','bin_arch_name', 'sum', 'study')], big4_asp)
big4_d14_asp$sum <- as.numeric(big4_d14_asp$sum)

protein_avgs <- ddply(big4_d14_asp, .(study,bin_arch_name), summarise, mean = mean(sum), CV = CV(sum))
names(protein_avgs)[2] <- 'protein_fn_group'
write_csv(protein_avgs, 'output/protein_avgs.csv')

x <- ggplot(big4_d14_asp, aes(x = sum, fill = bin_arch_name, color = bin_arch_name, ..count..))
x <- x + geom_density(alpha = 0.2)
x <- x + facet_grid(.~ study, scale = c('free'))
x <- x + xlab('[protein] (mg / m2)')
x <- x + theme_bw()
x <- x + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), 
               legend.title=element_blank(),
               legend.position="bottom",
               axis.line = element_line(colour = "black"))

tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)), base_size = 8)
tbl <- tableGrob(protein_avgs, rows=NULL, theme=tt)
table_size <- nrow(protein_avgs) * 0.3
# Plot chart and table into one object
grid.arrange(x, tbl,
             nrow=2,
             as.table=TRUE,
             bottom=TRUE,
             heights=c(4,table_size))




# total protein amount per sample

total_protein <- ddply(na.omit(protein_climate), .(sample), summarise, total_protein = sum(sum))
total_protein_clim <- merge(climate_locs, total_protein, by = 'sample')

plot(total_protein ~ tavg, data = total_protein_clim)
summary(lm(total_protein ~ tavg, data = total_protein_clim))
abline(lm(total_protein ~ tavg, data = total_protein_clim))
