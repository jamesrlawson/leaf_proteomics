### eucter and euccam ASP vs D14 comparisons (density plots)

source('../scripts/binProteins.R')

protein_samples_asp <- read_csv('../data/asp_proteins.csv')

## absolute [protein]

sample <- c('c','d','e','f','a','b','g','j','k','l','h','i')
euc_samples_asp <- data.frame(sample)
euc_samples_asp$species <- NA
euc_samples_asp[1:6,]$species <- 'euccam'
euc_samples_asp[7:12,]$species <- 'eucter'

eucsasp <- protein_asp[protein_asp$sample %in% euc_samples_asp$sample,]
eucsasp <- merge(eucsasp, euc_samples_asp)
eucsasp$source <- 'asp'

euc_samples_D14 <- read_csv('data/euccam_eucter_samples.csv')
eucsD14 <- protein_D14[protein_D14$sample %in% euc_samples_D14$sample,]
eucsD14 <- merge(eucsD14, euc_samples_D14)
eucsD14$source <- 'D14'

ASPvsD14 <- rbind(eucsD14,eucsasp)
ASPvsD14$sample <- NULL

ASPvsD14_ <- melt(ASPvsD14)
ASPvsD14_big4 <- ASPvsD14_[ASPvsD14_$variable %in% c('Rubisco','PSI','PSII','ATP_synthase_chloroplastic'),]

# density plot

x <- ggplot(ASPvsD14_big4, aes(x = value, fill = variable, color = variable, ..count..))
x <- x + geom_density(alpha = 0.2)
x <- x + facet_grid(species ~ source, scale = c('free'))
x <- x + xlab('[protein] (mg / m2)')
x <- x + theme_bw()
x <- x + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), 
               legend.title=element_blank(),
               legend.position="bottom",
               axis.line = element_line(colour = "black"))
x

## relative [protein]

sample <- c('c','d','e','f','a','b','g','j','k','l','h','i')
euc_samples_asp <- data.frame(sample)
euc_samples_asp$species <- NA
euc_samples_asp[1:6,]$species <- 'euccam'
euc_samples_asp[7:12,]$species <- 'eucter'

eucsasp <- protein_stand_asp[protein_stand_asp$sample %in% euc_samples_asp$sample,]
eucsasp <- merge(eucsasp, euc_samples_asp)
eucsasp$source <- 'asp'

euc_samples_D14 <- read_csv('data/euccam_eucter_samples.csv')
eucsD14 <- protein_stand_D14[protein_stand_D14$sample %in% euc_samples_D14$sample,]
eucsD14 <- merge(eucsD14, euc_samples_D14)
eucsD14$source <- 'D14'

ASPvsD14 <- rbind(eucsD14,eucsasp)
ASPvsD14$sample <- NULL

ASPvsD14_ <- melt(ASPvsD14)
ASPvsD14_big4 <- ASPvsD14_[ASPvsD14_$variable %in% c('Rubisco','PSI','PSII','ATP_synthase_chloroplastic'),]

# density plots

x <- ggplot(ASPvsD14_big4, aes(x = value, fill = variable, color = variable))
x <- x + geom_density(alpha = 0.2)
x <- x + facet_grid(species ~ source, scale = c('free'))
x <- x + xlab('[protein] (proportion of total)')
x <- x + theme_bw()
x <- x + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), 
               legend.title=element_blank(),
               legend.position="bottom",
               axis.line = element_line(colour = "black"))
x

summary(aov(PSII ~ species * source, ASPvsD14))

# summary table

blah <- ddply(ASPvsD14_big4, .(species,source,variable), summarise, mean = mean(value), CV = CV(value))
