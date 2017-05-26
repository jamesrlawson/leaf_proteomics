source('scripts/functions.R')

require(readr)
require(stringr)
library(dplyr)
library(tidyr)
require(reshape2)
require(plyr)

protein_samples <- read_csv('data/D14_protein_sites.csv')
mercator <- read_csv('data/asp_mercator.csv')
mercator_bins <- read_csv('data/asp_mercator_bins.csv')
sample_locations <- read_csv('data/sample_locations.csv')
climate <- read_csv('data/discovery_site_climate.csv')


protein_samples <- getProteinBins(protein_samples, mercator)

#protein_samples <- multiple_bins(protein_samples)

# sum values within overarching bins

#bin_arch <- c('_20.2_', '20.2.','_1.1.1_','_1.1.2_','_1.1_', '1.1.','_1.2_', '1.2.','_1.3_', '1.3.','_1.3.1_','_1.3.2_','_1.1.3_','_2_','_2.',
#              '_3_','_3.','_4_','_4.','_8.1_', '8.1.','_8.2_','_8.2.','_9_','_9.',
#                '_11_','_11.','_16_','_16.','_17_','_17.','_21_','_21.','_26.9_','_26.9.','_22_','_22.','_30_','_30.','_29.6_','_29.6.','_20.2.1_',
#              '_20.2.1.','_29_','_29.','_28_','_28.','_27_','_27.','_1.1.4_')

# note that for this to work, all bincodes in asp_mercator.csv and mercator_bins.csv must be prefixed and suffixed by '_'
# the following is a hierarchically ordered list of bincodes to be grepped over, where the first string in each list element defines the top level, 
# and therefore a grep using this element as the pattern will only find proteins in the top level of the hierarchy,
# while the second element allows the grep call to find all proteins below the top category
# it is essential that elements in this list are ordered sequentially in a hierarchical fashion (i.e. higher levels must come before lower levels)
# the order of elements which are at the same level is not relevant

bin_arch.list <- list(
              c('_20.2_', '_20.2.'),
              c('_20.2.1_','_20.2.1.'),
              c('_1.1_', '_1.1.'),
              c('_1.1.1_', '_1.1.1.'),
              c('_1.1.2_', '_1.1.2.'),
              c('_1.1.3_', '_1.1.3.'),
              c('_1.1.4_', '_1.1.4.'),
              c('_1.1.5_', '_1.1.5.'),
              c('_1.1.6_', '_1.1.6.'),
              c('_1.2_', '_1.2.'),
              c('_1.3_', '_1.3.'),
              c('_1.3.1_', '_1.3.1.'),
              c('_1.3.2_', '_1.3.2.'),
              c('_2_','_2.'),
              c('_3_','_3.'),
              c('_4_','_4.'),
              c('_8.1_', '_8.1.'),
              c('_8.2_','_8.2.'),
              c('_9_','_9.'),
              c('_11_','_11.'),
              c('_16_','_16.'),
              c('_17_','_17.'),
              c('_21_','_21.'),
              c('_26.9_','_26.9.'),
              c('_22_','_22.'),
              c('_30_','_30.'),
              c('_29_','_29.'),
              c('_29.6_','_29.6.'),
              c('_28_','_28.'),
              c('_27_','_27.'))

protein_samples$bin_arch <- NA

for(i in 1:length(bin_arch.list)) {
  
  x <-  bin_arch.list[[i]][1]
  rows1 <- nrow(protein_samples[grep(bin_arch.list[[i]][1], protein_samples$BINCODE, fixed=TRUE),])
  rows2 <- nrow(protein_samples[grep(bin_arch.list[[i]][2], protein_samples$BINCODE, fixed=TRUE),])
  
  if(nrow(protein_samples[grep(bin_arch.list[[i]][1], protein_samples$BINCODE, fixed=TRUE),]) > 0) {
   protein_samples[grep(bin_arch.list[[i]][1], protein_samples$BINCODE, fixed=TRUE),]$bin_arch <- bin_arch.list[[i]][1]
  }
  
  if(nrow(protein_samples[grep(bin_arch.list[[i]][2], protein_samples$BINCODE, fixed=TRUE),]) > 0) {
   protein_samples[grep(bin_arch.list[[i]][2], protein_samples$BINCODE, fixed=TRUE),]$bin_arch <- bin_arch.list[[i]][1]
   
  } 
   
}


protein_bins <- melt(protein_samples, id = c('Protein', 'bin_arch'))
names(protein_bins) <- c('Protein', 'bin_arch', 'sample', 'value')

protein_bins <- ddply(protein_bins, .(bin_arch, sample), summarise, sum = sum(as.numeric(value), na.rm=TRUE), nprot=length(value))
protein_bins <- na.omit(protein_bins)

# recombobulate arch categories (where lower levels need to be incorporated into an higher level)
 protein_bins$sample <- as.character(protein_bins$sample)

  protein_bins[protein_bins$bin_arch == '_1.1_',]$sum <-     protein_bins[protein_bins$bin_arch == '_1.1_',]$sum +
  protein_bins[protein_bins$bin_arch == '_1.1.1_',]$sum +
  protein_bins[protein_bins$bin_arch == '_1.1.2_',]$sum +
  protein_bins[protein_bins$bin_arch == '_1.1.3_',]$sum +
  protein_bins[protein_bins$bin_arch == '_1.1.4_',]$sum +
  protein_bins[protein_bins$bin_arch == '_1.1.5_',]$sum +
  protein_bins[protein_bins$bin_arch == '_1.1.6_',]$sum 

  protein_bins[protein_bins$bin_arch == '_1.3_',]$sum <-     protein_bins[protein_bins$bin_arch == '_1.3_',]$sum +
  protein_bins[protein_bins$bin_arch == '_1.3.1_',]$sum +
  protein_bins[protein_bins$bin_arch == '_1.3.2_',]$sum

#depends if we want to include heat stress with abiotic stress
  #protein_bins[protein_bins$bin_arch == '_20.2_',]$sum <- protein_bins[protein_bins$bin_arch == '_20.2_',]$sum +
  #protein_bins[protein_bins$bin_arch == '_20.2.1_',]$sum 

  protein_bins[protein_bins$bin_arch == '_29_',]$sum <- protein_bins[protein_bins$bin_arch == '_29_',]$sum +
  protein_bins[protein_bins$bin_arch == '_29.6_',]$sum

# add some summed bins for special cases

rubisco <- protein_bins[protein_bins$bin_arch %in% c('_1.3.1_', '_1.3.2_'),]
rubisco <- ddply(na.omit(rubisco), .(sample), summarise, sum = sum(sum))

rubisco$bin_arch <- '_1.3.1_+_1.3.2_'

photosystems <- protein_bins[protein_bins$bin_arch %in% c('_1.1.1_', '_1.1.2_'),]
photosystems <- ddply(na.omit(photosystems), .(sample), summarise, sum = sum(sum))
photosystems$bin_arch <- '_1.1.1_+_1.1.2_'

TCA_orgtrans <- protein_bins[protein_bins$bin_arch %in% c('_8.1_', '_8.2_'),]
TCA_orgtrans <- ddply(na.omit(TCA_orgtrans), .(sample), summarise, sum = sum(sum))
TCA_orgtrans$bin_arch <- '_8.1_+_8.2_'

redox <- protein_bins[protein_bins$bin_arch %in% c('_21_', '_26.9_'),]
redox <- ddply(na.omit(redox), .(sample), summarise, sum = sum(sum))
redox$bin_arch <- '_21_+_26.9_'

electron_transport_minATPsynth <- protein_bins[protein_bins$bin_arch %in% c('_1.1.3_', '_1.1.5_', '_1.1.6_'),]
electron_transport_minATPsynth <- ddply(na.omit(electron_transport_minATPsynth), .(sample), summarise, sum = sum(sum))
electron_transport_minATPsynth$bin_arch <- '_1.1.3_+_1.1.5_+_1.1.6_'

protein_bins <- rbind(protein_bins[,c('bin_arch','sample','sum')], rubisco, photosystems, TCA_orgtrans, redox, electron_transport_minATPsynth)
rm(rubisco,photosystems,TCA_orgtrans,redox, electron_transport_minATPsynth)
protein_bins$sample <- as.character(protein_bins$sample)

# merge in bin names from mercator_bins

protein_bins$bin_arch <- gsub(" ", "", protein_bins$bin_arch) # remove _'s

unique(protein_bins$bin_arch)[!unique(protein_bins$bin_arch) %in% mercator_bins$bin_arch]

protein_bins <- protein_bins[!protein_bins$bin_arch %in% c('_1.3.1_', '1.3.2_'),]

protein_bins$bin_arch_name <- NA


for(i in 1:length(unique(protein_bins$bin_arch))) {
  x <- unique(protein_bins$bin_arch)[i]
  y <- mercator_bins[grep(x, mercator_bins$bin_arch, fixed =TRUE),]$bin_name
  
  if (x %in% mercator_bins$bin_arch) {
    protein_bins[grep(x, paste(protein_bins$bin_arch, "_", sep = ""), fixed=TRUE),]$bin_arch_name <- y[1]
  } else {
    protein_bins[grep(x, paste(protein_bins$bin_arch, "_", sep = ""), fixed=TRUE),]$bin_arch_name <- NA
  }
}

# get climate data

#climate_locs <- merge(sample_locations, climate, by = c('Longitude', 'Latitude'))
#climate_locs <- climate_locs[!duplicated(climate_locs$sample),]
#climate_locs$Var.4 <- NULL
#rm(sample_locations,climate)
#climate_locs <- climate_locs[,c('tavg', 'prec', 'sample')]

#protein_climate <- merge(climate_locs, protein_bins, by = 'sample')
