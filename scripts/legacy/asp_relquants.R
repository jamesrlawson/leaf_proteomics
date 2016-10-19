# this script outputs protein functional category amounts for the Aspinwall dataset
# amounts are quantified relative to Rubisco

#protein bins ASPINWALL

source('../scripts/functions.R')

require(readr)
require(stringr)
library(dplyr)
library(tidyr)
require(reshape2)
require(plyr)

protein_samples <- read_csv('../data/asp_proteins.csv')
mercator <- read_csv('../data/asp_mercator.csv')
mercator_bins <- read_csv('../data/asp_mercator_bins.csv')
sample_locations <- read_csv('../data/sample_locations.csv')
climate <- read_csv('../data/discovery_site_climate.csv')


protein_samples <- getProteinBins(protein_samples, mercator)

#protein_samples <- multiple_bins(protein_samples)

# sum values within overarching bins

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


protein_bins <- populateProteinBins(protein_samples, bin_arch.list)


# get relative quants (standardised by Rubisco)



protein_bins1 <-protein_bins[,c('sample', 'bin_arch_name','sum')]
protein_bins1 <- na.omit(protein_bins1)
protein_bins1 <- spread(protein_bins1, key = bin_arch_name, value=sum)
protein_bins1[,2:26] <- protein_bins1[,2:26]/protein_bins1$Rubisco

blaz <- ddply(melt(protein_bins1[,c('Rubisco','Calvin_cycle','Cytochrome_b6f','PSI','PSII','ATP_synthase_chloroplastic', 'electron_transport_minATPsynth')]), 
              .(variable), 
              summarise,
              mean = mean(value, na.rm=TRUE),
              CV = CV(value))
write_csv(blaz, '../output/aspinwall_proteins_standRbc.csv')





