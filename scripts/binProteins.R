## PROCESS PROTEIN DATA INTO BINS AND STANDARDISE ##

source('scripts/functions.R')

require(readr)
require(stringr)
library(dplyr)
library(tidyr)
require(reshape2)
require(plyr)

## ASPINWALL DATA ##

#protein_samples_asp <- read_csv('data/asp_proteins.csv')
mercator <- read_csv('data/asp_mercator.csv')
mercator_bins <- read_csv('data/asp_mercator_bins.csv')

#total_protein_asp <- getTotalProtein(protein_samples_asp)

#protein_samples_asp <- getProteinBins(protein_samples_asp, mercator)

bin_arch.list <- list(

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
  c('_20_','_20.'),
  c('_20.1_','_20.1.'),
  c('_20.2_', '_20.2.'),
  c('_20.2.1_','_20.2.1.'),
  c('_21_','_21.'),
  c('_22_','_22.'),
  c('_26.9_','_26.9.'),
  c('_30_','_30.'),
  c('_29_','_29.'),
  c('_29.6_','_29.6.'),
  c('_28_','_28.'),
  c('_27_','_27.'))

#protein_bins_asp <- populateProteinBins(protein_samples_asp, bin_arch.list)

# absolute quants, unstandardised

#protein_asp <- na.omit(protein_bins_asp[,c('sample','bin_arch_name','sum')])
#protein_asp <- protein_asp[!protein_asp$sample %in% c('BINCODE','NAME'),]
#protein_asp <- protein_asp[!duplicated(protein_asp[,c('sample', 'bin_arch_name')]),]
#protein_asp <- spread(protein_asp, key = bin_arch_name, value=sum)
#protein_asp <- protein_asp[!protein_asp$sample %in% c('BINCODE', 'NAME'),]
#protein_asp$total_protein <- total_protein_asp$total_protein

# relative quants, standardised by total_protein

#protein_bins_asp_spread <- na.omit(protein_bins_asp[,c('sample','bin_arch_name','sum')])
#protein_bins_asp_spread <- spread(protein_bins_asp_spread, key = bin_arch_name, value=sum)
#protein_bins_asp_spread <- protein_bins_asp_spread[!protein_bins_asp_spread$sample %in% c('BINCODE', 'NAME'),]

#protein_stand_asp <- merge(protein_bins_asp_spread, total_protein_asp, by = 'sample')
#protein_stand_asp[,2:26] <- protein_stand_asp[,2:26]/protein_stand_asp$total_protein
#rm(protein_bins_asp_spread)

# relative quants, standardised by Rubisco

#protein_bins_asp <- protein_bins_asp[,c('sample', 'bin_arch_name','sum')]
#protein_bins_asp <- na.omit(protein_bins_asp)
#protein_bins_asp <- spread(protein_bins_asp, key = bin_arch_name, value=sum)
#protein_bins_asp[,2:26] <- protein_bins_asp[,2:26]/protein_bins_asp$Rubisco

#protein_RbcStand_asp <- ddply(melt(protein_bins_asp[,c('Rubisco','Calvin_cycle','Cytochrome_b6f','PSI','PSII','ATP_synthase_chloroplastic', 'electron_transport_minATPsynth')]), 
#                                  .(variable), 
#                                  summarise,
#                                  mean = mean(value, na.rm=TRUE),
#                                  CV = CV(value))

## DISCOVERY DATA ##

protein_samples_D14 <- read_csv('data/D14_protein_sites.csv')

total_protein_D14 <- getTotalProtein(protein_samples_D14)

protein_samples_D14 <- getProteinBins(protein_samples_D14, mercator)

protein_bins_D14 <- populateProteinBins(protein_samples_D14, bin_arch.list)

# relative quants, standardised by total_protein

protein_bins_D14_spread <- na.omit(protein_bins_D14[,c('sample','bin_arch_name','sum')])
protein_bins_D14_spread <- protein_bins_D14_spread[!protein_bins_D14_spread$sample %in% c('BINCODE','NAME'),]

protein_bins_D14_spread <- protein_bins_D14_spread[!duplicated(protein_bins_D14_spread[,c('sample', 'bin_arch_name')]),]

protein_bins_D14_spread <- spread(protein_bins_D14_spread, key = bin_arch_name, value=sum)
protein_bins_D14_spread <- protein_bins_D14_spread[!protein_bins_D14_spread$sample %in% c('BINCODE', 'NAME'),]

protein_stand_D14 <- merge(protein_bins_D14_spread, total_protein_D14, by = 'sample')
protein_stand_D14[,2:28] <- protein_stand_D14[,2:28]/protein_stand_D14$total_protein
rm(protein_bins_D14_spread)

# absolute quants, unstandardised

protein_D14 <- na.omit(protein_bins_D14[,c('sample','bin_arch_name','sum')])
protein_D14 <- protein_D14[!protein_D14$sample %in% c('BINCODE','NAME'),]
protein_D14 <- protein_D14[!duplicated(protein_D14[,c('sample', 'bin_arch_name')]),]
protein_D14 <- spread(protein_D14, key = bin_arch_name, value=sum)
protein_D14 <- protein_D14[!protein_D14$sample %in% c('BINCODE', 'NAME'),]
protein_D14$total_protein <- protein_stand_D14$total_protein
