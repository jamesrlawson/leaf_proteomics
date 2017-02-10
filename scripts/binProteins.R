## PROCESS PROTEIN DATA INTO BINS BASED ON MERCATOR FUNCTIONAL CATEGORIES ##

source('scripts/functions.R')

require(readr)
require(stringr)
library(dplyr)
library(tidyr)
require(reshape2)
require(plyr)

mercator <- read_csv('data/mercator/D14_mercator.csv')

mercator_bins <- read_csv('data/mercator/mercator_bins.csv')

# note that for this to work, all bincodes in asp_mercator.csv and mercator_bins.csv must be prefixed and suffixed by '_'
# the following is a hierarchically ordered list of bincodes to be grepped over, where the first string in each list element defines the top level, 
# and therefore a grep using this element as the pattern will only find proteins in the top level of the hierarchy,
# while the second element allows the grep call to find all proteins below the top category
# it is essential that elements in this list are ordered sequentially in a hierarchical fashion (i.e. higher levels must come before lower levels)
# the order of elements which are at the same level is not relevant

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

## DISCOVERY DATA ##

#protein_samples_D14 <- read_csv('data/D14_protein_sites.csv') # these are protein amounts calculated using ovalbumin equivalents
#protein_samples_D14 <- read_csv('data/protein_amounts_by_signal_fraction_perArea_D14.csv') # these are protein amounts calculated using signal intensity fraction
#protein_samples_D14 <- read_csv('data/D14_protein_GGLEP-DEDT.csv') # protein amounts calculated using D14 ion library, in avg(GGLEP/DEDT) equivalents
protein_samples_D14 <- read_csv('data/D14_protein_moles_GGLEP-DEDT.csv') # protein amounts as above but in moles (not multiplied by MW)

protein_samples_D14 <- protein_samples_D14[!names(protein_samples_D14) %in% c('YG029','YG031','YG030','YG028')]

total_protein_D14 <- getTotalProtein(protein_samples_D14)

protein_samples_D14 <- getProteinBins(protein_samples_D14, mercator)

protein_bins_D14 <- populateProteinBins(protein_samples_D14, bin_arch.list)
#protein_bins_D14 <- populateProteinBins_mean(protein_samples_D14, bin_arch.list)


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
