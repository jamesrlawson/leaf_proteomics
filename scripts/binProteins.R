## ASSIGN PROTEINS TO FUNCTIONAL CATEGORIES AND CALCULATE PROTEIN AMOUNTS IN EACH CATEGORY ##

# called by transformations.R

require(readr)
require(stringr)
require(dplyr)
require(tidyr)

source('scripts/functions.R')

mercator <- read_csv('data/proteomics_data/mercator/euc/mercator_beta/D14_mercator4_20180218.csv')

# 'mg_per_m2' and 'moles' switches are defined in transformations.R
if(mg_per_m2) {
  protein_samples_D14 <- read_csv('data/proteomics_data/proteomics/derived/euc/D14_protein_GGLEP-DEDT.csv') # protein amounts calculated using D14 ion library, in avg(GGLEP/DEDT) equivalents
}

if(moles) {
  protein_samples_D14 <- read_csv('data/proteomics_data/proteomics/derived/euc/D14_protein_moles_GGLEP-DEDT.csv') # protein amounts as above but in moles (not multiplied by MW)
}

# first add the mercator$NAME values for each protein in protein_samples_D14

protein_samples_D14 <- getProteinBins(protein_samples_D14, mercator)

# then import the names of categories we're interested in from mercator_names* and use grep to find all proteins associated with those categories
# put the results in instances of a list
# mercator_names.csv contains search terms for protein categories. These searches are run on mercator$NAMES. 
# search terms must be unique to the functional category to avoid non-target returns. 
# for example, a search for 'photosystem I' will also pick up proteins from 'photosystem II' - to avoid this we search for 'photosystem I\.'
# this works because all instances of proteins within 'photosystem I' are actually within subcategories. We'd miss some returns if there were proteins in the upper 'photosystem I' category.
  # N.B. the '\' is an 'escape' and must be used because .'s are special in regular expressions and mean 'anything'. By using the escape we will actually search for the character '.'
# search terms for top level categories can be made unique by using ' in front

mercator_names <- read.csv('data/proteomics_data/mercator/mercator4_names.csv', header=T, stringsAsFactors = F) 
mercator_names <- arrange(mercator_names, funccat)

func_assigned.list <- vector('list', length(mercator_names$funccat))

func_assigned <- data.frame()

for(i in 1:length(mercator_names$funccat)) {
  
  name <- mercator_names$funccat[i]
  
  proteins <- protein_samples_D14[grep(name, protein_samples_D14$NAME),]
  
  proteins$funccat <- name
  
  proteins <- distinct(proteins, Protein, .keep_all = TRUE)
  
  func_assigned.list[[i]] <- proteins
  
}

func_assigned <- rbind(func_assigned, do.call(rbind, func_assigned.list))
rm(func_assigned.list)

# generate df with summed protein amounts for each category

funccat_sums <- func_assigned %>% group_by(funccat) %>% summarise_at(vars(2:(ncol(protein_samples_D14)-2)), sum) %>% arrange(funccat)
funccat_sums$funccat <- mercator_names$funccat_rename

# transform df

funccat_sums_t <- t(funccat_sums[,2:ncol(funccat_sums)])
colnames(funccat_sums_t) <- funccat_sums$funccat
funccat_sums_t <- as.data.frame(funccat_sums_t)
funccat_sums_t$sample <- colnames(funccat_sums)[2:ncol(funccat_sums)]
funccat_sums <- funccat_sums_t
rm(funccat_sums_t)

# create some custom categories

  # funccat_sums$PSII_min_LHCII <- funccat_sums$photosystem_II - funccat_sums$LHC_I
  # funccat_sums$PSI_min_LHCI <- funccat_sums$photosystem_II - funccat_sums$LHC_II
  # funccat_sums$Photosystems <- funccat_sums$photosystem_I + funccat_sums$photosystem_II
  # funccat_sums$electron_transport_minATPsynth <- funccat_sums$other_electron_carrier + funccat_sums$cytochrome_b6f
  # funccat_sums$Rubisco <- funccat_sums$rubisco_large_subunit + funccat_sums$rubisco_small_subunit
  # funccat_sums$redox <- funccat_sums$redox + funccat_sums$glutathione_S_transferases
  # funccat_sums$stress <- funccat_sums$stress + funccat_sums$glutathione_S_transferases
  # 
  # funccat_sums$electron_transport <- funccat_sums$electron_transport_minATPsynth + funccat_sums$ATP_synthase_chloroplastic
  # funccat_sums$LHC <- funccat_sums$LHC_I + funccat_sums$LHC_II
  # funccat_sums$Photosystems_min_LHC <- funccat_sums$Photosystems - funccat_sums$LHC
  # funccat_sums$LHCI_per_PSI <- funccat_sums$LHC_I / funccat_sums$PSI_min_LHCI
  # funccat_sums$LHCII_per_PSII <- funccat_sums$LHC_II / funccat_sums$PSII_min_LHCII
  # funccat_sums$LHC_per_PS <- funccat_sums$LHC / funccat_sums$Photosystems

# add in total (detected) protein and get relative abundances

protein_D14 <- protein_samples_D14 %>% gather(key = sample, value = value, 2:(ncol(.)-2)) %>%
  group_by(sample) %>%
  dplyr::summarise(total_protein = sum(value, na.rm=TRUE)) %>%
  full_join(funccat_sums, by = 'sample')

protein_stand_D14 <- protein_D14 %>% 
  mutate_at(vars(3:(ncol(.)-3)), funs(. / total_protein)) # divide columns by total protein to get relative abundances

# cleanup

rm(func_assigned, mercator, mercator_names, proteins, funccat_sums, i, name)
gc(verbose = FALSE)

