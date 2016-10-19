# ASPINWALL PROTEIN AMOUNTS RECALCULATED USING SIGNAL AREA FRACTIONS

require(readr)
require(dplyr)
require(stringr)
require(tidyr)

FDR <- read_csv('data/aspinwall/aspinwall_FDR.csv')

FDR <- FDR[FDR$Decoy == FALSE,]

FDR <- FDR %>% mutate(true_id = as.numeric(any(.[8:ncol(FDR)] < 0.01))) %>% filter(true_id == 1)

ion_areas <- read_csv('data/aspinwall/aspinwall_ion_areas.csv')

top2 <- function(x) {
  sum(sort(x, decreasing = TRUE)[1:2], na.rm=TRUE)
}

top3 <- function(x) {
  mean(sort(x, decreasing = TRUE)[1:3], na.rm=TRUE)
}

protein_areas <-  ion_areas %>% 
                  group_by(Peptide, Protein) %>% 
                  summarise_at(vars(matches("sample")), top2) %>%
                  group_by(Protein) %>%
                  summarise_at(vars(matches("sample")), top3)

protein_MW <- read_csv('data/aspinwall/Protein_MW.csv')

protein_areas <- merge(protein_areas, protein_MW, by = 'Protein')

protein_fractions <- protein_areas

protein_fractions[,c(2:97)] <- protein_fractions[,c(2:97)] * protein_fractions$MW

protein_fractions[,2:97] <- protein_fractions[,2:97] / colSums(protein_fractions[,2:97])
                           

for(i in 2:97) {
  
  name_split <- unlist(str_split(names(protein_fractions)[i], " "))[1]
  names(protein_fractions)[i] <- name_split
  
}

asp_spreadsheet <- read.csv('data/aspinwall/aspinwall_spread.csv', header=TRUE, stringsAsFactors = TRUE)
asp_spreadsheet$ID <- gsub("[[:punct:]]", "", asp_spreadsheet$`sample.id..day.species.abbreviation.pot.`)

protein_amounts <- gather(protein_fractions, 'ID', 'protein_fraction', 2:97)
                
protein_amounts <- merge(protein_amounts, asp_spreadsheet, by = "ID")       
                
protein_amounts$protein_per_dry_weight <- protein_amounts$protein_fraction * protein_amounts$Protein.per.DW..mg.g..protein.assay.
                
protein_amounts$protein_per_leaf_area <- protein_amounts$protein_fraction * protein_amounts$Protein.per.Leaf.Area..mg.cm2..protein.assay. * 10000

##

protein_per_dw <- spread(protein_amounts[,c('Protein', 'ID', 'protein_per_dry_weight')], key = 'ID', value= 'protein_per_dry_weight')                
protein_per_leaf_area <- spread(protein_amounts[,c('Protein', 'ID', 'protein_per_leaf_area')], key = 'ID', value= 'protein_per_leaf_area')                

write_csv(protein_per_dw, 'output/protein_amounts_by_signal_fraction_perDW.csv')
write_csv(protein_per_leaf_area, 'output/protein_amounts_by_signal_fraction_perArea.csv')


##

protein_amounts$treatment <- paste(protein_amounts$day, protein_amounts$Temperature.treatment, protein_amounts$species.abbreviation)

CV <- function(x){
  sqrt(var(x))/mean(x)
}

fractions_CV <- protein_amounts %>% group_by(treatment, ID) %>%
                                    summarise(sum_protein_per_DW = sum(protein_per_dry_weight, na.rm=TRUE), sum_protein_per_area = sum(protein_per_leaf_area, na.rm=TRUE)) %>%
                                    group_by(treatment) %>%
                                    summarise(CV_protein_per_DW_fractions = CV(sum_protein_per_DW), CV_protein_per_area_fractions = CV(sum_protein_per_area), N_fractions = length(sum_protein_per_DW))
        
        

####### OVALB QUANT #######

ovalb_quant <- read.csv('data/aspinwall/asp_ovalbQuant.csv', header=TRUE, stringsAsFactors = FALSE)
names(ovalb_quant)[2:97] <- names(protein_areas)[2:97]

for(i in 2:97) {
  
  name_split <- unlist(str_split(names(ovalb_quant)[i], " "))[1]
  names(ovalb_quant)[i] <- name_split
  
}

ovalb_amounts <- gather(ovalb_quant, 'ID', 'protein_amount', 2:97)
ovalb_amounts <- merge(ovalb_amounts, asp_spreadsheet, by = "ID")   

ovalb_amounts$treatment <- paste(ovalb_amounts$day, ovalb_amounts$Temperature.treatment, ovalb_amounts$species.abbreviation)

ovalb_CV <- ovalb_amounts %>% group_by(treatment, ID) %>%
                      summarise(sum_protein_per_area = sum(protein_amount, na.rm=TRUE)) %>%
                      group_by(treatment) %>%
                      summarise(CV_protein_per_area_ovalb = CV(sum_protein_per_area), N_ovalb = length(sum_protein_per_area))


method_comparison <- merge(fractions_CV, ovalb_CV)

write_csv(method_comparison, 'output/protein_amount_calculation_method_comparison.csv')
