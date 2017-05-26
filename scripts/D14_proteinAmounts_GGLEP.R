require(dplyr)
require(stringr)
require(readr)
require(tidyr)

source('scripts/functions.R')

### calculate protein amounts in GGLEP_DEDT equivalents ###
 
  protein_areas <- read_csv('data/proteomics_data/proteomics/derived/euc/protein_areas_top2top2.csv')
  ion_areas <- read_csv('data/proteomics_data/proteomics/derived/euc/D14_ion_areas_new_ion_library.csv')
  
  # get GGLEP top2 areas for each sample
  
   # GGLEP <- ion_areas[ion_areas$Peptide == 'GGLEPINFQTAADQAR',]
   # GGLEP <- GGLEP %>% summarise_at(vars(10:323), top2)
  
  # find protein areas relative to GGLEP area
    # protein_areas[,2:ncol(protein_areas)] <- t(t(protein_areas[,2:ncol(protein_areas)])/as.vector(t(GGLEP))) 
  
    #ovalb <- protein_areas[protein_areas$Protein == "sp|OVAL_CHICK",] # or use ovalb top2top3
    #protein_areas[,2:ncol(protein_areas)] <- t(t(protein_areas[,2:ncol(protein_areas)])/as.vector(t(ovalb[,2:ncol(ovalb)]))) 
  
  # find protein areas relative to GGLEP/DEDT averaged area
  
    # get GGLEP/DEDT top2/top2avg
    
    GGLEP_DEDT <- ion_areas[ion_areas$Peptide %in% c('GGLEPINFQTAADQAR', 'DEDTQAMPFR'),]
    GGLEP_DEDT <- GGLEP_DEDT %>% group_by(Peptide) %>% summarise_at(vars(10:ncol(GGLEP_DEDT)), top2)
    
    GGLEP_DEDT <- data.frame(t(rowMeans(t(GGLEP_DEDT[,2:ncol(GGLEP_DEDT)]))))
    
    protein_areas[,2:ncol(protein_areas)] <- t(t(protein_areas[,2:ncol(protein_areas)])/as.vector(t(GGLEP_DEDT))) 
  
  # multiply by 5.64x10^-11 (2.5 * 10^-6 g/cm2 / MW of ovalbumin - 44287) to get moles per cm2
  
  protein_areas[,2:ncol(protein_areas)] <- protein_areas[,2:ncol(protein_areas)]*(5.64e-11)
  
  
  write_csv(protein_areas, "data/proteomics_data/proteomics/derived/euc/D14_protein_moles_GGLEP-DEDT.csv")
  
  # multiply by the molecular weight to get g/cm2
  
  source('scripts/protein_MW.R')
  
  protein_MW <- rbind(read_csv('data/proteomics_data/sequences/euc/D14_protein_MW.csv'), read_csv('data/proteomics_data/sequences/euc/extraMWs.csv'))
  
  protein_areas <- protein_areas[!protein_areas$Protein %in% setdiff(protein_areas$Protein, protein_MW$Protein),]
  
  protein_areas <- arrange(protein_areas, Protein)
  protein_MW <- arrange(protein_MW, Protein)
  
  setdiff(protein_areas$Protein, protein_MW$Protein)
  setdiff(protein_MW$Protein, protein_areas$Protein)
  
  protein_amounts <- protein_areas
  
  protein_amounts[,2:ncol(protein_amounts)] <- protein_amounts[,2:ncol(protein_amounts)] * protein_MW$MW # multiply by MW to get g/cm2
  
  protein_amounts[,2:ncol(protein_amounts)] <- protein_amounts[,2:ncol(protein_amounts)] * (10^7) # multiply by 10^07 to get mg/m2
  
  
  write_csv(protein_amounts,"data/proteomics_data/proteomics/derived/euc/D14_protein_GGLEP-DEDT.csv")


  
  
  
  
  
  qual_check <- function() {
    
    protein_assay <- read_csv('data/proteomics_data/sequences/euc/protein_assay.csv')
    
    quality_check <- protein_amounts[protein_amounts$Protein=="sp|OVAL_CHICK",]
    
    blah <- data.frame(t(quality_check[,2:ncol(quality_check)]))
    blah$sample <- rownames(blah)
    names(blah)[1] <- 'ovalb_amount_MS_mg_per_m2'
    
    blax <- merge(protein_assay, blah, by = 'sample')
    blax$ovalb_amount_assay_mg_per_m2 <- blax$assay_protein_per_area_mg_per_m2 * blax$percent_ovalb_of_total_protein / 100
    plot(blax$ovalb_amount_MS_mg_per_m2 ~ blax$ovalb_amount_assay_mg_per_m2, 
         xlab = "ovalbumin amount assay measured (mg/m2)",
         ylab = "ovalbumin amount MS measured (mg/m2)")

  }
  
  qual_check()
  
  
quality_check <- protein_amounts[protein_amounts$Protein=="sp|OVAL_CHICK",]
  
  
  







# find protein amounts according to target protein area fraction of total protein area

protein_fractions <- protein_areas

cols <- 2:(ncol(protein_fractions)-1) # indices denoting sample columns

protein_fractions[,cols] <- protein_fractions[,cols] * protein_fractions$MW
protein_fractions[,cols] <- t(t(protein_fractions[,cols])/colSums(protein_fractions[,cols]))

for(i in cols) {
  
  name_split <- unlist(str_split(names(protein_fractions)[i], " "))[1]
  names(protein_fractions)[i] <- name_split
  
}

protein_assay <- read_csv('data/proteomics_data/sequences/euc/protein_assay.csv')
protein_amounts <- gather(protein_fractions, 'sample', 'protein_fraction', cols)
protein_amounts <- merge(protein_amounts, protein_assay, by = "sample")       
protein_amounts$protein_per_dry_weight <- protein_amounts$protein_fraction * protein_amounts$assay_protein_per_area_mg_per_g
protein_amounts$protein_per_leaf_area <- protein_amounts$protein_fraction * protein_amounts$assay_protein_per_area_mg_per_m2
protein_per_dw <- spread(protein_amounts[,c('Protein', 'sample', 'protein_per_dry_weight')], key = 'sample', value= 'protein_per_dry_weight')                
protein_per_leaf_area <- spread(protein_amounts[,c('Protein', 'sample', 'protein_per_leaf_area')], key = 'sample', value= 'protein_per_leaf_area')                
write_csv(protein_per_dw, 'output/protein_amounts_by_signal_fraction_perDW_D14.csv')
write_csv(protein_per_leaf_area, 'output/protein_amounts_by_signal_fraction_perArea_D14.csv')

