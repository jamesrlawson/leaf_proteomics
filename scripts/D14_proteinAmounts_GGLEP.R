### calculate protein amounts in GGLEP equivalents ###
 
  protein_areas <- read_csv('data/large_files/protein_areas.csv')
  
  # get GGLEP top2 areas for each sample
  
  GGLEP <- ion_areas[ion_areas$Peptide == "GGLEPINFQTAADQAR",]
  GGLEP <- GGLEP %>% summarise_at(vars(10:323), top2)
  
  # find protein areas relative to GGLEP area
  protein_areas[,2:315] <- t(t(protein_areas[,2:315])/as.vector(t(GGLEP))) 
  
  # multiply by 5.64x10^-11 to get moles per cm2
  
  protein_areas[,2:315] <- protein_areas[,2:315]*(5.64e-11)
  
  # multiply by the molecular weight to get g/cm2
  
  source('scripts/protein_MW.R')
  
  protein_MW <- rbind(read_csv('data/D14_protein_MW.csv'), read_csv('data/extraMWs.csv'))
  
  protein_areas <- protein_areas[!protein_areas$Protein %in% setdiff(protein_areas$Protein, protein_MW$Protein),]
  
  protein_areas <- arrange(protein_areas, Protein)
  protein_MW <- arrange(protein_MW, Protein)
  
  setdiff(protein_areas$Protein, protein_MW$Protein)
  setdiff(protein_MW$Protein, protein_areas$Protein)
  
  protein_amounts <- protein_areas
  
  protein_amounts[,2:315] <- protein_amounts[,2:315] * protein_MW$MW # multiply by MW to get g/cm2
  
  protein_amounts[,2:315] <- protein_amounts[,2:315] * (1e07) # multiply by 10^07 to get mg/m2
  
  
  write_csv(protein_amounts,"data/D14_protein_GGLEP.csv")


  qual_check <- function() {
    
    protein_assay <- read_csv('data/protein_assay.csv')
    
    quality_check <- protein_amounts[protein_amounts$Protein=="sp|OVAL_CHICK",]
    
    blah <- data.frame(t(quality_check[,2:315]))
    blah$sample <- rownames(blah)
    names(blah)[1] <- 'ovalb_amount_MS_mg_per_m2'
    
    blax <- merge(protein_assay, blah, by = 'sample')
    blax$ovalb_amount_assay_mg_per_m2 <- blax$assay_protein_per_area_mg_per_m2 * blax$percent_ovalb_of_total_protein / 100
    plot(blax$ovalb_amount_MS_mg_per_m2 ~ blax$ovalb_amount_assay_mg_per_m2, 
         xlab = "ovalbumin amount assay measured (mg/m2)",
         ylab = "ovalbumin amount MS measured (mg/m2)")

  }
  
  qual_check()
  
# find protein amounts according to target protein area fraction of total protein area

protein_fractions <- protein_areas

cols <- 2:(ncol(protein_fractions)-1) # indices denoting sample columns

protein_fractions[,cols] <- protein_fractions[,cols] * protein_fractions$MW
protein_fractions[,cols] <- t(t(protein_fractions[,cols])/colSums(protein_fractions[,cols]))

for(i in cols) {
  
  name_split <- unlist(str_split(names(protein_fractions)[i], " "))[1]
  names(protein_fractions)[i] <- name_split
  
}

protein_assay <- read_csv('data/protein_assay.csv')
protein_amounts <- gather(protein_fractions, 'sample', 'protein_fraction', cols)
protein_amounts <- merge(protein_amounts, protein_assay, by = "sample")       
protein_amounts$protein_per_dry_weight <- protein_amounts$protein_fraction * protein_amounts$assay_protein_per_area_mg_per_g
protein_amounts$protein_per_leaf_area <- protein_amounts$protein_fraction * protein_amounts$assay_protein_per_area_mg_per_m2
protein_per_dw <- spread(protein_amounts[,c('Protein', 'sample', 'protein_per_dry_weight')], key = 'sample', value= 'protein_per_dry_weight')                
protein_per_leaf_area <- spread(protein_amounts[,c('Protein', 'sample', 'protein_per_leaf_area')], key = 'sample', value= 'protein_per_leaf_area')                
write_csv(protein_per_dw, 'output/protein_amounts_by_signal_fraction_perDW_D14.csv')
write_csv(protein_per_leaf_area, 'output/protein_amounts_by_signal_fraction_perArea_D14.csv')

