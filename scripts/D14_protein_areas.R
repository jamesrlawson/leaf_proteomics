require(readxl)
require(dplyr)
require(stringr)
require(readr)
require(tidyr)

load_excel_files <- function(path, sheet) { 
  
  load_excel_arrange <- function(path, sheet) {
    
    out <- read_excel(path, sheet)
    out <- arrange(out, Protein)
    return(out)
    
  }
  
  filenames = paste0(path, list.files(path))
  tables <- lapply(filenames, load_excel_arrange, sheet = sheet)
  do.call(cbind, tables)
}

# load ion areas, rename columns

euc_swath_reanalysed <- load_excel_files('data/large_files/swath/', 'Area - ions')
euc_swath_reanalysed <- cbind(euc_swath_reanalysed[,1:9], euc_swath_reanalysed[,10:377] %>% select(contains('sample')))
names(euc_swath_reanalysed)[10:323] <- do.call(rbind, str_split(names(euc_swath_reanalysed[,10:323]), pattern = " "))[,1]

# load FDR data, rename columns, remove decoys and discard proteins below identification probability threshold

euc_swath_reanalysed_FDR <- load_excel_files('data/large_files/swath/', 'FDR')
euc_swath_reanalysed_FDR <- cbind(euc_swath_reanalysed_FDR[,1:7], euc_swath_reanalysed_FDR[,10:363] %>% select(contains('sample')))
names(euc_swath_reanalysed_FDR)[8:319] <- do.call(rbind, str_split(names(euc_swath_reanalysed_FDR[,8:319]), pattern = " "))[,2]
euc_swath_reanalysed_FDR <- subset(euc_swath_reanalysed_FDR, Decoy == 0)

euc_swath_reanalysed_FDR$FDR <- apply(euc_swath_reanalysed_FDR[8:ncol(euc_swath_reanalysed_FDR)],1, function(x) as.numeric(any(x < 0.01)) )
nrow(euc_swath_reanalysed_FDR[euc_swath_reanalysed_FDR$FDR == 0,])

euc_swath_reanalysed_FDR <- subset(euc_swath_reanalysed_FDR, FDR == 1)

# ion areas constrained by FDR

ion_areas <- euc_swath_reanalysed[euc_swath_reanalysed$Peptide %in% euc_swath_reanalysed_FDR$Peptide,]
rm(euc_swath_reanalysed_FDR,euc_swath_reanalysed)

#write_csv(ion_areas, 'data/large_files/D14_ion_areas_new_ion_library.csv')

ion_areas <- read_csv('data/large_files/D14_ion_areas_new_ion_library.csv')

# get protein areas using top2top3 method

top2 <- function(x) {
  sum(sort(x, decreasing = TRUE)[1:2], na.rm=TRUE)
}

top3 <- function(x) {
  mean(sort(x, decreasing = TRUE)[1:3], na.rm=TRUE)
}

protein_areas <-  ion_areas %>% 
  group_by(Peptide, Protein) %>% 
  summarise_at(vars(10:323), top2) %>%
  group_by(Protein) %>%
  summarise_at(vars(3:316), top3)

write_csv(protein_areas, 'data/large_files/protein_areas.csv')


protein_areas <- read_csv('data/large_files/protein_areas.csv')


# get GGLEP top2 areas for each sample

GGLEP <- ion_areas[ion_areas$Peptide == "GGLEPINFQTAADQAR",]
GGLEP <- GGLEP %>% summarise_at(vars(10:323), top2)

# find protein areas relative to GGLEP area
protein_areas[,2:315] <- t(t(protein_areas[,2:315])/as.vector(t(GGLEP))) 

# multiply by 5.64x10^-11 to get moles per cm2

protein_areas[,2:315] <- protein_areas[,2:315]*(5.64e-11)

# multiply by the molecular weight to get g/cm2

protein_MW <- rbind(read_csv('data/D14_protein_MW.csv'), read_csv('data/extraMWs.csv'))

protein_areas <- protein_areas[!protein_areas$Protein %in% setdiff(protein_areas$Protein, protein_MW$Protein),] # a few MW's missing, not important?

protein_areas <- arrange(protein_areas, Protein)
protein_MW <- arrange(protein_MW, Protein)

#setdiff(protein_areas$Protein, protein_MW$Protein)
#setdiff(protein_MW$Protein, protein_areas$Protein)

protein_amounts <- protein_areas

protein_amounts[,2:315] <- protein_amounts[,2:315] * protein_MW$MW # multiply by MW to get g/cm2

protein_amounts[,2:315] <- protein_amounts[,2:315] * (1e07) # multiply by 10^07 to get mg/m2

quality_check<-protein_amounts[which(protein_amounts$Protein=="sp|OVAL_CHICK"),]

write.csv(data,"protein_areas_converted_per_leaf_area.csv")
































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



