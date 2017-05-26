require(readxl)
require(dplyr)
require(stringr)
require(readr)
require(tidyr)

source('scripts/functions.R')

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

euc_swath_reanalysed <- load_excel_files('data/proteomics_data/proteomics/raw/euc/swath', 'Area - ions')
euc_swath_reanalysed <- cbind(euc_swath_reanalysed[,1:9], euc_swath_reanalysed[,10:ncol(euc_swath_reanalysed)] %>% select(contains('sample')))
names(euc_swath_reanalysed)[10:ncol(euc_swath_reanalysed)] <- do.call(rbind, str_split(names(euc_swath_reanalysed[,10:ncol(euc_swath_reanalysed)]), pattern = " "))[,1]

euc_swath_reanalysed <- euc_swath_reanalysed %>% subset(., select=which(!duplicated(names(.)))) # remove duplicated columns (e.g KG022)


# load FDR data, rename columns, remove decoys and discard proteins below identification probability threshold

euc_swath_reanalysed_FDR <- load_excel_files('data/proteomics_data/proteomics/raw/euc/swath', 'FDR')
euc_swath_reanalysed_FDR <- cbind(euc_swath_reanalysed_FDR[,1:7], euc_swath_reanalysed_FDR[,10:ncol(euc_swath_reanalysed_FDR)] %>% dplyr::select(contains('sample')))
names(euc_swath_reanalysed_FDR)[8:ncol(euc_swath_reanalysed_FDR)] <- do.call(rbind, str_split(names(euc_swath_reanalysed_FDR[,8:ncol(euc_swath_reanalysed_FDR)]), pattern = " "))[,2]
euc_swath_reanalysed_FDR <- subset(euc_swath_reanalysed_FDR, Decoy == 0)

euc_swath_reanalysed_FDR$FDR <- apply(euc_swath_reanalysed_FDR[8:ncol(euc_swath_reanalysed_FDR)],1, function(x) as.numeric(sum(as.numeric(x < 0.01)) > 2) )

nrow(euc_swath_reanalysed_FDR[euc_swath_reanalysed_FDR$FDR == 0,])

euc_swath_reanalysed_FDR <- subset(euc_swath_reanalysed_FDR, FDR == 1)

# ion areas constrained by FDR

ion_areas <- euc_swath_reanalysed[euc_swath_reanalysed$Peptide %in% euc_swath_reanalysed_FDR$Peptide,]
rm(euc_swath_reanalysed_FDR,euc_swath_reanalysed)

write_csv(ion_areas, 'data/proteomics_data/proteomics/derived/euc/D14_ion_areas_new_ion_library.csv')

ion_areas <- read_csv('data/proteomics_data/proteomics/derived/euc/D14_ion_areas_new_ion_library.csv')

ion_areas <- ion_areas[-grep('RRRRR', ion_areas$Protein),]

# get protein areas using top2top3 method

#protein_areas <-  ion_areas %>% 
#  group_by(Peptide, Protein) %>% 
#  summarise_at(vars(10:323), top2) %>%
#  group_by(Protein) %>%
#  summarise_at(vars(3:316), top3)

#write_csv(protein_areas, 'data/protein_areas.csv')

# get protein areas using top2top2 method

protein_areas <-  ion_areas %>% 
  group_by(Peptide, Protein) %>% 
  summarise_at(vars(10:ncol(ion_areas)), top2) 

protein_areas <- protein_areas %>%
  group_by(Protein) %>%
  summarise_at(vars(3:ncol(protein_areas)), top2avg)

write_csv(protein_areas, 'data/proteomics_data/proteomics/derived/euc/protein_areas_top2top2.csv')




























