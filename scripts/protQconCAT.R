# QConCAT calculations for 8-sample proteaceae subset

require(readxl)
require(dplyr)
require(stringr)
require(readr)
require(tidyr)

source('scripts/functions.R')

# load ion areas, rename columns

protQconCAT_swath_reanalysed <- read_excel('data/proteomics_data/proteomics/raw/proteaceae/QconCAT test EDITED SWATH library 8 samples SVS 170610.xlsx', 'Area - ions')
protQconCAT_swath_reanalysed <- cbind(protQconCAT_swath_reanalysed[,1:9], protQconCAT_swath_reanalysed[,10:ncol(protQconCAT_swath_reanalysed)] %>% select(contains('sample')))
names(protQconCAT_swath_reanalysed)[10:ncol(protQconCAT_swath_reanalysed)] <- do.call(rbind, str_split(names(protQconCAT_swath_reanalysed[,10:ncol(protQconCAT_swath_reanalysed)]), pattern = " "))[,1]

protQconCAT_swath_reanalysed <- protQconCAT_swath_reanalysed %>% 
                                subset(., select=which(!duplicated(names(.)))) %>% # remove duplicated columns (e.g KG022) 
                                filter(`Ion Type` != 'b') # take out the b ions

# load FDR data, rename columns, remove decoys and discard proteins below identification probability threshold

protQconCAT_swath_reanalysed_FDR <- read_excel('data/proteomics_data/proteomics/raw/proteaceae/QconCAT test EDITED SWATH library 8 samples SVS 170610.xlsx', 'FDR')
protQconCAT_swath_reanalysed_FDR <- cbind(protQconCAT_swath_reanalysed_FDR[,1:7], protQconCAT_swath_reanalysed_FDR[,10:ncol(protQconCAT_swath_reanalysed_FDR)] %>% dplyr::select(contains('sample')))
names(protQconCAT_swath_reanalysed_FDR)[8:ncol(protQconCAT_swath_reanalysed_FDR)] <- do.call(rbind, str_split(names(protQconCAT_swath_reanalysed_FDR[,8:ncol(protQconCAT_swath_reanalysed_FDR)]), pattern = " "))[,2]
protQconCAT_swath_reanalysed_FDR <- subset(protQconCAT_swath_reanalysed_FDR, Decoy == 0)

protQconCAT_swath_reanalysed_FDR$FDR <- apply(protQconCAT_swath_reanalysed_FDR[8:ncol(protQconCAT_swath_reanalysed_FDR)],1, function(x) as.numeric(sum(as.numeric(x < 0.01)) > 2) )

nrow(protQconCAT_swath_reanalysed_FDR[protQconCAT_swath_reanalysed_FDR$FDR == 0,])

protQconCAT_swath_reanalysed_FDR <- subset(protQconCAT_swath_reanalysed_FDR, FDR == 1)

# ion areas constrained by FDR

ion_areas <- protQconCAT_swath_reanalysed[protQconCAT_swath_reanalysed$Peptide %in% protQconCAT_swath_reanalysed_FDR$Peptide,]
rm(protQconCAT_swath_reanalysed_FDR,protQconCAT_swath_reanalysed)

ion_areas <- ion_areas[-grep('RRRRR', ion_areas$Protein),]

write_csv(ion_areas, 'ignore/protQconCAT_ion_areas.csv')


#need to use the same y ions for the QconCAT peptides for every sample
#which y ions should they be?
#  for Lys labelled peptides the QconCAT peptide fragment MZ be close to x + (8/Fragment charge)
#  for Arg labelled peptides the QconCAT peptide fragment MZ be close to x + (10/Fragment charge)

#so, top2top should procede as usual, with the exception of proteins which have associated labelled peptides
#  make a df with these proteins excised, calculate top2top2 for non-labelled proteins
#    make a list of peptides with [+08] or [+10], then strip these out and filter ion_areas by peptide names

qconcat <- ion_areas[grep('(8|10)', ion_areas$Peptide),]
qconcat$Peptide_clean <- gsub('\\[\\+[0-9][0-9]\\]', '', qconcat$Peptide)

# create separate df's for proteins with and without qconcat peptides

other_ion_areas <- ion_areas[grep(paste(unique(qconcat$Peptide_clean), collapse='|'), ion_areas$Peptide, invert=TRUE),]

qconcat_ion_areas <- ion_areas[grep(paste(unique(qconcat$Peptide_clean), collapse='|'), ion_areas$Peptide),]

qconcat_ion_areas$ion_means <- rowMeans(qconcat_ion_areas[,10:ncol(qconcat_ion_areas)]) # get rowMeans across samples for each ion 

# strip out labels from peptide names so that we get averaged positions for both labelled and unlabelled ions
qconcat_ion_areas$Peptide_clean <- gsub('\\[\\+[0-9][0-9]\\]', '', qconcat_ion_areas$Peptide) 

# find ions associated with top2 ion_means, by peptide
qconcat_top2_ions <- dplyr::group_by(qconcat_ion_areas, Protein, Peptide) %>% top_n(2, ion_means)

# are the top 2 ions the same for qconcat peptides?

check_ions <- function(df) {

  df <- subset(qconcat_top2_ions, Protein == "AT1G42970.1")
  
  df <- subset(df, Peptide_clean == "VIITAPAK")
  
  df <- dplyr::group_by(df, Peptide) %>% dplyr::mutate(check = sum(`Fragment MZ`))
  
  label <- as.numeric(max(gsub('[^0-9]', "", df$Peptide)))
  
  max(df$check) - mean(df$check) < (label + 0.1) && max(df$check) - mean(df$check) > (label - 0.1)

}

x <- dplyr::group_by(qconcat_top2_ions, Protein, Peptide_clean) %>% dplyr::do(check_ions(.))

x <- ddply(qconcat_top2_ions, .(Protein, Peptide_clean), check_ions)

qconcat_ion_areas <- qconcat_ion_areas[grep(paste(unique(qconcat_top2_ions$`Fragment MZ`), collapse='|'), qconcat_ion_areas$`Fragment MZ`),]
qconcat_ion_areas$Peptide_clean <- NULL
qconcat_ion_areas$ion_means <- NULL

qconcat_ion_areas_labelled <- qconcat_ion_areas[grep('(8|10)', qconcat_ion_areas$Peptide),]
qconcat_ion_areas_unlabelled  <- qconcat_ion_areas[grep('(8|10)', qconcat_ion_areas$Peptide, invert=TRUE),]
rm(qconcat_ion_areas)

bla <- rbind(qconcat_ion_areas_labelled,qconcat_ion_areas_unlabelled)


#  then make a df with top2top2 done for non-qconcat, qconcat-labelled and qconcat-unlabelled proteins separately, 
#  where the ions used to calculate the area of qconcat peptides (labelled and unlabelled) are always the same

qconcat_protein_areas_labelled <- qconcat_ion_areas_labelled %>% 
  group_by(Peptide, Protein) %>% 
  summarise_at(vars(10:ncol(qconcat_ion_areas_labelled)), top2) 

qconcat_protein_areas_labelled <- qconcat_protein_areas_labelled %>%
  group_by(Protein) %>%
  summarise_at(vars(3:ncol(qconcat_protein_areas_labelled)), top2avg)


qconcat_protein_areas_unlabelled <- qconcat_ion_areas_unlabelled %>% 
  group_by(Peptide, Protein) %>% 
  summarise_at(vars(10:ncol(qconcat_ion_areas_unlabelled)), top2) 

qconcat_protein_areas_unlabelled <- qconcat_protein_areas_unlabelled %>%
  group_by(Protein) %>%
  summarise_at(vars(3:ncol(qconcat_protein_areas_unlabelled)), top2avg)


other_protein_areas <- other_ion_areas %>% 
  group_by(Peptide, Protein) %>% 
  summarise_at(vars(10:ncol(other_ion_areas)), top2) 

other_protein_areas <- other_protein_areas %>%
  group_by(Protein) %>%
  summarise_at(vars(3:ncol(other_protein_areas)), top2avg)

qconcat_prot_protein_areas <- rbind(qconcat_protein_areas_labelled, qconcat_protein_areas_unlabelled, other_protein_areas)
