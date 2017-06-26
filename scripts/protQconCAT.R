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

qconcat_ion_areas$ion_means <- rowMeans(qconcat_ion_areas[,10:ncol(qconcat_ion_areas)], na.rm=TRUE) # get rowMeans across samples for each ion 

# strip out labels from peptide names so that we get averaged positions for both labelled and unlabelled ions
qconcat_ion_areas$Peptide_clean <- gsub('\\[\\+[0-9][0-9]\\]', '', qconcat_ion_areas$Peptide) 

# find ions associated with top2 ion_means, by peptide
qconcat_top2_ions <- dplyr::group_by(qconcat_ion_areas, Protein, Peptide) %>% top_n(2, ion_means)

# are the top 2 ions the same for qconcat peptides?


# subset by Peptide_clean, then sort by fragment mz / fragment charge, 
# then take the top 2 pairs where the difference between them is the label mass

pairwise_avg_top2 <- function(bla) {
  
  bla <- subset(qconcat_ion_areas, Peptide_clean == "VIITAPAK")
  
  bla <- bla[order(bla$`Fragment MZ`),]
  
# if(length(unique(bla$Peptide_clean)) != length(unique(bla$Peptide))) { # i.e. if the df doesn't just contain either labelled or unlabelled peptides
    
    # could modify labelled fragment MZ's so that they're the same as unlabelled fragent MZ's... remove label amount and then round everything 
    
    #label <- as.numeric(max(gsub('[^0-9]', "", bla$Peptide)))
    
    bla$label <- as.numeric(gsub('[^0-9]', "", bla$Peptide))
    
    if(any(bla$label %in% '810')) {
      bla[bla$label %in% '810',]$label <- 18
    }
    
    bla[is.na(bla$label),]$label <- 0
    
    bla$Fragment_MZ_charge <- bla$`Fragment MZ` * bla$`Fragment Charge`
    
    bla$Fragment_MZ_charge_adj <- bla$Fragment_MZ_charge
    
    bla$Fragment_MZ_charge_adj <- bla$Fragment_MZ_charge_adj - bla$label
    
    bla$Fragment_MZ_charge_adj <- round(bla$Fragment_MZ_charge_adj, 0)
    
    # so now i can sort or run operations by fragment mz / charge
    
    # now need to subset dataset to only keep paired ions??
    
    # na.omit it first to get rid of rare ions
    
    bla <- na.omit(bla)
    
    # how do i get top2 pairs? get a pairwise average ion area?
    
    bla <- dplyr::group_by(bla, Fragment_MZ_charge_adj) %>% dplyr::summarise(pairwise_avg_ion_mean = mean(ion_means),
                                                                             labelled_unlabelled = length(ion_means)) %>%
      full_join(bla) %>% 
      filter(labelled_unlabelled == 2) %>%
   #   filter(ions_per_peptide == 2) %>%
      dplyr::group_by(Peptide) %>% top_n(2, pairwise_avg_ion_mean) %>%
      dplyr::group_by(Peptide) %>% dplyr::mutate(ions_per_peptide = length(Peptide))
      
    
    
    #browser()
    
    return(bla)
    
  }
  
}


x <- ddply(qconcat_ion_areas, .(Peptide_clean), pairwise_avg_top2)

# problem: for some peptides, there is only 1 ion that is present in all samples


VIITAPAK <- na.omit(ion_areas[grep('VIITAPAK', ion_areas$Peptide),])
VIITAPAK_x <- na.omit(x[grep('VIITAPAK', x$Peptide),])



# still not picking up ions that are doubly labelled
length(unique(x$Peptide_clean))
length(unique(qconcat_ion_areas$Peptide_clean))





check_ions <- function(df) {

  #df <- subset(qconcat_top2_ions, Protein == "ATCG00340.1")
  
 # df <- subset(df, Peptide_clean == "DKPVALSIVQAR")
  
  label <- as.numeric(max(gsub('[^0-9]', "", df$Peptide)))
  
  if(label %in% 810) {
   label <- 18 
  }
  
  df$check <- NA
  
  ion1_diff <- (sort(df$`Fragment MZ`)[2] - sort(df$`Fragment MZ`)[1]) * df[order(df$`Fragment MZ`),]$`Fragment Charge`[1]
  df[order(df$`Fragment MZ`),]$check[1:2] <- ion1_diff > (label - 0.1) && ion1_diff < (label + 0.1)
  
  if(length(unique(df$`Fragment MZ`)) > 2) {
    ion2_diff <- (sort(df$`Fragment MZ`)[4] - sort(df$`Fragment MZ`)[3]) * df[order(df$`Fragment MZ`),]$`Fragment Charge`[3]
    df[order(df$`Fragment MZ`),]$check[3:4] <- ion2_diff > (label - 0.1) && ion2_diff < (label + 0.1)
  }
}

x <- ddply(x, .(Protein, Peptide_clean), check_ions)



x <- ddply(qconcat_top2_ions, .(Protein, Peptide_clean), check_ions)
x_false <- subset(x, V1 == FALSE)

x1 <- qconcat_top2_ions[qconcat_top2_ions$Peptide_clean %in% x_false$Peptide_clean,]
x2 <- qconcat_ion_areas[qconcat_ion_areas$Peptide_clean %in% x_false$Peptide_clean,]


# appears that some of the top2 ions are not present in both sets. Need away to find ions that are present in both sets.

# could do top3, run the checks, subset out the FALSE's and then do top2?



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
