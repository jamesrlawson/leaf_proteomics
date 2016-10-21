# ASPINWALL PROTEIN AMOUNTS RECALCULATED USING SIGNAL AREA FRACTIONS

require(readr)
require(dplyr)
require(stringr)
require(tidyr)

FDR <- read.csv('data/aspinwall/aspinwall_FDR.csv', header=TRUE, stringsAsFactors = FALSE)

FDR <- FDR[FDR$Decoy == FALSE,]

FDR$FDR <- apply(FDR[8:ncol(FDR)],1, function(x) as.numeric(any(x < 0.01)) )

FDR <- subset(FDR, FDR == 1)

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

cols <- 2:(ncol(protein_fractions)-1) # indices denoting sample columns

protein_fractions[,cols] <- protein_fractions[,cols] * protein_fractions$MW

protein_fractions[,cols] <- t(t(protein_fractions[,cols])/colSums(protein_fractions[,cols]))
                           

for(i in cols) {
  
  name_split <- unlist(str_split(names(protein_fractions)[i], " "))[1]
  names(protein_fractions)[i] <- name_split
  
}

asp_spreadsheet <- read.csv('data/aspinwall/aspinwall_spread.csv', header=TRUE, stringsAsFactors = TRUE)
asp_spreadsheet$ID <- gsub("[[:punct:]]", "", asp_spreadsheet$`sample.id..day.species.abbreviation.pot.`)

protein_amounts <- gather(protein_fractions, 'ID', 'protein_fraction', cols)
                
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
names(ovalb_quant)[cols] <- names(protein_areas)[cols]

for(i in cols) {
  
  name_split <- unlist(str_split(names(ovalb_quant)[i], " "))[1]
  names(ovalb_quant)[i] <- name_split
  
}

ovalb_amounts <- gather(ovalb_quant, 'ID', 'protein_amount', cols)
ovalb_amounts <- merge(ovalb_amounts, asp_spreadsheet, by = "ID")   

ovalb_amounts$treatment <- paste(ovalb_amounts$day, ovalb_amounts$Temperature.treatment, ovalb_amounts$species.abbreviation)

ovalb_CV <- ovalb_amounts %>% group_by(treatment, ID) %>%
                      summarise(sum_protein_per_area = sum(protein_amount, na.rm=TRUE)) %>%
                      group_by(treatment) %>%
                      summarise(CV_protein_per_area_ovalb = CV(sum_protein_per_area), N_ovalb = length(sum_protein_per_area))


method_comparison <- merge(fractions_CV, ovalb_CV)

write_csv(method_comparison, 'output/protein_amount_calculation_method_comparison.csv')

mean(method_comparison$CV_protein_per_area_fractions)
mean(method_comparison$CV_protein_per_area_ovalb)

## ovalb correction factors

fraction_oval <- subset(protein_fractions, Protein == "sp|OVAL_CHICK")
x <- fraction_oval[1,cols]


blah <- asp_spreadsheet %>% select(sample.id..day.species.abbreviation.pot., X.Ovalbumin.of.Total.Protein) %>% mutate(X.Ovalbumin.of.Total.Protein = X.Ovalbumin.of.Total.Protein / 100) %>%
                            arrange(sample.id..day.species.abbreviation.pot.)
                        
blah$sample <- gsub("[[:punct:]]", "", blah$`sample.id..day.species.abbreviation.pot.`)


y <- data.frame(cbind(fraction_oval_MS = t(x), sample = names(x)))

j <- merge(blah, y, by = 'sample')

names(j)[4] <- 'fraction_oval_MS'
j$fraction_oval_MS <- as.numeric(as.character(j$fraction_oval_MS))

plot(j$fraction_oval_MS ~ j$X.Ovalbumin.of.Total.Protein)
model <- lm(j$fraction_oval_MS ~ j$X.Ovalbumin.of.Total.Protein)
abline(model)
summary(model)

# ovalb MS amounts calculated from individual peptides 

ov_peps <-  c('GGLEPINFQTAADQAR',
              'HIATNAVLFFGR',
              'ISQAVHAAHAEINEAGR',
              'YPILPEYLQC[Pye]VK',
              'LTEWTSSNVMEER',
              'DILNQITKPNDVYSFSLASR',
              'DEDTQAMPFR')

ion_areas <- read_csv('data/aspinwall/aspinwall_ion_areas.csv')

top2 <- function(x) {
  sum(sort(x, decreasing = TRUE)[1:2], na.rm=TRUE)
}

ion_areas_ov <- subset(ion_areas, Protein == "sp|OVAL_CHICK" & Peptide %in% ov_peps)

peptide_areas_ov <-  ion_areas_ov %>% 
  group_by(Peptide) %>% 
  summarise_at(vars(matches("sample")), top2) %>%
  mutate(protein_MW = protein_MW[protein_MW$Protein == "sp|OVAL_CHICK",]$MW)

#peptide_areas_ov[,cols] <- peptide_areas_ov[,cols] * peptide_areas_ov$protein_MW
peptide_areas_ov[,cols] <- t(t(peptide_areas_ov[,cols])/colSums(protein_areas[,cols]))

for(i in cols) {
  
  name_split <- unlist(str_split(names(peptide_areas_ov)[i], " "))[1]
  names(peptide_areas_ov)[i] <- name_split
  
}

for(i in 1:nrow(peptide_areas_ov)) {
  
 # browser()
 
  fraction_oval <- peptide_areas_ov[i,]
  x <- fraction_oval[1,cols]
  
  blah <- asp_spreadsheet %>% select(sample.id..day.species.abbreviation.pot., X.Ovalbumin.of.Total.Protein) %>% mutate(X.Ovalbumin.of.Total.Protein = X.Ovalbumin.of.Total.Protein / 100) %>%
    arrange(sample.id..day.species.abbreviation.pot.)
  
  blah$sample <- gsub("[[:punct:]]", "", blah$`sample.id..day.species.abbreviation.pot.`)
  
  
  y <- data.frame(cbind(t(x),names(x)))
  names(y) <- c('fraction_oval_MS', 'sample')
  y$fraction_oval_MS <- as.numeric(as.character(y$fraction_oval_MS))
  y$sample <- as.character(y$sample)
  
  j <- merge(blah, y, by = 'sample')
  
  names(j)[4] <- 'fraction_oval_MS'
  j$fraction_oval_MS <- as.numeric(as.character(j$fraction_oval_MS))
  
  plot(j$fraction_oval_MS ~ j$X.Ovalbumin.of.Total.Protein, main = peptide_areas_ov$Peptide[i], ylab = 'Fraction of ovalbumin (MS)', xlab = 'Ovalbumin % of assay protein')
  model <- lm(j$fraction_oval_MS ~ j$X.Ovalbumin.of.Total.Protein)
  abline(model)
  summary(model)
  
  
}

# myoglobin MS amounts calculated from 

ion_areas_myo <- subset(ion_areas, Protein == "sp|MYG_HORSE")

ion_areas_myo_ <-  ion_areas_myo %>% 
  group_by(Peptide) %>% 
  summarise_at(vars(matches("sample")), top3)

ion_areas_myo_[,cols] <-  t(t(ion_areas_myo_[,cols])/colSums(protein_areas[,cols]))

ion_areas_myo_$CV <- 2
ion_areas_myo_$mean <- 2

for(i in 1:nrow(ion_areas_myo_)) {
  
  x <- as.numeric(as.vector(ion_areas_myo_[i,cols]))
  ion_areas_myo_$mean[i] <- mean(x, na.rm=TRUE)
  ion_areas_myo_$CV[i] <- CV(x)
  
}

ion_areas_myo_top3 <- ion_areas_myo_[,c(1,98,99)]


# top2 / top 3


ion_areas_myo <- subset(ion_areas, Protein == "sp|MYG_HORSE")

ion_areas_myo_ <-  ion_areas_myo %>% 
  group_by(Peptide, Protein) %>% 
  summarise_at(vars(matches("sample")), top2) %>%
  ungroup() %>%
  summarise_at(vars(matches("sample")), top3)

ion_areas_myo_top2top3 <-  as.vector(as.numeric(ion_areas_myo_[1,]))/colSums(protein_areas[,cols])

CV(ion_areas_myo_top2top3)
mean(ion_areas_myo_top2top3)


# recreate plot with top2top3 values

