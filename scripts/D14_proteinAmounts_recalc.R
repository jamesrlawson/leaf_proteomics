# D14 PROTEIN AMOUNTS RECALCULATED USING SIGNAL AREA FRACTIONS

require(readr)
require(dplyr)
require(stringr)
require(tidyr)

FDR <- read_csv('data/D14_FDR.csv')

FDR <- FDR[FDR$Decoy == FALSE,]

FDR$FDR <- apply(FDR[8:ncol(FDR)],1, function(x) as.numeric(any(x < 0.01)) )

FDR <- subset(FDR, FDR == 1)

ion_areas <- read_csv('data/D14_ion_areas_all_samples.csv')

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

protein_assay <- read_csv('data/protein_assay.csv')

protein_amounts <- gather(protein_fractions, 'sample', 'protein_fraction', cols)

protein_amounts <- merge(protein_amounts, protein_assay, by = "sample")       

protein_amounts$protein_per_dry_weight <- protein_amounts$protein_fraction * protein_amounts$assay_protein_per_area_mg_per_g

protein_amounts$protein_per_leaf_area <- protein_amounts$protein_fraction * protein_amounts$assay_protein_per_area_mg_per_m2

protein_per_dw <- spread(protein_amounts[,c('Protein', 'sample', 'protein_per_dry_weight')], key = 'sample', value= 'protein_per_dry_weight')                
protein_per_leaf_area <- spread(protein_amounts[,c('Protein', 'sample', 'protein_per_leaf_area')], key = 'sample', value= 'protein_per_leaf_area')                

write_csv(protein_per_dw, 'output/protein_amounts_by_signal_fraction_perDW_D14.csv')
write_csv(protein_per_leaf_area, 'output/protein_amounts_by_signal_fraction_perArea_D14.csv')















## ovalb correction factors

fraction_oval <- subset(protein_fractions, Protein == "sp|OVAL_CHICK")
x <- fraction_oval[1,cols]

blah <- protein_assay %>% select(sample, percent_ovalb_of_total_protein) %>% mutate(percent_ovalb_of_total_protein = percent_ovalb_of_total_protein / 100) %>%
  arrange(sample)

y <- data.frame(cbind(fraction_oval_MS = t(x), sample = names(x)))

j <- merge(blah, y, by = 'sample')

names(j)[3] <- 'fraction_oval_MS'
j$fraction_oval_MS <- as.numeric(as.character(j$fraction_oval_MS))

plot(j$fraction_oval_MS ~ j$percent_ovalb_of_total_protein)
model <- lm(j$fraction_oval_MS ~ j$percent_ovalb_of_total_protein)
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

ion_areas <- read_csv('data/D14_ion_areas_all_samples.csv')

top2 <- function(x) {
  sum(sort(x, decreasing = TRUE)[1:2], na.rm=TRUE)
}

ion_areas_ov <- subset(ion_areas, Protein == "sp|OVAL_CHICK" & Peptide %in% ov_peps)

peptide_areas_ov <-  ion_areas_ov %>% 
  group_by(Peptide) %>% 
  summarise_at(vars(matches("sample")), top2) %>%
  mutate(protein_MW = protein_MW[protein_MW$Protein == "sp|OVAL_CHICK",]$MW)

peptide_areas_ov[,cols] <- t(t(peptide_areas_ov[,cols])/colSums(protein_areas[,cols]))

for(i in cols) {
  
  name_split <- unlist(str_split(names(peptide_areas_ov)[i], " "))[1]
  names(peptide_areas_ov)[i] <- name_split
  
}

# produce plots and output 'x' as 

my.list <- vector("list", nrow(peptide_areas_ov))

for(i in 1:nrow(peptide_areas_ov)) {
  
  # browser()
  
  fraction_oval <- peptide_areas_ov[i,]
  x <- fraction_oval[1,cols]
  
  blah <- protein_assay %>% select(sample, percent_ovalb_of_total_protein) %>% mutate(percent_ovalb_of_total_protein = percent_ovalb_of_total_protein / 100) %>%
    arrange(sample)
  
  y <- data.frame(cbind(fraction_oval_MS = t(x), sample = names(x)))
  
  j <- merge(blah, y, by = 'sample')
  
  names(j)[3] <- 'fraction_oval_MS'
  j$fraction_oval_MS <- as.numeric(as.character(j$fraction_oval_MS))
  
  plot(j$fraction_oval_MS ~ j$percent_ovalb_of_total_protein, main = peptide_areas_ov$Peptide[i], ylab = 'Fraction of ovalbumin (MS)', xlab = 'Ovalbumin % of assay protein')
  model <- lm(j$fraction_oval_MS ~ j$percent_ovalb_of_total_protein)
  abline(model)
  
  pep <- peptide_areas_ov$Peptide[i]
  
  model.stats <-  data.frame(cbind(pep,
                                   round(summary(model)$r.squared,3),
                                   round(summary(model)$coefficients[,4][2],2),
                                   model$coefficients[1],
                                   model$coefficients[2]))
  
  names(model.stats) <- c('submodel','R2','pval', 'intercept', 'slope') 
  rownames(model.stats) <- NULL
  
  my.list[[i]] <- model.stats
  
 
}

x <- data.frame()

output <- rbind(x, do.call(rbind, my.list))
output <- output %>% mutate(R2 = as.numeric(as.character(R2))) %>% arrange(desc(R2))


# avg of ovalb peptides with top 4 R2 # 

ovalb_top4 <- c(#'DEDTQAMPFR',
                'GGLEPINFQTAADQAR',
                'LTEWTSSNVMEER',
                'YPILPEYLQC[Pye]VK')

ovalb_top4avg <- peptide_areas_ov %>% filter(Peptide %in% ovalb_top4) %>%
                                      summarise_at(cols, mean)
                            
blah <- protein_assay %>% select(sample, percent_ovalb_of_total_protein) %>% mutate(percent_ovalb_of_total_protein = percent_ovalb_of_total_protein / 100) %>%
  arrange(sample)

y <- data.frame(cbind(fraction_oval_MS = t(ovalb_top4avg), sample = names(ovalb_top4avg)))

j <- merge(blah, y, by = 'sample')

names(j)[3] <- 'fraction_oval_MS'
j$fraction_oval_MS <- as.numeric(as.character(j$fraction_oval_MS))

plot(j$fraction_oval_MS ~ j$percent_ovalb_of_total_protein, main = "avg of 4 ovalb peptides with strongest model fits", ylab = 'Fraction of ovalbumin (MS)', xlab = 'Ovalbumin % of assay protein')
model <- lm(j$fraction_oval_MS ~ j$percent_ovalb_of_total_protein)
abline(model)
summary(model)



## some bias here (and spurious values of 'fraction_oval_MS' may be due to peptides (?) not being identified & counted 
## because they have no matches in the ion library
## how well are individual samples matched to the ion library?
## a good diagnostic is: samplewise, ion area used in top2/top3 calculations vs total ion area
## wait I don't understand this properly
## but here's the code..

sumtop3 <- function(x) {
  sum(sort(x, decreasing = TRUE)[1:3], na.rm=TRUE)
}

sum_ion_areas_in_top2top3   <-  ion_areas %>% 
                                group_by(Peptide, Protein) %>% 
                                summarise_at(vars(matches("sample")), top2) %>%
                                group_by(Protein) %>%
                                summarise_at(vars(matches("sample")), sumtop3) %>%
                                summarise_at(vars(matches("sample")), sum)

sum_ion_areas_all <- ion_areas %>% 
                     summarise_at(vars(matches("sample")), sum, na.rm=TRUE)


x <-  t(sum_ion_areas_in_top2top3) / t(sum_ion_areas_all)
rownames(x) <- names(ion_areas)[10:323]
hist(x, main = 'sum of ion areas used in top2top3 / sum of all ion areas')
