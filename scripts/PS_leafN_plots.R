source('scripts/transformations.R')

## plots showing relationships between Asat/Amax and various traits / variables; leaf N as well

# Amax (A1500)

blah <- subset(climate_locs, Cond > 0.05)

blah <- blah[!is.na(blah$photo_max),]
blah$leaf_age <- as.factor(blah$leaf_age)

levels(blah$leaf_age)<-c('new', 'mid', 'old')

boxplot(photo_max ~ leaf_age, blah, ylab = 'Photosynthetic rate (A1500) per leaf area', xlab = 'leaf age')

blax <- subset(protein_climate_D14, bin_arch_name == "total_protein")
plot(photo_max ~ sum, blax, ylab = 'Photosynthetic rate (A1500) per leaf area', xlab = 'total protein per leaf area (mg/mm2)')

plot(photo_max ~ tavg, climate_locs, xlab = 'tavg (oC)', ylab = 'Photosynthetic rate (A1500) per leaf area')
abline(lm(photo_max ~ tavg, climate_locs))
summary(lm(photo_max ~ tavg, climate_locs))

plot(photo_max ~ prec, climate_locs, xlab = 'prec (mm)', ylab = 'Photosynthetic rate (A1500) per leaf area')
abline(lm(photo_max ~ prec, climate_locs))

plot(photo_max ~ LMA_g_per_m2, climate_locs, xlab = "LMA (g/m2)", ylab = 'Photosynthetic rate (A1500) per leaf area')
summary(lm(photo_max ~ LMA_g_per_m2, climate_locs))


# leaf_N

plot(total_protein ~ N, merge(protein_D14, climate_locs, by = 'sample'))
abline(lm(total_protein ~ N, merge(protein_D14, climate_locs, by = 'sample')))
summary(lm(total_protein ~ N, merge(protein_D14, climate_locs, by = 'sample')))

climate_locs$N_per_area <- climate_locs$N * climate_locs$LMA_g_per_m2

plot(total_protein ~ N_per_area, merge(protein_D14, climate_locs, by = 'sample'), xlab = "N_per_area (mg/m2)", ylab = "total_protein (mg/m2)")
abline(lm(total_protein ~ N_per_area, merge(protein_D14, climate_locs, by = 'sample')))
summary(lm(total_protein ~ N_per_area, merge(protein_D14, climate_locs, by = 'sample')))

plot(photo_max ~ N_per_area, climate_locs)
abline(lm(photo_max ~ N_per_area, climate_locs))
summary(lm(photo_max ~ N_per_area, climate_locs))

#plot(photo_amb ~ N_per_area, climate_locs)
#abline(lm(photo_amb ~ N_per_area, climate_locs))
#summary(lm(photo_amb ~ N_per_area, climate_locs))


# leaf N vs leaf age

climate_locs$leaf_age <- factor(climate_locs$leaf_age, levels = c('new', 'mid', 'old'))
plot(N ~ as.factor(leaf_age), climate_locs)

# protein assay protein amounts

protein_assays <- read_csv('data/protein_assay.csv')
protein_assays <- protein_assays[protein_assays$sample %in% climate_locs$sample,]
climate_locs <- merge(protein_assays, climate_locs, all.y=TRUE, by = 'sample')

bla <- merge(protein_D14, climate_locs, by = 'sample')

plot(bla$total_protein / bla$LMA_g_per_m2 ~ bla$N, 
     ylab = 'mass spec protein (mg/g)', xlab = "N %") 

plot(bla$assay_protein_per_area_mg_per_g ~ bla$N, 
     ylab = 'assay protein (mg/g)', xlab = "N %") 

plot(total_protein ~ N_per_area, bla, ylab = 'mass spec protein (mg/m2)', xlab = "N per area (mg/m2)") 
plot(assay_protein_per_area_mg_per_m2 ~ N_per_area, bla, ylab = 'assay protein (mg/m2)', xlab = "N per area (mg/m2)") 

plot(total_protein ~ assay_protein_per_area_mg_per_m2, bla)

# is there a lower cutoff of Cond which renders the distribution of photomax normal?

plot(climate_locs$photo_max ~ climate_locs$Ci)
plot(climate_locs$photo_max ~ climate_locs$Cond)

climate_locs$photo_max <- as.numeric(climate_locs$photo_max)
blah <- subset(climate_locs, Cond > 0.35) 
shapiro.test(blah$photo_max)# yes but that's silly
plot(blah$photo_max ~ blah$Cond)

# a better cutoff is 0.05, which was used by Steve (actually he used 0.04) in the field as a 'discard sample' cutoff
blah <- subset(climate_locs, Cond > 0.05)
plot(blah$photo_max ~ blah$Cond)

