## plots showing relationships between Asat/Amax and various traits / variables; leaf N as well

# Asat (A400)

blah <- merge(climate_locs, leaf_age, by = 'sample')

blah <- blah[!is.na(blah$photo_amb),]
blah$leaf_age <- as.factor(blah$leaf_age)

levels(blah$leaf_age)<-c('new', 'mid', 'old')

boxplot(photo_amb ~ leaf_age, blah, ylab = 'Photosynthetic rate (A400) per leaf area', xlab = 'leaf age')

blax <- subset(protein_climate_D14, bin_arch_name == "total_protein")
plot(photo_amb ~ sum, blax, xlab = 'Photosynthetic rate (A1500) per leaf area', ylab = 'total protein per leaf area (mg/mm2)')


plot(photo_amb ~ tavg, climate_locs, xlab = 'tavg (oC)', ylab = 'Photosynthetic rate (A400) per leaf area')
abline(lm(photo_amb ~ tavg, climate_locs))
summary(lm(photo_amb ~ tavg, climate_locs))

plot(photo_amb ~ prec, climate_locs, xlab = 'prec (mm)', ylab = 'Photosynthetic rate (A400) per leaf area')
abline(lm(photo_amb ~ prec, climate_locs))
summary(lm(photo_amb ~ prec, climate_locs))

plot(photo_amb ~ LMA_g_per_m2, climate_locs, xlab = "LMA (g/m2)", ylab = 'Photosynthetic rate (A400) per leaf area')
abline(lm(photo_amb ~ LMA_g_per_m2, climate_locs))
abline(lm(photo_amb ~ LMA_g_per_m2, climate_locs))
summary(lm(photo_amb ~ LMA_g_per_m2, climate_locs))

plot(photo_amb ~ N, climate_locs)
abline(lm(photo_amb ~ N, climate_locs))
summary(lm(photo_amb ~ N, climate_locs))



# Amax (A1500)

blah <- merge(climate_locs, leaf_age, by = 'sample')

blah <- blah[!is.na(blah$photo_max),]
blah$leaf_age <- as.factor(blah$leaf_age)

levels(blah$leaf_age)<-c('new', 'mid', 'old')

boxplot(photo_max ~ leaf_age, blah, ylab = 'Photosynthetic rate (A1500) per leaf area', xlab = 'leaf age')

blax <- subset(protein_climate_D14, bin_arch_name == "total_protein")
plot(photo_max ~ sum, blax, xlab = 'Photosynthetic rate (A1500) per leaf area', ylab = 'total protein per leaf area (mg/mm2)')


plot(photo_max ~ tavg, climate_locs, xlab = 'tavg (oC)', ylab = 'Photosynthetic rate (A1500) per leaf area')
abline(lm(photo_max ~ tavg, climate_locs))
summary(lm(photo_max ~ tavg, climate_locs))

plot(photo_max ~ prec, climate_locs, xlab = 'prec (mm)', ylab = 'Photosynthetic rate (A1500) per leaf area')
abline(lm(photo_max ~ prec, climate_locs))
summary(lm(photo_max ~ prec, climate_locs))

plot(photo_max ~ LMA_g_per_m2, climate_locs, xlab = "LMA (g/m2)", ylab = 'Photosynthetic rate (A1500) per leaf area')
abline(lm(photo_max ~ LMA_g_per_m2, climate_locs))
abline(lm(photo_max ~ LMA_g_per_m2, climate_locs))
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

plot(photo_amb ~ N_per_area, climate_locs)
abline(lm(photo_amb ~ N_per_area, climate_locs))
summary(lm(photo_amb ~ N_per_area, climate_locs))


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

plot(bla$photo_max ~ bla$Ci)
plot(bla$photo_max ~ bla$Cond)







