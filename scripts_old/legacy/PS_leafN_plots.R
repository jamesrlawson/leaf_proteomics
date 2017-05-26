source('scripts/transformations.R')

## plots showing relationships between Asat/Amax and various traits / variables; leaf N as well

# Amax (A1500)

blah <- subset(subset(climate_locs, Cond > 0.05), Cond > 0.05)

blah <- blah[!is.na(blah$photo_max),]
blah$leaf_age <- as.factor(blah$leaf_age)

levels(blah$leaf_age)<-c('new', 'mid', 'old')

boxplot(photo_max ~ leaf_age, blah, ylab = 'Photosynthetic rate (A1500) per leaf area', xlab = 'leaf age')

blax <- subset(protein_climate_D14, bin_arch_name == "total_protein")
plot(photo_max ~ sum, blax, ylab = 'Photosynthetic rate (A1500) per leaf area', xlab = 'total protein per leaf area (mg/mm2)')

plot(photo_max ~ tavg, subset(climate_locs, Cond > 0.05), xlab = 'tavg (oC)', ylab = 'Photosynthetic rate (A1500) per leaf area')
abline(lm(photo_max ~ tavg, subset(climate_locs, Cond > 0.05)))
summary(lm(photo_max ~ tavg, subset(climate_locs, Cond > 0.05)))

plot(photo_max ~ prec, subset(climate_locs, Cond > 0.05), xlab = 'prec (mm)', ylab = 'Photosynthetic rate (A1500) per leaf area')
abline(lm(photo_max ~ prec, subset(climate_locs, Cond > 0.05)))

plot(photo_max ~ LMA_g_per_m2, subset(climate_locs, Cond > 0.05), xlab = "LMA (g/m2)", ylab = 'Photosynthetic rate (A1500) per leaf area')
summary(lm(photo_max ~ LMA_g_per_m2, subset(climate_locs, Cond > 0.05)))


# leaf_N

plot(total_protein ~ N, merge(protein_D14, subset(climate_locs, Cond > 0.05), by = 'sample'))
abline(lm(total_protein ~ N, merge(protein_D14, subset(climate_locs, Cond > 0.05), by = 'sample')))
summary(lm(total_protein ~ N, merge(protein_D14, subset(climate_locs, Cond > 0.05), by = 'sample')))

subset(climate_locs, Cond > 0.05)$N_per_area <- subset(climate_locs, Cond > 0.05)$N * subset(climate_locs, Cond > 0.05)$LMA_g_per_m2

plot(total_protein ~ N_per_area, merge(protein_D14, subset(climate_locs, Cond > 0.05), by = 'sample'), xlab = "N_per_area (mg/m2)", ylab = "total_protein (mg/m2)")
abline(lm(total_protein ~ N_per_area, merge(protein_D14, subset(climate_locs, Cond > 0.05), by = 'sample')))
summary(lm(total_protein ~ N_per_area, merge(protein_D14, subset(climate_locs, Cond > 0.05), by = 'sample')))

plot(photo_max ~ N_per_area, subset(climate_locs, Cond > 0.05))
abline(lm(photo_max ~ N_per_area, subset(climate_locs, Cond > 0.05)))
summary(lm(photo_max ~ N_per_area, subset(climate_locs, Cond > 0.05)))

plot(photo_max ~ N, subset(climate_locs, Cond > 0.05))
abline(lm(photo_max ~ N, subset(climate_locs, Cond > 0.05)))
summary(lm(photo_max ~ N, subset(climate_locs, Cond > 0.05)))


# leaf N vs leaf age

subset(climate_locs, Cond > 0.05)$leaf_age <- factor(subset(climate_locs, Cond > 0.05)$leaf_age, levels = c('new', 'mid', 'old'))
plot(N ~ as.factor(leaf_age), subset(climate_locs, Cond > 0.05))




