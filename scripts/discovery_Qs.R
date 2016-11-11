# Discovery Questions

# Q1: As leaves age (but while still alive and photosynthetically active), 
# does the amount of leaf N in refractory stress- and defense related proteins increase, 
# at the expense of photosynthetic proteins, and potentially account for the decline in photosynthetic capacity with age?

source('scripts/transformations.R')

protein_D14_age$stress_vs_photosynth <- protein_D14_age$stress / (protein_D14_age$Photorespiration + protein_D14_age$Light_reactions + protein_D14_age$Calvin_cycle)
protein_D14_age$leaf_age <- factor(protein_D14_age$leaf_age, levels = c("new", "mid", "old"))
plot(stress ~ leaf_age, protein_D14_age, ylab = "stress protein amount (mg/m2)")
plot((protein_D14_age$Photorespiration + protein_D14_age$Light_reactions + protein_D14_age$Calvin_cycle) ~ leaf_age, protein_D14_age, ylab = "photosynthesis protein amount (mg/m2)")
plot(stress_vs_photosynth ~ leaf_age, protein_D14_age, ylab = "stress / photosynthesis protein ratio")

# Q2: Do differences among species in fraction of leaf N resorbed 
# reflect differences in the proportion of N that is in refractory proteins towards the end of leaf life?

source('scripts/resorption.R')

ba <- subset(ba, leaf_age == 'old')

plot(ba$stress ~ ba$resorption_N, ylab = 'stress protein (proportion) [old leaves]', xlab = 'N resorption (senecent leaf N / old leaf N)')
abline(lm(ba$stress ~ ba$resorption_N))
summary(lm(ba$stress ~ ba$resorption_N))



# Q3: Do the Npool changes associated with leaf aging (as hypothesized under Q1) come about faster in short-lived leaves?

LMA_LWC <- read_csv('data/LMA_LWC.csv')
LMA_LWC <- LMA_LWC[LMA_LWC$sample %in% climate_locs$sample,]
LMA_LWC$LMA_g_per_m2  <- as.numeric(LMA_LWC$LMA_g_per_m2)
LMA_LWC$LWC_percent  <- as.numeric(LMA_LWC$LWC_percent)

protein_D14_age$leaf_age <- factor(protein_D14_age$leaf_age, levels = c("new", "mid", "old")

summary(lm(stress ~ leaf_age * LMA_g_per_m2, merge(protein_D14_age, LMA_LWC, by = 'sample')))
anova(lm(stress ~ leaf_age * LMA_g_per_m2, merge(protein_D14_age, LMA_LWC, by = 'sample')))

# Q4: Do differences among species in N resorption efficiency and proficiency 
# reflect different fractions of N placed in constitutive defenses from the outset of the leaf's life?

# Q5: Are ratios of stress- and defense-related proteins to photosynthesis related proteins higher at lower latitudes?

plot(stress_vs_photosynth ~ Latitude, merge(climate_locs, protein_D14_age), ylab = "stress / photosynthesis protein ratio")
abline(lm(stress_vs_photosynth ~ Latitude, merge(climate_locs, protein_D14_age)))
summary(lm(stress_vs_photosynth ~ Latitude, merge(climate_locs, protein_D14_age)))
  
  # no, but...
plot(stress ~ Latitude, merge(climate_locs, protein_D14_age), ylab = "stress protein amount (mg/m2)")
abline(lm(stress ~ Latitude, merge(climate_locs, protein_D14_age)))
summary(lm(stress ~ Latitude, merge(climate_locs, protein_D14_age)))

# Q6: Do the ratios of stress- and defense-related proteins, and photosynthesis-related proteins differ with aridity?

plot(stress_vs_photosynth ~ alpha, merge(climate_locs, protein_D14_age), xlab = 'Cramer Prentice Alpha (inverse of aridity)', ylab = "stress / photosynthesis protein ratio")
abline(lm(stress_vs_photosynth ~ alpha, merge(climate_locs, protein_D14_age)))
summary(lm(stress_vs_photosynth ~ alpha, merge(climate_locs, protein_D14_age)))

plot((Photorespiration + Light_reactions + Calvin_cycle) ~ alpha, merge(climate_locs, protein_D14_age), xlab = 'Cramer Prentice Alpha (inverse of aridity)', ylab = "photosynth proteins (mg/m2)")
abline(lm((Photorespiration + Light_reactions + Calvin_cycle) ~ alpha, merge(climate_locs, protein_D14_age)))
summary(lm((Photorespiration + Light_reactions + Calvin_cycle) ~ alpha, merge(climate_locs, protein_D14_age)))

plot(stress ~ alpha, merge(climate_locs, protein_D14_age), xlab = 'Cramer Prentice Alpha (inverse of aridity)', ylab = "stress (mg/m2)")
abline(lm(stress ~ alpha, merge(climate_locs, protein_D14_age)))
summary(lm(stress ~ alpha, merge(climate_locs, protein_D14_age)))

