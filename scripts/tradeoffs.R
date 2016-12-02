# TRADEOFFS

source('scripts/transformations.R')

protein_D14$PS_ratio <- protein_D14$PSI/protein_D14$PSII

protein_clim_D14 <- merge(protein_D14, climate_locs, by = 'sample')

# PSI vs PSII

plot(PSI ~ PSII, protein_D14, xlab = "PSII (umol / m2)", ylab = "PSI (umol/m2)", main = "PSI/PSII stoichiometry in wild Eucalyptus")
summary(lm(PSI ~ PSII, protein_D14))
abline(lm(PSI ~ PSII, protein_D14))


require(ggplot2)

p <- ggplot(protein_D14, aes(x = PSII, y = PSI))
p <- p + geom_point()
p <- p + geom_smooth()
print(p)



plot(PSI ~ PSII, protein_stand_D14)
summary(lm(PSI ~ PSII, protein_stand_D14))


plot(PSI ~ prec, protein_clim_D14)
plot(PSII ~ prec, protein_clim_D14)
plot(PS_ratio ~ prec, protein_clim_D14)
summary(lm(PS_ratio ~ prec, protein_clim_D14))
abline(lm(PS_ratio ~ prec, protein_clim_D14))


plot(PS_ratio ~ alpha, protein_clim_D14)
summary(lm(PS_ratio ~ alpha, protein_clim_D14))

plot(PS_ratio ~ gap, protein_clim_D14)
summary(lm(PS_ratio ~ gap, protein_clim_D14))
abline(lm(PS_ratio ~ gap, protein_clim_D14))
plot(PSI ~ gap, protein_clim_D14)
summary(lm(PSI ~ gap, protein_clim_D14))

# can variation in PSI ~ PSII can be accounted for by env variables?

require(MuMIn)

blah <- lm(PSI ~ PSII, protein_clim_D14)
blax <- lm(PSI ~ PSII + prec, protein_clim_D14)
blaz <- lm(PSI ~ PSII * prec, protein_clim_D14)

AIC(blah, blax, blaz)

summary(blax)
anova(blax)

protein_clim_D14.naomit <- protein_clim_D14[!is.na(protein_clim_D14$gap),]

blah <- lm(PSI ~ PSII, protein_clim_D14.naomit)
blax <- lm(PSI ~ PSII + gap, protein_clim_D14.naomit)
blaz <- lm(PSI ~ PSII * gap, protein_clim_D14.naomit)

AIC(blah, blax, blaz)

summary(blax)
anova(blax)


require(vegan)

PS_ratio.varpart <- varpart(protein_clim_D14.naomit$PSI,
                            ~ PSII,
                            ~ gap,
                            data = protein_clim_D14.naomit)
PS_ratio.varpart
plot(PS_ratio.varpart)


# light rxns vs calvin cycle

plot(Light_reactions ~ Calvin_cycle, protein_clim_D14)

protein_clim_D14_stand <- merge(protein_stand_D14, climate_locs, by = 'sample')

plot(Light_reactions ~ Calvin_cycle, protein_clim_D14_stand)
summary(lm(Light_reactions ~ Calvin_cycle, protein_clim_D14_stand))

plot(PSI ~ Calvin_cycle, protein_clim_D14_stand) ##
summary(lm(PSI ~ Calvin_cycle, protein_clim_D14_stand)) ##
abline(lm(PSI ~ Calvin_cycle, protein_clim_D14_stand)) ##


plot(PSII ~ Calvin_cycle, protein_clim_D14_stand) 
summary(lm(PSII ~ Calvin_cycle, protein_clim_D14_stand))

plot(Photosystems ~ Calvin_cycle, protein_clim_D14_stand)
summary(lm(Photosystems ~ Calvin_cycle, protein_clim_D14_stand))

plot(PSI ~ Rubisco, protein_clim_D14_stand) ##
summary(lm(PSI ~ Rubisco, protein_clim_D14_stand)) ##

plot(PSII ~ Rubisco, protein_clim_D14_stand)
summary(lm(PSII ~ Rubisco, protein_clim_D14_stand))

plot(Photosystems ~ Rubisco, protein_clim_D14_stand)
summary(lm(Photosystems ~ Rubisco, protein_clim_D14_stand))

plot(Cytochrome_b6f ~ Rubisco, protein_clim_D14_stand)
summary(lm(Cytochrome_b6f ~ Rubisco, protein_clim_D14_stand))

plot(Cytochrome_b6f ~ Rubisco, protein_clim_D14_stand)

summary(lm(Cytochrome_b6f ~ Rubisco, protein_clim_D14_stand))

# light rxns / calvin cycle ratio vs env's

protein_stand_D14$lightDark_ratio <- protein_stand_D14$Light_reactions / protein_stand_D14$Calvin_cycle
protein_clim_D14_stand <- merge(protein_stand_D14, climate_locs, by = 'sample')

plot(lightDark_ratio ~ prec, protein_clim_D14_stand)
plot(lightDark_ratio ~ gap, protein_clim_D14_stand)

summary(lm(lightDark_ratio ~ gap, protein_clim_D14_stand))
abline(lm(lightDark_ratio ~ gap, protein_clim_D14_stand))

protein_stand_D14$PSIRubisco_ratio <- protein_stand_D14$PSI / protein_stand_D14$Rubisco
protein_clim_D14_stand <- merge(protein_stand_D14, climate_locs, by = 'sample')

plot(PSIRubisco_ratio ~ prec, protein_clim_D14_stand)
plot(PSIRubisco_ratio ~ gap, protein_clim_D14_stand)

summary(lm(PSIRubisco_ratio ~ prec, protein_clim_D14_stand))
abline(lm(PSIRubisco_ratio ~ prec, protein_clim_D14_stand))
