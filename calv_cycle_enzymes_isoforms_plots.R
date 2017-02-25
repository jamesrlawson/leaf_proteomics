# Relationships between Calvin cycle enzymes and rubisco small subunit / activase isoforms

source('scripts/transformations.R')

# pull out particular proteins from protein_samples_D14 and arrange

rbact <- protein_samples_D14[protein_samples_D14$Protein %in% c('eucgr.j02030.1.p', 'eucgr.j01234.1.p','eucgr.l03031.1.p','eucgr.b02310.1.p','eucgr.b02532.1.p'),] # rubisco activase, chaperonin
rb <- protein_samples_D14[protein_samples_D14$Protein %in% c('eucgr.c03525.1.p', # rubisco large subunit
                                                             'tr|t1qkk4|t1qkk4_eucgl'),]
rb <- rbind(rb, rbact)
rb <- rbind(rb, protein_samples_D14[protein_samples_D14$Protein %in% c('eucgr.e01261.1.p', 'eucgr.b01439.1.p', 'eucgr.f01476.1.p'),]) # PRK, PGLP, PGK

rb <- rbind(rb, protein_samples_D14[protein_samples_D14$Protein %in% c('eucgr.j01502.1.p', 'eucgr.j01502.2.p', 'eucgr.c00150.1.p', 'eucgr.b03013.1.p', 'eucgr.k02223.1.p'),]) # small subunit of rubisco


rb <- arrange(rb, Protein)
names <- rb$Protein
rb <- as.data.frame(t(rb[,2:314]))

names(rb) <- names
rb$sample <- rownames(rb)
rownames(rb) <- NULL

rb <- merge(rb, total_protein_D14)

rb[,2:(ncol(rb)-1)] <- lapply(rb[,2:(ncol(rb)-1)],function(x) as.numeric(as.character(x)))

rb[,2:(ncol(rb)-1)] <- rb[,2:(ncol(rb)-1)]/rb$total_protein # standardise to total protein amount
str(rb)

# sum components of multi-protein complexes

rb$rb_largesub <- rowSums(rb[,c('eucgr.c03525.1.p',
                                'tr|t1qkk4|t1qkk4_eucgl'),])

rb$rb_act <- rowSums(rb[,c('eucgr.j02030.1.p', 
                           'eucgr.j01234.1.p',
                           'eucgr.l03031.1.p',
                           'eucgr.b02310.1.p')])


# merge with full protein and env data sets

bla <- merge(rb, protein_D14)
bla <- merge(bla, climate_locs)
bla <- merge(bla, replicates)

# aggregate to ID-wise means

data_means <- bla %>% group_by(ID) %>% summarise(Rubisco_mean = mean(Rubisco, na.rm=TRUE),
                                                 Calvin_cycle_mean = mean(Calvin_cycle, na.rm=TRUE),
                                                 PRK_mean = mean(eucgr.e01261.1.p, na.rm=TRUE),
                                                 PGLP_mean = mean(eucgr.b01439.1.p, na.rm=TRUE),
                                                 PGK_mean = mean(eucgr.f01476.1.p, na.rm=TRUE),
                                                 Photorespiration_mean = mean(Photorespiration, na.rm=TRUE),
                                                 rbact_mean = mean(rb_act, na.rm=TRUE),
                                                 rbL_mean = mean(rb_largesub, na.rm=TRUE),
                                                 chaperonin_mean = mean(eucgr.b02532.1.p, na.rm=TRUE),
                                                 total_protein_mean = mean(total_protein),
                                                 Photosystems_mean = mean(Photosystems, na.rm=TRUE),
                                                 Light_rxns_mean = mean(Photosystems, na.rm=TRUE),
                                                 rb_small1 = mean(eucgr.j01502.1.p, n.rm=TRUE),
                                                 rb_small2 = mean(eucgr.j01502.2.p, na.rm=TRUE),
                                                 rb_small3 = mean(eucgr.c00150.1.p, na.rm=TRUE),
                                                 rb_small4 = mean(eucgr.b03013.1.p, na.rm=TRUE),
                                                 rb_small5 = mean(eucgr.k02223.1.p, na.rm=TRUE)) %>%
  full_join(bla, by = 'ID') %>%
  mutate(Calv_ex_rubisco = Calvin_cycle_mean - Rubisco_mean)



# small subunit of rubisco plots

dir.create('output/figures/20170217/rbsmall_isoforms/', showWarnings = FALSE)

png(filename = paste('output/figures/20170217/rbsmall_isoforms/', 'rb_small1-tavg_', 'R2-', round(summary(lm(rb_small1 ~ tavg, data_means))$r.squared,3), '.png', sep=""))
plot(rb_small1 ~ tavg, data_means)
abline(lm(rb_small1 ~ tavg, data_means))
dev.off()
png(filename = paste('output/figures/20170217/rbsmall_isoforms/', 'rb_small2-tavg_', 'R2-', round(summary(lm(rb_small2 ~ tavg, data_means))$r.squared,3), '.png', sep=""))
plot(rb_small2 ~ tavg, data_means)
abline(lm(rb_small2 ~ tavg, data_means))
dev.off()
png(filename = paste('output/figures/20170217/rbsmall_isoforms/', 'rb_small3-tavg_', 'R2-', round(summary(lm(rb_small3 ~ tavg, data_means))$r.squared,3), '.png', sep=""))
plot(rb_small3 ~ tavg, data_means)
abline(lm(rb_small3 ~ tavg, data_means))
dev.off()
png(filename = paste('output/figures/20170217/rbsmall_isoforms/', 'rb_small4-tavg_', 'R2-', round(summary(lm(rb_small4 ~ tavg, data_means))$r.squared,3), '.png', sep=""))
plot(rb_small4 ~ tavg, data_means)
abline(lm(rb_small4 ~ tavg, data_means))
dev.off()
png(filename = paste('output/figures/20170217/rbsmall_isoforms/', 'rb_small5-tavg_', 'R2-', round(summary(lm(rb_small5 ~ tavg, data_means))$r.squared,3), '.png', sep=""))
plot(rb_small5 ~ tavg, data_means)
abline(lm(rb_small5 ~ tavg, data_means))
dev.off()


# calvin cycle ex. rubisco vs tavg, prec, light, d13C, total protein

dir.create('output/figures/20170217/calv_exRubisco_vs_envvars_otherProteins/', showWarnings = FALSE)

png(filename = paste('output/figures/20170217/calv_exRubisco_vs_envvars_otherProteins/', 'Calv_ex_rubisco-tavg_', 'R2-', round(summary(lm(Calv_ex_rubisco ~ tavg, data_means))$r.squared,3), '.png', sep=""))
plot(Calv_ex_rubisco ~ tavg, data_means)
abline(lm(Calv_ex_rubisco ~ tavg, data_means))
dev.off()

png(filename = paste('output/figures/20170217/calv_exRubisco_vs_envvars_otherProteins/', 'Calv_ex_rubisco-log10prec_', 'R2-', round(summary(lm(Calv_ex_rubisco ~ log10(prec), data_means))$r.squared,3), '.png', sep=""))
plot(Calv_ex_rubisco ~ log10(prec), data_means)
abline(lm(Calv_ex_rubisco ~ log10(prec), data_means))
dev.off()

png(filename = paste('output/figures/20170217/calv_exRubisco_vs_envvars_otherProteins/', 'Calv_ex_rubisco-gap_mean_', 'R2-', round(summary(lm(Calv_ex_rubisco ~ gap_mean, data_means))$r.squared,3), '.png', sep=""))
plot(Calv_ex_rubisco ~ gap_mean, data_means)
abline(lm(Calv_ex_rubisco ~ gap_mean, data_means))
dev.off()

png(filename = paste('output/figures/20170217/calv_exRubisco_vs_envvars_otherProteins/', 'Calv_ex_rubisco-tavg_', 'R2-', round(summary(lm(Calv_ex_rubisco ~ leafrad_mean, data_means))$r.squared,3), '.png', sep=""))
plot(Calv_ex_rubisco ~ leafrad_mean, data_means)
abline(lm(Calv_ex_rubisco ~ leafrad_mean, data_means))
dev.off()

#png(filename = paste('output/figures/20170217/calv_exRubisco_vs_envvars_otherProteins/', 'Calv_ex_rubisco-d13C_', 'R2-', round(summary(lm(Calv_ex_rubisco ~ d13C, data_means))$r.squared,3), '.png', sep=""))
#  plot(Calv_ex_rubisco ~ d13C, data_means)
#  abline(lm(Calv_ex_rubisco ~ d13C, data_means))
#dev.off()

png(filename = paste('output/figures/20170217/calv_exRubisco_vs_envvars_otherProteins/', 'Calv_ex_rubisco-total_protein_mean_', 'R2-', round(summary(lm(Calv_ex_rubisco ~ total_protein_mean, data_means))$r.squared,3), '.png', sep=""))
plot(Calv_ex_rubisco ~ total_protein_mean, data_means)
abline(lm(Calv_ex_rubisco ~ total_protein_mean, data_means))
dev.off()

png(filename = paste('output/figures/20170217/calv_exRubisco_vs_envvars_otherProteins/', 'Calv_ex_rubisco-Photosystems_mean_', 'R2-', round(summary(lm(Calv_ex_rubisco ~ Photosystems_mean, data_means))$r.squared,3), '.png', sep=""))
plot(Calv_ex_rubisco ~ Photosystems_mean, data_means)
abline(lm(Calv_ex_rubisco ~ Photosystems_mean, data_means))
dev.off()

png(filename = paste('output/figures/20170217/calv_exRubisco_vs_envvars_otherProteins/', 'Calv_ex_rubisco-total_protein_mean_', 'R2-', round(summary(lm(Calv_ex_rubisco ~ total_protein_mean, data_means))$r.squared,3), '.png', sep=""))
plot(Calv_ex_rubisco ~ total_protein_mean, data_means)
abline(lm(Calv_ex_rubisco ~ total_protein_mean, data_means))
dev.off()

png(filename = paste('output/figures/20170217/calv_exRubisco_vs_envvars_otherProteins/', 'Calv_ex_rubisco-Light_rxns_mean_', 'R2-', round(summary(lm(Calv_ex_rubisco ~ Light_rxns_mean, data_means))$r.squared,3), '.png', sep=""))
plot(Calv_ex_rubisco ~ Light_rxns_mean, data_means)
abline(lm(Calv_ex_rubisco ~ Light_rxns_mean, data_means))
dev.off()



# rubisco activase vs PGK, PGLP, PRK, calv ex rb, RbcL, tavg, prec, light, d13C, total_protein

dir.create('output/figures/20170217/rbact_vs_envvars_calvEnzymes/', showWarnings = FALSE)

png(filename = paste('output/figures/20170217/rbact_vs_envvars_calvEnzymes/', 'rbact_mean-PGK_mean_', 'R2-', round(summary(lm(rbact_mean ~ PGK_mean, data_means))$r.squared,3), '.png', sep=""))
plot(rbact_mean ~ PGK_mean, data_means)
abline(lm(rbact_mean ~ PGK_mean, data_means))
dev.off()

png(filename = paste('output/figures/20170217/rbact_vs_envvars_calvEnzymes/', 'rbact_mean-PRK_mean_', 'R2-', round(summary(lm(rbact_mean ~ PRK_mean, data_means))$r.squared,3), '.png', sep=""))
plot(rbact_mean ~ PRK_mean, data_means)
abline(lm(rbact_mean ~ PRK_mean, data_means))
dev.off()

png(filename = paste('output/figures/20170217/rbact_vs_envvars_calvEnzymes/', 'rbact_mean-PGLP_mean_', 'R2-', round(summary(lm(rbact_mean ~ PGLP_mean, data_means))$r.squared,3), '.png', sep=""))
plot(rbact_mean ~ PGLP_mean, data_means)
abline(lm(rbact_mean ~ PGLP_mean, data_means))
dev.off()

png(filename = paste('output/figures/20170217/rbact_vs_envvars_calvEnzymes/', 'rbact_mean-Calv_ex_rubisco_', 'R2-', round(summary(lm(rbact_mean ~ Calv_ex_rubisco, data_means))$r.squared,3), '.png', sep=""))
plot(rbact_mean ~ Calv_ex_rubisco, data_means)
dev.off()

png(filename = paste('output/figures/20170217/rbact_vs_envvars_calvEnzymes/', 'rbact_mean-rbL_mean_', 'R2-', round(summary(lm(rbact_mean ~ rbL_mean, data_means))$r.squared,3), '.png', sep=""))
plot(rbact_mean ~ rbL_mean, data_means)
abline(lm(rbact_mean ~ rbL_mean, data_means))
dev.off()

png(filename = paste('output/figures/20170217/rbact_vs_envvars_calvEnzymes/', 'rbact_mean-tavg_', 'R2-', round(summary(lm(rbact_mean ~ tavg, data_means))$r.squared,3), '.png', sep=""))
plot(rbact_mean ~ tavg, data_means)
abline(lm(rbact_mean ~ tavg, data_means))
dev.off()

png(filename = paste('output/figures/20170217/rbact_vs_envvars_calvEnzymes/', 'rbact_mean-log10prec_', 'R2-', round(summary(lm(rbact_mean ~ log10(prec), data_means))$r.squared,3), '.png', sep=""))
plot(rbact_mean ~ log10(prec), data_means)
abline(lm(rbact_mean ~ log10(prec), data_means))
dev.off()

png(filename = paste('output/figures/20170217/rbact_vs_envvars_calvEnzymes/', 'rbact_mean-gap_mean_', 'R2-', round(summary(lm(rbact_mean ~ gap_mean, data_means))$r.squared,3), '.png', sep=""))
plot(rbact_mean ~ gap_mean, data_means)
abline(lm(rbact_mean ~ gap_mean, data_means))
dev.off()

png(filename = paste('output/figures/20170217/rbact_vs_envvars_calvEnzymes/', 'rbact_mean-leafrad_mean_', 'R2-', round(summary(lm(rbact_mean ~ leafrad_mean, data_means))$r.squared,3), '.png', sep=""))
plot(rbact_mean ~ leafrad_mean, data_means)
abline(lm(rbact_mean ~ leafrad_mean, data_means))
dev.off()

#png(filename = paste('output/figures/20170217/', 'rbact_mean-PGK_mean_', 'R2-', round(summary(lm(rbact_mean ~ PGK_mean, data_means))$r.squared,3), '.png', sep=""))
#  plot(rbact_mean ~ d13C, data_means)
#  summary(lm(rbact_mean ~ d13C, data_means))
#dev.off()

png(filename = paste('output/figures/20170217/rbact_vs_envvars_calvEnzymes/', 'rbact_mean-total_protein_mean_', 'R2-', round(summary(lm(rbact_mean ~ total_protein_mean, data_means))$r.squared,3), '.png', sep=""))
plot(rbact_mean ~ total_protein_mean, data_means)
abline(lm(rbact_mean ~ total_protein_mean, data_means))
dev.off()



# activase / PGK vs tavg, prec, light, d13c

dir.create('output/figures/20170217/rbactPerPGK_vs_envvars/', showWarnings = FALSE)

png(filename = paste('output/figures/20170217/rbactPerPGK_vs_envvars/', 'rbactPerPGK_mean-tavg_', 'R2-', round(summary(lm((rbact_mean/PGK_mean) ~ tavg, data_means))$r.squared,3), '.png', sep=""))
plot((rbact_mean/PGK_mean) ~ tavg, data_means)
abline(lm((rbact_mean/PGK_mean) ~ tavg, data_means))
dev.off()



blax <- lm(scale(rbact_mean) ~ scale(PGK_mean)  + scale(tavg), data_means)
summary(blax)
blaz <- lm(scale(data_means$rbact_mean/data_means$total_protein_mean) ~  scale(data_means$PGK_mean/data_means$total_protein_mean) + scale(tavg) , data_means)
summary(blaz)

blax <- lm(scale(data_means$rbact_mean/data_means$total_protein_mean) ~  scale(data_means$PGK_mean/data_means$total_protein_mean) + tavg  , data_means)

plot((rbact_mean/total_protein_mean) ~ tavg, data_means)

require(vegan)
require(MuMIn)

dredge(blax)

up <- varpart(data_means$rbact_mean/data_means$total_protein_mean,
              ~ data_means$PGK_mean/data_means$total_protein_mean,
              ~ data_means$tavg)
up
plot(up)

plot(data_means$rbact_mean/data_means$total_protein_mean ~ data_means$tavg)
summary(lm(data_means$rbact_mean/data_means$total_protein_mean ~ data_means$tavg))

plot(PGK_mean/total_protein_mean ~ tavg, data_means)


png(filename = paste('output/figures/20170217/rbactPerPGK_vs_envvars/', 'rbactPerPGK_mean-log10prec_', 'R2-', round(summary(lm((rbact_mean/PGK_mean) ~ log10(prec), data_means))$r.squared,3), '.png', sep=""))
plot((rbact_mean/PGK_mean) ~ log10(prec), data_means)
abline((rbact_mean/PGK_mean) ~ log10(prec), data_means)
dev.off()

png(filename = paste('output/figures/20170217/rbactPerPGK_vs_envvars/', 'rbactPerPGK_mean-gap_mean_', 'R2-', round(summary(lm((rbact_mean/PGK_mean) ~ gap_mean, data_means))$r.squared,3), '.png', sep=""))
plot((rbact_mean/PGK_mean) ~ gap_mean, data_means)
abline(lm((rbact_mean/PGK_mean) ~ gap_mean, data_means))
dev.off()

png(filename = paste('output/figures/20170217/rbactPerPGK_vs_envvars/', 'rbactPerPGK_mean-leafrad_mean_', 'R2-', round(summary(lm((rbact_mean/PGK_mean) ~ leafrad_mean, data_means))$r.squared,3), '.png', sep=""))
plot((rbact_mean/PGK_mean) ~ leafrad_mean, data_means)
abline(lm((rbact_mean/PGK_mean) ~ leafrad_mean, data_means))
dev.off()
#png(filename = paste('output/figures/20170217/rbactPerPGK_vs_envvars/', 'rbactPerPGK_mean-log10prec_', 'R2-', round(summary(lm((rbact_mean/PGK_mean) ~ log10(prec), data_means))$r.squared,3), '.png', sep=""))
#  plot((rbact_mean/PGK_mean) ~ d13C, data_means)
#  abline(lm((rbact_mean/PGK_mean) ~ d13C, data_means))
#dev.off()

png(filename = paste('output/figures/20170217/rbactPerPGK_vs_envvars/', 'rbactPerPGK_mean-total_protein_mean_', 'R2-', round(summary(lm((rbact_mean/PGK_mean) ~ total_protein_mean, data_means))$r.squared,3), '.png', sep=""))
  plot((rbact_mean/PGK_mean) ~ total_protein_mean, data_means)
  abline(lm((rbact_mean/PGK_mean) ~ total_protein_mean, data_means))
dev.off()


# oxidation : carboxylation ratios

dir.create('output/figures/20170217/oxidation_carboxylation/', showWarnings = FALSE)

png(filename = paste('output/figures/20170217/oxidation_carboxylation/', 'PGLP_mean-PRK_mean_', 'R2-', round(summary(lm(PGLP_mean ~ PRK_mean, data_means))$r.squared,3), '.png', sep=""))
plot(PGLP_mean ~ PRK_mean, data_means)
plot(PGLP_mean ~ PRK_mean, data_means)
dev.off()

png(filename = paste('output/figures/20170217/oxidation_carboxylation/', 'PGLP_per_PRK_mean-tavg_', 'R2-', round(summary(lm(PGLP_mean ~ tavg, data_means))$r.squared,3), '.png', sep=""))
plot((PGLP_mean / PRK_mean) ~ tavg, data_means)
abline(lm((PGLP_mean / PRK_mean) ~ tavg, data_means))
dev.off()

png(filename = paste('output/figures/20170217/oxidation_carboxylation/', 'PGLP_per_PRK_mean-tavg_', 'R2-', round(summary(lm(PGLP_mean ~ tavg, data_means))$r.squared,3), '.png', sep=""))
plot((PGLP_mean / PRK_mean) ~ tavg, data_means)
abline(lm((PGLP_mean / PRK_mean) ~ tavg, data_means))
dev.off()

png(filename = paste('output/figures/20170217/oxidation_carboxylation/', 'PGLP_per_PRK_mean-gap_mean_', 'R2-', round(summary(lm(PGLP_mean ~ gap_mean, data_means))$r.squared,3), '.png', sep=""))
plot((PGLP_mean / PRK_mean) ~ gap_mean, data_means)
abline(lm((PGLP_mean / PRK_mean) ~ leafrad_mean, data_means))
dev.off()

# misc

dir.create('output/figures/20170217/misc/', showWarnings = FALSE)

data_means$Calv_ex_act <- (data_means$Calvin_cycle_mean - data_means$rbact_mean)

png(filename = paste('output/figures/20170217/misc/', 'rbact_mean-Calv_ex_act_', 'R2-', round(summary(lm(rbact_mean ~ Calv_ex_act, data_means))$r.squared,3), '.png', sep=""))
  plot(rbact_mean ~ Calv_ex_act, data_means)
  abline(lm(rbact_mean ~ Calv_ex_act, data_means))
dev.off()

png(filename = paste('output/figures/20170217/misc/', 'rbactPerRBL_mean-Calv_ex_act_', 'R2-', round(summary(lm((rbact_mean / rbL_mean) ~ tavg, data_means))$r.squared,3), '.png', sep=""))
  plot((rbact_mean / rbL_mean) ~ tavg, data_means)
  abline(lm((rbact_mean / rbL_mean) ~ tavg, data_means))
dev.off()


png(filename = paste('output/figures/20170217/misc/', 'rbactPerRBL_mean-Calv_ex_act_', 'R2-', round(summary(lm((rbact_mean / rbL_mean) ~ tavg, data_means))$r.squared,3), '.png', sep=""))
  plot(chaperonin_mean ~ Calvin_cycle_mean, data_means)
  abline(lm(chaperonin_mean ~ Calvin_cycle_mean, data_means))
dev.off()