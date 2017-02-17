# rubisco activase

# add in these: 'eucgr.l03031.1.p' 'eucgr.b02310.1.p'


rb <- protein_samples_D14[protein_samples_D14$Protein %in% c('eucgr.j02030.1.p', 'eucgr.j01234.1.p','eucgr.l03031.1.p','eucgr.b02310.1.p',
                                                             'eucgr.j01502.1.p',
                                                             'eucgr.j01502.2.p',
                                                             'eucgr.c00150.1.p',
                                                             'eucgr.b03013.1.p',
                                                             'eucgr.k02223.1.p'),]
rb <- arrange(rb, Protein)

names <- rb$Protein

rb <- as.data.frame(t(rb[,2:314]))

names(rb) <- names
rb$sample <- rownames(rb)
rownames(rb) <- NULL

rb <- merge(rb, total_protein_D14)

rb[,2:10] <- lapply(rb[,2:10],function(x) as.numeric(as.character(x)))

rb[,2:10] <- rb[,2:10]/rb$total_protein
str(rb)

rb$rb_smallsub <- rowSums(rb[,c('eucgr.j01502.1.p', 'eucgr.j01502.2.p', 'eucgr.c00150.1.p', 'eucgr.b03013.1.p', 'eucgr.k02223.1.p')])

rb$rb_act <- rowSums(rb[,c('eucgr.j02030.1.p', 'eucgr.j01234.1.p','eucgr.l03031.1.p','eucgr.b02310.1.p')])

bla <- merge(rb, climate_locs)

plot(bla$eucgr.j01234.1.p ~ bla$tavg, ylab = 'eucgr.j01234.1.p', xlab = 'tavg (oC)')
abline(lm(bla$eucgr.j01234.1.p ~ bla$tavg))
summary(lm(bla$eucgr.j01234.1.p ~ bla$tavg))

blax <- merge(bla, replicates)

blax <- blax %>% group_by(ID) %>% summarise(rbact_mean = mean(eucgr.j01234.1.p, na.rm=TRUE),
                                            rbact_SE = SE(eucgr.j01234.1.p))

blax <- merge(blax, merge(climate_locs, replicates))

plot(blax$rbact_mean ~ blax$tavg)

rbact.lm <- summary(lm(blax$rbact_mean ~ blax$tavg))

round(rbact.lm$coefficients[,4][2],2)


h <- ggplot(blax, aes(x = tavg, y = rbact_mean)) + geom_point() 
h <- h + geom_errorbar(aes(ymin = rbact_mean - rbact_SE, ymax = rbact_mean + rbact_SE), width = 0.2, alpha = 0.8)
h <- h + geom_smooth(method = 'lm')
h <- h + expand_limits(y=0)
h <- h + annotate("text" , x = max(blax$tavg)*0.3, y = max(blax$rbact_mean), label = paste('R2 =', round(rbact.lm$r.squared,3), sep = " "))
h <- h + annotate("text" , x = max(blax$tavg)*0.3, y = max(blax$rbact_mean)*0.9, label = paste('pval =', round(rbact.lm$coefficients[,4][2],2), sep = " "))
h <- h + ylab('mean rubisco activase')
h <- h + theme_bw()
h <- h + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), 
               legend.title=element_blank(),
               legend.position="bottom",
               axis.line = element_line(colour = "black"))
#h



blax <- merge(bla, replicates)

blax <- blax %>% group_by(ID) %>% summarise(rb_small_mean = mean(eucgr.k02223.1.p, na.rm=TRUE),
                                            rb_small_SE = SE(eucgr.k02223.1.p))

blax <- merge(blax, merge(climate_locs, replicates))

plot(blax$rb_small_mean ~ blax$tavg)

rb_small.lm <- summary(lm(blax$rb_small_mean ~ blax$tavg))

round(rb_small.lm$coefficients[,4][2],2)


h <- ggplot(blax, aes(x = tavg, y = rb_small_mean)) + geom_point() 
h <- h + geom_errorbar(aes(ymin = rb_small_mean - rb_small_SE, ymax = rb_small_mean + rb_small_SE), width = 0.2, alpha = 0.8)
h <- h + geom_smooth(method = 'lm')
h <- h + annotate("text" , x = max(blax$tavg)*0.3, y = max(blax$rb_small_mean), label = paste('R2 =', round(rb_small.lm$r.squared,3), sep = " "))
h <- h + annotate("text" , x = max(blax$tavg)*0.3, y = max(blax$rb_small_mean)*0.9, label = paste('pval =', round(rb_small.lm$coefficients[,4][2],2), sep = " "))
h <- h + ylab('mean eucgr.k02223.1.p')
h <- h + theme_bw()
h <- h + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), 
               legend.title=element_blank(),
               legend.position="bottom",
               axis.line = element_line(colour = "black"))
h

plot(blax$rb_small_SE ~ blax$tavg)

plot(bla$rb_smallsub_mean ~ bla$tavg)

plot(bla$rb_act_mean/bla$rb_smallsub_mean ~ bla$tavg)

plot(eucgr.j01502.1.p ~ tavg, bla)
summary(lm(eucgr.j01502.1.p ~ tavg, bla))

plot(eucgr.j01502.2.p ~ tavg, bla)
summary(lm(eucgr.j01502.2.p ~ tavg, bla))

plot(eucgr.c00150.1.p ~ tavg, bla)
summary(lm(eucgr.c00150.1.p ~ tavg, bla))

plot(eucgr.b03013.1.p ~ tavg, bla)
summary(lm(eucgr.b03013.1.p ~ tavg, bla))

plot(as.numeric(as.character(eucgr.k02223.1.p)) ~ tavg, bla)
summary(lm(as.numeric(as.character(eucgr.k02223.1.p)) ~ tavg, bla))

 
ip <- merge(bla,blax)

c('eucgr.j01502.1.p', 'eucgr.j01502.2.p', 'eucgr.c00150.1.p', 'eucgr.b03013.1.p', 'eucgr.k02223.1.p')



# large subunit of rubisco

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

plot(rb_small1 ~ tavg, data_means)
plot(rb_small2 ~ tavg, data_means)
plot(rb_small3 ~ tavg, data_means)
plot(rb_small4 ~ tavg, data_means)
plot(rb_small5 ~ tavg, data_means)




### plots ###
# calvin cycle ex. rubisco vs tavg, prec, light, d13C, total protein

plot(Calv_ex_rubisco ~ tavg, data_means)
abline(lm(Calv_ex_rubisco ~ tavg, data_means))
summary(lm(Calv_ex_rubisco ~ tavg, data_means))

plot(Calv_ex_rubisco ~ log10(prec), data_means)
abline(lm(Calv_ex_rubisco ~ log10(prec), data_means))
summary(lm(Calv_ex_rubisco ~ log10(prec), data_means))

plot(Calv_ex_rubisco ~ gap_mean, data_means)

plot(Calv_ex_rubisco ~ log10(leafrad_mean), data_means)

plot(Calv_ex_rubisco ~ d13C, data_means)

plot(Calv_ex_rubisco ~ total_protein_mean, data_means)


# rubisco activase vs PGK, PGLP, PRK, calv ex rb, RbcL, tavg, prec, light, d13C, total_protein

plot(rbact_mean ~ PGK_mean, data_means)
summary(lm(rbact_mean ~ PGK_mean, data_means))

plot(Rubisco_mean ~ PGK_mean, data_means)
summary(lm(Rubisco_mean ~ PGK_mean, data_means))

plot(rbact_mean ~ PRK_mean, data_means)
summary(lm(rbact_mean ~ PRK_mean, data_means))

plot(rbact_mean ~ PGLP_mean, data_means)
summary(lm(rbact_mean ~ PGLP_mean, data_means))


plot(rbact_mean ~ Calv_ex_rubisco, data_means)
plot(rbact_mean ~ rbL_mean, data_means)
summary(lm(rbact_mean ~ rbL_mean, data_means))
abline(lm(rbact_mean ~ rbL_mean, data_means))


plot(rbact_mean ~ tavg, data_means)
plot(rbact_mean ~ log10(prec), data_means)
plot(rbact_mean ~ gap_mean, data_means)
plot(rbact_mean ~ leafrad_mean, data_means)
plot(rbact_mean ~ d13C, data_means)
plot(rbact_mean ~ total_protein_mean, data_means)

# activase / PGK vs tavg, prec, light, d13c

plot((rbact_mean/PGK_mean) ~ tavg, data_means)
summary(lm((rbact_mean/PGK_mean) ~ tavg, data_means))
abline(lm((rbact_mean/PGK_mean) ~ tavg, data_means))


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

plot((rbact_mean/PGK_mean) ~ log10(prec), data_means)
plot((rbact_mean/PGK_mean) ~ gap_mean, data_means)
plot((rbact_mean/PGK_mean) ~ leafrad_mean, data_means)
plot((rbact_mean/PGK_mean) ~ d13C, data_means)
plot((rbact_mean/PGK_mean) ~ total_protein_mean, data_means)

# calvin cycle ex rubisco vs total protein, photosystems, light reactions, 

plot(Calvin_cycle_mean ~ Photosystems_mean, data_means)
plot(Calv_ex_rubisco ~ Photosystems_mean, data_means)

plot(Calv_ex_rubisco ~ total_protein_mean, data_means)
plot(Calvin_cycle_mean ~ total_protein_mean, data_means)

plot(Calv_ex_rubisco ~ Light_rxns_mean, data_means)
plot(Calvin_cycle_mean ~ Light_rxns_mean, data_means)
  

# oxidation : carboxylation ratios

plot(PGLP_mean ~ PRK_mean, data_means)

plot((PGLP_mean / PRK_mean) ~ tavg, data_means)
abline(lm((PGLP_mean / PRK_mean) ~ tavg, data_means))
summary(lm((PGLP_mean / PRK_mean) ~ tavg, data_means))

plot((PGLP_mean / PRK_mean) ~ leafrad_mean, data_means)
abline(lm((PGLP_mean / PRK_mean) ~ leafrad_mean, data_means))
summary(lm((PGLP_mean / PRK_mean) ~ leafrad_mean, data_means))

plot((PGLP_mean / PRK_mean) ~ gap_mean, data_means)
abline(lm((PGLP_mean / PRK_mean) ~ leafrad_mean, data_means))
summary(lm((PGLP_mean / PRK_mean) ~ leafrad_mean, data_means))



data_means$Calv_ex_act <- (data_means$Calvin_cycle_mean - data_means$rbact_mean)

plot(rbact_mean ~ Calv_ex_act, data_means)
summary(lm(rbact_mean ~ Calv_ex_act, data_means))
abline(lm(rbact_mean ~ Calv_ex_act, data_means))

plot((rbact_mean / rbL_mean) ~ tavg, data_means)
summary(lm((rbact_mean / rbL_mean) ~ tavg, data_means))


plot(chaperonin_mean ~ Calvin_cycle_mean, data_means)

## PRK and calvin cycle enzymes

source('scripts/transformations.R')

Currently 1.3.3 (Calvin cycle): eucgr.f04463.1.p, eucgr.e01261.1.p

Currently 4.1.11 (glycolysis): eucgr.f04463.1.p


1.2.1, PS.photorespiration.phosphoglycolate phosphatase, eucgr.b01439.1.p


source('scripts/transformations.R')

rb <- protein_samples_D14[protein_samples_D14$Protein %in% c('eucgr.e01261.1.p', 'eucgr.b01439.1.p', 'eucgr.f01476.1.p'),]
rb <- arrange(rb, Protein)

names <- rb$Protein

rb <- as.data.frame(t(rb[,2:314]))

names(rb) <- names
rb$sample <- rownames(rb)
rownames(rb) <- NULL

rb <- merge(rb, total_protein_D14)

rb[,2:4] <- lapply(rb[,2:4],function(x) as.numeric(as.character(x)))

#rb[,2:4] <- rb[,2:4]/rb$total_protein
str(rb)

bla <- merge(rb, protein_D14)
bla <- merge(bla, climate_locs)
bla <- merge(bla, replicates)

data_means <- bla %>% group_by(ID) %>% summarise(Rubisco_mean = mean(Rubisco, na.rm=TRUE),
                                                 PRK_mean = mean(eucgr.e01261.1.p, na.rm=TRUE),
                                                 PGLP_mean = mean(eucgr.b01439.1.p, na.rm=TRUE),
                                                 PGK_mean = mean(eucgr.f01476.1.p, na.rm=TRUE),
                                                 Photorespiration_mean = mean(Photorespiration, na.rm=TRUE)) %>%
  full_join(bla, by = 'ID')
  
plot((PGK_mean / Rubisco_mean) ~ tavg, data_means)
abline(lm((PGK_mean / Rubisco_mean) ~ tavg, data_means))
summary(lm((PGK_mean / Rubisco_mean) ~ tavg, data_means))

plot((PGK_mean / Rubisco_mean) ~ PGLP_mean, data_means)
abline(lm((PGK_mean / Rubisco_mean) ~ PGLP_mean, data_means))
summary(lm((PGK_mean / Rubisco_mean) ~ PGLP_mean, data_means))

 

plot((PGK_mean / Rubisco_mean) ~ leafrad_mean, data_means)
plot((PGK_mean / Rubisco_mean) ~ gap_mean, data_means)
plot((PGK_mean / Rubisco_mean) ~ log10(prec), data_means)
plot((PGK_mean / Rubisco_mean) ~ d13C, data_means)
plot((PGK_mean / Rubisco_mean) ~ total_protein, data_means)


plot(eucgr.e01261.1.p ~ Rubisco, bla)




















plot(bla$eucgr.j02030.1.p ~ bla$tavg)
abline(lm(bla$eucgr.j02030.1.p ~ bla$tavg))
summary(lm(bla$eucgr.j02030.1.p ~ bla$tavg))



blax <- merge(bla, replicates)

blax <- blax %>% group_by(ID) %>% summarise(rbact_mean = mean(eucgr.j02030.1.p, na.rm=TRUE),
                                            rbact_SE = SE(eucgr.j02030.1.p))

blax <- merge(blax, merge(climate_locs, replicates))

plot(blax$rbact_mean ~ blax$tavg)

rbact.lm <- summary(lm(blax$rbact_mean ~ blax$tavg))

round(rbact.lm$coefficients[,4][2],2)


h <- ggplot(blax, aes(x = tavg, y = rbact_mean)) + geom_point() 
h <- h + geom_errorbar(aes(ymin = rbact_mean - rbact_SE, ymax = rbact_mean + rbact_SE), width = 0.2, alpha = 0.8)
h <- h + geom_smooth(method = 'lm')
h <- h + expand_limits(y=0)
h <- h + annotate("text" , x = max(blax$tavg)*0.3, y = max(blax$rbact_mean), label = paste('R2 =', round(rbact.lm$r.squared,3), sep = " "))
h <- h + annotate("text" , x = max(blax$tavg)*0.3, y = max(blax$rbact_mean)*0.9, label = paste('pval =', round(rbact.lm$coefficients[,4][2],2), sep = " "))
h <- h + ylab('mean eucgr.j02030.1.p')
h <- h + theme_bw()
h <- h + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), 
               legend.title=element_blank(),
               legend.position="bottom",
               axis.line = element_line(colour = "black"))
h




