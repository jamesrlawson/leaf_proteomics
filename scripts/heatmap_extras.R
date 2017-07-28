# other things to add to the heatmap

rbact <- protein_samples_D14[protein_samples_D14$Protein %in% c('eucgr.j02030.1.p', 'eucgr.j01234.1.p'),]
rbact_ <- t(rbact[,2:327])

rbact_ <- as.data.frame(rbact_)
colnames(rbact_) <- rbact$Protein
rbact_$sample <- rownames(rbact_)
rbact_$rbact_sum <- rowSums(rbact_[,1:2])

# https://phytozome.jgi.doe.gov/pz/portal.html#!results?search=0&crown=1&star=1&method=4908&searchText=eucgr.k00881.1.p&offset=0
iso <- protein_samples_D14[protein_samples_D14$Protein %in% c('eucgr.k00881.1.p'),]
iso_ <- t(iso[,2:327])
iso_ <- as.data.frame(iso_)
colnames(iso_) <- iso$Protein
iso_$sample <- rownames(iso_)
names(iso_)[1] <- 'isoprene_synthase'


require(broom)
require(Hmisc)

# heatmaps looking at correlations between protein funccats

# uses 291 observations / 224

inc_photosynthesis = TRUE  # different from include_photosynthesis so it doesn't immediately affect transformations.R
include_photosynthesis = TRUE
include_d13C = TRUE
include_leaf_N = TRUE
include_leaf_P = TRUE
include_soil_N = TRUE
include_soil_P = TRUE
include_chlorophyll = TRUE

env_vars <- c('prec',
              'tavg',
              'leafrad_mean',
              'gap_mean')

if(include_soil_P) {
  env_vars <- c(env_vars, 'soil_P') 
}

if(include_soil_N) {
  env_vars <- c(env_vars, 'soil_N') 
}


trait_vars <- c()

if(include_d13C) {
  trait_vars <- c(trait_vars, 'd13C') 
}

if(include_leaf_P) {
  trait_vars <- c(trait_vars, 'Parea_mean') 
}

trait_vars <- c(trait_vars, 'LMA_mean')


if(include_leaf_N) {
  trait_vars <- c(trait_vars, 'Narea_mean') 
}

if(include_chlorophyll) {
  trait_vars <- c(trait_vars, 'mg_Cl_total_per_m2') 
}

if(inc_photosynthesis) {
  trait_vars <- c(trait_vars, 'photo_max', 'Cond') 
}

prot_vars <- c('total_protein_mean',
               'rubisco_mean', 
               'rbact_mean',
               'calvin_cycle_mean', 
               'photosystems_mean', 
               'ATP_synthase_chloroplastic_mean',
               'electron_transport_mean',
               'photorespiration_mean',
               'protein_mean',
               'heatstress_mean',
               'isoprene_synthase_mean')

cor_vars <- c(env_vars, trait_vars, prot_vars)

source('scripts/transformations_photomax.R')

source('scripts/prep_data.R')

data <- merge(data, rbact_, by = 'sample')
data$rbact_sum <- data$rbact_sum / data$total_protein

data <- merge(data, iso_, by = 'sample')
data$isoprene_synthase <- data$isoprene_synthase / data$total_protein

means <- data %>% dplyr::group_by(ID) %>% dplyr::summarise(calvin_cycle_mean = mean(calvin_cycle, na.rm=TRUE),
                                                           rubisco_mean = mean(Rubisco, na.rm=TRUE),
                                                           photosystems_mean = mean(Photosystems, na.rm=TRUE),
                                                           ATP_synthase_chloroplastic_mean = mean(ATP_synthase_chloroplastic, na.rm=TRUE),
                                                           electron_transport_mean = mean(electron_transport_minATPsynth, na.rm=TRUE),
                                                           photorespiration_mean = mean(photorespiration,na.rm=TRUE),
                                                           protein_mean = mean(protein, na.rm=TRUE),
                                                           heatstress_mean = mean(stress.abiotic.heat, na.rm=TRUE),
                                                           rbact_mean = mean(rbact_sum, na.rm=TRUE),
                                                           isoprene_synthase_mean = mean(isoprene_synthase, na.rm=TRUE))
data <-  full_join(data, means)
data <- distinct(data, calvin_cycle_mean, rubisco_mean, photosystems_mean, ATP_synthase_chloroplastic_mean, electron_transport_mean, photorespiration_mean, protein_mean, 
                 heatstress_mean, rbact_mean, isoprene_synthase_mean, .keep_all=TRUE)

prot <- data
prot$calvin_cycle_mean <- prot$calvin_cycle_mean - prot$rubisco_mean
prot$prec <- log10(prot$prec)
cormat.prot_rel <- Hmisc::rcorr(as.matrix(prot[,cor_vars]), type = 'pearson')


source('scripts/prep_data_mg_per_mm2.R')

data <- merge(data, rbact_, by = 'sample')
data <- merge(data, iso_, by = 'sample')

means <- data %>% dplyr::group_by(ID) %>% dplyr::summarise(calvin_cycle_mean = mean(calvin_cycle, na.rm=TRUE),
                                                           rubisco_mean = mean(Rubisco, na.rm=TRUE),
                                                           photosystems_mean = mean(Photosystems, na.rm=TRUE),
                                                           ATP_synthase_chloroplastic_mean = mean(ATP_synthase_chloroplastic, na.rm=TRUE),
                                                           electron_transport_mean = mean(electron_transport_minATPsynth, na.rm=TRUE),
                                                           photorespiration_mean = mean(photorespiration,na.rm=TRUE),
                                                           protein_mean = mean(protein, na.rm=TRUE),
                                                           heatstress_mean = mean(stress.abiotic.heat, na.rm=TRUE),
                                                           rbact_mean = mean(rbact_sum, na.rm=TRUE),
                                                           isoprene_synthase_mean = mean(isoprene_synthase, na.rm=TRUE))
data <-  full_join(data, means)
data <- distinct(data, calvin_cycle_mean, rubisco_mean, photosystems_mean, ATP_synthase_chloroplastic_mean, electron_transport_mean, photorespiration_mean, protein_mean, 
                 heatstress_mean, rbact_mean, isoprene_synthase_mean, .keep_all=TRUE)

prot <- data
prot$calvin_cycle_mean <- prot$calvin_cycle_mean - prot$rubisco_mean
prot$prec <- log10(prot$prec)
cormat.prot_abs <- Hmisc::rcorr(as.matrix(prot[,cor_vars]), type = 'pearson')


cormat_all <- cormat.prot_abs$r
cormat_all[upper.tri(cormat_all)] <- cormat.prot_rel$r[upper.tri(cormat.prot_rel$r)]
cormat_all <- melt(cormat_all, na.rm=TRUE)

cormat_all.P <- cormat.prot_abs$P
cormat_all.P[upper.tri(cormat_all.P)] <- cormat.prot_rel$P[upper.tri(cormat.prot_rel$P)]
cormat_all.P <- melt(cormat_all.P, na.rm=FALSE)
cormat_all.P[is.na(cormat_all.P$value),]$value <- 0

bla <- full_join(cormat_all, cormat_all.P, by = c('Var1', 'Var2'))


p <- ggplot(bla, aes(x = Var1, y = Var2))
p <- p  + geom_raster(data = subset(bla, value.y < 0.05), aes(fill = value.x))
p <- p + ggtitle('Correlation heatmap (lower = abs, upper = rel)')
p <- p + scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                              midpoint = 0, limit = c(-1,1), space = "Lab", 
                              name="Pearson\nCorrelation")
p <- p + theme_minimal()  + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
                                  axis.title=element_blank()) + coord_fixed()
p

rm(include_photosynthesis,
   include_d13C,
   include_leaf_N,
   include_leaf_P,
   include_soil_N,
   include_soil_P,
   include_chlorophyll)
