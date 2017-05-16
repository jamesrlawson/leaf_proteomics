require(broom)
require(Hmisc)

# heatmaps looking at correlations between protein funccats

inc_photosynthesis = TRUE  # different from include_photosynthesis so it doesn't immediately affect transformations.R
include_d13C = TRUE
include_leaf_N = TRUE
include_leaf_P = TRUE
include_soil_N = TRUE
include_soil_P = TRUE
include_chlorophyll = TRUE

env_vars <- c('prec',
              'tavg',
              'leaf_rad',
              'gap')

if(include_soil_N) {
  env_vars <- c(env_vars, 'soil_N') 
}

if(include_soil_P) {
  env_vars <- c(env_vars, 'soil_P') 
}

trait_vars <- c('LMA_g_per_m2')

if(include_d13C) {
  trait_vars <- c(trait_vars, 'd13C') 
}

if(include_leaf_N) {
  trait_vars <- c(trait_vars, 'N_per_area') 
}

if(include_leaf_P) {
  trait_vars <- c(trait_vars, 'P_per_area') 
}

if(include_chlorophyll) {
  trait_vars <- c(trait_vars, 'mg_Cl_total_per_m2') 
}

source('scripts/transformations.R')

source('scripts/prep_data.R')
prot <- data
prot$calvin_cycle <- prot$calvin_cycle - prot$Rubisco
prot$prec <- log10(prot$prec)
cormat.prot_rel <- Hmisc::rcorr(as.matrix(prot[,cor_vars]), type = 'pearson')


source('scripts/prep_data_mg_per_mm2.R')
prot <- data
prot$calvin_cycle <- prot$calvin_cycle - prot$Rubisco
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
p <- p + scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", 
                              midpoint = 0, limit = c(-1,1), space = "Lab", 
                              name="Pearson\nCorrelation")
p <- p + theme_minimal()  + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + coord_fixed()
p

rm(include_photosynthesis,
   include_d13C,
   include_leaf_N,
   include_leaf_P,
   include_soil_N,
   include_soil_P,
   include_chlorophyll)
