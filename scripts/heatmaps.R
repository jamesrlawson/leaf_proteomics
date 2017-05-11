# heatmaps looking at correlations between protein funccats

include_photosynthesis = TRUE
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

if(include_photosynthesis) {
  trait_vars <- c(trait_vars, 'photo_max', 'Cond') 
}

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

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

reorder_cormat_env <- function(cormat, trait = trait_vars, env = env_vars) {
  # Use correlation between variables as distance, forces env_vars to front
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
#  env_vars <- c('LMA_g_per_m2',
#                'prec',
#                'tavg',
#                'leaf_rad',
#                'gap')
  env_vars_ <- which(colnames(cormat) %in% env)
  trait_vars_ <- which(colnames(cormat) %in% trait)
  hc$order <- c(hc$order[hc$order %in% env_vars_],
                hc$order[hc$order %in% trait_vars_],
                hc$order[!hc$order %in% env_vars_])
  cormat <-cormat[hc$order, hc$order]
}
  
protein_cats <- c('Rubisco', 
                  'calvin_cycle', 
                  'Photosystems', 
                  'ATP_synthase_chloroplastic',
                  'photorespiration',
                  'protein',
                  'stress',
                  'TCA_org_transformation',
                  'total_protein')

cor_vars <- c(protein_cats, trait_vars, env_vars)

source('scripts/prep_data.R')

prot <- data
prot$calvin_cycle <- prot$calvin_cycle - prot$Rubisco
prot$prec <- log10(prot$prec)
cormat.prot <- cor(na.omit(prot[,cor_vars]))
cormat.prot <- reorder_cormat_env(cormat.prot)
blah_rows <- rownames(cormat.prot)
blah_cols <- colnames(cormat.prot)
cormat.prot <- melt(cormat.prot, na.rm=TRUE)

p <- ggplot(cormat.prot, aes(x = Var1, y = Var2, fill = value))
p <- p  + geom_tile()
p <- p + ggtitle('Correlation heatmap (relative protein amounts)')
p <- p + scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", 
                              midpoint = 0, limit = c(-1,1), space = "Lab", 
                              name="Pearson\nCorrelation")
p <- p + theme_minimal()  + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + coord_fixed()
p


mg_per_m2 = FALSE
moles = TRUE

source('scripts/transformations.R')
source('scripts/prep_data_mg_per_mm2.R')

prot <- data
prot$prec <- log10(prot$prec)
cormat.prot <- cor(na.omit(prot[,cor_vars]))
cormat.prot <- reorder_cormat_env(cormat.prot)
blah_rows <- rownames(cormat.prot)
blah_cols <- colnames(cormat.prot)
cormat.prot <- melt(cormat.prot, na.rm=TRUE)

p <- ggplot(cormat.prot, aes(x = Var1, y = Var2, fill = value))
p <- p  + geom_tile()
p <- p + ggtitle('Correlation heatmap (protein amounts in mol/area)')
p <- p + scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", 
                              midpoint = 0, limit = c(-1,1), space = "Lab", 
                              name="Pearson\nCorrelation")
p <- p + theme_minimal()  + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + coord_fixed()
p

rm(mg_per_m2, moles)



rm(include_photosynthesis,
include_d13C,
include_leaf_N,
include_leaf_P,
include_soil_N,
include_soil_P,
include_chlorophyll)