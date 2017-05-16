require(broom)

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

### relative amounts ###

source('scripts/prep_data.R')

prot <- data
prot$calvin_cycle <- prot$calvin_cycle - prot$Rubisco
prot$prec <- log10(prot$prec)

#zap <- na.omit(prot[,cor_vars])
zap <- prot[,cor_vars]

cor.list <- list(NULL)
length(cor.list) <- length(zap)^2

for(i in seq_along(zap)){
  for(j in seq_along(zap)){
    cor.list[[(i-1)*length(zap) + j]] <- 
      glance(cor.test(zap[, i], zap[, j])) %>%
      mutate(x = names(zap)[i],
             y = names(zap)[j])
  }
}

paz <- bind_rows(cor.list)

names(paz)[names(paz) %in% 'x'] <- 'Var1'
names(paz)[names(paz) %in% 'y'] <- 'Var2'



cormat.prot_rel <- cor(prot[,cor_vars], use = 'pairwise.complete.obs')

cormat.prot <- reorder_cormat_env(cormat.prot_rel)
cormat.prot <- melt(cormat.prot, na.rm=TRUE)

paz$Var1 <- factor(paz$Var1, levels = levels(cormat.prot$Var1))
paz$Var2 <- factor(paz$Var2, levels = levels(cormat.prot$Var2))

blax <- full_join(cormat.prot, paz, by = c('Var1', 'Var2')) %>%
        distinct(Var1, Var2, .keep_all=TRUE)

blax$p.value <- p.adjust(blax$p.value , method = 'BH')


p <- ggplot(blax, aes(x = Var1, y = Var2))
p <- p  + geom_tile(data = subset(blax, p.value < 0.05), aes(fill = estimate))
p <- p + ggtitle('Correlation heatmap (relative protein amounts)')
p <- p + scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", 
                              midpoint = 0, limit = c(-1,1), space = "Lab", 
                              name="Pearson\nCorrelation")
p <- p + theme_minimal()  + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + coord_fixed()
p


### absolute molar amounts ###



mg_per_m2 = FALSE
moles = TRUE

source('scripts/transformations.R')
source('scripts/prep_data_mg_per_mm2.R')

prot <- data
prot$calvin_cycle <- prot$calvin_cycle - prot$Rubisco
prot$prec <- log10(prot$prec)

#zap <- na.omit(prot[,cor_vars])
zap <- prot[,cor_vars]

cor.list <- list(NULL)
length(cor.list) <- length(zap)^2

for(i in seq_along(zap)){
  for(j in seq_along(zap)){
    cor.list[[(i-1)*length(zap) + j]] <- 
      glance(cor.test(zap[, i], zap[, j])) %>%
      mutate(x = names(zap)[i],
             y = names(zap)[j])
  }
}

paz <- bind_rows(cor.list)

names(paz)[names(paz) %in% 'x'] <- 'Var1'
names(paz)[names(paz) %in% 'y'] <- 'Var2'


cormat.prot_abs <- cor(prot[,cor_vars], use = 'pairwise.complete.obs')

cormat.prot <- reorder_cormat_env(cormat.prot_abs)
cormat.prot <- melt(cormat.prot, na.rm=TRUE)

paz$Var1 <- factor(paz$Var1, levels = levels(cormat.prot$Var1))
paz$Var2 <- factor(paz$Var2, levels = levels(cormat.prot$Var2))

blax <- full_join(cormat.prot, paz, by = c('Var1', 'Var2')) %>%
  distinct(Var1, Var2, .keep_all=TRUE)

blax$p.value <- p.adjust(blax$p.value , method = 'BH')


p <- ggplot(blax, aes(x = Var1, y = Var2))
p <- p  + geom_tile(data = subset(blax, p.value < 0.05), aes(fill = estimate))
p <- p + ggtitle('Correlation heatmap (protein amounts in mol/area)')
p <- p + scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", 
                              midpoint = 0, limit = c(-1,1), space = "Lab", 
                              name="Pearson\nCorrelation")
p <- p + theme_minimal()  + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + coord_fixed()
p

rm(mg_per_m2, moles)




### both ###


source('scripts/prep_data.R')

prot <- data
prot$calvin_cycle <- prot$calvin_cycle - prot$Rubisco
prot$prec <- log10(prot$prec)


if(inc_photosynthesis) {
  
  # alter trait_vars and cor_vars to include gas exchange vars
  if(!any(trait_vars %in% 'photo_max')) {
    trait_vars <- c(trait_vars, 'photo_max', 'Cond') 
  }
  cor_vars <- c(protein_cats, trait_vars, env_vars)
  
  # add place holder columns for Cond and photo_max
  
  prot$Cond <- runif(nrow(prot))
  prot$photo_max <- runif(nrow(prot))

  cormat.prot_rel <- Hmisc::rcorr(as.matrix(prot[,cor_vars]), type = 'pearson')
  
  include_photosynthesis = TRUE
  
  source('scripts/transformations.R')
  source('scripts/prep_data.R')
  
  rm(include_photosynthesis)
  
  prot_photo <- data
  prot_photo$calvin_cycle <- prot_photo$calvin_cycle - prot_photo$Rubisco
  prot_photo$prec <- log10(prot_photo$prec)
  
  cormat.prot_photo <- Hmisc::rcorr(as.matrix(prot_photo[,cor_vars]), type = 'pearson')
  cormat.prot_rel$r[,which(rownames(cormat.prot_rel$r) %in% c('photo_max', 'Cond'))] <- cormat.prot_photo$r[,which(rownames(cormat.prot_photo$r) %in% c('photo_max', 'Cond'))]
  cormat.prot_rel$P[,which(rownames(cormat.prot_rel$P) %in% c('photo_max', 'Cond'))] <- cormat.prot_photo$P[,which(rownames(cormat.prot_photo$P) %in% c('photo_max', 'Cond'))]
  
}



source('scripts/prep_data_mg_per_mm2.R')
prot <- data
prot$calvin_cycle <- prot$calvin_cycle - prot$Rubisco
prot$prec <- log10(prot$prec)


if(inc_photosynthesis) {
  
  # alter trait_vars and cor_vars to include gas exchange vars
  if(!any(trait_vars %in% 'photo_max')) {
    trait_vars <- c(trait_vars, 'photo_max', 'Cond') 
  }
  cor_vars <- c(protein_cats, trait_vars, env_vars)
  
  # add place holder columns for Cond and photo_max
  
  prot$Cond <- runif(nrow(prot))
  prot$photo_max <- runif(nrow(prot))
  
  cormat.prot_abs <- Hmisc::rcorr(as.matrix(prot[,cor_vars]), type = 'pearson')
  
  include_photosynthesis = TRUE
  
  source('scripts/transformations.R')
  source('scripts/prep_data_mg_per_mm2.R')
  
  rm(include_photosynthesis)
  
  prot_photo <- data
  prot_photo$calvin_cycle <- prot_photo$calvin_cycle - prot_photo$Rubisco
  prot_photo$prec <- log10(prot_photo$prec)
  
  cormat.prot_photo <- Hmisc::rcorr(as.matrix(prot_photo[,cor_vars]), type = 'pearson')
  cormat.prot_abs$r[,which(rownames(cormat.prot_abs$r) %in% c('photo_max', 'Cond'))] <- cormat.prot_photo$r[,which(rownames(cormat.prot_photo$r) %in% c('photo_max', 'Cond'))]
  cormat.prot_abs$P[,which(rownames(cormat.prot_abs$P) %in% c('photo_max', 'Cond'))] <- cormat.prot_photo$P[,which(rownames(cormat.prot_photo$P) %in% c('photo_max', 'Cond'))]
  
}



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

