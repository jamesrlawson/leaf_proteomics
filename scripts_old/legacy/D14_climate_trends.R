# generates heatmaps to look at correlations between environmental variables, and between relative and absolute protein functional category amounts

require(ggplot2)
require(vegan)

source('../scripts/binProteins.R')

## D14 CLIMATE TRENDS ##

sample_locations <- read_csv('../data/sample_locations.csv')
climate <- read_csv('../data/discovery_site_climate.csv')

climate_locs <- merge(sample_locations, climate, by = c('Longitude', 'Latitude'))
climate_locs <- climate_locs[!duplicated(climate_locs$sample),]
climate_locs$Var.4 <- NULL
rm(sample_locations,climate)

climate_locs <- climate_locs[climate_locs$sample %in% unique(protein_stand_D14$sample),]


## correlations ##

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

# climate

cormat.clim <- cor(climate_locs[,c(4:23)])
cormat.clim <- reorder_cormat(cormat.clim)

c <- ggplot(melt(cormat.clim, na.rm=TRUE), aes(x = Var1, y = Var2, fill = value))
c <- c + geom_tile()
c <- c + scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", 
                              midpoint = 0, limit = c(-1,1), space = "Lab", 
                              name="Pearson\nCorrelation")
c <- c + theme_minimal()  + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + coord_fixed()
c

clim.pca <- prcomp(climate_locs[,c(4:22)], scale = TRUE, retx=TRUE, center=TRUE)
clim.pca
plot(clim.pca)
clim.pca$rotation[,1:2]


PC1 <- clim.pca$x[,1]
PC2 <- clim.pca$x[,2]

plot(PC1 ~ protein_stand_D14$PSII)


# relative [protein]

cormat.protRel <- cor(protein_stand_D14[,c(2:28)])
cormat.protRel <- reorder_cormat(cormat.protRel)
cormat.protRel <- melt(cormat.protRel, na.rm=TRUE)

p <- ggplot(cormat.protRel, aes(x = Var1, y = Var2, fill = value))
p <- p  + geom_tile()
p <- p + scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", 
                              midpoint = 0, limit = c(-1,1), space = "Lab", 
                              name="Pearson\nCorrelation")
p <- p + theme_minimal()  + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + coord_fixed()
p


# absolute [protein]

cormat.protAbs <- cor(protein_D14[,c(2:28)])
cormat.protAbs <- reorder_cormat(cormat.protAbs)
cormat.protAbs <- melt(cormat.protAbs, na.rm=TRUE)

a <- ggplot(cormat.protAbs, aes(x = Var1, y = Var2, fill = value))
a <- a  + geom_tile()
a <- a + scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", 
                              midpoint = 0, limit = c(-1,1), space = "Lab", 
                              name="Pearson\nCorrelation")
a <- a + theme_minimal()  + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + coord_fixed()
a




# RDA rel. [protein] ~ climate

climate_locs <- climate_locs[order(climate_locs$sample),]
protein_stand_D14 <- protein_stand_D14[order(protein_stand_D14$sample),]

protein_stand_clim.rda <- rda(protein_stand_D14[,c(2:27)], climate_locs[,c(4:22)], scale=TRUE)
#plot(protein_stand_clim.rda)            
protein_stand_clim.rda
summary(protein_stand_clim.rda)

# RDA abs. [protein] ~ climate

climate_locs <- climate_locs[order(climate_locs$sample),]
protein_D14 <- protein_D14[order(protein_D14$sample),]

protein_clim.rda <- rda(protein_D14[,c(2:28)], climate_locs[,c(4:22)], scale=TRUE)
plot(protein_clim.rda)            
protein_clim.rda
summary(protein_clim.rda)


## LINEAR REGRESSIONS

#relative [protein]
protein_stand_D14_long <- melt(protein_stand_D14)
names(protein_stand_D14_long) <- c('sample', 'bin_arch_name', 'sum')

protein_climate_D14_stand <- merge(climate_locs, protein_stand_D14_long, by = 'sample')
protein_climate_D14_stand$bin_arch_name <- as.character(protein_climate_D14_stand$bin_arch_name)

#absolute [protein]
protein_D14_long <- melt(protein_D14)
names(protein_D14_long) <- c('sample', 'bin_arch_name', 'sum')

protein_climate_D14 <- merge(climate_locs, protein_D14_long, by = 'sample')
protein_climate_D14$bin_arch_name <- as.character(protein_climate_D14$bin_arch_name)


# linear regressions over rel. [protein]

tavg_lms <- regression_output(na.omit(protein_climate_D14_stand), tavg, sum)
tavg_lms <- tavg_lms[order(tavg_lms$R2, decreasing =TRUE),]
tavg_lms
write.csv(tavg_lms, '../output/MAT_models_rel.csv')

prec_lms <- regression_output(na.omit(protein_climate_D14_stand), prec, sum)
prec_lms <- prec_lms[order(prec_lms$R2, decreasing =TRUE),]
prec_lms
write.csv(tavg_lms, '../output/MAP_models_rel.csv')

soilN_lms <- regression_output(na.omit(protein_climate_D14_stand), soilN, sum)
soilN_lms <- soilN_lms[order(soilN_lms$R2, decreasing =TRUE),]
soilN_lms
write.csv(soilN_lms, '../output/soilN_models_rel.csv')

# linear regressions over abs [protein]

tavg_lms <- regression_output(na.omit(protein_climate_D14), tavg, sum)
tavg_lms <- tavg_lms[order(tavg_lms$R2, decreasing =TRUE),]
tavg_lms
write.csv(tavg_lms, '../output/MAT_models_abs.csv')

prec_lms <- regression_output(na.omit(protein_climate_D14), prec, sum)
prec_lms <- prec_lms[order(prec_lms$R2, decreasing =TRUE),]
prec_lms
write.csv(tavg_lms, '../output/MAT_models_abs.csv')



# [protein] separated by leaf age

leaf_age <- read_csv('../data/leaf_age.csv')
leaf_age <- leaf_age[unique(leaf_age$sample) %in% protein_stand_D14$sample,]
protein_stand_D14_age <- merge(leaf_age, protein_stand_D14, by = 'sample')

protein_stand_D14_new <- subset(protein_stand_D14_age, leaf_age == "new")
protein_stand_D14_mid <- subset(protein_stand_D14_age, leaf_age == "mid")
protein_stand_D14_old <- subset(protein_stand_D14_age, leaf_age == "old")


#NEW
cormat.prot_new <- cor(protein_stand_D14_new[,c(3:27)])
cormat.prot_new <- reorder_cormat(cormat.prot_new)
blah_rows <- rownames(cormat.prot_new)
blah_cols <- colnames(cormat.prot_new)
cormat.prot_new <- melt(cormat.prot_new, na.rm=TRUE)

p <- ggplot(cormat.prot_new, aes(x = Var1, y = Var2, fill = value))
p <- p  + geom_tile()
p <- p + ggtitle('New leaves')
p <- p + scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", 
                              midpoint = 0, limit = c(-1,1), space = "Lab", 
                              name="Pearson\nCorrelation")
p <- p + theme_minimal()  + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + coord_fixed()
p

protein_climate_D14_new_stand <- merge(climate_locs, melt(protein_stand_D14_new), by = 'sample')
names(protein_climate_D14_new_stand)[24:25] <- c('bin_arch_name', 'sum')
protein_climate_D14_new_stand$bin_arch_name <- as.character(protein_climate_D14_new_stand$bin_arch_name)


tavg_lms <- regression_output(na.omit(protein_climate_D14_new_stand), tavg, sum)
tavg_lms <- tavg_lms[order(tavg_lms$R2, decreasing =TRUE),]
tavg_lms
write.csv(tavg_lms, '../output/MAT_models_rel_new.csv')

prec_lms <- regression_output(na.omit(protein_climate_D14_new_stand), prec, sum)
prec_lms <- prec_lms[order(prec_lms$R2, decreasing =TRUE),]
prec_lms
write.csv(prec_lms, '../output/MAP_models_rel_new.csv')


#MID
cormat.prot_mid <- cor(protein_stand_D14_mid[,c(3:27)])
#cormat.prot_mid <- reorder_cormat(cormat.prot_mid)
cormat.prot_mid <- cormat.prot_mid[blah_rows,]
cormat.prot_mid <- cormat.prot_mid[,blah_cols]

cormat.prot_mid <- melt(cormat.prot_mid, na.rm=TRUE)

p <- ggplot(cormat.prot_mid, aes(x = Var1, y = Var2, fill = value))
p <- p  + geom_tile()
p <- p + ggtitle('Middle aged leaves')
p <- p + scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", 
                              midpoint = 0, limit = c(-1,1), space = "Lab", 
                              name="Pearson\nCorrelation")
p <- p + theme_minimal()  + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + coord_fixed()
p


protein_climate_D14_mid_stand <- merge(climate_locs, melt(protein_stand_D14_mid), by = 'sample')
names(protein_climate_D14_mid_stand)[24:25] <- c('bin_arch_name', 'sum')
protein_climate_D14_mid_stand$bin_arch_name <- as.character(protein_climate_D14_mid_stand$bin_arch_name)


tavg_lms <- regression_output(na.omit(protein_climate_D14_mid_stand), tavg, sum)
tavg_lms <- tavg_lms[order(tavg_lms$R2, decreasing =TRUE),]
tavg_lms
write.csv(tavg_lms, '../output/MAT_models_rel_mid.csv')

prec_lms <- regression_output(na.omit(protein_climate_D14_mid_stand), prec, sum)
prec_lms <- prec_lms[order(prec_lms$R2, decreasing =TRUE),]
prec_lms
write.csv(prec_lms, '../output/MAP_models_rel_mid.csv')


#OLD
cormat.prot_old <- cor(protein_stand_D14_old[,c(3:27)])
#cormat.prot_old <- reorder_cormat(cormat.prot_old)
cormat.prot_old <- cormat.prot_old[blah_rows,]
cormat.prot_old <- cormat.prot_old[,blah_cols]

cormat.prot_old <- melt(cormat.prot_old, na.rm=TRUE)

p <- ggplot(cormat.prot_old, aes(x = Var1, y = Var2, fill = value))
p <- p  + geom_tile()
p <- p + ggtitle('Old leaves')
p <- p + scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", 
                              midpoint = 0, limit = c(-1,1), space = "Lab", 
                              name="Pearson\nCorrelation")
p <- p + theme_minimal()  + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + coord_fixed()
p

protein_climate_D14_old_stand <- merge(climate_locs, melt(protein_stand_D14_old), by = 'sample')
names(protein_climate_D14_old_stand)[24:25] <- c('bin_arch_name', 'sum')
protein_climate_D14_old_stand$bin_arch_name <- as.character(protein_climate_D14_old_stand$bin_arch_name)


tavg_lms <- regression_output(na.omit(protein_climate_D14_old_stand), tavg, sum)
tavg_lms <- tavg_lms[order(tavg_lms$R2, decreasing =TRUE),]
tavg_lms
write.csv(tavg_lms, '../output/MAT_models_rel_old.csv')

prec_lms <- regression_output(na.omit(protein_climate_D14_old_stand), prec, sum)
prec_lms <- prec_lms[order(prec_lms$R2, decreasing =TRUE),]
prec_lms
write.csv(prec_lms, '../output/MAP_models_rel_old.csv')





