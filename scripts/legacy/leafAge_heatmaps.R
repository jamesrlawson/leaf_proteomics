# heatmaps looking at correlations between protein amounts separated by leaf age

source('../scripts/transformations.R')

# rel. [protein] separated by leaf age

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


# abs. [protein] separated by leaf age

leaf_age <- read_csv('../data/leaf_age.csv')
leaf_age <- leaf_age[unique(leaf_age$sample) %in% protein_D14$sample,]
protein_D14_age <- merge(leaf_age, protein_D14, by = 'sample')

protein_D14_new <- subset(protein_D14_age, leaf_age == "new")
protein_D14_mid <- subset(protein_D14_age, leaf_age == "mid")
protein_D14_old <- subset(protein_D14_age, leaf_age == "old")

#NEW
cormat.prot_new <- cor(protein_D14_new[,c(3:27)])
cormat.prot_new <- reorder_cormat(cormat.prot_new)
blah_rows <- rownames(cormat.prot_new)
blah_cols <- colnames(cormat.prot_new)
cormat.prot_new <- melt(cormat.prot_new, na.rm=TRUE)

p <- ggplot(cormat.prot_new, aes(x = Var1, y = Var2, fill = value))
p <- p  + geom_tile()
p <- p + ggtitle('New leaves')
p <- p + scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", 
                              midpoint = 0, limit = c(-1,1), space = "Lab", 
                              name="Pearson\nCorrelsation")
p <- p + theme_minimal()  + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + coord_fixed()
p

#MID
cormat.prot_mid <- cor(protein_D14_mid[,c(3:27)])
#cormat.prot_mid <- reorder_cormat(cormat.prot_mid)
cormat.prot_mid <- cormat.prot_mid[blah_rows,]
cormat.prot_mid <- cormat.prot_mid[,blah_cols]

cormat.prot_mid <- melt(cormat.prot_mid, na.rm=TRUE)

p <- ggplot(cormat.prot_mid, aes(x = Var1, y = Var2, fill = value))
p <- p  + geom_tile()
p <- p + ggtitle('Middle aged leaves')
p <- p + scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", 
                              midpoint = 0, limit = c(-1,1), space = "Lab", 
                              name="Pearson\nCorrelsation")
p <- p + theme_minimal()  + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + coord_fixed()
p


#OLD
cormat.prot_old <- cor(protein_D14_old[,c(3:27)])
#cormat.prot_old <- reorder_cormat(cormat.prot_old)
cormat.prot_old <- cormat.prot_old[blah_rows,]
cormat.prot_old <- cormat.prot_old[,blah_cols]

cormat.prot_old <- melt(cormat.prot_old, na.rm=TRUE)

p <- ggplot(cormat.prot_old, aes(x = Var1, y = Var2, fill = value))
p <- p  + geom_tile()
p <- p + ggtitle('Old leaves')
p <- p + scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", 
                              midpoint = 0, limit = c(-1,1), space = "Lab", 
                              name="Pearson\nCorrelsation")
p <- p + theme_minimal()  + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + coord_fixed()
p


