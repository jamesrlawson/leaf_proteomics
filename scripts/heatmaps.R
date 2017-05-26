# heatmaps looking at correlations between protein funccats

source('scripts/transformations.R')

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}



# for imposing mapman structure...

mercator <- read_csv('data/proteomics_data/mercator/euc/D14_mercator_20170217.csv')

protein_samples_D14 <- getProteinBins(protein_samples_D14, mercator)

# then import the names of categories we're interested in from mercator_names* and use grep to find all proteins associated with those categories
# put the results in instances of a list
# mercator_names.csv contains search terms for protein categories. These searches are run on mercator$NAMES. 
# search terms must be unique to the functional category to avoid non-target returns. 
# for example, a search for 'photosystem I' will also pick up proteins from 'photosystem II' - to avoid this we search for 'photosystem I\.'
# this works because all instances of proteins within 'photosystem I' are actually within subcategories. We'd miss some returns if there were proteins in the upper 'photosystem I' category.
# N.B. the '\' is an 'escape' and must be used because .'s are special in regular expressions and mean 'anything'. By using the escape we will actually search for the character '.'
# search terms for top level categories can be made unique by using ' in front

mercator_names <- read.csv('data/proteomics_data/mercator/mercator_names.csv', header=T, stringsAsFactors = F) 
mercator_names <- arrange(mercator_names, funccat)

func_assigned.list <- vector('list', length(mercator_names$funccat))

func_assigned <- data.frame()

for(i in 1:length(mercator_names$funccat)) {
  
  name <- mercator_names$funccat[i]
  
  proteins <- protein_samples_D14[grep(name, protein_samples_D14$NAME),]
  
  proteins$funccat <- mercator_names$funccat[i]
  
  proteins <- distinct(proteins, Protein, .keep_all = TRUE)
  
  func_assigned.list[[i]] <- proteins
  
}

func_assigned <- rbind(func_assigned, do.call(rbind, func_assigned.list))
rm(func_assigned.list)







source('scripts/prep_data.R')

prot <- data
prot$calvin_cycle <- prot$calvin_cycle - prot$Rubisco

cormat.prot <- cor(prot[,c('Rubisco', 
                           'calvin_cycle', 
                           'Photosystems', 
                           'ATP_synthase_chloroplastic',
                           'photorespiration',
                           'protein',
                           'stress',
                           'TCA_org_transformation')])
cormat.prot <- reorder_cormat(cormat.prot)
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
prot$calvin_cycle <- prot$calvin_cycle - prot$Rubisco

cormat.prot <- cor(prot[,c('Rubisco', 
                           'calvin_cycle', 
                           'Photosystems', 
                           'ATP_synthase_chloroplastic',
                           'photorespiration',
                           'protein',
                           'stress',
                           'TCA_org_transformation')])
cormat.prot <- reorder_cormat(cormat.prot)
blah_rows <- rownames(cormat.prot)
blah_cols <- colnames(cormat.prot)
cormat.prot <- melt(cormat.prot, na.rm=TRUE)

p <- ggplot(cormat.prot, aes(x = Var1, y = Var2, fill = value))
p <- p  + geom_tile()
p <- p + ggtitle('Correlation heatmap (absolute protein amounts)')
p <- p + scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", 
                              midpoint = 0, limit = c(-1,1), space = "Lab", 
                              name="Pearson\nCorrelation")
p <- p + theme_minimal()  + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + coord_fixed()
p

rm(mg_per_m2, moles)


protein_samples_D14 <- read_csv('data/proteomics_data/proteomics/derived/euc/D14_protein_GGLEP-DEDT.csv')
allprot_prots <- protein_samples_D14$Protein
allprot <- t(protein_samples_D14[,2:ncol(protein_samples_D14)])
colnames(allprot) <- allprot_prots
cormat.allprot <- cor(allprot)
cormat.allprot <- reorder_cormat(cormat.allprot)
cormat.allprot <- melt(cormat.allprot, na.rm=TRUE)

p <- ggplot(cormat.allprot, aes(x = Var1, y = Var2, fill = value))
p <- p + geom_tile()
p <- p + ggtitle('Correlation heatmap (absolute protein amounts)')
p <- p + scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", 
                              midpoint = 0, limit = c(-1,1), space = "Lab", 
                              name="Pearson\nCorrelation")
p <- p + theme_minimal()  + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + coord_fixed()
p