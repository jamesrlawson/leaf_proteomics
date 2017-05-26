source('scripts/functions.R')
source('scripts/transformations.R')

# stacked column graphs

# get list of top 50, 100, 200, 500, 2000 proteins
  # not so simple as ranks will have to be generated separately for each sample
# then get protein amounts for funccats using those proteins
# then bind df's together and create stacked column plots

require(dplyr)

protein_samples_D14 <- read_csv('data/D14_protein_GGLEP-DEDT.csv') # protein amounts calculated using D14 ion library, in avg(GGLEP/DEDT) equivalents

mercator <- read_csv('data/mercator/D14_mercator_20170217.csv')

mercator_names <- read.csv('data/mercator/mercator_names.csv', header=T, stringsAsFactors = F) 
mercator_names <- dplyr::arrange(mercator_names, funccat)

bins <- c(1,5,10,25,50,75,100,200,500,1000,2000)

#bins <- seq(1,2000,by=2)
#bins <- c(1,5,10)

#bins <- c(1,5,10,20,30,40,50)


protein_amounts_rankBins <- vector('list', length(bins)) 
protein_amounts_rankBins_total_protein <- vector('list', length(bins)) 

bins_df <- data.frame()

for(l in 1:length(bins)) {

  n <- bins[l]
  
  my.list <- vector('list', length(ncol(protein_samples_D14)-2))
  my.list_protein <- vector('list', length(ncol(protein_samples_D14)-2))
  
  for(i in 2:(ncol(protein_samples_D14)-2)) {
   
    sample <-  names(protein_samples_D14)[i]  
  
    sub <- cbind(protein_samples_D14[1],protein_samples_D14[i])
    
    sub <- sub[order(-sub[,2]),][1:n,] # order by protein amount and take top 1:n rows
  
    total_protein <- data.frame(cbind(sample,total_protein=sum(sub[2]), n=n))
    
    sub <- getProteinBins(sub, mercator)
    
    func_assigned.list <- vector('list', length(mercator_names$funccat))
    
    func_assigned <- data.frame()
    
    for(j in 1:length(mercator_names$funccat)) {
      
        name <- mercator_names$funccat[j]
        
        proteins <- sub[grep(name, sub$NAME),]
        
        if(nrow(proteins) > 0) {
        proteins$funccat <- mercator_names$funccat[j]
        
        proteins <- distinct(proteins, Protein, .keep_all = TRUE)
        
        } else {
         
          proteins[1,] <- c(NA,NA,NA,name)
          proteins$funccat <- mercator_names$funccat[j]
          
        }
        
        proteins$BINCODE <- NULL
        func_assigned.list[[j]] <- proteins
      
        }
  
    func_assigned <- rbind(func_assigned, do.call(rbind, func_assigned.list))
  
    my.list[[i]] <- func_assigned
    my.list_protein[[i]] <- total_protein 

  }
  
  func_assigned <- my.list[[2]]

  for(k in 3:length(my.list)) {
  
    func_assigned <- full_join(func_assigned, my.list[[k]], all=TRUE)
  
  }

  
  func_assigned[,c(2,5:ncol(func_assigned))] <- lapply(func_assigned[,c(2,5:ncol(func_assigned))], function(x) as.numeric(x))
  
  
    # generate df with summed protein amounts for each category
    
    funccat_sums <- func_assigned %>% group_by(funccat) %>% summarise_at(vars(c(2,5:ncol(func_assigned))), sum, na.rm=TRUE) %>% arrange(funccat)
    funccat_sums <- dplyr::arrange(funccat_sums, funccat)
    funccat_sums$funccat <- dplyr::arrange(mercator_names[mercator_names$funccat %in% funccat_sums$funccat,], funccat)$funccat_rename
    
    # transform df
    
    funccat_sums_t <- t(funccat_sums[,2:ncol(funccat_sums)])
    colnames(funccat_sums_t) <- funccat_sums$funccat
    funccat_sums_t <- as.data.frame(funccat_sums_t)
    funccat_sums_t$sample <- colnames(funccat_sums)[2:ncol(funccat_sums)]
    funccat_sums <- funccat_sums_t
    rm(funccat_sums_t)
    
    # create some custom categories
    
    funccat_sums$PSII_min_LHCII <- funccat_sums$photosystem_II - funccat_sums$LHC_I
    funccat_sums$PSI_min_LHCI <- funccat_sums$photosystem_II - funccat_sums$LHC_II
    funccat_sums$Photosystems <- funccat_sums$photosystem_I + funccat_sums$photosystem_II
    funccat_sums$electron_transport_minATPsynth <- funccat_sums$other_electron_carrier + funccat_sums$cytochrome_b6f
    funccat_sums$Rubisco <- funccat_sums$rubisco_large_subunit + funccat_sums$rubisco_small_subunit
    funccat_sums$redox <- funccat_sums$redox + funccat_sums$glutathione_S_transferases
    funccat_sums$stress <- funccat_sums$stress + funccat_sums$glutathione_S_transferases
    
    funccat_sums$electron_transport <- funccat_sums$electron_transport_minATPsynth + funccat_sums$ATP_synthase_chloroplastic
    funccat_sums$LHC <- funccat_sums$LHC_I + funccat_sums$LHC_II
    funccat_sums$Photosystems_min_LHC <- funccat_sums$Photosystems - funccat_sums$LHC
    funccat_sums$LHCI_per_PSI <- funccat_sums$LHC_I / funccat_sums$PSI_min_LHCI
    funccat_sums$LHCII_per_PSII <- funccat_sums$LHC_II / funccat_sums$PSII_min_LHCII
    funccat_sums$LHC_per_PS <- funccat_sums$LHC / funccat_sums$Photosystems

    funccat_sums$n <- n
    
    protein_amounts_rankBins[[l]] <- funccat_sums
    names(protein_amounts_rankBins)[l] <- n
    
    protein_amounts_rankBins_total_protein[[l]] <- do.call(rbind, my.list_protein)
    
}
  
protein_amounts_rankBins_all <- do.call(rbind, protein_amounts_rankBins)
protein_amounts_rankBins_total_protein_all <- do.call(rbind, protein_amounts_rankBins_total_protein)
protein_amounts_rankBins_total_protein_all[] <- lapply(protein_amounts_rankBins_total_protein_all[], function(x) as.character(x))
protein_amounts_rankBins_total_protein_all[2:3] <- lapply(protein_amounts_rankBins_total_protein_all[2:3], function(x) as.numeric(x))


samp <- protein_amounts_rankBins_all$sample 
protein_amounts_rankBins_all$sample <- NULL
protein_amounts_rankBins_all$sample <- samp

hid <- full_join(protein_amounts_rankBins_all, protein_amounts_rankBins_total_protein_all, by = c('n', 'sample')) %>% 
  select(- LHC_per_PS, - LHCI_per_PSI, - LHCII_per_PSII) %>%
  mutate_at(vars(1:(ncol(.)-3)), funs(. / total_protein))

hid$calvin_cycle <- hid$calvin_cycle - hid$Rubisco

hid <- full_join(select(protein_D14, total_protein, sample), select(hid, total_protein, sample, n), by = 'sample') %>% 
  mutate(total_protein = total_protein.y / total_protein.x) %>%
  select(total_protein, sample, n) %>%
 # full_join(select(hid, -total_protein), by = c('sample', 'n'))
  full_join(hid, by = c('sample', 'n'))

po <- hid$total_protein.y
pa <- hid$total_protein.x
hid$total_protein <- pa
#hid$total_protein_mgPerm2 <- po
rm(po)
rm(pa)
hid$total_protein.x <- NULL
hid$total_protein.y <- NULL

hid <- na.omit(hid)

### kruft ##

hid_other <- bla[names(bla) %in%c("n",                                              
                                  "cell_wall",                                     
                                  "DNA",                                    
                                  "glycolysis",                                    
                                  "hormone_metabolism",                          
                                  "lipid_metabolism",                              
                                  "mitochondrial_electron_transport_ATP_synthesis", 
                                  "redox",                                         
                                  "RNA",                                            
                                  "secondary_metabolism",                          
                                  "signalling",                                     
                                  "stress",                                        
                                  "TCA_org_transformation",                         
                                  "CHO_metabolism",     
                                  "polyamine_metabolism", 
                                  "electron_transport_minATPsynth",
                                  "total_protein")]
hid_other$other <- rowSums(hid_other[,c(2:(ncol(hid_other)-1))])
#hid$other <- hid_other$other

### kruft ###

hid <- hid[,c('n', 'total_protein', 'photorespiration','calvin_cycle', 'Rubisco', 'Photosystems', 'ATP_synthase_chloroplastic', 'protein', 'stress')]
hid$other <- 1 - rowSums(hid[,3:ncol(hid)])



write_csv(hid, 'output/protein_amounts_rankBins.csv')

bla <- group_by(hid, n) %>% summarise_at(vars(2:ncol(hid)), mean, na.rm=TRUE)

bla_SE <- group_by(hid, n) %>% summarise_at(vars(2:ncol(hid)), SE)

#blai <- bla[,c('n', 'total_protein', 'other', 'photorespiration','calvin_cycle', 'Rubisco', 'Photosystems_min_LHC', 'LHC', 'electron_transport_minATPsynth', 'ATP_synthase_chloroplastic', 'protein', 'stress')]

#blai_SE <- bla_SE[,c('n', 'total_protein',  'other', 'photorespiration', 'calvin_cycle', 'Rubisco', 'Photosystems_min_LHC', 'LHC', 'electron_transport_minATPsynth', 'ATP_synthase_chloroplastic', 'protein', 'stress')]

blai <- bla[,c('n', 'total_protein', 'other', 'photorespiration','calvin_cycle', 'Rubisco', 'Photosystems', 'ATP_synthase_chloroplastic', 'protein', 'stress')]

blai_SE <- bla_SE[,c('n', 'total_protein',  'other', 'photorespiration', 'calvin_cycle', 'Rubisco', 'Photosystems', 'ATP_synthase_chloroplastic', 'protein', 'stress')]


bax <- gather(blai, key = 'funccat', value = 'protein_amount', -n) %>%
  #filter(funccat %in% c('total_protein','photorespiration','calvin_cycle', 'Rubisco', 'Photosystems_min_LHC', 'LHC', 'electron_transport_minATP_synthesis', 'ATP_synthase_chloroplastic', 'protein', 'stress')) 
  filter(!funccat %in% c('n'))

bax_SE <- gather(blai_SE, key = 'funccat', value = 'protein_SE', -n) %>%
  #filter(funccat %in% c('total_protein', 'photorespiration', 'calvin_cycle', 'Rubisco', 'Photosystems_min_LHC', 'LHC', 'electron_transport_minATP_synthesis', 'ATP_synthase_chloroplastic', 'protein', 'stress')) 
  filter(!funccat %in% c('n'))

bob <- full_join(bax, bax_SE, by = c('n', 'funccat'))


bob$funccat <- factor(bob$funccat, levels = c('total_protein', 
                                              'Rubisco', 
                                              'calvin_cycle', 
                                              'Photosystems', 
                                              'ATP_synthase_chloroplastic',
                                              'photorespiration',
                                              'protein',
                                              'stress',
                                              'other'))
levels(bob$funccat)

require(ggplot2)
cumsumplot <- ggplot(bob, aes(y = protein_amount, x =n, colour=funccat, shape = funccat, order=as.numeric(funccat))) + geom_point() + 
  geom_smooth(method='loess', formula = y ~ log(x), se=FALSE, aes(linetype = funccat)) 
#  geom_smooth(method='loess', se=FALSE) +
#  geom_ribbon(aes(ymin = protein_amount - protein_SE, ymax = protein_amount + protein_SE, fill=funccat), alpha = 0.5) 
cumsumplot <- cumsumplot + theme_classic()
cumsumplot <- cumsumplot + theme(text = element_text(size = 17))
cumsumplot <- cumsumplot + xlim(0,500) + xlab('Protein rank abundance') + ylab('Proportion of protein accounted for')
cumsumplot

require(ggplot2)
cumsumplot <- ggplot(bob[!bob$funccat %in% 'total_protein',], aes(y = protein_amount, x =n, colour = funccat, fill=funccat, shape = funccat, order=as.numeric(funccat))) + geom_bar(stat='identity', alpha=0.5) + 
  geom_smooth(method='loess', formula = y ~ log(x), se=FALSE, aes(linetype = funccat)) 
#  geom_smooth(method='loess', se=FALSE) +
#  geom_ribbon(aes(ymin = protein_amount - protein_SE, ymax = protein_amount + protein_SE, fill=funccat), alpha = 0.5) 
cumsumplot <- cumsumplot + theme_classic()
cumsumplot <- cumsumplot + theme(text = element_text(size = 17))
cumsumplot <- cumsumplot + xlim(0,500) +xlab('Protein rank abundance') + ylab('Proportion of protein accounted for')
cumsumplot <- cumsumplot + scale_x_continuous(breaks = unique(bob$n))
cumsumplot


write_csv(bla, 'output/cumsum_data.csv')

group_by(filter(bob, funccat != 'total_protein'), n) %>% dplyr::summarise(check = sum(protein_amount))

