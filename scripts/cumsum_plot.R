# cumulative protein abundance plots

source('scripts/binProteins.R')


all_protein <- data.frame()
my.list <- vector("list", length(2:(ncol(protein_samples_D14)-2)))

for(i in 2:(ncol(protein_samples_D14)-2)) {

sample <- names(protein_samples_D14)[i]  
  
out <- select(protein_samples_D14, Protein, get(sample)) %>% 
         arrange(desc(get(sample))) %>%
         mutate(cumsum = cumsum(get(sample)), rank = seq(1, nrow(.), by = 1)) 
names(out)[2] <- 'amount'
out$sample <- sample

my.list[[i-1]] <- out

}

all_protein <- rbind(all_protein, do.call(rbind, my.list))

blax <- all_protein %>% dplyr::group_by(rank) %>% dplyr::summarise(mean_cumsum = mean(cumsum), SE_cumsum = SE(cumsum))
blax$maxcumsum <- max(blax$mean_cumsum)

ggplot(blax, aes(x = rank, y = mean_cumsum/maxcumsum)) + geom_point(size = 0.1, alpha = 0.1) + 
  geom_ribbon(aes(ymin=mean_cumsum/maxcumsum - SE_cumsum/maxcumsum, ymax = mean_cumsum/maxcumsum + SE_cumsum/maxcumsum), alpha = 0.3) 




## cumsums for functional categories

protein_samples <- protein_samples_D14

protein_samples$bin_arch <- NA

for(i in 1:length(bin_arch.list)) {
  
  x <-  bin_arch.list[[i]][1]
  rows1 <- nrow(protein_samples[grep(bin_arch.list[[i]][1], protein_samples$BINCODE, fixed=TRUE),])
  rows2 <- nrow(protein_samples[grep(bin_arch.list[[i]][2], protein_samples$BINCODE, fixed=TRUE),])
  
  if(nrow(protein_samples[grep(bin_arch.list[[i]][1], protein_samples$BINCODE, fixed=TRUE),]) > 0) {
    protein_samples[grep(bin_arch.list[[i]][1], protein_samples$BINCODE, fixed=TRUE),]$bin_arch <- bin_arch.list[[i]][1]
  }
  
  if(nrow(protein_samples[grep(bin_arch.list[[i]][2], protein_samples$BINCODE, fixed=TRUE),]) > 0) {
    protein_samples[grep(bin_arch.list[[i]][2], protein_samples$BINCODE, fixed=TRUE),]$bin_arch <- bin_arch.list[[i]][1]
    
  } 
  
}



for(i in 1:length(unique(all_protein$sample))) {
  
  sample_subset <- all_protein[all_protein$sample %in% all_protein$sample[i],]
  
  func_group_prots <- protein_samples[grepl('_1.1', protein_samples$bin_arch),]$Protein

  for(i in 1:max(all_protein$rank))
  
    rank_group <- c(1:i) # I don't think I even need the rank group, I can just redo the cumsums for the intersect(sample_subset$Protein , func_group_prots)
                         # check out the rank_subset output df
    rank_subset_proteins <- intersect(sample_subset[sample_subset$rank %in% rank_group,]$Protein, func_group_prots)
    
    rank_subset <- sample_subset[sample_subset$Protein %in% rank_subset_proteins,]
  
    rank_subset$cumsum <- cumsum(rank_subset$amount)
    
  my.list[[i-1]] <- rank_subset
  
}




# what are the most abundant functional groups?

jah <- gather(protein_D14, sample) 
names(jah)[2] <- 'func_cat'
jahx <- jah %>% dplyr::group_by(func_cat) %>% dplyr::summarise(func_sum = sum(value)) %>% arrange(desc(func_sum))



