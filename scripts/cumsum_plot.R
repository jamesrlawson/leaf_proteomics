# cumulative protein abundance plots

source('scripts/binProteins.R')

require(ggplot2)

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
blax$proportion <- blax$mean_cumsum/blax$maxcumsum

cumsumplot <- ggplot(blax, aes(x = rank, y = mean_cumsum/maxcumsum)) + geom_point(size = 0.1, alpha = 0.1) + 
 # geom_vline(xintercept = 469, size = 0.5, alpha = 0.5) +
#  geom_hline(yintercept = 0.9, size = 0.5, alpha = 0.5) +
  geom_vline(xintercept = 100, size = 0.5, alpha = 0.5) +
  geom_hline(yintercept = 0.725, size = 0.5, alpha = 0.5) +
  
  
  geom_ribbon(aes(ymin=mean_cumsum/maxcumsum - SE_cumsum/maxcumsum, ymax = mean_cumsum/maxcumsum + SE_cumsum/maxcumsum), alpha = 0.3, fill = 'blue') +
  xlab('Protein abundance rank') + ylab('Cumulative proportion of protein') +
  theme_classic() +
  theme(text = element_text(size = 17))

cumsumplot



## cumsums for functional categories

# find the ranks of proteins for a given functional group?
# then for each rank represented by the func group of choice, find the proportion of maxcumsum accounted for by the cumsum of yer functional group

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

# get traces of proportions accounted for by defined functional categories vs rank

get_funccat_proportions <- function(func_group_prots) {

my.list <- vector("list", length(unique(all_protein$sample)))
func_cat_proportions <- data.frame()

for(i in 1:length(unique(all_protein$sample))) {
  
  sample_subset <- all_protein[all_protein$sample %in% unique(all_protein$sample)[i],]
  
#  func_group_prots <- c(protein_samples[grepl('_1.1.1', protein_samples$bin_arch),]$Protein, protein_samples[grepl('_1.1.2', protein_samples$bin_arch),]$Protein) # photosystems
  
#  func_group_prots <- protein_samples[grepl('_1.1', protein_samples$bin_arch),]$Protein # Light reactions
  
#  func_group_prots <- protein_samples[grepl('_1.3', protein_samples$bin_arch),]$Protein # Calvin_cycle
  
  subset_proteins <- intersect(sample_subset$Protein, func_group_prots)
  
  bla <- sample_subset[sample_subset$Protein %in% subset_proteins,]
  
  bla$cumsum <- cumsum(bla$amount)
  
 # bla$total_cumsum <- blax[blax$rank %in% bla$rank,]$mean_cumsum
  
  all_protein_samplesub <- all_protein[all_protein$sample %in% unique(all_protein$sample)[i],]
  
  bla$total_cumsum <- all_protein_samplesub[all_protein_samplesub$rank %in% bla$rank,]$cumsum
  
  bla$proportion_of_total <- bla$cumsum / bla$total_cumsum
  
 # plot(proportion_of_total ~ rank, bla)
  
  #my.list[[i-1]] <- rank_subset
  my.list[[i]] <- bla
  
}

func_cat_proportions <- rbind(func_cat_proportions, do.call(rbind, my.list))
return(func_cat_proportions)

}


func_group_prots_photosytems <- c(protein_samples[grepl('_1.1.1', protein_samples$bin_arch),]$Protein, protein_samples[grepl('_1.1.2', protein_samples$bin_arch),]$Protein) # photosystems
func_cat_proportions_photosystems <- get_funccat_proportions(func_group_prots_photosytems)

func_group_prots_LRs <- protein_samples[grepl('_1.1', protein_samples$bin_arch),]$Protein # Light reactions
func_cat_proportions_LRs <- get_funccat_proportions(func_group_prots_photosytems)

func_group_prots_Calvin <- protein_samples[grepl('_1.3', protein_samples$bin_arch),]$Protein # Calvin_cycle
func_cat_proportions_Calvin <- get_funccat_proportions(func_group_prots_Calvin)

func_group_prots_photoresp <- protein_samples[grepl('_1.2', protein_samples$bin_arch),]$Protein # photorespiration
func_cat_proportions_photoresp <- get_funccat_proportions(func_group_prots_photoresp)

func_group_prots_electron_transport <- c(protein_samples[grepl('_1.1.3', protein_samples$bin_arch),]$Protein, protein_samples[grepl('_1.1.5', protein_samples$bin_arch),]$Protein, protein_samples[grepl('_1.1.6', protein_samples$bin_arch),]$Protein) # photosynthetic electron transport
func_cat_proportions_electron_transport <- get_funccat_proportions(func_group_prots_electron_transport)

func_group_prots_chloroplasticATPsynth <- protein_samples[grepl('_1.1.4', protein_samples$bin_arch),]$Protein # chloroplastic ATP synthase
func_cat_proportions_chloroplasticATPsynth <- get_funccat_proportions(func_group_prots_chloroplasticATPsynth)

func_group_prots_protein <- protein_samples[grepl('_29', protein_samples$bin_arch),]$Protein # protein synthesis & folding
func_cat_proportions_protein <- get_funccat_proportions(func_group_prots_protein)


cumsumplot <- ggplot(blax, aes(x = rank, y = mean_cumsum/maxcumsum)) + geom_point(size = 0.1, alpha = 0.1) + 
  geom_ribbon(aes(ymin=mean_cumsum/maxcumsum - SE_cumsum/maxcumsum, ymax = mean_cumsum/maxcumsum + SE_cumsum/maxcumsum), alpha = 0.3) 
cumsumplot <- cumsumplot + geom_point(data = func_cat_proportions_LRs, aes(x = rank, y = proportion_of_total), size = 0.1, alpha = 0.1, colour = 'green')
cumsumplot <- cumsumplot + geom_point(data = func_cat_proportions_Calvin, aes(x = rank, y = proportion_of_total), size = 0.1, alpha = 0.1, colour = 'red')
cumsumplot <- cumsumplot + geom_point(data = func_cat_proportions_photoresp, aes(x = rank, y = proportion_of_total), size = 0.1, alpha = 0.1, colour = 'blue')

cumsumplot <- cumsumplot + xlim(0,200)

cumsumplot


LRs_proportions_summary <- func_cat_proportions_LRs %>% dplyr::group_by(rank) %>% dplyr::summarise(mean_proportion = mean(proportion_of_total), SE_proportion = SE(proportion_of_total))
LRs_proportions_summary$maxcumsum <- max(LRs_proportions_summary$mean_proportion)
LRs_proportions_summary$ID = 'LR'

Calvin_proportions_summary <- func_cat_proportions_Calvin %>% dplyr::group_by(rank) %>% dplyr::summarise(mean_proportion = mean(proportion_of_total), SE_proportion = SE(proportion_of_total))
Calvin_proportions_summary$maxcumsum <- max(Calvin_proportions_summary$mean_proportion)
Calvin_proportions_summary$ID = 'calv'

photoresp_proportions_summary <- func_cat_proportions_photoresp %>% dplyr::group_by(rank) %>% dplyr::summarise(mean_proportion = mean(proportion_of_total), SE_proportion = SE(proportion_of_total))
photoresp_proportions_summary$maxcumsum <- max(photoresp_proportions_summary$mean_proportion)
photoresp_proportions_summary$ID <- 'photoresp'

photosystems_proportions_summary <- func_cat_proportions_photosystems %>% dplyr::group_by(rank) %>% dplyr::summarise(mean_proportion = mean(proportion_of_total), SE_proportion = SE(proportion_of_total))
photosystems_proportions_summary$maxcumsum <- max(photosystems_proportions_summary$mean_proportion)
photosystems_proportions_summary$ID = 'photosys'

electron_transport_proportions_summary <- func_cat_proportions_electron_transport %>% dplyr::group_by(rank) %>% dplyr::summarise(mean_proportion = mean(proportion_of_total), SE_proportion = SE(proportion_of_total))
electron_transport_proportions_summary$maxcumsum <- max(electron_transport_proportions_summary$mean_proportion)
electron_transport_proportions_summary$ID = 'Etrans'

chloroplasticATPsynth_proportions_summary <- func_cat_proportions_chloroplasticATPsynth %>% dplyr::group_by(rank) %>% dplyr::summarise(mean_proportion = mean(proportion_of_total), SE_proportion = SE(proportion_of_total))
chloroplasticATPsynth_proportions_summary$maxcumsum <- max(chloroplasticATPsynth_proportions_summary$mean_proportion)
chloroplasticATPsynth_proportions_summary$ID = 'ATPsynth'



blox <- blax 

blox$mean_cumsum <- blax$mean_cumsum/blax$maxcumsum
blox$SE_cumsum <- blax$SE_cumsum/blax$maxcumsum

blox$ID <- 'total'

names(blox) <- names(photoresp_proportions_summary)

proportions_summary <- rbind(photosystems_proportions_summary, chloroplasticATPsynth_proportions_summary, electron_transport_proportions_summary, Calvin_proportions_summary, photoresp_proportions_summary, blox)
#proportions_summary <- rbind(electron_transport_proportions_summary,photoresp_proportions_summary,blox)

#proportions_summary[is.na(proportions_summary$SE_proportion),]$SE_proportion <- 0.001


cumsumplot <- ggplot(proportions_summary, aes(x = rank, y = mean_proportion, fill = ID, colour = ID)) + geom_point(size = 0.1, alpha = 0.2) + 
  geom_ribbon(aes(ymin=mean_proportion - SE_proportion, ymax = mean_proportion + SE_proportion), alpha = 0.5)
cumsumplot <- cumsumplot + xlim(0,200) 
cumsumplot <- cumsumplot + scale_colour_manual(values = c("blue", "red", "orange", "purple", "green", "gray"))
cumsumplot <- cumsumplot + scale_fill_manual(values = c("blue", "red", "orange", "purple", "green", "gray"))
cumsumplot <- cumsumplot + theme_classic()
cumsumplot <- cumsumplot + theme(text = element_text(size = 17))


cumsumplot





# mean and CV protein abudnances
require(dplyr) 
require(readr)
meanCVprotein_abunds <- protein_D14_long %>% group_by(bin_arch_name) %>%
                        dplyr::summarise(mean_abund = mean(sum, na.rm=TRUE), CV_abund = CV(sum), sd_abund = sd(sum, na.rm=TRUE))
write_csv(meanCVprotein_abunds, 'output/protein_abunds.csv')

# what are the most abundant functional groups?

jah <- gather(protein_D14, sample) 
names(jah)[2] <- 'func_cat'
jahx <- jah %>% dplyr::group_by(func_cat) %>% dplyr::summarise(func_sum = sum(value)) %>% arrange(desc(func_sum))



