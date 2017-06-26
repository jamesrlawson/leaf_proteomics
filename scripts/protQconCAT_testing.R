pairwise_avg_top2 <- function(bla) {

#bla <- subset(qconcat_ion_areas, Peptide_clean == "AVDSLVPIGR")

bla <- bla[order(bla$`Fragment MZ`),]

if(length(unique(bla$Peptide_clean)) != length(unique(bla$Peptide))) {

# could modify labelled fragment MZ's so that they're the same as unlabelled fragent MZ's... remove label amount and then round everything 

#label <- as.numeric(max(gsub('[^0-9]', "", bla$Peptide)))

bla$label <- as.numeric(gsub('[^0-9]', "", bla$Peptide))

if(any(bla$label %in% '810')) {
  bla[bla$label %in% '810',]$label <- 18
}

bla[is.na(bla$label),]$label <- 0

bla$Fragment_MZ_charge <- bla$`Fragment MZ` * bla$`Fragment Charge`

bla$Fragment_MZ_charge_adj <- bla$Fragment_MZ_charge

bla$Fragment_MZ_charge_adj <- bla$Fragment_MZ_charge_adj - bla$label

bla$Fragment_MZ_charge_adj <- round(bla$Fragment_MZ_charge_adj, 0)

# so now i can sort or run operations by fragment mz / charge

# now need to subset dataset to only keep paired ions??

# na.omit it first to get rid of rare ions

bla <- na.omit(bla)

# how do i get top2 pairs? get a pairwise average ion area?

bla <- dplyr::group_by(bla, Fragment_MZ_charge_adj) %>% dplyr::summarise(pairwise_avg_ion_mean = mean(ion_means),
                                                                                           labelled_unlabelled = length(ion_means)) %>%
full_join(bla) %>% 
  filter(labelled_unlabelled == 2) %>%
  dplyr::group_by(Peptide) %>% top_n(2, pairwise_avg_ion_mean) %>%
  dplyr::group_by(Peptide) %>% dplyr::mutate(ions_per_peptide = length(Peptide))


#browser()

return(bla)

}

}


x <- ddply(qconcat_top2_ions, .(Peptide_clean), pairwise_avg_top2)

