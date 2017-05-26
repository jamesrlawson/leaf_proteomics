
require(seqinr)
require(readr)
require(dplyr)
require(stringr)

blah <- read.fasta('data/proteomics_data/sequences/euc/Egrandis_Eglobulus_chloroplast_Myrtales_At_mt_cRAP_160405.fasta', seqtype='AA', as.string=TRUE, strip.desc=TRUE)

EGrandis_AA_seqs <- data.frame(do.call(rbind, blah))

write.csv(EGrandis_AA_seqs, 'data/proteomics_data/sequences/euc/Egrandis_Eglobulus_chloroplast_Myrtales_At_mt_cRAP_160405.csv')

# had to manually remove proteins with | at the end (e.g. human, bov etc.)

EGrandis_AA_seqs <- read_csv('data/proteomics_data/sequences/euc/Egrandis_Eglobulus_chloroplast_Myrtales_At_mt_cRAP_160405.csv')
names(EGrandis_AA_seqs) <- c('Protein', 'AA_sequence')


if(!exists("ion_areas")) {
  ion_areas <- read_csv('data/proteomics_data/proteomics/derived/euc/D14_ion_areas_new_ion_library.csv')
}

D14_euc_AA_seqs <- EGrandis_AA_seqs[EGrandis_AA_seqs$Protein %in% unique(ion_areas$Protein),]

write_csv(D14_euc_AA_seqs, 'data/proteomics_data/sequences/euc/Egrandis_Eglobulus_chloroplast_Myrtales_At_mt_cRAP_160405_D14eucSubset.csv')

