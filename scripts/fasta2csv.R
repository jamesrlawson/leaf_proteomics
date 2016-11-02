install.packages('seqinr')

require(seqinr)
require(readr)

blah <- read.fasta('data/large_files/Egrandis_Eglobulus_chloroplast_Myrtales_At_mt_cRAP_160405.fasta', seqtype='AA', as.string=TRUE, strip.desc=TRUE)

EGrandis_AA_seqs <- data.frame(do.call(rbind, blah))

write.csv(EGrandis_AA_seqs, 'data/large_files/Egrandis_Eglobulus_chloroplast_Myrtales_At_mt_cRAP_160405.csv')

EGrandis_AA_seqs <- read_csv('data/large_files/Egrandis_Eglobulus_chloroplast_Myrtales_At_mt_cRAP_160405.csv')

ion_areas <- read_csv('data/large_files/D14_ion_areas_new_ion_library.csv')


D14_euc_AA_seqs <- EGrandis_AA_seqs[EGrandis_AA_seqs$Protein %in% ion_areas$Protein,]

write_csv(D14_euc_AA_seqs, 'data/large_files/Egrandis_Eglobulus_chloroplast_Myrtales_At_mt_cRAP_160405_D14eucSubset.csv')
