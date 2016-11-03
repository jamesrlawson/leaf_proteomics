# protein molecular weights #

require(Peptides)
require(readr)

source('scripts/fasta2csv.R')

D14_euc_AA_seqs <- read_csv('data/large_files/Egrandis_Eglobulus_chloroplast_Myrtales_At_mt_cRAP_160405_D14eucSubset.csv')

D14_euc_AA_seqs$MW <- apply(D14_euc_AA_seqs[,2],1, mw)

write_csv(D14_euc_AA_seqs, 'data/D14_protein_MW_seqs.csv')

write_csv(select(D14_euc_AA_seqs, - AA_sequence), 'data/D14_protein_MW.csv')


