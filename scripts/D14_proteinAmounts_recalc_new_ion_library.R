require(readxl)
require(dplyr)
require(stringr)

load_excel_files <- function(path, sheet) { 
                      filenames = paste0(path, list.files(path))
                      tables <- lapply(filenames, read_excel, sheet = sheet)
                      do.call(cbind, tables)
                    }

# load ion areas, rename columns

euc_swath_reanalysed <- load_excel_files('data/large_files/', 'Area - ions')
euc_swath_reanalysed <- cbind(euc_swath_reanalysed[,1:9], euc_swath_reanalysed[,10:377] %>% select(contains('sample')))
names(euc_swath_reanalysed)[10:323] <- do.call(rbind, str_split(names(euc_swath_reanalysed[,10:323]), pattern = " "))[,1]

# load FDR data, rename columns, remove decoys and discard proteins below identification probability threshold

euc_swath_reanalysed_FDR <- load_excel_files('data/large_files/', 'FDR')
euc_swath_reanalysed_FDR <- cbind(euc_swath_reanalysed_FDR[,1:7], euc_swath_reanalysed_FDR[,10:363] %>% select(contains('sample')))
names(euc_swath_reanalysed_FDR)[8:319] <- do.call(rbind, str_split(names(euc_swath_reanalysed_FDR[,8:319]), pattern = " "))[,2]
euc_swath_reanalysed_FDR <- subset(euc_swath_reanalysed_FDR, Decoy == 0)

euc_swath_reanalysed_FDR$FDR <- apply(euc_swath_reanalysed_FDR[8:ncol(euc_swath_reanalysed_FDR)],1, function(x) as.numeric(any(x < 0.01)) )
nrow(euc_swath_reanalysed_FDR[euc_swath_reanalysed_FDR$FDR == 0,])

euc_swath_reanalysed_FDR <- subset(euc_swath_reanalysed_FDR, FDR == 1)

# ion areas constrained by FDR

ion_areas <- euc_swath_reanalysed[euc_swath_reanalysed$Peptide %in% euc_swath_reanalysed_FDR$Peptide,]
rm(euc_swath_reanalysed_FDR,euc_swath_reanalysed)
