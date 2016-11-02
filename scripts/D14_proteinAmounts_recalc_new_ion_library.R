require(readxl)
require(dplyr)
require(stringr)
require(readr)

load_excel_files <- function(path, sheet) { 
  
                        load_excel_arrange <- function(path, sheet) {
                          
                          out <- read_excel(path, sheet)
                          out <- arrange(out, Protein)
                          return(out)
                          
                        }
  
                      filenames = paste0(path, list.files(path))
                      tables <- lapply(filenames, load_excel_arrange, sheet = sheet)
                      do.call(cbind, tables)
}

# load ion areas, rename columns

euc_swath_reanalysed <- load_excel_files('data/large_files/swath/', 'Area - ions')
euc_swath_reanalysed <- cbind(euc_swath_reanalysed[,1:9], euc_swath_reanalysed[,10:377] %>% select(contains('sample')))
names(euc_swath_reanalysed)[10:323] <- do.call(rbind, str_split(names(euc_swath_reanalysed[,10:323]), pattern = " "))[,1]

# load FDR data, rename columns, remove decoys and discard proteins below identification probability threshold

euc_swath_reanalysed_FDR <- load_excel_files('data/large_files/swath/', 'FDR')
euc_swath_reanalysed_FDR <- cbind(euc_swath_reanalysed_FDR[,1:7], euc_swath_reanalysed_FDR[,10:363] %>% select(contains('sample')))
names(euc_swath_reanalysed_FDR)[8:319] <- do.call(rbind, str_split(names(euc_swath_reanalysed_FDR[,8:319]), pattern = " "))[,2]
euc_swath_reanalysed_FDR <- subset(euc_swath_reanalysed_FDR, Decoy == 0)

euc_swath_reanalysed_FDR$FDR <- apply(euc_swath_reanalysed_FDR[8:ncol(euc_swath_reanalysed_FDR)],1, function(x) as.numeric(any(x < 0.01)) )
nrow(euc_swath_reanalysed_FDR[euc_swath_reanalysed_FDR$FDR == 0,])

euc_swath_reanalysed_FDR <- subset(euc_swath_reanalysed_FDR, FDR == 1)

# ion areas constrained by FDR

ion_areas <- euc_swath_reanalysed[euc_swath_reanalysed$Peptide %in% euc_swath_reanalysed_FDR$Peptide,]
rm(euc_swath_reanalysed_FDR,euc_swath_reanalysed)

write_csv(ion_areas, 'data/large_files/D14_ion_areas_new_ion_library.csv')

# calculate peptide top2 for ovalb peptides

top2 <- function(x) {
  sum(sort(x, decreasing = TRUE)[1:2], na.rm=TRUE)
}

ion_areas_ovalb <- ion_areas[ion_areas$Peptide %in% c('GGLEPINFQTAADQAR',
                                                      'HIATNAVLFFGR',
                                                      'ISQAVHAAHAEINEAGR',
                                                      'YPILPEYLQC[Pye]VK',
                                                      'LTEWTSSNVMEER',
                                                      'DILNQITKPNDVYSFSLASR',
                                                      'DEDTQAMPFR'),]

peptide_areas_ov <-  ion_areas_ovalb %>% 
  group_by(Peptide) %>% 
  summarise_at(vars(10:323), top2)

# check distributions

blah <- data.frame(t(peptide_areas_ov[,2:315]))
names(blah) <- peptide_areas_ov$Peptide

hist(blah$DEDTQAMPFR)
hist(blah$DILNQITKPNDVYSFSLASR)
hist(blah$GGLEPINFQTAADQAR)
hist(blah$HIATNAVLFFGR)
hist(blah$ISQAVHAAHAEINEAGR)
hist(blah$`YPILPEYLQC[Pye]VK`)

  # summary stats
  
  peptide_areas_ov$mean <- apply(peptide_areas_ov[2:315],1, mean, na.rm=TRUE)
  peptide_areas_ov$median <- apply(peptide_areas_ov[2:315],1, median, na.rm=TRUE)
  peptide_areas_ov$sd <- apply(peptide_areas_ov[2:315],1, sd, na.rm=TRUE)
  peptide_areas_ov$min <- apply(peptide_areas_ov[2:315],1, min, na.rm=TRUE)
  peptide_areas_ov$max <- apply(peptide_areas_ov[2:315],1, max, na.rm=TRUE)
  
  
  peptide_areas_ov_stats <- select(peptide_areas_ov, Peptide, mean, median, sd, min, max)
  peptide_areas_ov_stats$med_minus_min_over_sd <- (peptide_areas_ov_stats$median - peptide_areas_ov_stats$min) / peptide_areas_ov_stats$sd
  peptide_areas_ov_stats$max_minus_med_over_sd <- (peptide_areas_ov_stats$max - peptide_areas_ov_stats$median) / peptide_areas_ov_stats$sd
  
  
  write_csv(peptide_areas_ov_stats, 'output/peptide_areas_ov_stats.csv')


# ovalbumin standard plots: ovalb peptide fractional area vs assayed ovalb concentration 
  
  # do top2top3 analysis for all proteins to get total protein area
  
  ion_areas <- read_csv('data/large_files/D14_ion_areas_new_ion_library.csv')
  
  
  top3 <- function(x) {
    mean(sort(x, decreasing = TRUE)[1:3], na.rm=TRUE)
  }
  
  protein_areas <-  ion_areas %>% 
    group_by(Peptide, Protein) %>% 
    summarise_at(vars(10:323), top2) %>%
    group_by(Protein) %>%
    summarise_at(vars(3:316), top3)
  
  # find area of ovalb peptides / total protein area
  
  total_protein_area <- colSums(protein_areas[,2:315])
  peptide_areas_ovalb_fract <- peptide_areas_ov
  peptide_areas_ovalb_fract[,2:315] <- t(t(peptide_areas_ovalb_fract[,2:315])/colSums(protein_areas[,2:315]))

  # ovalb peptide area / total area vs fraction ovalb / total protein by protein assay
  
  protein_assay <- read_csv('data/protein_assay.csv')
  
    # produce plots and output df with model coeffs for each peptide
  
    my.list <- vector("list", nrow(peptide_areas_ovalb_fract))
    
    for(i in 1:nrow(peptide_areas_ovalb_fract)) {
      
      # browser()
      
      fraction_oval <- peptide_areas_ovalb_fract[i,]
      x <- fraction_oval[1,cols]
      
      blah <- protein_assay %>% select(sample, percent_ovalb_of_total_protein) %>% mutate(percent_ovalb_of_total_protein = percent_ovalb_of_total_protein / 100) %>%
        arrange(sample)
      
      y <- data.frame(cbind(fraction_oval_MS = t(x), sample = names(x)))
      
      j <- merge(blah, y, by = 'sample')
      
      names(j)[3] <- 'fraction_oval_MS'
      j$fraction_oval_MS <- as.numeric(as.character(j$fraction_oval_MS))
      
      plot(j$fraction_oval_MS ~ j$percent_ovalb_of_total_protein, main = peptide_areas_ovalb_fract$Peptide[i], ylab = 'Fraction of ovalbumin (MS)', xlab = 'Ovalbumin % of assay protein')
      model <- lm(j$fraction_oval_MS ~ j$percent_ovalb_of_total_protein)
      abline(model)
      
      pep <- peptide_areas_ovalb_fract$Peptide[i]
      
      model.stats <-  data.frame(cbind(pep,
                                       round(summary(model)$r.squared,3),
                                       round(summary(model)$coefficients[,4][2],2),
                                       model$coefficients[1],
                                       model$coefficients[2]))
      
      names(model.stats) <- c('submodel','R2','pval', 'intercept', 'slope') 
      rownames(model.stats) <- NULL
      
      my.list[[i]] <- model.stats
      
    }
    
    x <- data.frame()
    
    output <- rbind(x, do.call(rbind, my.list))
    output <- output %>% mutate(R2 = as.numeric(as.character(R2))) %>% arrange(desc(R2))
    
    write_csv(output, 'output/ovalb_peptide_standards.csv')
    