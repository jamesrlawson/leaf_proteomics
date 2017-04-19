# finds sample numbers in main database using the fields in the licor data remarks as search terms
# also adds in area adjustments

require(plyr)
require(readr)
require(stringr)

licor <- read_csv('data/licor/licor_data2.csv')
samples <- read_csv('data/licor/sample_data.csv')

licor <- licor[!is.na(licor$licor_remarks),]

getSample <- function(licor_df, samples) {
  
  #licor_df <- licor_df[!is.na(licor$licor_remarks),]
  
#  browser()
  
  licor_df$sample <- NA
  
  datesub <- unique(licor_df$licor_date)
  samples_df <- subset(samples, samples$date == datesub)

  for(i in 1:nrow(licor_df)) {
    
    if(any(grep(licor_df$licor_remarks[i], paste(samples_df$species_confirmed, samples_df$biological_rep, " ", samples_df$leaf_age, sep = "")))) {
      
      samp <- samples_df[grep(licor_df$licor_remarks[i], paste(samples_df$species_confirmed, samples_df$biological_rep, " ", samples_df$leaf_age, sep = "")),]$sample_number
      licor_df$sample[i] <- samp
      
    } else {
      
      if(any(grep(licor_df$licor_remarks[i], paste(samples_df$species_confirmed, samples_df$leaf_age, sep = " ")))) {
        
        samp <- samples_df[grep(licor_df$licor_remarks[i], paste(samples_df$species_confirmed, samples_df$leaf_age, sep = " ")),]$sample_number
        licor_df$sample[i] <- samp
        
      }
      
    }
    
  }
  
  return(licor_df)
  
}



licor_samples <- ddply(licor, .(licor_date), getSample, samples)

licor_samples$HHMMSS <- NULL

write_csv(licor_samples, 'data/licor/licor_samples.csv')






licor <- read_csv('data/licor/licor_samples.csv')
samples <- read_csv('data/licor/area_adjustments.csv')

licor <- licor[!is.na(licor$sample),]


for(i in 1:nrow(licor)) {

  #browser()
  
  if(any(grep(licor$sample[i], samples$sample_number, fixed=TRUE))) {
      
    new_area <-  samples[grep(licor$sample[i], samples$sample_number, fixed=TRUE),]$leaf_area
    
    licor$Area[i] <- new_area
  
  }

}

write_csv(licor, 'data/licor/licor_area_adj.csv')













