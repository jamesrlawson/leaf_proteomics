mg_per_m2 = TRUE
moles = FALSE

## ASSIGN PROTEINS TO FUNCTIONAL CATEGORIES AND CALCULATE PROTEIN AMOUNTS IN EACH CATEGORY ##

# called by transformations.R

require(readr)
require(stringr)
require(dplyr)
require(tidyr)

source('scripts/functions.R')

mercator <- read_csv('data/proteomics_data/mercator/euc/D14_mercator_20170217.csv')

# 'mg_per_m2' and 'moles' switches are defined in transformations.R
if(mg_per_m2) {
  protein_samples_D14 <- read_csv('data/proteomics_data/proteomics/derived/euc/D14_protein_GGLEP-DEDT.csv') # protein amounts calculated using D14 ion library, in avg(GGLEP/DEDT) equivalents
}

if(moles) {
  protein_samples_D14 <- read_csv('data/proteomics_data/proteomics/derived/euc/D14_protein_moles_GGLEP-DEDT.csv') # protein amounts as above but in moles (not multiplied by MW)
}

# first add the mercator$NAME values for each protein in protein_samples_D14

protein_samples_D14 <- getProteinBins(protein_samples_D14, mercator)

protein_samples_D14$mean <- rowMeans(protein_samples_D14[,c(2:(ncol(protein_samples_D14)-3))]) # calculate mean values across all samples

# this will create data for a mean across all samples

get_sunburstData_mean <- function(column) {
  
  to_sunburst <- dplyr::select(protein_samples_D14, Protein, NAME, mean)
  
  # some gsubs to get names into the right format for string splitting
  
  to_sunburst$NAME <- sapply(to_sunburst$NAME, 
                             function(x){gsub(pattern = "\\'",
                                              replacement = "", x)})
  
  to_sunburst$NAME <- sapply(to_sunburst$NAME, 
                             function(x){gsub(pattern = "-",
                                              replacement = "_", x)})
  
  to_sunburst$NAME <- sapply(to_sunburst$NAME, 
                             function(x){gsub(pattern = " ",
                                              replacement = "_", x)})
  
  to_sunburst$NAME <- sapply(to_sunburst$NAME, 
                             function(x){gsub(pattern = ",_",
                                              replacement = ",", x)})
  
  to_sunburst$NAME <- sapply(to_sunburst$NAME, 
                             function(x){gsub(pattern = "([0-9])(,)([0-9])", # create 3 capturing groups (first number)(comma)(second number)
                                              replacement = "\\1\\3", x)}) # replacement is first capture and third capture
  
  to_sunburst$NAME <- sapply(to_sunburst$NAME, 
                             function(x){gsub(pattern = "\\.",
                                              replacement = "-", x)})
  
  to_sunburst$NAME <- sapply(to_sunburst$NAME, 
                             function(x){gsub(pattern = " ",
                                              replacement = "", x)})
  
  x <- stringr::str_split(to_sunburst$NAME, ',') # split by ',' which should separate multiple functional category assignments
  
  y <- rbind.fill(lapply(x,function(y){as.data.frame(t(y),stringsAsFactors=FALSE)})) # build a df from x which fills blanks with NA's
  
  to_sunburst$NAME <- y[,column] # to_sunburst's name is now equal to just the first functional assignment
  
  to_sunburst <- select(to_sunburst, NAME, mean)
  
  bla <- stringr::str_split(to_sunburst$NAME, '-') # now split out the name into multiple columns
  
  bla <- rbind.fill(lapply(bla,function(y){as.data.frame(t(y),stringsAsFactors=FALSE)}))
  
  to_sunburst <- cbind(bla, mean = to_sunburst$mean) # and add the calculated mean protein amount to the df 

  to_sunburst <- to_sunburst %>% group_by(V1,V2,V3,V4,V5,V6,V7) %>% dplyr::summarise(mean = sum(mean)) %>% ungroup(.) %>% filter(!V1 %in% "") # combine values for proteins in same V7 level
    
  return(to_sunburst)
  
}

bla <- get_sunburstData_mean(1)



# use a loop to ask, for a given row, does the last column have a parent column? if not then create. 
# then does second to last column have a parent? if not then create...

blanames <- names(bla[,1:7])

for(i in rev(2:length(blanames))) {
  
parent_col <- blanames[i-1]

child_col <- blanames[i]

for(j in 1:nrow(bla)) {

  vec_j <- bla[j,1:7]

  parent <- as.character(vec_j[(i-1)])
  

  is.na(bla[j,child_col])
  
  
 # if(!any(bla[[parent_col]] %in% parent)) {
  
  if(!is.na(bla[j,child_col])) { # if there isn't an NA value in the child column, create an extra row. this will make dupes but we can handle that...
  
    vec_parent <- vec_j
    vec_parent[(i-1):7] <- NA
    vec_parent <- cbind(vec_parent, mean = NA)
    
    bla <- rbind(bla, vec_parent)
      
    bla <- distinct(bla, V1,V2,V3,V4,V5,V6,V7,mean, .keep_all=TRUE)
    
  }
  
#  browser()
  
}

}

write_csv(bla, 'output/sunburst.csv')


# this will create sunburst data for every sample

get_sunburstData_all <- function(column) {
  
  to_sunburst <- dplyr::select(protein_samples_D14, -BINCODE)
  
  # some gsubs to get names into the right format for string splitting
  
  to_sunburst$NAME <- sapply(to_sunburst$NAME, 
                             function(x){gsub(pattern = "\\'",
                                              replacement = "", x)})
  
  to_sunburst$NAME <- sapply(to_sunburst$NAME, 
                             function(x){gsub(pattern = "-",
                                              replacement = "_", x)})
  
  to_sunburst$NAME <- sapply(to_sunburst$NAME, 
                             function(x){gsub(pattern = " ",
                                              replacement = "_", x)})
  
  to_sunburst$NAME <- sapply(to_sunburst$NAME, 
                             function(x){gsub(pattern = ",_",
                                              replacement = ",", x)})
  
  to_sunburst$NAME <- sapply(to_sunburst$NAME, 
                             function(x){gsub(pattern = "([0-9])(,)([0-9])", # create 3 capturing groups (first number)(comma)(second number)
                                              replacement = "\\1\\3", x)}) # replacement is first capture and third capture
  
  to_sunburst$NAME <- sapply(to_sunburst$NAME, 
                             function(x){gsub(pattern = "\\.",
                                              replacement = "-", x)})
  
  to_sunburst$NAME <- sapply(to_sunburst$NAME, 
                             function(x){gsub(pattern = " ",
                                              replacement = "", x)})
  
  x <- stringr::str_split(to_sunburst$NAME, ',') # split by ',' which should separate multiple functional category assignments
  
  y <- rbind.fill(lapply(x,function(y){as.data.frame(t(y),stringsAsFactors=FALSE)})) # build a df from x which fills blanks with NA's
  
  to_sunburst$NAME <- y[,column] # to_sunburst's name is now equal to just the first functional assignment
  
  bla <- stringr::str_split(to_sunburst$NAME, '-') # now split out the name into multiple columns
  
  bla <- rbind.fill(lapply(bla,function(y){as.data.frame(t(y),stringsAsFactors=FALSE)}))
  
  to_sunburst <- cbind(bla, to_sunburst) # and create a df with the split name and all the protein amount data
  
  # combine values for proteins in same V7 level, for every sample
  to_sunburst <- group_by(to_sunburst, V1,V2,V3,V4,V5,V6,V7) %>% 
    dplyr::summarise_at(vars(9:(ncol(to_sunburst)-1)), sum) %>% ungroup(.) %>% filter(!V1 %in% "") 
  
  replicates <- read_csv('data/misc_data/replicates.csv') 
  replicates <- replicates[replicates$sample %in% names(protein_samples_D14)[2:327],]
  
  to_sunburst <- gather(to_sunburst, key = 'sample', value = 'value', 8:333) %>% 
    full_join(replicates, by = 'sample') %>% 
    dplyr::group_by(V1,V2,V3,V4,V5,V6,V7,ID) %>%
    dplyr::summarise(ID_mean = mean(value, na.rm=TRUE)) %>%
    spread(key=ID, value=ID_mean)
  
  return(to_sunburst)
  
}

mercator <- read_csv('data/proteomics_data/mercator/euc/D14_mercator_20170217.csv')

# 'mg_per_m2' and 'moles' switches are defined in transformations.R
if(mg_per_m2) {
  protein_samples_D14 <- read_csv('data/proteomics_data/proteomics/derived/euc/D14_protein_GGLEP-DEDT.csv') # protein amounts calculated using D14 ion library, in avg(GGLEP/DEDT) equivalents
}

if(moles) {
  protein_samples_D14 <- read_csv('data/proteomics_data/proteomics/derived/euc/D14_protein_moles_GGLEP-DEDT.csv') # protein amounts as above but in moles (not multiplied by MW)
}

# first add the mercator$NAME values for each protein in protein_samples_D14

protein_samples_D14 <- getProteinBins(protein_samples_D14, mercator)

bla <- get_sunburstData_all(1)





















# use nested loop

# for i in unique(bla$V6), subset bla by unique(bla$V6) in bla$V6, 
# then sub[sub$V7 == "",]$mean = sub[sub$V7 == "",]$mean + sum(sub[!sub$V6 == "",]$mean)

blanames <- names(bla[,1:7])

for(i in rev(2:length(blanames))) {
  
 child_col <- blanames[i]
 parent_col <- blanames[i-1]
 
 #parent_col_uniques <- as.vector(na.omit(unique(bla[[parent_col]])))
 
 parent_col_uniques <- as.vector(na.omit(bla[[parent_col]]))
 
 
  for(j in 1:length(parent_col_uniques)) {
    
    target_row <- parent_col_uniques[j]
    
    parent_col_sub <- cbind(bla[,c(parent_col, child_col)],mean = bla$mean)
  
    blah <- subset(parent_col_sub, get(parent_col) == target_row)
    
    #if(nrow(blah > 1)) {
    
  #  if(any(is.na(bla[child_col]))) { 
    if(any(is.na(blah[child_col]))) { 
    
      x <-  which(bla[parent_col] == target_row)
      
      y <-  which(is.na(bla[child_col]))
      
      x <- x[x%in%y]
      
      check <- is.na(blah[which(is.na(blah[child_col])),]$mean) 
      
      if(length(check) > 1) {
       
       # browser() 
        
      }
      
      if(!check) {
        
      bla[x,'mean'] <- blah[which(is.na(blah[child_col])),]$mean + sum(blah[which(!is.na(blah[child_col])),]$mean, na.rm=TRUE)
  
      } else {
       
        bla[x,'mean'] <- sum(blah[which(!is.na(blah[child_col])),]$mean, na.rm=TRUE)
        
      }
      
      }
    
  }

}

#bla <- bla[!is.na(bla$mean),]

write_csv(bla, 'output/sunburst.csv')




















