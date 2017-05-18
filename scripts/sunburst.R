## ASSIGN PROTEINS TO FUNCTIONAL CATEGORIES AND CALCULATE PROTEIN AMOUNTS IN EACH CATEGORY ##

# called by transformations.R

require(readr)
require(stringr)
require(dplyr)
require(tidyr)

source('scripts/functions.R')

mercator <- read_csv('data/mercator/D14_mercator_20170217.csv')

# 'mg_per_m2' and 'moles' switches are defined in transformations.R
if(mg_per_m2) {
  protein_samples_D14 <- read_csv('data/D14_protein_GGLEP-DEDT.csv') # protein amounts calculated using D14 ion library, in avg(GGLEP/DEDT) equivalents
}

if(moles) {
  protein_samples_D14 <- read_csv('data/D14_protein_moles_GGLEP-DEDT.csv') # protein amounts as above but in moles (not multiplied by MW)
}

# first add the mercator$NAME values for each protein in protein_samples_D14

protein_samples_D14 <- getProteinBins(protein_samples_D14, mercator)

# then import the names of categories we're interested in from mercator_names* and use grep to find all proteins associated with those categories
# put the results in instances of a list
# mercator_names.csv contains search terms for protein categories. These searches are run on mercator$NAMES. 
# search terms must be unique to the functional category to avoid non-target returns. 
# for example, a search for 'photosystem I' will also pick up proteins from 'photosystem II' - to avoid this we search for 'photosystem I\.'
# this works because all instances of proteins within 'photosystem I' are actually within subcategories. We'd miss some returns if there were proteins in the upper 'photosystem I' category.
# N.B. the '\' is an 'escape' and must be used because .'s are special in regular expressions and mean 'anything'. By using the escape we will actually search for the character '.'
# search terms for top level categories can be made unique by using ' in front

mercator_names <- read.csv('data/mercator/mercator_names.csv', header=T, stringsAsFactors = F) 
mercator_names <- arrange(mercator_names, funccat)

func_assigned.list <- vector('list', length(mercator_names$funccat))

func_assigned <- data.frame()

for(i in 1:length(mercator_names$funccat)) {
  
  name <- mercator_names$funccat[i]
  
  proteins <- protein_samples_D14[grep(name, protein_samples_D14$NAME),]
  
  proteins$funccat <- mercator_names$funccat[i]
  
  proteins <- distinct(proteins, Protein, .keep_all = TRUE)
  
  func_assigned.list[[i]] <- proteins
  
}

func_assigned <- rbind(func_assigned, do.call(rbind, func_assigned.list))
rm(func_assigned.list)

func_assigned$mean <- rowMeans(func_assigned[,c(2:(ncol(func_assigned)-4))])



to_sunburst <- dplyr::select(func_assigned, Protein, NAME, mean) %>% dplyr::group_by(NAME) %>% dplyr::summarise(funccat_sum = sum(mean))



to_sunburst$NAME <- sapply(to_sunburst$NAME, 
       function(x){gsub(pattern = "\\'",
                        replacement = "", x)})

to_sunburst$NAME <- sapply(to_sunburst$NAME, 
                           function(x){gsub(pattern = "\\.",
                                            replacement = "-", x)})

to_sunburst$NAME <- sapply(to_sunburst$NAME, 
                           function(x){gsub(pattern = " ",
                                            replacement = "", x)})


require(stringr)

x <- str_split(to_sunburst$NAME, ',')
y <- do.call(rbind, x)

to_sunburst$NAME1 <- y[,1]

to_sunburst <- select(to_sunburst, NAME1, funccat_sum)

to_sunburst$funccat_sum <- as.integer(to_sunburst$funccat_sum)


bla <- str_split(to_sunburst$NAME1, '-')
bla <- do.call(rbind, bla)

for(i in 1:nrow(bla)) {
  branch_length <- max(which(!duplicated(bla[i,])==TRUE))
    if(branch_length < 9) {
      #bla[i,c((branch_length+1):ncol(bla))] <- bla[i,branch_length]
     bla[i,c((branch_length+1):ncol(bla))] <- 'end'
    }
}

for(i in 1:nrow(bla)) {
  branch_length <- max(which(!duplicated(bla[i,])==TRUE))
  if(branch_length < 9) {
    #bla[i,c((branch_length+1):ncol(bla))] <- bla[i,branch_length]
    bla[i,c((branch_length+1):ncol(bla))] <- ""
  }
}


to_sunburst$NAME1 <- unite(as.data.frame(bla), all, 1:9, sep = '-')

to_sunburst$NAME1 <- sapply(to_sunburst$NAME1, 
                           function(x){gsub(pattern = "----",
                                            replacement = "", x)})


to_sunburst$NAME1 <- sapply(to_sunburst$NAME1, 
                          function(x){gsub(pattern = "end-",
                                           replacement = "end", x)})
sunburst(to_sunburst)

