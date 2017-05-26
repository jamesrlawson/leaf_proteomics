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

protein_samples_D14$mean <- rowMeans(protein_samples_D14[,c(2:(ncol(protein_samples_D14)-3))])
  
get_sunburstData <- function(column) {

  to_sunburst <- dplyr::select(protein_samples_D14, Protein, NAME, mean) %>% dplyr::group_by(NAME) %>% dplyr::summarise(funccat_sum = sum(mean))
  
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
  
  x <- stringr::str_split(to_sunburst$NAME, ',')
  
  y <- rbind.fill(lapply(x,function(y){as.data.frame(t(y),stringsAsFactors=FALSE)}))
  
  to_sunburst$NAME <- y[,column]
  
  to_sunburst <- select(to_sunburst, NAME, funccat_sum)
  
  
  bla <- stringr::str_split(to_sunburst$NAME, '-')
  
  bla <- rbind.fill(lapply(bla,function(y){as.data.frame(t(y),stringsAsFactors=FALSE)}))
  
  
  for(i in 1:nrow(bla)) {

    if(any(is.na(bla[i,]))) {
      branch_end <- min(which(as.numeric(is.na(bla[i,]))==1))
      bla[i,branch_end] <- 'end'
    }
    
  }
  
  
  
  to_sunburst$NAME <- unite(as.data.frame(bla), all, 1:ncol(bla), sep = '-')
  
  to_sunburst$NAME <- sapply(to_sunburst$NAME, 
                             function(x){gsub(pattern = '<NA>',
                                              replacement = "", x, fixed=TRUE)})
  
  
  to_sunburst$NAME <- sapply(to_sunburst$NAME, 
                              function(x){gsub(pattern = "end.+",
                                               replacement = "end", x)})
  
  
  return(to_sunburst)

}

bla <- get_sunburstData(1)

bla <- merge(bla, filter(get_sunburstData(2), NAME != 'end'), all.x= TRUE, by = 'NAME')
bla$funccat_sum.x[!is.na(bla$funccat_sum.y)] <- bla$funccat_sum.x[!is.na(bla$funccat_sum.y)]+bla$funccat_sum.y[!is.na(bla$funccat_sum.y)]
bla <- select(bla, -funccat_sum.y) 
names(bla)[2] <- 'funccat_sum'

bla <- merge(bla, filter(get_sunburstData(3), NAME != 'end'), all.x= TRUE, by = 'NAME')
bla$funccat_sum.x[!is.na(bla$funccat_sum.y)] <- bla$funccat_sum.x[!is.na(bla$funccat_sum.y)]+bla$funccat_sum.y[!is.na(bla$funccat_sum.y)]
bla <- select(bla, -funccat_sum.y) 
names(bla)[2] <- 'funccat_sum'

bla <- merge(bla, filter(get_sunburstData(4), NAME != 'end'), all.x= TRUE, by = 'NAME')
bla$funccat_sum.x[!is.na(bla$funccat_sum.y)] <- bla$funccat_sum.x[!is.na(bla$funccat_sum.y)]+bla$funccat_sum.y[!is.na(bla$funccat_sum.y)]
bla <- select(bla, -funccat_sum.y) 
names(bla)[2] <- 'funccat_sum'

bla <- merge(bla, filter(get_sunburstData(5), NAME != 'end'), all.x= TRUE, by = 'NAME')
bla$funccat_sum.x[!is.na(bla$funccat_sum.y)] <- bla$funccat_sum.x[!is.na(bla$funccat_sum.y)]+bla$funccat_sum.y[!is.na(bla$funccat_sum.y)]
bla <- select(bla, -funccat_sum.y) 
names(bla)[2] <- 'funccat_sum'

bla <- merge(bla, filter(get_sunburstData(6), NAME != 'end'), all.x= TRUE, by = 'NAME')
bla$funccat_sum.x[!is.na(bla$funccat_sum.y)] <- bla$funccat_sum.x[!is.na(bla$funccat_sum.y)]+bla$funccat_sum.y[!is.na(bla$funccat_sum.y)]
bla <- select(bla, -funccat_sum.y) 
names(bla)[2] <- 'funccat_sum'

bla <- merge(bla, filter(get_sunburstData(7), NAME != 'end'), all.x= TRUE, by = 'NAME')
bla$funccat_sum.x[!is.na(bla$funccat_sum.y)] <- bla$funccat_sum.x[!is.na(bla$funccat_sum.y)]+bla$funccat_sum.y[!is.na(bla$funccat_sum.y)]
bla <- select(bla, -funccat_sum.y) 
names(bla)[2] <- 'funccat_sum'

bla <- group_by(bla, NAME) %>% dplyr::summarise(funccat_sum = sum(funccat_sum)) %>% ungroup(.)


require(sunburstR)
sunburst(bla)





#### convert to json

bla <- get_sunburstData(1)
blax <- str_split(bla$NAME, pattern = '-') 

y <- rbind.fill(lapply(blax,function(y){as.data.frame(t(y),stringsAsFactors=FALSE)}))[,1:4]
bla <- cbind(y,bla[,-1])


bla <- bla %>% group_by(V1,V2,V3,V4) %>% dplyr::summarise(funccat_sum = sum(funccat_sum)) %>% ungroup(.)

write_csv(bla, 'output/sunburst.csv')

require(RJSONIO)

makeList<-function(x){
  if(ncol(x)>2){
    listSplit<-split(x[-1],x[1],drop=T)
    lapply(names(listSplit),function(y){list(name=y,children=makeList(listSplit[[y]]))})
  }else{
    lapply(seq(nrow(x[1])),function(y){list(name=x[,1][y],funccat_sum=x[,2][y])})
  }
}

mylist <- makeList(bla)


jsonOut<-toJSON(list(name="MyData",children=makeList(bla)))
cat(jsonOut)

write(jsonOut, 'output/jsonOut.json')


blo <- bla[grep('lightreaction', bla$NAME),]

x <- str_split(blo$NAME, pattern = ',')
y <- rbind.fill(lapply(x,function(y){as.data.frame(t(y),stringsAsFactors=FALSE)}))[1]

View(select(blo, NAME, mean))