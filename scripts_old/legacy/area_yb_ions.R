# this code calculates the number of y ions which have a greater mean area than the b ion with the greatest mean area
# this number is output as area_ions$numArea_YgreaterthanB
# a binary values defining 

peptides <- read.csv('data/peptides.csv', header=T, stringsAsFactors = F)
area_ions <- read.csv('data/12spp_recal_Arabidopsis_library_recal_160722_SVS_160831_areaIons.csv', header=T, stringsAsFactors = F)

area_ions$numArea_YgreaterthanB <- NA
area_ions$avgIntensity <- NA

for(i in 1:nrow(peptides)) {
  
   thisPep <- peptides$Peptide[i]
    
   if(any(area_ions$Peptide %in% thisPep)) {
    
     peptide_areas <- area_ions[area_ions$Peptide %in% thisPep,] # subset peptide_areas df by current peptide
      
     peptide_areas$area_sums <- rowMeans(peptide_areas[,10:45])
     
     blah <- data.frame(cbind(peptide_areas$Ion.Type, peptide_areas$area_sums))
    
     blah$X2 <- as.numeric(as.character(blah$X2))
     blah$X1 <- as.character(blah$X1)
     
     blah <- na.omit(blah)
     
     area_ions[area_ions$Peptide %in% thisPep,]$avgIntensity <- rowMeans(peptide_areas[,10:45])
     
    # browser()
     
     # do the top two ions contain missing values?
     area_ions$missingvals <- 0
     peptide_areas <- peptide_areas[order(peptide_areas$avgIntensity, decreasing=TRUE),] 
     area_ions[area_ions$Peptide %in% thisPep,]$missingvals <- any(is.na(subset(peptide_areas[c(1:2),], Ion.Type == "y"))) # logical indicating presence of missing vals for top 2 y ions
       
     if(any(blah$X1 %in% 'b')) { # subset has to contain a 'b' ion for this line to work
      area_ions[area_ions$Peptide %in% thisPep,]$numArea_YgreaterthanB <- length(blah$X2[blah$X2 > max(blah[blah$X1 %in% 'b',]$X2)]) 
     } else {
     
     area_ions[area_ions$Peptide %in% thisPep,]$numArea_YgreaterthanB <- length(blah$X2) 
     
   }
   
  }
 
}

write.csv(area_ions, 'output/area_ions.csv')

blah <- area_ions[area_ions$Peptide %in % peptides$Peptide,]

write.csv(blah, 'output/area_ions_trunc.csv')
 

# avg of top 3 peptides per protein, summed top 2 ions per peptide

peptides <- read.csv('data/peptides.csv', header=T, stringsAsFactors = F)
area_ions <- read.csv('data/12spp_recal_Arabidopsis_library_recal_160722_SVS_160831_areaIons.csv', header=T, stringsAsFactors = F)

area_ions <- area_ions[area_ions$Peptide %in% peptides$Peptide,]

require(plyr)

sumTop2ions <- ddply(area_ions, .(Protein, Peptide), summarise, sumTop2ions = sum(Fragment.MZ[c(1:2)]))
avgTop3Peptides <- ddply(sumTop2ions, .(Protein), summarise, avgTop3Peptides = mean(sumTop2ions[c(1:3)], na.rm=TRUE))

avgTop3Peptides$ovbStand <- avgTop3Peptides$avgTop3Peptides / avgTop3Peptides[avgTop3Peptides$Protein %in% "sp|OVAL_CHICK",]$avgTop3Peptides
avgTop3Peptides <- avgTop3Peptides[order(avgTop3Peptides$ovbStand, decreasing=TRUE),]

