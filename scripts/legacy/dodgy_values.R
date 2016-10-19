# finds suspicious values in the main database (according to CV and tech rep additions)

require(plyr)

DP14_data <- read.csv('../data/kruft/DP14_Database_160914.csv', header=T, stringsAsFactors = F, na.strings = c("","#DIV/0!"))
DP14_data <- DP14_data[!DP14_data$leaf_age %in% 'sen',]
  
DP14_data$leaf_FW_check <- as.numeric(DP14_data$leaf_FW_check)
DP14_data$LWC <- as.numeric(DP14_data$LWC)
DP14_data$LMA <- as.numeric(DP14_data$LMA)
names(DP14_data)[names(DP14_data) == 'whole_leaf_area_per_whole_leaf_FW.cm2.g.'] <- 'freshSLA'
DP14_data$freshSLA <- as.numeric(DP14_data$freshSLA)
DP14_data$freshSLA[DP14_data$freshSLA == 0] <- NA

DP14_data_dodgy_FW <- rbind(subset(DP14_data, leaf_FW_check > 0.02),
                    subset(DP14_data, leaf_FW_check < -0.1))
DP14_data_dodgy_LWC <- ddply(DP14_data[!is.na(DP14_data$LWC),], .(biological_rep, species_confirmed), summarise, LWC_CV = CV(LWC))
DP14_data_dodgy_LMA <- ddply(DP14_data[!is.na(DP14_data$LMA),], .(biological_rep, species_confirmed), summarise, LMA_CV = CV(LMA))
DP14_data_dodgy_freshSLA <- ddply(DP14_data[!is.na(DP14_data$freshSLA),], .(biological_rep, species_confirmed), summarise, freshSLA_CV = CV(freshSLA))

DP14_data_dodgy_CV <- merge(DP14_data_dodgy_LWC, DP14_data)
DP14_data_dodgy_CV <- merge(DP14_data_dodgy_LMA, DP14_data_dodgy_CV, all=TRUE)
DP14_data_dodgy_CV <- merge(DP14_data_dodgy_freshSLA, DP14_data_dodgy_CV, all=TRUE)


write.csv(DP14_data_dodgy_FW, '../output/DP14_data_dodgy_FW.csv')
write.csv(DP14_data_dodgy_CV, '../output/DP14_data_dodgy_CV.csv')


DP14_data_dodgy_LWC <- ddply(DP14_data[!is.na(DP14_data$LWC),], .(leaf_age, species_confirmed), summarise, LWC_CV = CV(LWC))
DP14_data_dodgy_LMA <- ddply(DP14_data[!is.na(DP14_data$LMA),], .(leaf_age, species_confirmed), summarise, LMA_CV = CV(LMA))
DP14_data_dodgy_freshSLA <- ddply(DP14_data[!is.na(DP14_data$freshSLA),], .(biological_rep, species_confirmed), summarise, freshSLA_CV = CV(freshSLA))

DP14_data_dodgy_CV <- merge(DP14_data_dodgy_LWC, DP14_data)
DP14_data_dodgy_CV <- merge(DP14_data_dodgy_LMA, DP14_data_dodgy_CV, all=TRUE)
DP14_data_dodgy_CV <- merge(DP14_data_dodgy_freshSLA, DP14_data_dodgy_CV, all=TRUE)
