# temperature dependences of Gamma and Km
.Rgas <- function()8.314
Tk <- function(x)x+273.15

# Arrhenius
arrh <- function(Tleaf, Ea){
  exp((Ea * (Tk(Tleaf) - 298.15)) / (298.15 * .Rgas() * Tk(Tleaf))) 
}

#Gamma Star
TGammaStar <- function(Tleaf,  Egamma=37830.0, value25=42.75){  
  value25*arrh(Tleaf,Egamma)
}

# Km 
Tkm <- function(Tleaf) {
  Kc <- 404.9*arrh(Tleaf,79430)
  Ko <- 278400*arrh(Tleaf,36380)
  Oi <- 205000
  Km <- Kc * (1+Oi/Ko)
  return(Km)
}

# function to implement one-point method
Vcmax <- function(Photo, Ci, Tleaf) {
  
  Km <- Tkm(Tleaf)
  Gamstar <- TGammaStar(Tleaf)
  Vm <- Photo/((Ci-Gamstar)/(Ci+Km) - 0.015)
  Vm25 <- Vm / (arrh(Tleaf,51560))
  return(Vm25)
}


Jmax <- function(Photo, Ci, Tleaf) {
  
  Gamstar <- TGammaStar(Tleaf)
  Jm <- 4*Photo/((Ci-Gamstar)/(Ci+2*Gamstar) - 0.0075)
  return(Jm)
}

plot(blah$Tleaf.x)


include_photosynthesis=TRUE
include_leaf_N = TRUE

source('scripts/transformations.R')

source('scripts/prep_data_mg_per_mm2.R')

bla <- select(data, sample, photo_max, Ci, Tleaf)# %>%
    #   filter(Ci > 150) %>%
    #   filter(Ci < 350)
bla <- na.omit(bla)

bla$Jmax <- Jmax(bla$photo_max, bla$Ci, bla$Tleaf)

blah <- merge(data, bla, by = c('sample','Ci'))

blah <- filter(blah, Jmax > 30)

plot(data$photo_max ~ data$N_per_area)

blah$etrans <- blah$electron_transport_minATPsynth + blah$ATP_synthase_chloroplastic
cor(blah$Jmax,blah$etrans)

plot(blah$Jmax ~ blah$electron_transport_minATPsynth)

bla.lm <- lm(Jmax ~ protein + stress + Rubisco + TCA_org_transformation + cytochrome_b6f + other_electron_carrier + Photosystems + calvin_cycle + photorespiration + ATP_synthase_chloroplastic + electron_transport + glycolysis + hormone_metabolism, blah)

bla.dredge <- dredge(bla.lm, m.max = 3, extra = c('R^2', 'adjR^2'))

plot(blah$Jmax ~ blah$protein)
plot(blah$Jmax ~ blah$cytochrome_b6f)

plot(protein ~ cytochrome_b6f, blah)
