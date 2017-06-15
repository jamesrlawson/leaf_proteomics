# temperature dependences of Gamma and Km
.Rgas <- function()8.314
Tk <- function(x)x+273.15

# Arrhenius
arrh <- function(Tleaf, Ea){
  exp((Ea * (Tk(Tleaf) - 298.15)) / (298.15 * .Rgas() * Tk(Tleaf))) 
}

# modified Arrhenius function
marrh <- function(Tleaf, Ea, delS, Ed = 200000, Tref = 298.15) {
  t1 <- arrh(Tleaf, Ea)
  t2 <- 1 + exp((delS * Tref - Ed) / .Rgas() / Tref) 
  t3 <- 1 + exp((delS * Tk(Tleaf) - Ed) / .Rgas() / Tk(Tleaf))
  return(t1*t2/t3)
}

#Gamma Star
TGammaStar <- function(Tleaf,  Egamma=37830.0, value25=42.75){  
  value25*arrh(Tleaf,Egamma)
}

# Km 
TKm <- function(Tleaf) {
  Kc <- 404.9*arrh(Tleaf,79430)
  Ko <- 278400*arrh(Tleaf,36380)
  Oi <- 205000
  Km <- Kc * (1+Oi/Ko)
  return(Km)
}

# function to implement one-point method
# T parameters for Km, G* from Bernacchi et al. (2001)
# T parameters for Vcmax from Kumarathunge (fitted to many Euc species)
Vcmax <- function(Photo, Ci, Tleaf) {
  
  Km <- TKm(Tleaf)
  Gamstar <- TGammaStar(Tleaf)
  Vm <- Photo/((Ci-Gamstar)/(Ci+Km) - 0.015)
  Vm25 <- Vm/marrh(Tleaf, Ea=64930, delS=636)
  
  return(Vm25)
}

# function to implement one-point method
# T parameters for G* from Bernacchi et al. (2001)
# T parameters for Jmax from Kumarathunge (fitted to many Euc species)
Jmax <- function(Photo, Ci, Tleaf) {
  
  Gamstar <- TGammaStar(Tleaf)
  Jm <- 4*Photo/((Ci-Gamstar)/(Ci+2*Gamstar) - 0.015/2)
  Jm25 <- Jm/marrh(Tleaf, Ea=28900, delS=631)
  
  return(Jm25)
}




include_photosynthesis=TRUE
include_leaf_N = TRUE

source('scripts/transformations.R')

source('scripts/prep_data_mg_per_mm2.R')

bla <- select(data, sample, photo_amb, Ci, Tleaf)# %>%
#   filter(Ci > 150) %>%
#   filter(Ci < 350)
bla <- na.omit(bla)

bla$Vcmax <- Vcmax(bla$photo_amb, bla$Ci, bla$Tleaf)

blah <- merge(data, bla, by = c('sample','Ci', 'photo_amb'))

#blah <- filter(blah, Jmax > 30)
blah <- filter(blah, photo_amb > 1)

plot(blah$photo_amb ~ blah$N_per_area)

blah$etrans <- blah$electron_transport_minATPsynth + blah$ATP_synthase_chloroplastic
cor(blah$Vcmax,blah$etrans)

plot(blah$Vcmax ~ blah$electron_transport_minATPsynth)

blah <- filter(blah, Vcmax < 150)

bla.lm <- lm(Vcmax ~ protein + stress + Rubisco + TCA_org_transformation + cytochrome_b6f + other_electron_carrier + Photosystems + calvin_cycle + photorespiration + ATP_synthase_chloroplastic + electron_transport + glycolysis + hormone_metabolism, blah)

require(MuMIn)
options(na.action = "na.fail")
bla.dredge <- dredge(bla.lm, m.max = 1, extra = c('R^2', 'adjR^2'))

plot(blah$Jmax ~ blah$protein)
plot(blah$Jmax ~ blah$cytochrome_b6f)

plot(protein ~ cytochrome_b6f, blah)
