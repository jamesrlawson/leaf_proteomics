

library(testthat)


test_that("nrow in climate_locs is consistent after each merge / duplicate removal", {
 
  for(n in 2:length(climate_locs_test[climate_locs_test != 0])) {
  
   expect_that(climate_locs_test[n], equals(climate_locs_test[n-1]), info = n)

  }
  
})


source('scripts/prep_data.R')

bla <- data %>% group_by(ID) %>% dplyr::summarise(n_per_ID = length(sample)) %>% full_join(select(data, sample, ID))

test_that("n=9 for all species-site combinations", {
  
  expect_that(unique(bla$n_per_ID), equals(9))
  
})

source('scripts/transformations.R')
