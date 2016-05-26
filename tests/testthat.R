library(testthat)
library(banocc)

knitr::knit('testthat/testthat_objects/generate_test_samples.Rnw')
test_check("banocc")
