context("Step 0 - compiling bayesStanModel")

load("testthat_objects/bayesModel_test.RData")

test_that("bayesStanModel compiles to give an object of class 'stanmodel'", {
    expect_is(bayes_model, 'stanmodel')
})
