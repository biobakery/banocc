context("Step 0 - compiling bayesStanModel")

test_that("bayesStanModel compiles to give an object of class 'stanmodel'", {
    bayes_model <- rstan::stan_model(model_code=bayesStanModel)
    expect_is(bayes_model, 'stanmodel')
    save(bayes_model, file="bayesModel_test.RData")
})
