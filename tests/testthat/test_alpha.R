context("Test alpha screening") 

test_that("Fake test", { 
  T = 100
  N = 25
  set.seed(123)
  rets  = matrix(rnorm(T * N), nrow = T, ncol = N)
  Sigma = cov(rets)
  expect_equal(dim(Sigma), c(25, 25)) 
}) 
