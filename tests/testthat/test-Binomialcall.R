test_that("Function  works", {
  expect_no_error(Binomial_call(1,100,100,0.15,0.02,10))
})
test_that("Negative time value does not work", {
  expect_error(Binomial_call(1,100,100,0.15,0.02,-1))
})
test_that("Negative steps does not work", {
  expect_error(Binomial_call(-1,100,100,0.15,0.02,10))
})
