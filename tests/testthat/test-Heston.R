
test_that("function works", {
  expect_no_error( HESTON(rho=-0.5, sigma=0.1, lambda=2, t=10, r=0.02, theta=0.1^2, V0=0.1^2, x0=log(100), n_steps=5, paths=10))
})

test_that('Negative values breaks',{
  expect_error( HESTON(rho=-0.5, sigma=0.1, lambda=2, t=10, r=0.02, theta=0.1^2, V0=0.1^2, x0=log(100), n_steps=5, paths=-1))
  expect_error( HESTON(rho=-0.5, sigma=0.1, lambda=2, t=10, r=0.02, theta=0.1^2, V0=0.1^2, x0=log(100), n_steps=-1, paths=10))
  expect_error(HESTON(rho=-0.5, sigma=0.1, lambda=2, t=-1, r=0.02, theta=0.1^2, V0=0.1^2, x0=log(100), n_steps=5, paths=10))
})

test_that('Missing sigma breaks',{
  expect_error( mcLIBOR(rep(0.5,5),-1,matrix(NA,5,5),10))
})
test_that('Negative sigma breaks',{
  expect_error( mcLIBOR(rep(0.5,5),-1,matrix(-1,5,5),10))
})

test_that('Incorrect dimensions break',{
  expect_error( mcLIBOR(rep(0.5,4),-1,matrix(-1,5,5),10))
})
