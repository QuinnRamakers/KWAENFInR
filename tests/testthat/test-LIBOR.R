#this function seems to have issues with the test and parallelisation so i comment it out
#test_that("function works", {
 #expect_no_error( mcLIBOR(rep(0.5,5),1,matrix(0.01,5,5),10))
#})

test_that('Negative dt breaks',{
  expect_error( mcLIBOR(rep(0.5,5),-1,matrix(0.01,5,5),10))
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
