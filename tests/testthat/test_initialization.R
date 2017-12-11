context("Test initialization()")

# 1
test_that("there are no NA's in initialization", {
  expect_equal(sum(is.na(initialization(30,30))),0)
  #expect_equal(sum(is.na(initialization(30,30))),1) # this obviously fails (but I did it to check)
})

# 2
test_that("only contain logical values in initialization",{
  inital_pop <- initialization(10, 15)
  inital_unique <- unique(sort(inital_pop))

  expect_equal(class(inital_pop), "matrix")
  expect_equal(class(inital_pop[,1]), "logical")

  expect_equal(inital_unique, c(F,T))
  expect_less_than(sum(inital_pop), length(inital_pop))
  expect_more_than(sum(inital_pop), 0)
})

# 3
test_that("dimension of initialization is C by P",{
  expect_equal(dim(initialization(10, 15)), c(10,15))
})
