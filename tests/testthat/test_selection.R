context("Test selection()")

#1
test_that("selection returns a list", {

  dat <- mtcars
  X <- as.matrix(mtcars[-1])
  y <- as.matrix(mtcars[1])
  C <- ncol(dat) - 1
  P <- as.integer(1.5 * C)
  pop <- initialization(C, P)
  model <- lm
  fitnessFunction <- AIC
  fitScores <- fitness(pop, y, X, fitnessFunction, model, dat)
  offspringNum <- 2
  selResult <- selection(pop, fitnessFunction, model, fitScores, offspringNum, method = 1, dat, K = 2)

  expect_equal(is(selResult, 'vector'), TRUE)
  expect_equal(length(selResult$indexParent1), 1)

})

# 2
test_that("number of parent1 and parent2 is half of the size of offspringNum", {

  dat <- mtcars
  X <- as.matrix(mtcars[-1])
  y <- as.matrix(mtcars[1])
  C <- ncol(dat) - 1
  P <- as.integer(1.5 * C)
  pop <- initialization(C, P)
  model <- lm
  fitnessFunction <- AIC
  fitScores <- fitness(pop, y, X, fitnessFunction, model, dat)
  offspringNum <- 2
  selResult <- selection(pop, fitnessFunction, model, fitScores, offspringNum, method = 1, dat, K = 2)


  expect_equal(length(selResult$indexParent1), 1)
  expect_equal(length(selResult$indexParent2), 1)

})
