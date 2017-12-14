context("Test select()")

#1
test_that("select returns a list", {

  dat <- mtcars
  C <- ncol(dat) - 1
  P <- as.integer(1.5 * C)
  pop <- initialization(C, P)
  y <- as.matrix(dat[1])
  X <- as.matrix(dat[-1])
  fitnessFunction <- AIC
  model <- lm
  fitScores <- fitness(pop, y, X, fitnessFunction, model, dat)
  offspringNum <- 4
  selResult <- selection(pop, fitnessFunction, model, fitScores, offspringNum, method = 1, dat, K = 2)
  nextGeneration(pop, selResult, offspringNum)

  expect_equal(is(select(dat), 'list'), TRUE)

})
