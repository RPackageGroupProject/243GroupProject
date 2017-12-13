
test_that("nextGeneration() returns a matrix with the same column number of the previous matirx", {

  library(MASS)
  X<-as.matrix(Boston$medv)
  y<-as.matrix(Boston[1:13])
  pop <- initialization(ncol(Boston) - 1, 30)

  fitScore <- fitness(pop, X, y, fitnessFunction = AIC, model = lm)
  selResult<- selection(pop, fitScore, offspringNum = 6, dat=Boston, method = 3, K=3)
  popNext <- nextGeneration(pop, selResult, offspringNum = 6)

  expect_equal(ncol(popNext), ncol(pop))

})




test_that("nextGeneration() returns a matrix that uodate the fitness value", {

  library(MASS)
  X<-as.matrix(Boston$medv)
  y<-as.matrix(Boston[1:13])
  pop <- initialization(ncol(Boston) - 1, 30)

  fitScore <- fitness(pop, X, y, fitnessFunction = AIC, model = lm)
  selResult<- selection(pop, fitScore, offspringNum = 6, dat=Boston, method = 3, K=3)
  popNext <- nextGeneration(pop, selResult, offspringNum = 6)
  fitScoreNext <- fitness(popNext, X, y, fitnessFunction = AIC, model = lm)

  expect_equal(fitScoreNext <= fitScore, TRUE)

})




test_that("nextGeneration() returns a matrix that is no more than offspringNum different columns compared with previous matrix", {

  library(MASS)
  X<-as.matrix(Boston$medv)
  y<-as.matrix(Boston[1:13])
  pop <- initialization(ncol(Boston) - 1, 30)

  fitScore <- fitness(pop, X, y, fitnessFunction = AIC, model = lm)
  selResult<- selection(pop, fitScore, offspringNum = 6, dat=Boston, method = 3, K=3)
  popNext <- nextGeneration(pop, selResult, offspringNum = 6)
  fitScoreNext <- fitness(popNext, X, y, fitnessFunction = AIC, model = lm)

  numChange <- length(which(colSums(abs(pop - popNext))!=0))

  expect_equal(numChange <= offspringNum, TRUE)

})
