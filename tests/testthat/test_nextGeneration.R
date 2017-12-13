context("Test nextGeneration()")

# 1
test_that("nextGeneration() returns a matrix with the same column number of the previous matirx", {
  
  require(MASS)
  dat <- Boston
  fitnessFunction <- AIC
  model <- lm
  y<-as.matrix(Boston$medv)
  X<-as.matrix(Boston[1:13])
  pop <- initialization(ncol(Boston) - 1, 30)
  
  fitScore <- fitness(pop, y, X, fitnessFunction, model, dat)
  selResult<- selection(pop, fitnessFunction, model, fitScore, offspringNum = 6, method = 3, dat, K=3)
  popNext <- nextGeneration(pop, selResult, offspringNum = 6)
  
  expect_equal(ncol(popNext), ncol(pop))
  
})



# 2
test_that("nextGeneration() returns a matrix not a vector", {
  
  require(MASS)
  dat <- Boston
  fitnessFunction <- AIC
  model <- lm
  y<-as.matrix(Boston$medv)
  X<-as.matrix(Boston[1:13])
  pop <- initialization(ncol(Boston) - 1, 30)
  
  fitScore <- fitness(pop, y, X, fitnessFunction, model, dat)
  selResult<- selection(pop, fitnessFunction, model, fitScore, offspringNum = 6, method = 3, dat, K=3)
  popNext <- nextGeneration(pop, selResult, offspringNum = 6)
  
  expect_equal(is(popNext, 'matrix'), TRUE)
  
})



# 3
test_that("nextGeneration() returns a matrix that is no more than offspringNum different columns compared with previous matrix", {
  
  require(MASS)
  dat <- Boston
  fitnessFunction <- AIC
  model <- lm
  offspringNum <- 6
  y<-as.matrix(Boston$medv)
  X<-as.matrix(Boston[1:13])
  pop <- initialization(ncol(Boston) - 1, 30)
  
  fitScore <- fitness(pop, y, X, fitnessFunction, model, dat)
  selResult<- selection(pop, fitnessFunction, model, fitScore, offspringNum, method = 3, dat, K=3)
  popNext <- nextGeneration(pop, selResult, offspringNum)
  fitScoreNext <- fitness(popNext, y, X, fitnessFunction, model, dat)
  
  numChange <- length(which(colSums(abs(pop - popNext))!=0))
  
  expect_equal(numChange <= offspringNum, TRUE)
  
})




# 4
test_that("nextGeneration() returns a matrix that update the fitness value", {
  
  require(MASS)
  dat <- Boston
  fitnessFunction <- AIC
  model <- lm
  y<-as.matrix(Boston$medv)
  X<-as.matrix(Boston[1:13])
  pop <- initialization(ncol(Boston) - 1, 30)
  
  fitScore <- fitness(pop, y, X, fitnessFunction, model, dat)
  selResult<- selection(pop, fitnessFunction, model, fitScore, offspringNum = 6, method = 3, dat, K=3)
  popNext <- nextGeneration(pop, selResult, offspringNum = 6)
  fitScoreNext <- fitness(popNext, y, X, fitnessFunction, model, dat)
  
  expect_equal(sum(fitScoreNext) <= sum(fitScore), TRUE)
  
})
