context("Test fitness()")

#1
test_that("fitness() returns a lower value for a known better model than a known worse model", {

  # two known models from mtcars dataset
  mod_best <- c(F, F, F, F, T, T, F, T, F, F)
  mod_other <- c(F, T, T, F, F, T, F, T, F, F)
  known_pop <- cbind(mod_best, mod_other)

  # make sure when I tweak the fitness function, I get mod_best lower AIC than mod_other.
  X <- as.matrix(mtcars[-1])
  y <- as.matrix(mtcars[1])
  fit <- fitness(known_pop, y, X, AIC, lm, mtcars)

  expect_equal(fit[1] <= fit[2], TRUE)

})

# 2
test_that("fitness() returns a vector and not a matrix", {

  C <- ncol(mtcars) - 1
  P <- as.integer(1.5 * C)
  pop <- initialization(C, P)
  X <- as.matrix(mtcars[-1])
  y <- as.matrix(mtcars[1])
  fitScores <- fitness(pop, y, X, AIC, lm, mtcars)
  expect_equal(is(fitScores, 'vector'), TRUE)

})
