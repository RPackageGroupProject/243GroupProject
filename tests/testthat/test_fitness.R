

test_that("fitness() returns a lower value for a known better model than a known worse model", {

  # two known models from mtcars dataset
  mod_best<-c(FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,FALSE,TRUE,FALSE,FALSE)
  mod_other<-c(FALSE,TRUE,TRUE,FALSE,FALSE,TRUE,FALSE,TRUE,FALSE,FALSE)
  known_pop<-cbind(mod_best,mod_other)

  # make sure when I tweak the fitness function, I get mod_best lower AIC than mod_other.
  fit<-fitness(known_pop,mtcars,AIC,lm)
  expect_equal(fit[1]<=fit[2],TRUE)

})

test_that("fitness() returns a vector and not a matrix"){
  pop <- initialization(ncol(dat) - 1)
  model<-lm
  fitnessFunction<-AIC
  data<-mtcars
  X<-as.matrix(mtcars[-1])
  y<-as.matrix(mtcars[1])
  fit<-fitness(pop,y,X,fitnessFunction,model)
  expect_equal(is(fit,'vector'),TRUE)
}
