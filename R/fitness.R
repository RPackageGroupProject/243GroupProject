#######################################
#### Calculate Fitness
#######################################

#' Calculate Fitness based on Fitness Function
#'
#' @description Calculate the fitness value for each individuals in the
#' generation based on the passed in fitness function.
#'
#' @usage fitness(pop, dat, fitnessFunction, model)
#'
#' @param pop boleans matrix determined by \code{GA::initialization()}
#'
#' @param X data matrix with rows as observations and columns as predictors in the model.
#' @param y response variable vector
#'
#' @param fitnessFunction fitness function that takes in an lm or glm model and
#' returns a numerical fitness of that model. Users can choose AIC or BIC or
#' even define by themselves. If the argument is missing, the default is AIC.
#'
#' @param model the linear model that user wants to use to fit in the data,
#' can be either \code{lm} or \code{glm}.
#'
#' @return Returns a matrix containing one row with \code{ncol(pop)}
#' observations of the fitness scores of each chromosomes.

fitness <- function(pop,y,X,fitnessFunction,model) {

  # Number of chromosomes
  P <- ncol(pop)

  # Number of genes must equal number of predictors
  stopifnot(nrow(pop) == ncol(X))

  # place holder vector for score on fitness function
  fit <- c()

  for (i in 1:P) {

    # Select columns of data matrix X, based on chromosome
    chosen <- pop[,i]
    Xsel<-X[,chosen]

    # check for null model
    if (ncol(Xsel)==0){Xsel=1}

    # Calculate the fitness using lm() or glm()
    score <- fitnessFunction(model(y~Xsel, data = dat))
    stopifnot(length(score)==1)

    # update the place holder vector
    fit <- c(fit, score)
  }

  # convert the result to a matrix
  return(fit)
}
