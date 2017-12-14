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
#'
#' @param y response variable vector
#'
#' @param fitnessFunction fitness function that takes in an lm or glm model and
#' returns a numerical fitness of that model. Users can choose AIC or BIC or
#' even define by themselves.
#'
#' @param model the linear model that user wants to use to fit in the data,
#' can be either \code{lm} or \code{glm}.
#'
#' @param dat data frame containing the predictors in the model.
#' First column should be the response variable.
#'
#' @param ... additional arguments to pass to regression model
#'
#' @return Returns a matrix containing one row with \code{ncol(pop)}
#' observations of the fitness scores of each chromosomes.

fitness <- function(pop, y, X, fitnessFunction, model, dat,...) {

  # Number of chromosomes
  P <- ncol(pop)

  # Number of genes must equal number of predictors
  stopifnot(nrow(pop) == ncol(X))

  # place holder vector for score on fitness function
  fit <- c()

  for (i in 1 : P) {

    # Select columns of data matrix X, based on chromosome
    chosen <- pop[, i]
    XSel <- as.matrix(X[, chosen])


    # Check for null model and calculate the fitness
    if (ncol(XSel) == 0) {
      score <- fitnessFunction(model(y ~ 1, data = dat,...))
    } else {
      score <- fitnessFunction(model(y ~ XSel, data = dat,...))
    }

    # Check
    stopifnot(length(score) == 1)

    # update the place holder vector
    fit <- c(fit, score)

  }

  # return the result
  return(fit)
}



