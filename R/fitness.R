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
#' @param dat data frame containing the predictors in the model.
#' First column should be the response variable.
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
#' 
#' @examples 
#' dat <- mtcars
#' pop <- initialization(ncol(dat) - 1)
#' model <- lm
#' fitnessFunction <- AIC
#' fitness(pop, dat, fitnessFunction, model)
fitness <- function(pop, dat, fitnessFunction, model) {
  
  # Number of chromosomes
  P <- ncol(pop)
  
  # Variables are columns of the dataset
  response <- colnames(dat)[1]
  predictors <- colnames(dat)[-1]
  
  # Input parameters must be named matrices
  stopifnot(!is.null(P) && !is.null(response) && !is.null(predictors))
  
  # Number of genes must equal number of predictors
  stopifnot(nrow(pop) == length(predictors))
  
  # place holder vector for score on fitness function
  fit <- c()
  
  # Process each chromosome 
  for (i in 1:P) {
    
    # Build a formula using the chosen predictors
    chosen <- pop[,i]
    
    # When at least one predictor is active
    if(sum(chosen)) {
      form <- as.formula(paste(response, "~",
                               paste(predictors[chosen], collapse = "+")))
      
      # Calculate the fitness using the provided fitness function
      score <- fitnessFunction(model(formula = form, data = dat))
      
      # update the place holder vector
      fit <- c(fit, score)
    }
    
    # When all predictors are inactive
    else {
      form <- as.formula(paste(response, "~", 1))
      
      # Calculate the fitness using the provided fitness function
      fitnessFunction(model(formula = form, data = dat))
      
      # update the place holder vector
      fit <- c(fit, score)
    }
  }
  
  # convert the result to a matrix
  return(matrix(fit, nrow = 1))
  
}




