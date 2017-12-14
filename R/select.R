#######################################
####  Genetic Algorithm
#######################################

#' Uses a genetic algorithm for variable selection in either lm or glm models
#'
#' @usage
#' select(dat)
#'
#' @param dat data frame containing the predictors in the model.
#' First column should be the response variable.
#'
#' @param P number of chromosomes, same as the size of generation. If not
#' specified, the default is set to be 1.5 * C where C is chromosome length.
#'
#' @param numGens total number of generations, default to be 100.
#'
#' @param G the proportion of the current generation to be replaced by the
#' offspring to construct the next generation, should be a numeric number in
#' the range of (0, 1], default to be 0.2 which is 20 percent.
#'
#' @param fitnessFunction fitness function that takes in an lm or glm model and
#' returns a numerical fitness of that model. Users can choose AIC or BIC or
#' even define by themselves, but need to make sure a lower fitness scores
#' indicates the corresponding model is better.
#'
#' @param method the selection mechanism that user wants to apply to select
#' parents, can be choosen from 1 to 3; 1 indicates selecting both parents with
#' probability proportional to ranking; 2 indicates selecting one parent with
#' probability proportional to ranking and one parent randomly, and 3 indicates
#' selecting with method of tournament selection; default is method 1.
#'
#' @param model the linear model that user wants to use to fit in the data,
#' can be either \code{lm} or \code{glm}; default to be \code{lm}.
#'
#' @param K number of groups to partition the population into; default is 2.
#'
#' @param verbose logical; if TRUE (default) prints the progress of algorithm
#'
#' @param ... additional arguments to pass to regression model
#'
#' @details The algorithm:
#' (1) First initializes population,
#' For g generations; do:
#' (2) calculates fitness of models and selects parent pairs to breed
#' (3) breeds the parent pairs, obtain the children
#' (4) replaces the least fit individuals in current generation with the
#'     children to obtain the next generation
#'
#' @return Returns a list with the fittest model and the corresponding fitness
#' score, together with a matrix of the population fitness across generations
#' (useful for plotting)
#'
#' @export
#'
#' @examples
#' select(mtcars)

select <- function(dat,
                   P = NULL, numGens = NULL, G = NULL, fitnessFunction = NULL,
                   method = NULL, model = NULL, K = NULL, verbose=NULL,...){

  # Number of variables, same as chromosome length
  C <- ncol(dat) - 1

  # Make default inputs
  if (is.null(P)) {P <- as.integer(1.5 * C)}
  if (is.null(numGens)) {numGens <- 100}
  if (is.null(G)) {G <- 0.5}
  if (is.null(fitnessFunction)) {fitnessFunction <- AIC}
  if (is.null(method)) {method <- 1}
  if (is.null(model)) {model <- lm}
  if (is.null(K)) {K <- 2}
  if (is.null(verbose)){verbose <- TRUE}

  # Check if inputs are valid
  if (is.null(dat)) {stop("Need to pass in data.")}
  if (P <= 0 || P != as.integer(P)) {stop("P needs to be a whole number greater than 0.")}
  if (numGens <= 0 || numGens != as.integer(numGens)) {
    stop("numGens needs to be a whole number greater than 0.")
    }
  if (numGens < 30) {warning("Number of generations is too small.")}
  if (G <= 0 || G > 1) {
    stop("Generation gap should be a numeric number in the range of (0, 1]")
  }
  if (!is.function(fitnessFunction)) {
    stop("Need to pass in a function into the argument fitnessFunction.")
  }
  if (method != 1 & method != 2 & method != 3 ) {
    stop("The method argument only takes in 1, 2, or 3.")
  }
  if (!identical(model, lm) && !identical(model, glm)) {
    stop("The regression model should be either lm or glm.")
  }


  # Initialize population
  pop <- initialization(C = C, P = P) # generate random starting population

  # Obtain the number of offspring (offspringNum) needed to be generated
  selectPop <- ceiling(P * G)
  if (selectPop %% 2 == 0) {
    offspringNum <- selectPop
  }
  else if (selectPop + 1 >= P) {
    offspringNum <- selectPop - 1
  }
  else {
    offspringNum <- selectPop + 1
  }

  # Obtain the inputs X and y for fitness()
  X <- as.matrix(dat[-1]) # predictors matrix
  y <- as.matrix(dat[1])  # response variable vector


  # Loop over the genetic algorithm
  for (gen in seq(numGens)) {

    # Obtain fitness scores for each model
    fitScores <- fitness(pop = pop, y= y, X = X,
                         fitnessFunction = fitnessFunction,
                         model = model, dat = dat,...)

    # Selection of parents
    selResult <- selection(pop = pop, fitnessFunction = fitnessFunction,
                           model = model, fitScores = fitScores,
                           offspringNum = offspringNum,
                           method = method, dat = dat, K = K)

    # store fitness for each generation to check algorithm is improving
    if (gen == 1) {
      fitnessScores <- matrix(selResult$fitnessScores, ncol = 1)
      }
    else {
      fitnessScores <- cbind(fitnessScores,
                             matrix(selResult$fitnessScores, ncol = 1))
      }

    if((gen %% 10==0) & (verbose)){
      print(paste('Generation: ', gen, ' Fitness: ', min(selResult$fitnessScores)))
    }

    # select next generation by updating the current generation
    # do not update if it is the last generation
    if (gen < numGens) {
      pop <- nextGeneration(pop, selResult, offspringNum)
    }
    else {
      print("ALgorithm Ends")
    }

  }

  # fit the best model
  response <- colnames(dat)[1]
  predictors <- colnames(dat)[-1]
  chosen <- pop[, which.min(fitnessScores[, numGens])]

  # fit the best model
  form <- as.formula(paste(response, "~",
                           paste(predictors[chosen], collapse = "+")))
  fittestMod <- model(formula = form, data = dat,...)
  fittestScore <- fitnessFunction(model(formula = form, data = dat,...))

  # return result as a list
  genResult <- list()
  genResult$fittestModel <- fittestMod
  genResult$fittestScore <- fittestScore
  genResult$fitnessScores <- fitnessScores
  return(genResult)
}




