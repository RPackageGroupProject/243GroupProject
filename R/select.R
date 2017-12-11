#######################################
####  Genetic Algorithm
#######################################

#' Uses a genetic algorithm for variable selection in either lm or glm models
#'
#' @usage
#' select(mod,dat,f)
#'
#' @param model the linear model that user wants to use to fit in the data,
#' can be either \code{lm} or \code{glm}.
#'
#' @param dat data frame containing the predictors in the model.
#' First column should be the response variable.
#'
#' @param f fitness function that takes in an lm or glm model and returns a
#' a numerical 'qualification' of that model. Users can choose AIC or BIC or
#' even define by themselves. If the f argument is missing, the default is AIC.
#'
#' @param ... additional arguments to pass to the model function. DID NOT IMPLEMENT YET
#'
#' @param G numeric number in the range of (0, 1].
#'
#' @details The algorithm:
#' (1) First initializes population,
#' For g generations; do:
#' (2) calculates fitness of models and selects parent pairs to breed
#' (3) breeds the parent pairs, choosing the first child
#' (4) replaces the parents with the children
#'
#' @return Returns a list with the fittest model and
#' a matrix of the population fitness across generations (useful for plotting)
#'
#' @export
#'
#' @examples
#' select(mtcars)

select <- function(dat,
                   P = NULL, numGens = NULL, G = NULL, fitnessFunction = NULL,
                   method = NULL, model = NULL, family = NULL){

  # Number of variables, same as chromosome length
  C <- ncol(dat) - 1

  # Make default inputs
  if (is.null(P)) {P <- as.integer(1.5 * C)}
  if (is.null(numGens)) {numGens <- 100}
  if (is.null(G)) {G <- 0.5}
  if (is.null(fitnessFunction)) {fitnessFunction <- AIC}
                                                  # function(fit,...){return(extractAIC(fit,...)[2])}
  if (is.null(method)) {method <- 1}
  if (is.null(model)) {model <- lm}
  if (is.null(family)) {family <- gaussian}

  # Check if inputs are valid


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
                         model = model, dat = dat)

    # Selection of parents
    selResult <- selection(pop = pop, fitScores = fitScores,
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

    print(paste('Generation: ', gen, ' Fitness: ', min(selResult$fitnessScores)))

    # select next generation by updating the current generation
    pop <- nextGeneration(pop, selResult, offspringNum)

  }

  # fit the best model
  response <- colnames(dat)[1]
  predictors <- colnames(dat)[-1]
  chosen <- pop[, which.min(fitnessScores[, numGens])]

  # fit the best model
  form <- as.formula(paste(response, "~",
                           paste(predictors[chosen], collapse = "+")))
  fittestMod <- model(formula = form, data = dat)
  fittestScore <- fitnessFunction(model(formula = form, data = dat))

  # return
  genResult <- list()
  genResult$fittestModel <- fittestMod
  genResult$fittestScore <- fittestScore
  genResult$fitnessScores <- fitnessScores
  return(genResult)
}




