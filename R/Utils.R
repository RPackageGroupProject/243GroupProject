#######################################
####  Initialization
#######################################

#' Randomly Generate the Initial Generation
#'
#' @description Bernoulli sampling to form the initial generation of the
#' Genetic Algorithm. The size of the initial generation is set to be the
#' integer closest to P = 1.5*C.
#'
#' @usage initialization(C)
#'
#' @param C chromosome length (number of predictor variables in the data)
#'
#' @details This function produces initial generation given the chromosome length
#' for Genetic Algorithm. The row of the generated boolean matrix represents the
#' locus of a gene(variable), and the column of the matrix represents different
#' members of the first generation(different models). The number of first
#' generation is defined to be integer closest to 1.5C.
#'
#' @return A Bolean Matrix with dimension C by 1.5*C where each column
#' representing a chromosome, in which T marks that the gene (variable)
#' as active and F as inactive.
#'
#' @examples
#' initialization(10)
initialization <- function(C, P = as.integer(1.5 * C)){

  # Bernoulli sampling T or F to obtain C of them
  # A T in a locus of a gene (variable) means the perticular variable is included
  replicate(P, sample(c(F,T), size = C, replace = T))

}


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


#######################################
####  Parents Selection
#######################################

#' Select Parents Generated From \code{GA::initialization()}
#' Based on Given Fitness Function
#'
#' @description Selecting parents based on fitness ranks, with AIC as the
#' default fitness criteria. Alternatively, we can use tournament selection.
#'
#' @usage
#' selection(pop, f = AIC, dat)
#' selection(pop, f = BIC, dat)
#'
#' @param pop boleans matrix determined by \code{GA::initialization()}
#'
#' @param f fitness function that takes in an lm or glm model and returns a
#' a numerical 'qualification' of that model. Users can choose AIC or BIC or
#' even define by themselves. If the f argument is missing, the default is AIC.
#'
#' @param dat data frame containing the predictors in the model.
#' First column should be the response variable.
#'
#' @param model the linear model that user wants to use to fit in the data,
#' can be either \code{lm} or \code{glm}.
#'
#' @param ... additional arguments to pass to the model function.
#'
#' @details This function selects parents to breed by fitting linear regression
#' model to each possible parent generated by \code{GA::initialization()}
#' based on fitness rank with AIC as the default fitness function.
#'
#' @return Return a list, where containing chosen parents, fittest parent
#' and fitting score for each parent based on fitting function.
#'
#' @examples
#' dat <- mtcars # use the built-in mtcars data
#' C <- dim(dat)[2]-1 #Number of variables
#' pop <- initialization(C) #produce boleans matrix
#' selection(pop, f = AIC, dat)
selection <- function(pop, f, dat, model, ...){

  # number of chromosomes
  P <- ncol(pop)

  # Variables are columns of the dataset
  response <- colnames(dat)[1]
  predictors <- colnames(dat)[-1]

  # Input parameters must be named matrices
  stopifnot(!is.null(P) && !is.null(response) && !is.null(predictors))

  # Number of genes must equal number of predictors
  stopifnot(nrow(pop)==length(predictors))

  # Process each chromosome in parallel
  fit <- foreach::`%dopar%`(foreach::foreach(i = 1:P, .combine = c),{

    # Build a formula using the chosen predictors
    chosen <- pop[,i]

    # Account for situation when all predictors are F
    if(sum(chosen)){
      form <- as.formula(paste(response, "~",
                               paste(predictors[chosen], collapse = "+")))

      # Calculate the fitness using the provided fitness function
      #f(model(formula = form, data = dat,...))
      f(model(formula = form, data = dat))
    }else{
      Inf
    }
  })

  # Compute a vector of probability weights
  # Since lowest fitness is the best, take the reverse rank
  fitness <- 2 * rank(-fit) / (P * (P + 1))

  # Selecting Parents Method 1: Select both parents proportional to their fitness
  # Sample from the P chromosomes, with weights specified in fitness,
  # with replacement, to generate a parenting population of size P.
  # Note there are duplicates within the parenting population.
  parent1_ind <- sample(x = 1:P, size = P, replace = T, prob = fitness)
  parent2_ind <- sample(x = 1:P, size = P, replace = T, prob = fitness)

  # Selecting parents Method 1: Select both parents proportional to their fitness
  #### The second method being using selections probabilities proportional to ranking
  #### to select one parent and randonly select one parent

  #parent1_ind <- sample(x = 1:P, size = P, replace = T, prob = fitness)
  #parent2_ind <- sample(x = 1:P, size = P, replace = T)

  # And maybe we should implement the tournament selection?


  # Build the selection result as a list
  sel.result <- list()
  #sel.result$parents <- parents  # we don't need to return the parent chromosomes here. They are in pop. We can just return the indices
  sel.result$parent1_ind<-parent1_ind
  sel.result$parent2_ind<-parent2_ind
  sel.result$fittest <- pop[,which.max(fitness)] # Keep a copy of the fittest individual.
  sel.result$fit <- fit # keep for plotting

  # Return the result
  return(sel.result)
}


#######################################
####  Mutation
#######################################

#' Chromosome Mutation
#'
#' @description Make a chromosome with fixed probability 0.01 to mutate
#'
#' @usage mutation(chr)
#'
#' @param chr a logical vector representing an individual chromosome.
#'
#' @details This function makes a chromosome have a 1% fixed chance to mutate
#' in each locale. If mutation happens at one locale, it will make the value
#' in that locale from T to F(or F to T).
#'
#' @return Return a mutated chromosome vector with the same length as input one.
#'
#' @examples
#' ind<-initialization(10)[,1]
#' mutation(ind)
mutation <- function(chr){

  # For each element of chr determine whether to mutate with 1% prob
  mutate <- sample(c(T,F), length(chr), prob = c(0.01, 0.99), replace = T)

  # 'exclusive or' will toggle F to T and T to F
  xor(chr, mutate)

}


#######################################
####  Crossover
#######################################

#' Genetic Operatior: Chromosome Crossover and Mutation
#'
#' @description Make 2 individual parent chromosomes crossover and
#' mutate when breeding offsprings.
#'
#' @usage crossover(chr1,chr2)
#'
#' @param chr1,chr2 a numeric vectors represents individual parent chromosome.
#'
#' @details This function makes two individual parent chromosomes
#' perfrom crossover and mutation when breeding next generation.
#' Note that the crossover is simply one-point crossover
#' and the mutation is based on \code{GA::mutation()}.
#'
#' @return A matrix with each column representing the
#' offsprings from a process of crossover and mutation.
#'
#' @examples
#' ind1<-initialization(10)[,1]
#' ind2<-initialization(10)[,2]
#' crossover(ind1,ind2)
crossover <- function(chr1, chr2){

  C1 <- length(chr1)
  C2 <- length(chr2)

  # Make sure we're working with equal-sized vectors
  stopifnot(C1 == C2)

  #randomly select a point as the point of crossover.
  k <- sample(1:(C1-1), 1)

  # Return the crossed chromosomes together
  cbind(
    mutation(c(chr1[1:k], chr2[(k + 1):C1])),
    mutation(c(chr2[1:k], chr1[(k + 1):C1]))
  )
}


#######################################
####  Choose Next Generation
#######################################

#' Breed the Selected Parents Generated From \code{GA::selection()}
#'
#' @usage
#' nextGeneration(pop,selectedParents)
#'
#' @param pop boleans matrix determined by \code{GA::initialization()}
#'
#' @param selectedParents list returnd by \code{GA::selection()}
#' At the minimum, it needs to contain indices for the selected parents 1 and 2.
#'
#' @param ... additional arguments to pass to the model function. DID NOT IMPLEMENT
#'
#' @details Breeding uses \code{GA::crossover()} for each pair of parents
#' The first child is chosen to replace that pair of parents
#' Children chromosomes replace all parent chromosomes in the population
#'
#' @return Returns a bolean matrix to replace the origin pop matrix
#'
#' @examples
#' dat <- mtcars # use the built-in mtcars data
#' C <- dim(dat)[2] - 1 #Number of variables
#' pop <- initialization(C) #produce boleans matrix
#' selectedParents <- selection(pop, f = AIC, dat)
#' pop <- nextGeneration(pop, selectedParents)
nextGeneration <- function(pop, selectedParents){

  chromes1 <- pop[,selectedParents$parent1_ind] # ??? - this should be making a copy
  chromes2 <- pop[,selectedParents$parent2_ind]

  # Breed to create new generation
  for (pair_ind in seq(ncol(pop))){
    # XXX - can we vectorize cross-over somehow.
    # XXX - it'd be nice to just use C instead of ncol(pop)

    # breed parents
    children_chromes <- crossover(chromes1[,pair_ind], chromes2[,pair_ind])

    # select first child, mutate, and put into pop
    # we should consider keeping the fittest X% of parents rather than replacing them all
    pop[,pair_ind] <- mutation(children_chromes[,1])

  }
  return(pop)

}
