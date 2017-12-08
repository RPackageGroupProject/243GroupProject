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

#' Select Parents to Breed Based on Given Selection Mechanism
#'
#' @description Selects parents based on specified selection mechanism. The
#' possibile selection mechanism includes (1)selecting both parents with
#' probability proportional to ranking, (2)selecting one parent with probability
#' proportional to ranking and one parent randomly, and (4)selecting with
#' method of tournament selection.
#'
#' @usage
#' selection(pop, fit, dat, method = 1)
#'
#' @param pop boleans matrix determined by \code{GA::initialization()}.
#'
#' @param fitScore fitness scores calculated by \code{GA::fitness()}. It must
#' be the case that a lower fitness score means a better model.
#'
#' @param offspringNum number of offspring generated to update the generation.
#'
#' @param dat data frame containing the predictors in the model.
#' First column should be the response variable.
#'
#' @param method the selection mechanism that user wants to apply to select
#' parents, can be choosen from 1 to 3; 1 indicates selecting both parents with
#' probability proportional to ranking; 2 indicates selecting one parent with
#' probability proportional to ranking and one parent randomly, and 3 indicates
#' selecting with method of tournament selection.
#'
#' @details This function selects parents to breed based on the passed in
#' selection mechanism. Select offspringNum / 2 number of parents from the
#' passed in population pop, by specific selection mechanism and the fitness
#' scores fitScore obtained from \code{GA::fitness()}.
#'
#' @return Returns a list containing the index for each of parent 1, the index
#' for each of parent 2, the fittest individual, and the fitness scores for the
#' whole population.
#'
#' @examples
#' dat <- mtcars
#' pop <- initialization(ncol(dat) - 1)
#' model <- lm
#' fitnessFunction <- AIC
#' fitScore <- fitness(pop, dat, fitnessFunction, model)
#' offspringNum <- 10
#' selection(pop, fitScore, offspringNum, dat, method = 1)
selection <- function(pop, fitScore, offspringNum, dat, method){

  # Size of the generation
  P <- ncol(pop)

  # Sample from the P chromosomes, with weights specified in fitnessProb,
  # with replacement, to obtain a parent population of size offspringNum /2.
  # Note there are duplicates within the parent population.

  # Method 1: Select both parents with probability proportional to ranking
  if (method == 1) {

    # Compute a vector of probability weights proportional to ranking
    # Since lowest score is the best, take the reverse rank
    fitnessProb <- 2 * rank(-fitScore) / (P * (P + 1))

    # Index for parent 1
    indexParent1 <- sample(x = 1:P, size = offspringNum / 2, replace = T,
                           prob = fitnessProb)
    # Index for parent 2
    indexParent2 <- sample(x = 1:P, size = offspringNum / 2, replace = T,
                           prob = fitnessProb)

  }

  # Method 2: Select one parent with probability proportional to ranking, one
  # parent randomly
  else if (method == 2) {

    # Compute a vector of probability weights proportional to ranking
    # Since lowest score is the best, take the reverse rank
    fitnessProb <- 2 * rank(-fitScore) / (P * (P + 1))

    # Index for parent 1
    indexParent1 <- sample(x = 1:P, size = offspringNum / 2, replace = T,
                           prob = fitnessProb)

    # Index for parent 2
    indexParent2 <- sample(x = 1:P, size = offspringNum / 2, replace = T)

  }

  else {
    stop("Need to pass in a selection mechanism.")
  }

  # Build the selection result as a list
  selResult <- list()
  selResult$indexParent1<-indexParent1
  selResult$indexParent2<-indexParent2
  selResult$fittest <- pop[,which.max(fitnessProb)] # Keep a copy of the fittest individual.
  selResult$fitnessScores <- fitScore # keep for plotting

  # Return the result
  return(selResult)
}


#######################################
####  tournament selection (need to be edited)
#######################################

#' Select Parents Generated From \code{GA::initialization()}
#' Based on Given Fitness Function
#'
#' @description Selecting parents based on fitness ranks, with AIC as the
#' default fitness criteria. Alternatively, we can use tournament selection.
#'
#' @usage
#' tournament(K)
#'
#' @param K integer that smaller than P, to determine the number of groups
#'
#' @param
#'
#' @param dat data frame containing the predictors in the model.
#' First column should be the response variable.
#'
#' @param model the linear model that user wants to use to fit in the data,
#' can be either \code{lm} or \code{glm}.
#'
#' @param ... additional arguments to pass to the model function.
#'
#' @details This function selects parents to breed based on tournament selection.
#' model to each possible parent generated by \code{GA::initialization()}
#' based on fitness rank with AIC as the default fitness function.
#' Randomly partition the last generation population into K disjoint subset of equal size.
#' The best individual in each group is chosen as a parent. Addtional partitioning are
#' carryed out until sufficient parent have been generated.
#'
#' @return Returns a list that contain two bolean matrix to replace the origin pop matrix
#'
tournament <- function(K, C, P, ...){

  #intialize 5*P number population
  pop <- initialization(C, 5P)

  P <- length(ncol(pop))

  #define function to find the index of the best fit in the assigned group
  myFit <- function(k){
    ind <- order(fitness(pop[,GroupInd[,k]], dat))[1]
    return(ind)
  }

  #initialize some index needed in while loop
  goodFitTol <- rep(0, P)
  m <- 1:K
  choose <- 0
  numRest <- P

  while(choose < P){

    #sample the index of the population to implement the randomly choosing
    sampInd <- sample(numRest- numRest %% k)

    #Partition into K disjoint subsets, each column represents a group
    GroupInd <- matrix(sampInd, ncol=k)

    #Find the index of best individual in each group
    goodFit <- sapply(m, myFit)+floor(numRest/k)*(m-1)

    #count the total number of good fit indivduals
    choose <- choose + K

    #remove the indivdual that was chosen
    sample <- sampInd[-goodFit]

    #combine the good fit individuals' index to a vector
    goodFitTol <- cbind(goodFitTol, goodFit)

    #number of rest poplation will be used in the next iteration
    numRest <- numRest-K
  }

  #if the number of choose larger than the number of population(P) we need, choose the first P
  if(choose > P ){
    goodFitTol <- goodFitTol[1:p]
  }
  #Randomly pair the parents
  randomPairInd <- sample(1:P, P/2)

  allParent1 <- pop[,randomPairInd]
  allParent2 <- pop[,-randomPairInd]

  sel.result <- list()
  sel.result$allParent1<-allParent1
  sel.result$allParent2<-allParent2

  return(sel.result)
}


#######################################
####  Mutation
#######################################

#' Chromosome Mutation
#'
#' @description Make a chromosome with fixed probability 0.01 to mutate.
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

  # Randomly select a point as the point of crossover.
  k <- sample(1 : (C1 - 1), 1)

  # Return the crossed chromosomes together
  cbind(
    mutation(c(chr1[1 : k], chr2[(k + 1) : C1])),
    mutation(c(chr2[1 : k], chr1[(k + 1) : C1]))
  )
}


#######################################
####  Choose Next Generation *need to check
#######################################

#' Breed the Selected Parents Generated From \code{GA::selection()}
#'
#' @usage
#' nextGeneration(pop,selResult)
#'
#' @param pop boleans matrix determined by \code{GA::initialization()}
#'
#' @param selResult list returnd by \code{GA::selection()}
#' At the minimum, it needs to contain indices for the selected parents 1 and 2.
#'
#' @param G numeric number in the range of (0, 1].
#'
#' @param ... additional arguments to pass to the model function. DID NOT IMPLEMENT
#'
#' @details Breeding uses \code{GA::crossover()} for each pair of parents
#' The Generation Gap G is the proportion of the generation to be replaced by generated offspring.
#' If the number of the selected population produced by generation gap is not an even number, choose the
#' nearest larger even integer.
#'
#' Producing the G*P number of offspring to replace the least fit individual in parent generation
#'
#' @return Returns a list that contain two bolean matrix to replace the origin pop matrix
#'
#' @examples
#' dat <- mtcars # use the built-in mtcars data
#' C <- dim(dat)[2] - 1 #Number of variables
#' pop <- initialization(C) #produce boleans matrix
#' selResult <- selection(pop, f = AIC, dat)
#' pop <- nextGeneration(pop, selResult)
nextGeneration <- function(pop, selResult, G){

  # P <- length(ncol(pop))
  # selectPop <- P * G
  # offspringNum <- if ((selectPop) %% 2 == 0) selectPop else if(selectPop +1  >= P) selectPop else selectPop +1

  # Size of the generation
  P <- ncol(pop)

  # Obtain the number of offspring needed to be generated
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

  # Extract the parents
  allParent1 <- pop[,selResult$indexParent1]
  allParent2 <- pop[,selResult$indexParent2]

  # ???
  # lessFitInd1 <- order(fitness(allParent1), decreasing = TRUE)[1 : offspringNum / 2]
  # lessFitInd2 <- order(fitness(allParent2), decreasing = TRUE)[1 : offspringNum / 2]

  # Exttract the ones that needed to be replaced by offspring
  lessFitInd <- order(selResult$fitnessScores, decreasing = T)[1 : offspringNum]


  # Breed to create new generation
  for (i in offspringNum / 2){

    # Each time generate two offspring, stored in columns
    childrenChromes <- crossover(allParent1[, i], allParent2[, i])

    # Select first child and put into pop
    # we should consider keeping the fittest X% of parents rather than replacing them all
    # allParent1[,lessFitIndex1[i]] <- mutation(childrenChromes[,1])
    # allParent2[,lessFitIndex2[i]] <- mutation(childrenChromes[,2])

    # Select first child and put into pop
    pop[, lessFitInd[i]] <- childrenChromes[, 1]
    pop[, lessFitInd[i + offspringNum / 2]] <- childrenChromes[, 2]

  }

  # sel.result <- list()
  # sel.result$allParent1 <- allParent1
  # sel.result$allParent2 <- allParent2

  # return(sel.result)

  return(pop)
}


