#######################################
####  Parents Selection
#######################################

#' Select Parents to Breed Based on Given Selection Mechanism
#'
#' @description Selects parents based on specified selection mechanism. The
#' possibile selection mechanism includes (1)selecting both parents with
#' probability proportional to ranking, (2)selecting one parent with probability
#' proportional to ranking and one parent randomly, and (3)selecting with
#' method of tournament selection.
#'
#' @usage
#' selection(pop, fitScores, offspringNum, method = 1, K)
#'
#' @param pop boleans matrix determined by \code{GA::initialization()}.
#'
#' @param fitnessFunction fitness function that takes in an lm or glm model and
#' returns a numerical fitness of that model. Users can choose AIC or BIC or
#' even define by themselves.
#'
#' @param model the linear model that user wants to use to fit in the data,
#' can be either \code{lm} or \code{glm}.
#'
#' @param fitScores fitness scores calculated by \code{GA::fitness()}. It must
#' be the case that a lower fitness score means a better model.
#'
#' @param offspringNum number of offspring generated to update the generation.
#'
#' @param method the selection mechanism that user wants to apply to select
#' parents, can be choosen from 1 to 3; 1 indicates selecting both parents with
#' probability proportional to ranking; 2 indicates selecting one parent with
#' probability proportional to ranking and one parent randomly, and 3 indicates
#' selecting with method of tournament selection.
#'
#' @param dat data frame containing the predictors in the model.
#' First column should be the response variable.
#'
#' @param K number of groups to partition the population into for tournament selection
#'
#' @details This function selects parents to breed based on the passed in
#' selection mechanism. Select (offspringNum / 2) pairs of parents from the
#' passed in population pop, by specific selection mechanism and the fitness
#' scores fitScores obtained from \code{GA::fitness()}.
#'
#' @return Returns a list containing the index for each of parent 1, the index
#' for each of parent 2, the fittest individual, and the fitness scores for the
#' whole population.


selection <- function(pop, fitnessFunction, model, fitScores, offspringNum, method, dat, K){

  # Size of the generation
  P <- ncol(pop)

  # Sample from the P chromosomes, with weights specified in fitnessProb,
  # with replacement, to obtain a parent population of size offspringNum /2.
  # Note there are duplicates within the parent population.

  # Method 1: Select both parents with probability proportional to ranking
  if (method == 1) {

    # Compute a vector of probability weights proportional to ranking
    # Since lowest score is the best, take the reverse rank
    fitnessProb <- 2 * rank(-fitScores) / (P * (P + 1))

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
    fitnessProb <- 2 * rank(-fitScores) / (P * (P + 1))

    # Index for parent 1
    indexParent1 <- sample(x = 1:P, size = offspringNum / 2, replace = T,
                           prob = fitnessProb)

    # Index for parent 2
    indexParent2 <- sample(x = 1:P, size = offspringNum / 2, replace = T)

  }

  else if (method == 3){

    #define function to find the index of the best fit in the assigned group
    myFit <- function(m){
      X <- as.matrix(dat[-1])
      y <- as.matrix(dat[1])
      ind <- order(fitness(pop[, GroupInd[, m]], y, X, fitnessFunction, model, dat))[1]
      return(ind)
    }

    #initialize some index needed in while loop
    m <- 1 : K

    #index of how many number of we already choose
    choose <- 0

    #Number of the population after the select iteration
    numRest <- P
    popInd <- 1 : P
    goodFitTol <- NULL

    while(choose < offspringNum){

      #sample the index of the population to implement the randomly choosing
      sampInd <- sample(popInd, numRest- numRest %% K)

      #Partition into K disjoint subsets, each column represents a group
      GroupInd <- matrix(sampInd, ncol = K)

      #Find the index of best individual in each group
      goodFit <- sapply(m, myFit) + floor(numRest / K) * (m - 1)

      #count the total number of good fit indivduals
      choose <- choose + K

      #remove the indivdual that was chosen
      popInd <- popInd[-goodFit]

      #combine the good fit individuals' index to a vector
      goodFitTol <- c(goodFitTol, goodFit)

      #number of rest poplation will be used in the next iteration
      numRest <- numRest - K
    }

    #if the number of choose larger than the number of population(P) we need, choose the first P
    if(choose > offspringNum){
      goodFitTol <- goodFitTol[1 : offspringNum]
    }
    #Randomly pair the parents
    randomPairInd <- sample(1 : offspringNum, offspringNum / 2)

    indexParent1 <- goodFitTol[randomPairInd]
    indexParent2 <- goodFitTol[-randomPairInd]


  }

  else {
    stop("Need to choose a selection mechanism.")
  }

  # Build the selection result as a list
  selResult <- list()
  selResult$indexParent1<-indexParent1
  selResult$indexParent2<-indexParent2
  # selResult$fittest <- pop[,which.max(fitnessProb)] # Keep a copy of the fittest individual.
  selResult$fitnessScores <- fitScores # keep for plotting

  # Return the result
  return(selResult)
}






