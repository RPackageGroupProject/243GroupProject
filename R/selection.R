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






