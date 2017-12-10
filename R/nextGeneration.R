#######################################
####  Choose Next Generation
#######################################

#' Breed the Selected Parents Generated From \code{GA::selection()}
#'
#' @usage nextGeneration(pop, selResult, G)
#'
#' @param pop boleans matrix determined by \code{GA::initialization()}
#'
#' @param selResult list returnd by \code{GA::selection()}
#' At the minimum, it needs to contain indices for the selected parents 1 and 2.
#'
#' @param offspringNum number of offspring generated to update the generation.
#'
#' @details Breeding uses \code{GA::crossover()} for each pair of parents. The
#' Generation Gap G is the proportion of the generation to be replaced by
#' generated offspring. If the number of the selected population produced by
#' generation gap is not an even number, choose the nearest bigger even number.
#' Uses the offspringNum number of offspring to replace the least fit
#' individuals in parent generation
#'
#' @return Returns a C by P maatrix containing the population for the next
#' generation.
#'
#' @examples
#' dat <- mtcars # use the built-in mtcars data
#' C <- dim(dat)[2] - 1 #Number of variables
#' pop <- initialization(C) #produce boleans matrix
#' selResult <- selection(pop, f = AIC, dat)
#' pop <- nextGeneration(pop, selResult)
nextGeneration <- function(pop, selResult, offspringNum){

  # Size of the generation
  P <- ncol(pop)

  # Extract the parents
  allParent1 <- pop[,selResult$indexParent1]
  allParent2 <- pop[,selResult$indexParent2]

  # Extract the index of the individuals in the parenting generation that
  # needed to be replaced by offspring
  lessFitInd <- order(selResult$fitnessScores, decreasing = T)[1 : offspringNum]

  # Breed to create new generation
  for (i in offspringNum / 2){

    # Each time generate two offspring, stored in columns
    childrenChromes <- crossover(allParent1[, i], allParent2[, i])

    # Select first child and put into pop
    pop[, lessFitInd[i]] <- childrenChromes[, 1]
    pop[, lessFitInd[i + offspringNum / 2]] <- childrenChromes[, 2]

  }

  # Output the population for the next generation
  return(pop)
}


