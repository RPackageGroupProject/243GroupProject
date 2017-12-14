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
#' dat <- mtcars
#' C <- ncol(dat) - 1
#' P <- as.integer(1.5 * C)
#' pop <- initialization(C, P)
#' fitnessFunction <- AIC
#' model <- lm
#' y <- as.matrix(dat[1])
#' X <- as.matrix(dat[-1])
#' fitScores <- fitness(pop, y, X, fitnessFunction, model, dat)
#' offspringNum <- 4
#' selResult <- selection(pop, fitnessFunction, model, fitScores, offspringNum, method = 1, dat, K = 2)
#' nextGeneration(pop, selResult, offspringNum)

nextGeneration <- function(pop, selResult, offspringNum){

  # Size of the generation
  P <- ncol(pop)

  # Extract the parents
  allParent1 <- as.matrix(pop[, selResult$indexParent1],
                          ncol = length(selResult$indexParent1))
  allParent2 <- as.matrix(pop[, selResult$indexParent2],
                          ncol = length(selResult$indexParent2))

  # Extract the index of the individuals in the parenting generation that
  # needed to be replaced by offspring
  lessFitInd <- order(selResult$fitnessScores, decreasing = T)[1 : offspringNum]

  # Breed to create new generation
  for (i in offspringNum / 2){

    # Each time generate two offspring, stored in columns
    childrenChromes <- crossover(allParent1[, i], allParent2[, i])

    # Put both children into the population
    pop[, lessFitInd[i]] <- childrenChromes[, 1]
    pop[, lessFitInd[i + offspringNum / 2]] <- childrenChromes[, 2]

  }

  # Output the population for the next generation
  return(pop)
}


