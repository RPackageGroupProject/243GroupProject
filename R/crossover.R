#######################################
####  Crossover
#######################################

#' Genetic Operator: Chromosome Crossover and Mutation
#'
#' @description Make 2 individual parent chromosomes crossover and
#' mutate when breeding offsprings.
#'
#' @usage crossover(chr1, chr2)
#'
#' @param chr1,chr2 a numeric vectors represents individual parent chromosome.
#'
#' @details This function makes two individual parent chromosomes perfrom
#' crossover and mutation when breeding offspring.
#' Note that the crossover is simply one-point crossover and the mutation
#' is based on \code{GA::mutation()}.
#'
#' @return A C by 2 matrix with each column representing the offspring from a genetic
#' operator of crossover and mutation, C is the chromosome length.
#'
#' @examples
#' ind1 <- initialization(10, 15)[, 1]
#' ind2 <- initialization(10, 15)[, 2]
#' crossover(ind1, ind2)

crossover <- function(chr1, chr2){

  # length of each chromosome
  C1 <- length(chr1)
  C2 <- length(chr2)

  # Make sure we're working with equal-sized vectors
  stopifnot(C1 == C2)

  # Randomly select a point as the point of crossover.
  k <- sample(1 : (C1 - 1), 1)

  # Return the crossed chromosomes together
  cbind(mutation(c(chr1[1 : k], chr2[(k + 1) : C1])),
        mutation(c(chr2[1 : k], chr1[(k + 1) : C1]))
  )
}


