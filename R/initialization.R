#######################################
####  Initialization
#######################################

#' Generate the Initial Generation Randomly
#'
#' @description Bernoulli sampling to form the initial generation of the
#' Genetic Algorithm.
#'
#' @usage initialization(C, P)
#'
#' @param C chromosome length (number of predictor variables)
#'
#' @details This function produces initial generation given the chromosome
#' length for Genetic Algorithm and the population size. The row of the
#' generated boolean matrix represents the locus of a gene(variable), and the
#' column of the matrix represents different members of the first generation.
#'
#' @return A Bolean Matrix with dimension C by P where each column
#' representing a chromosome, in which T marks that the gene(variable)
#' as active and F as inactive.


initialization <- function(C, P){

  # Bernoulli sampling T or F to obtain C of them
  # T in a locus of a gene (variable) means the perticular variable is included
  replicate(P, sample(c(F,T), size = C, replace = T))

}



