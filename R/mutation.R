#######################################
####  Mutation
#######################################

#' Chromosome Mutation
#'
#' @description Make a chromosome with fixed probability 0.01 to mutate.
#'
#' @usage mutation(chr)
#'
#' @param chr a boolean vector representing an individual chromosome.
#'
#' @details This function makes a chromosome have a 1% fixed chance to mutate
#' in each locale. If mutation happens at one locale, it will make the value
#' in that locale from T to F or from F to T.
#'
#' @return Return a mutated chromosome vector with the same length as input.


mutation <- function(chr){

  # For each element of chr determine whether to mutate with 1% prob
  mutate <- sample(c(T,F), length(chr), prob = c(0.01, 0.99), replace = T)

  # 'exclusive or' will toggle F to T and T to F
  xor(chr, mutate)

}



