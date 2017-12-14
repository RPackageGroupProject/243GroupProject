#######################################
####  Plotting Results
#######################################

#' Plotting the Progress of the Genetic Algorithm
#'
#' @description Make a chromosome with fixed probability 0.01 to mutate.
#'
#' @usage plotGA(results)
#'
#' @details This function plots the output list from \code{GA::select()}. The plot shows fitness
#' scores of each model in each generation as well as the best fitting model (in green).
#'
#' @export
#'
#' @examples
#' results <- select(mtcars)
#' plotGA(results)

plotGA <- function(results){

  P<-nrow(results$fitnessScores)
  numGens<-ncol(results$fitnessScores)
  X<-matrix(rep(seq(numGens),P,numGens),nrow=P,ncol=numGens)
  Y<-t(results$fitnessScores)
  plot(X,Y,xlab='generation',ylab='AIC',pch=19,cex=0.5)
  lines(seq(numGens),apply(results$fitnessScores,FUN=min,MARGIN=2),lty=1,col='green')

}

