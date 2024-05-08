#' @title Compute Bernoulli numbers up to a given index.
#' @description Compute Bernoulli numbers up to a given index using the recursive algorithm.
#'
#' @param n Integer specifying the index up to which Bernoulli numbers should be computed.
#' @param lista_entera Logical indicating whether to return the sequence of Bernoulli numbers as a character vector (default: TRUE). If FALSE, only the last Bernoulli number in the sequence is returned as a character.
#'
#' @return If \code{lista_entera} is \code{TRUE}, returns a character vector containing the computed Bernoulli numbers up to index \code{n}. If \code{lista_entera} is \code{FALSE}, returns the last computed Bernoulli number as a character.
#' @export BernoulliNos
#' @examples
#' BernoulliNos(25, F)
BernoulliNos <- function(n, lista_entera = TRUE) {

  my_packages=c("MASS", "gmp")
  not_installed = my_packages[!(my_packages %in% installed.packages()[, "Package"])]
  if(length(not_installed)) install.packages(not_installed)
  library(gmp)
  library(MASS)

    if (n == 0) {
    return("1")
  } else if (n == 1) {
    return("-1/2")
  }
  n <- n + 1
  bernoullis <- as.bigq(c())
  fila <- as.bigq(c())
  for (i in 1:n) {
    fila <- c(fila, as.bigq(1, i))
  }

  bernoullis <- c(bernoullis, fila[1])
  fila <- fila[-1]
  bernoullis <- c(bernoullis, -fila[1])

  for (i in 1:(n - 2)) {
    fila_anterior <- fila
    fila <- as.bigq(c())
    for (j in 1:(length(fila_anterior) - 1)) {
      fila <- c(fila, (fila_anterior[j] - fila_anterior[j + 1]) * j)
    }
    bernoullis <- c(bernoullis, fila[1])
  }

  if (lista_entera) {
    return(as.character(bernoullis))
  } else {
    return(as.character(bernoullis[length(bernoullis)]))
  }
}

