#' @title Taylor Series Approximation
#' @description
#' This function computes the Taylor series approximation of a given mathematical expression around a specified point.
#' The Taylor series is expanded up to a specified number of iterations.
#' The function returns either the sum of the Taylor series terms or displays a step-by-step matrix with detailed information.
#'
#' @param expr Mathematical expression represented as an R expression.
#' @param x0 Point around which the Taylor series is expanded (default: 0).
#' @param iterations Number of iterations for the Taylor series expansion (default: 2).
#' @param xs Value at which to evaluate the Taylor series approximation.
#' @param stepByStep Logical indicating whether to display the computation step-by-step (default: FALSE).
#'
#' @return The sum of the Taylor series terms evaluated at the point \code{xs}.
#'
#' @export TaylorP
#'
#' @examples
#' # Example: Compute Taylor series approximation
#' fn <- expression(sin(x) * cos(x))
#' result <- TaylorP(expr = fn, x0 = pi/3, iterations = 5, xs = 0.3 * pi, stepByStep = TRUE)

TaylorP = function(expr, x0=0, iterations=2, xs, stepByStep = F){
  if (iterations <= 0){
    print("Iterations must be greater than 0")
    return(xs)
  }

  x = x0
  terms = eval(expr, x0)
  matrix = matrix(c(0, as.character(expr), terms, terms, sum(terms)), ncol = 5)
  colnames(matrix) = list("Degree", "Dfunction", "F(x0)", "Terms", "P(xs)")

  for (i in 1:iterations){
    expr = D(expr, "x")
    eval_expr = eval(expr)
    terms = c(terms, (eval_expr*(xs-x0)^i)/(factorial(i)))
    matrix = rbind(matrix, c(i, expr, eval_expr, terms[length(terms)], sum(terms)))
  }

  if (stepByStep){
    View(matrix)
    return(sum(terms))
  } else{
    return(sum(terms))
  }
}

