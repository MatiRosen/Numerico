#' Perform the Regula Falsi method for root finding.
#'
#' This function implements the Regula Falsi (False Position) method to find a root of the function represented by the expression \code{expr} within the interval [\code{lowerLim}, \code{upperLim}].
#' It iteratively refines the estimate of the root based on the linear interpolation between the function values at the interval endpoints.
#'
#' @title Regula Falsi Root Finding
#' @description
#' This function uses the Regula Falsi method to find a root of the function represented by the expression \code{expr} within the interval [\code{lowerLim}, \code{upperLim}].
#' It iteratively refines the estimate of the root based on linear interpolation until convergence or the specified number of iterations.
#'
#' @param expr Expression representing the function for which the root is to be found.
#' @param lowerLim Lower limit of the interval for root searching.
#' @param upperLim Upper limit of the interval for root searching.
#' @param iterations Maximum number of iterations for the Regula Falsi method.
#' @param digits Number of significant digits to use in computations (default: 15).
#' @param stepByStep Logical indicating whether to display the iteration steps (default: FALSE).
#'
#' @return The estimated root of the function within the specified interval.
#'
#' @export regulaFalsi
#'
#' @examples
#' # Example: Regula Falsi method for root finding
#' regulaFalsi(expression(x^2 - 2), 0, 2, 6, digits = 15, stepByStep = TRUE)
regulaFalsi = function(expr, lowerLim, upperLim, iterations, digits = 15, stepByStep = F){
  options(digits=digits)

  if (lowerLim == upperLim)
    return("Limits should not be the same.")

  if (iterations < 1)
    return("Iterations should be greater than 0.")

  f = function(x){
    eval(expr)
  }

  if (f(lowerLim) * f(upperLim) >= 0)
    return("Bolzano's theorem is not verified.")

  separatorLine = paste(rep("-", 10), collapse = "")

  if(lowerLim >= upperLim){
    print("The interval limits are reversed! Fixing it...")
    aux = lowerLim
    lowerLim = upperLim
    upperLim = aux
  }

  counter = 0
  fx = 1
  result = matrix(0, ncol=2)

  while (counter < iterations && fx != 0){
    fUpperLim = f(upperLim)
    fLowerLim = f(lowerLim)
    x = upperLim - fUpperLim * (upperLim - lowerLim) / (fUpperLim - fLowerLim)
    fx = f(x)
    result = rbind(result, c(x, fx))

    if (f(x) * f(lowerLim) < 0){
      upperLim = x
    } else {
      lowerLim = x
    }

    counter = counter + 1
  }

  if (stepByStep){
    colnames(result) = c("xi", "f(xi)")
    rownames(result) = c(0 : (nrow(result) - 1))
    print(result[-1, ])
    print(separatorLine)
  }

  return(result[nrow(result), 1])
}
