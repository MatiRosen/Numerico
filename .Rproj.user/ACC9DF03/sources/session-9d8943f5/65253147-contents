#' @title Newton-Raphson Root Finding
#' @description
#' This function uses the Newton-Raphson method to find a root of the function represented by the expression \code{expr} within the interval [\code{lowerLim}, \code{upperLim}].
#' It iteratively refines the estimate of the root based on the function's derivative until convergence or the specified number of iterations.
#'
#' @param expr Expression representing the function for which the root is to be found.
#' @param lowerLim Lower limit of the interval for root searching.
#' @param upperLim Upper limit of the interval for root searching.
#' @param iterations Maximum number of iterations for the Newton-Raphson method.
#' @param digits Number of significant digits to use in computations (default: 15).
#' @param stepByStep Logical indicating whether to display the iteration steps (default: FALSE).
#'
#' @return The estimated root of the function within the specified interval.
#'
#' @export newtonR
#'
#' @examples
#' # Example: Newton-Raphson method for root finding
#' newtonR(expression(x^2 - 2), 0, 2, 6, digits = 15, stepByStep = TRUE)

newtonR = function(expr, lowerLim, upperLim, iterations, digits = 15, stepByStep = F){
  options(digits=digits)
  f = function(x, expr){
    eval(expr)
  }

  if (lowerLim == upperLim)
    return("Limits should not be the same.")

  if (iterations < 1)
    return("Iterations should be greater than 0.")


  if (f(lowerLim, expr) * f(upperLim, expr) >= 0)
    return("Bolzano's theorem is not verified.")

  separatorLine = paste(rep("-", 10), collapse = "")

  if(lowerLim >= upperLim){
    print("The interval limits are reversed! Fixing it...")
    aux = lowerLim
    lowerLim = upperLim
    upperLim = aux
  }

  derivate = D(expr, "x")
  secondD = D(derivate, "x")
  x = lowerLim

  if (f(upperLim, expr) * f(upperLim, secondD) > 0){
    x = upperLim
  }

  counter = 1
  fx = f(x, expr)
  result = matrix(c(x, fx), ncol=2)

  while (counter <= iterations && fx != 0){
    xDerivate = D(expr, "x")

    x = x - (fx / f(x, xDerivate))
    fx = f(x, expr)
    result = rbind(result, c(x, fx))

    counter = counter + 1
  }

  if (stepByStep){
    colnames(result) = c("xi", "f(xi)")
    rownames(result) = c(0 : (nrow(result) - 1))
    print(result)
    print(separatorLine)
  }


  return(result[nrow(result), 1])
}

