#' @title Perform the Bisection Method to find a root of a given function within a specified interval.
#' @description This function applies the Bisection Method to find a root of a given function within a specified interval.
#'
#' @param expr An expression representing the function for which a root is to be found. The function should be defined in terms of `x`.
#' @param lowerLim The lower limit of the interval within which the root is sought.
#' @param upperLim The upper limit of the interval within which the root is sought.
#' @param iterations The maximum number of iterations to perform.
#' @param digits Number of significant digits to be displayed in the result (default: 15).
#' @param stepByStep Logical indicating whether to display intermediate results at each iteration (default: FALSE).
#'
#' @return If `stepByStep` is `TRUE`, prints a table of intermediate results and returns the approximate root.
#'         If `stepByStep` is `FALSE`, returns the approximate root directly.
#' @export bisection
#'
#' @examples
#' bisection(expression(x^2-2), 0, 2, 6, 15, T)
bisection = function(expr, lowerLim, upperLim, iterations, digits = 15, stepByStep = F){
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

  counter = 1
  result = matrix(0, ncol=2)
  fx = 1

  while(counter <= iterations &&  fx != 0){
    x = (lowerLim + upperLim) / 2
    fx = f(x)
    result = rbind(result, c(x, fx))

    if(f(lowerLim) * f(x) < 0){
      upperLim = x
    } else{
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
