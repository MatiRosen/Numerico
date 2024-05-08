#' @title Trapezoidal Rule Numerical Integration
#' @description
#' This function applies the Trapezoidal rule to compute the integral of a function \code{func} from \code{a} to \code{b} with \code{n} subintervals.
#' The Trapezoidal rule is a numerical method for approximating definite integrals using trapezoids to approximate the area under the curve.
#'
#' @param func The function to be integrated.
#' @param a The lower limit of integration.
#' @param b The upper limit of integration.
#' @param n The number of subintervals (must be greater than 1).
#'
#' @return The approximate value of the integral of the function \code{func} from \code{a} to \code{b} using the Trapezoidal rule.
#'
#' @export print
#'
#' @examples
#' # Define the function to be integrated (e.g., x^2 - 2)
#' func <- function(x) {
#'   y <- x^2 - 2
#'   return(y)
#' }
#'
#' # Perform numerical integration using the Trapezoidal rule with n = 12 subintervals
#' result <- IntegTrap(func, 0, 2, 12)
#' print(result)
IntegTrap<-function(func,a,b,n){
  h =(b-a)/n
  if(n>1){c=c(1:(n-1))} else{
    c=0
  }
  integral=(h/2)*(func(a)+func(b)+(2*sum(func(a+(c*h)))))
  return(integral)
}

