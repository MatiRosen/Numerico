#' @title Simpson's Rule Numerical Integration
#' @description
#' This function applies Simpson's rule to compute the integral of a function \code{func} from \code{a} to \code{b} with \code{n} subintervals.
#' Simpson's rule is a numerical method for approximating definite integrals using quadratic interpolating polynomials.
#'
#' @param func The function to be integrated.
#' @param a The lower limit of integration.
#' @param b The upper limit of integration.
#' @param n The number of subintervals (must be even).
#'
#' @return The approximate value of the integral of the function \code{func} from \code{a} to \code{b} using Simpson's rule.
#'
#' @export IntegSimpson
#'
#' @examples
#' # Define the function to be integrated (e.g., x^2 - 2)
#' func <- function(x) {
#'   y <- x^2 - 2
#'   return(y)
#' }
#'
#' # Perform numerical integration using Simpson's rule with n = 12 subintervals
#' result <- IntegSimpson(func, 0, 2, 12)
#' print(result)
IntegSimpson<-function(func,a,b,n) {
  if ((n/2)!=round((n/2))) {
    return("n debe ser par")
  } else {
    h = (b-a)/n
  }
  if(n>2) {odds=seq(1,n-1,by=2)} else {odds=1}
  if(n>2) {evens=seq(2,n-2,by =2)} else {evens=0}
  s.odds=sum(func((a+h*odds)))
  if(n>2) {s.evens=sum(func((a+h*evens)))} else {s.evens=0}
  integral=(h/3)*(func(a)+func(b)+(4*s.odds)+(2*s.evens))
  return(integral)
}

