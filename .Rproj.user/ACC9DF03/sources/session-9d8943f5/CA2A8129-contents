#' @title Perform numerical integration using Simpson's 3/8 rule.
#'
#' @description This function performs numerical integration using Simpson's 3/8 rule to approximate the integral of a given function over a specified interval.
#'
#' @param fn The function to be integrated.
#' @param a The lower limit of integration.
#' @param b The upper limit of integration.
#' @param n The number of subintervals (must be a multiple of 3).
#'
#' @return The approximate value of the integral of the function \code{fn} from \code{a} to \code{b} using Simpson's 3/8 rule.
#'
#' @export IntegSimp3.8
#'
#' @examples
#' # Define the function to be integrated (e.g., x^2 - 2)
#' fn <- function(x) {
#'   y <- x^2 - 2
#'   return(y)
#' }
#'
#' # Perform numerical integration using Simpson's 3/8 rule with n = 12 subintervals
#' result <- IntegSimp3.8(fn, 0, 2, 12)
#' print(result)
#'

IntegSimp3.8<-function(fn,a,b,n){
  if ((n/3)!=round((n/3))) {
    return("n must be multiple of 3")
  } else{
      h = (b-a)/n
  }
  if(n>3) {m1=seq(1,n-1,3)} else {m1=1}
  if(n>3) {m2=seq(2,n-1,3)} else {m2=2}
  if(n>3) {m3=seq(3,n-3,3)} else {m3=0}
  M1=sum(fn(a+(h*m1)))
  M2=sum(fn(a+(h*m2)))
  if (n>3){M3=sum(fn(a+(h*m3)))} else {M3=0}
  integral= ((3/8)*h)*(fn(a)+fn(b)+(3*M1)+(2*M3)+(3*M2))
  return(integral)
}



