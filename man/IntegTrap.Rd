% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IntegTrap.R
\name{IntegTrap}
\alias{IntegTrap}
\title{Trapezoidal Rule Numerical Integration}
\usage{
IntegTrap(func, a, b, n)
}
\arguments{
\item{func}{The function to be integrated.}

\item{a}{The lower limit of integration.}

\item{b}{The upper limit of integration.}

\item{n}{The number of subintervals (must be greater than 1).}
}
\value{
The approximate value of the integral of the function \code{func} from \code{a} to \code{b} using the Trapezoidal rule.
}
\description{
This function applies the Trapezoidal rule to compute the integral of a function \code{func} from \code{a} to \code{b} with \code{n} subintervals.
The Trapezoidal rule is a numerical method for approximating definite integrals using trapezoids to approximate the area under the curve.
}
\examples{
# Define the function to be integrated (e.g., x^2 - 2)
func <- function(x) {
  y <- x^2 - 2
  return(y)
}

# Perform numerical integration using the Trapezoidal rule with n = 12 subintervals
result <- IntegTrap(func, 0, 2, 12)
print(result)
}
