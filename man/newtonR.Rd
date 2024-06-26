% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NewtonRaphson.R
\name{newtonR}
\alias{newtonR}
\title{Newton-Raphson Root Finding}
\usage{
newtonR(expr, lowerLim, upperLim, iterations, digits = 15, stepByStep = F)
}
\arguments{
\item{expr}{Expression representing the function for which the root is to be found.}

\item{lowerLim}{Lower limit of the interval for root searching.}

\item{upperLim}{Upper limit of the interval for root searching.}

\item{iterations}{Maximum number of iterations for the Newton-Raphson method.}

\item{digits}{Number of significant digits to use in computations (default: 15).}

\item{stepByStep}{Logical indicating whether to display the iteration steps (default: FALSE).}
}
\value{
The estimated root of the function within the specified interval.
}
\description{
This function uses the Newton-Raphson method to find a root of the function represented by the expression \code{expr} within the interval [\code{lowerLim}, \code{upperLim}].
It iteratively refines the estimate of the root based on the function's derivative until convergence or the specified number of iterations.
}
\examples{
# Example: Newton-Raphson method for root finding
newtonR(expression(x^2 - 2), 0, 2, 6, digits = 15, stepByStep = TRUE)
}
