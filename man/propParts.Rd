% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ProportionalParts.R
\name{propParts}
\alias{propParts}
\title{Secant Root Finding}
\usage{
propParts(expr, lowerLim, upperLim, iterations, digits = 15, stepByStep = F)
}
\arguments{
\item{expr}{Expression representing the function for which the root is to be found.}

\item{lowerLim}{Lower limit of the interval for root searching.}

\item{upperLim}{Upper limit of the interval for root searching.}

\item{iterations}{Maximum number of iterations for the Secant method.}

\item{digits}{Number of significant digits to use in computations (default: 15).}

\item{stepByStep}{Logical indicating whether to display the iteration steps (default: FALSE).}
}
\value{
The estimated root of the function within the specified interval.
}
\description{
This function uses the Secant method to find a root of the function represented by the expression \code{expr} within the interval [\code{lowerLim}, \code{upperLim}].
It iteratively refines the estimate of the root based on the secant line until convergence or the specified number of iterations.
}
\examples{
# Example: Secant method for root finding
propParts(expression(x^2 - 2), 0, 2, 6, digits = 15, stepByStep = TRUE)
}
