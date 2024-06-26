% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TaylorPolynomial.R
\name{TaylorP}
\alias{TaylorP}
\title{Taylor Series Approximation}
\usage{
TaylorP(expr, x0 = 0, iterations = 2, xs, stepByStep = F)
}
\arguments{
\item{expr}{Mathematical expression represented as an R expression.}

\item{x0}{Point around which the Taylor series is expanded (default: 0).}

\item{iterations}{Number of iterations for the Taylor series expansion (default: 2).}

\item{xs}{Value at which to evaluate the Taylor series approximation.}

\item{stepByStep}{Logical indicating whether to display the computation step-by-step (default: FALSE).}
}
\value{
The sum of the Taylor series terms evaluated at the point \code{xs}.
}
\description{
This function computes the Taylor series approximation of a given mathematical expression around a specified point.
The Taylor series is expanded up to a specified number of iterations.
The function returns either the sum of the Taylor series terms or displays a step-by-step matrix with detailed information.
}
\examples{
# Example: Compute Taylor series approximation
fn <- expression(sin(x) * cos(x))
result <- TaylorP(expr = fn, x0 = pi/3, iterations = 5, xs = 0.3 * pi, stepByStep = TRUE)
}
