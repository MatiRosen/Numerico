% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SplinesNaK.R
\name{SplinesNaK}
\alias{SplinesNaK}
\title{Cubic Spline with Known Second Derivatives}
\usage{
SplinesNaK(
  A,
  xs,
  showMatrix = FALSE,
  showPoly = FALSE,
  showGraph = FALSE,
  enlarge = FALSE,
  limits.y = c(-100, 100)
)
}
\arguments{
\item{A}{Matrix containing the data points. The first column corresponds to the x-values and the second column to the f(x) values.}

\item{xs}{The point at which the spline will be evaluated.}

\item{showMatrix}{Logical indicating whether to display the matrix used in spline computation (default: FALSE).}

\item{showPoly}{Logical indicating whether to display the polynomial coefficients (default: FALSE).}

\item{showGraph}{Logical indicating whether to display a graphical representation of the spline (default: FALSE).}

\item{enlarge}{Logical indicating whether to enlarge the graph window (default: FALSE).}

\item{limits.y}{Limits for the y-axis in the graphical representation (default: c(-100, 100)).}
}
\value{
The value of the cubic spline evaluated at \code{xs}.
}
\description{
This function computes the value of a cubic spline at a specified point \code{xs} using the given data points \code{A},
with known second derivatives at the endpoints (natural spline). It provides options to display the matrix used in the spline computation,
the polynomial coefficients, and/or a graphical representation of the spline.
}
\examples{
# Example: Compute the cubic spline with known second derivatives at a point
x <- c(1, 2, 4, 5, 6, 7, 9)
fx <- c(144, 56, 35, 22, 78, 3, 17)
A <- cbind(x, fx)
SplinesNaK(A, 3.5)

}
