#' @title Algebraic Interpolation
#' @description
#' This function uses the Vandermonde matrix method to perform algebraic interpolation based on given data points (x, fx) in matrix \code{A}.
#' It computes the coefficients of the algebraic polynomial interpolant and evaluates the polynomial at specified values \code{xs}.
#'
#' @param A Matrix containing the data points (x, fx), where the first column represents x-values and the second column represents corresponding function values fx.
#' @param xs Values at which to evaluate the interpolated polynomial.
#' @param showExpr Logical indicating whether to display the interpolated polynomial expression (default: FALSE).
#' @param showMatrix Logical indicating whether to display the coefficients matrix (default: FALSE).
#' @param dec Number of decimal places to round coefficients when displaying the polynomial expression (default: 2).
#' @param showGraph Logical indicating whether to plot the interpolated polynomial (default: FALSE).
#'
#' @return A vector containing the interpolated function values at specified \code{xs}.
#'
#' @export interpAlgebraic
#'
#' @examples
#' # Define the data points (x, fx)
#' x <- c(1, 2, 3, 4, 5, 6, 7)
#' fx <- c(144, 56, 35, 22, 78, 3, 17)
#' A <- cbind(x, fx)
#'
#' # Perform algebraic interpolation and evaluate at x = 5.5
#' interpAlgebraic(A, 5.5, showExpr = TRUE, showMatrix = TRUE, showGraph = TRUE)
interpAlgebraic=function(A,xs,showExpr=F,showMatrix=F,dec=2,showGraph=F){
  separatorLine = paste(rep("-", 10), collapse = "")
  x=A[,1]
  fx=A[,2]
  ln=length(x)
  M=matrix(rep(0,ln*ln),ln,ln)
  for(i in 1:ln){
    for(j in 1:ln){
      M[i,j]=x[i]^(ln-j)
    }
  }
  invCoef=solve(M,as.matrix(fx))
  aux=xs
  fxs=0
  coef=rev(invCoef)
  for(i in 1:length(coef)){
    fxs=fxs+aux^(i-1)*(coef[i])
  }
  #*******************************************************
  if(showGraph==T){
    for (i in 1:length(coef)-1) {
      expr = as.character(coef[1])
      for (j in 1:i) {
        expr = paste(expr," + ",as.character(coef[j+1]), " * ", "(x^",as.character(j),")",sep="")
      }
    }
    poly=function(x) {eval(parse(text=expr))}

    curve(poly, from = min(x), to = max(x), main = "Algebraic Polynomial",
          col = "black", lwd = 2,type = "l",xlab = "X",ylab = "P(x)")
    points(x,fx,pch = 21,cex=1.5,lwd=0.5,bg="red")
    points(xs, poly(xs), pch = 21, cex = 1.8, lwd = 0.9, bg = "yellow")
  }
#*******************************************************
  if(showMatrix==T){
    matrix = matrix(c(seq(0,ln-1),coef), byrow=F,nrow = ln, ncol = 2)
    colnames(matrix) = list("Degree", "Coeff")
    print("Coefficients Matrix: ")
    print(matrix)

    print(separatorLine)
      }
#*******************************************************
  if(showExpr==T){
    print(" Expression : ")
    rCoef=round(coef,dec)

  for (i in 1:length(rCoef)-1) {
    expr = as.character(rCoef[1])
    for (j in 1:i) {
      expr = paste(expr," + ",as.character(rCoef[j+1]), " * ", "(x^",as.character(j),")",sep="")
    }
  }
  cat("A(x) = ", expr, "\n")

  print(separatorLine)
  }

  #*******************************************************
  names(fxs)="f(xs)"
  return(fxs)
}
