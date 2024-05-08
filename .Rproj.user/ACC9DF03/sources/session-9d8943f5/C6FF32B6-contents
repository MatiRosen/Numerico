#' @title Lagrange Polynomial Interpolation
#' @description
#' This function computes the Lagrange polynomial interpolation based on the data matrix \code{A}, which contains data points (x, fx).
#' It evaluates the Lagrange polynomial at the specified values \code{xs}.
#'
#' @param A Matrix containing the data points (x, fx).
#' @param xs Values at which to evaluate the Lagrange polynomial.
#' @param showExpr Logical indicating whether to display the interpolated Lagrange polynomial expression (default: FALSE).
#' @param showMatrix Logical indicating whether to display the coefficient matrix (default: FALSE).
#' @param showGraph Logical indicating whether to plot the interpolated Lagrange polynomial (default: FALSE).
#'
#' @return A vector containing the interpolated function values at specified \code{xs}.
#'
#' @export interpLagrange
#'
#' @examples
#' # Example: Lagrange interpolation
#' x <- c(1, 2, 3, 4, 5, 6, 7)
#' fx <- c(144, 56, 35, 22, 78, 3, 17)
#' A <- cbind(x, fx)
#' interpLagrange(A, 5.5, showExpr = FALSE, showMatrix = FALSE, showGraph = FALSE)

interpLagrange=function(A,xs,showExpr=F,showMatrix=F,showGraph=F){
  x=A[,1]
  fx=A[,2]
  ln=length(x)
  terms=c()
  den=c()
  separatorLine = paste(rep("-", 10), collapse = "")

  for(i in 1:ln){
    den[i] = 1
    for(y in 1:ln){
      if(y != i){
        den[i]=den[i]*(x[i]-x[y])
      }
    }
    num = 1
    for(y in 1:ln){
      if(y != i){
        num=num*(xs-x[y])
      }
    }
    terms[i]=fx[i]*num/den[i]
  }

  fxs=(sum(terms))
#******************************************************
  if(showGraph==T | showExpr==T){
    expr=c()

    for (i in 1:ln) {
      expr[i]=paste("(",as.character(fx[i]))

      for(j in 1:ln){
        if(j != i){
          expr[i] = paste(expr[i],paste(")*(x-",as.character(x[j],")")))
        }
      }

      expr[i] = paste(expr[i],")","/",as.character(den[i]))
    }

    expr=paste(expr,collapse = " + ")

    poly=function(x) {
      eval(parse(text=paste(expr,collapse = " + ")))
    }
    #******************************************************
    if(showGraph==T){
      curve(poly, from = min(x), to = max(x), main = "Lagrange Polynomial",
        col = "black", lwd = 2,type = "l",xlab = "X",ylab = "P(x)")
      points(x,fx,pch = 21,cex=1.5,lwd=0.5,bg="red")
      points(xs, poly(xs), pch = 21, cex = 1.8, lwd = 0.9, bg = "yellow")
    }
    #******************************************************
    if(showExpr==T){
      print(" Expression : ")
      cat("A(x) = ", expr, "\n")
      print(separatorLine)
    }
  }
   #******************************************************
  if(showMatrix==T){
     matrix = matrix(c(seq(1,ln),terms), byrow=F,nrow = ln, ncol = 2)
     colnames(matrix) = list("L(i)", "Term")
     print("Coefficients Matrix: ")
     print(matrix)
     print(separatorLine)
   }
  #******************************************************
  names(fxs)="f(xs)"
  return(fxs)
}
