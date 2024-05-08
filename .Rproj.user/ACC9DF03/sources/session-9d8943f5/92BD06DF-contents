#' @title Cubic Spline with Specified Derivatives at Endpoints
#'
#' @description
#' This function computes the value of a cubic spline at a specified point \code{xs} using the given data points \code{A},
#' where the derivatives at the endpoints are specified. It provides options to display the matrix used in the spline computation,
#' the polynomial coefficients, and/or a graphical representation of the spline.
#'
#' @param A Matrix containing the data points. The first column corresponds to the x-values and the second column to the f(x) values.
#' @param deriv Vector of length 2 specifying the derivatives at the endpoints. The first element corresponds to the derivative at the first endpoint,
#' and the second element corresponds to the derivative at the last endpoint.
#' @param xs The point at which the spline will be evaluated.
#' @param showMatrix Logical indicating whether to display the matrix used in spline computation (default: FALSE).
#' @param showPoly Logical indicating whether to display the polynomial coefficients (default: FALSE).
#' @param showGraph Logical indicating whether to display a graphical representation of the spline (default: FALSE).
#' @param enlarge Logical indicating whether to enlarge the graph window (default: FALSE).
#' @param limits.y Limits for the y-axis in the graphical representation (default: c(-100, 100)).
#'
#' @return The value of the cubic spline evaluated at \code{xs}.
#'
#' @examples
#' # Example: Compute the cubic spline with specified derivatives at the endpoints at a point
#' x <- c(1, 2, 3, 4, 5, 6, 7)
#' fx <- c(144, 56, 35, 22, 78, 3, 17)
#' A <- cbind(x, fx)
#' deriv <- c(1, 0) # Derivative at first endpoint = 1, derivative at last endpoint = 0
#' SplinesClamp(A, deriv, 3.5)
#'
#' @export SplinesClamp
SplinesClamp=function(A,deriv,xs,showMatrix=FALSE,showPoly=FALSE,showGraph=FALSE,enlarge=FALSE,limits.y=c(-100,100)){
  separatorLine = paste(rep("-", 10), collapse = "")
  x=A[,1]
  fx=A[,2]
  nData=length(x)
  nPol=nData-1
  a=fx
  h=matrix(0,nPol,1)
  for (i in 1:nPol) {
    h[i]=(x[i+1]-x[i])}
  A=matrix(0,nData,nData)
  A[1,1:2]=c(2*h[1],h[1])
  A[nData,nPol:nData]=c(h[nPol],2*h[nPol])
  F=(matrix(0,nData,1))
  F[1]=3*((a[2]-a[1])/h[1])-3*deriv[1]
  F[nData]=deriv[2]-3*((a[nData]-a[nPol])/h[nPol])
  for (i in 2:nPol) {
    A[i,i-1]=h[i-1]
    A[i,i]=2*(h[i-1]+h[i])
    A[i,i+1]=h[i]
    F[i]=3*((a[i+1]-a[i])/h[i])-3*((a[i]-a[i-1])/h[i-1])
    F[nData]=3*deriv[2]-3*(a[nData]-a[nPol])/h[nPol]
  }
#-----------------------------------------------------------
  c=solve(A)%*%F
  d=c()
  b=c()
  for (i in 1:nPol) {
    d[i]=(c[i+1]-c[i])/(3*h[i])
    b[i]=(a[i+1]-a[i])/h[i]-h[i]*(2*c[i]+c[i+1])/3
  }
#-----------------------------------------------------------
  for(i in 1:nPol){
    if(xs>=x[i] & xs<=x[i+1]){
      fxs=a[i]+b[i]*(xs-x[i])+c[i]*(xs-x[i])^2+d[i]*(xs-x[i])^3
    }
  }
#*******************************************************
  if(showMatrix==TRUE){
    print(separatorLine)
    cat("\n")
    sep=rep("|",length(F))
    Frame=as.data.frame(cbind(sep,A,sep,F,sep))
    colnames(Frame)=c()
    rownames(Frame)=c()
    print("Matrix")
    print(Frame)
  }
#*****************************************************
  if(showPoly==TRUE){
    print(separatorLine)
    cat("\n")
    print("Polynomials  Si(x)")
    cat("\n")
    orden=seq(0,(nPol-1))
    M=cbind(orden,a[1:(length(a)-1)],b,c[1:(length(c)-1)],d)
    colnames(M)=c("i", "ai","bi","ci", "di")
    print(M)
  }
#*******************************************************
 if(showGraph==TRUE){
    if(enlarge==TRUE){
      plot(x,fx,pch = 21,cex=1.5,lwd=0.5,bg="red",
           main = "Spline Natural",xlab = "X",ylab = "SplNat(x)",
           ylim=limits.y)

    }else{
      plot(x,fx,pch = 21,cex=1.5,lwd=0.5,bg="red",
               main = "Spline Natural",xlab = "X",ylab = "SplNat(x)")
    }
    for(p in 1:length(b)){
      expr=as.character(a[p])
      expr=paste(expr,"+",as.character(b[p])," * (x-",as.character(x[p]),")")
      expr=paste(expr,"+",as.character(c[p]),"*(x-",as.character(x[p]),")^2")
      expr=paste(expr,"+",as.character(d[p]),"*(x-",as.character(x[p]),")^3")
      poly=function(x) {eval(parse(text=expr))}
      curve(poly, from = x[p], to = x[p+1], col = "black", lwd = 2,type = "l",add=T)
      points(xs, fxs, pch = 21, cex = 1.8, lwd = 0.9, bg = "yellow")
    }
  }
#*******************************************************
  names(fxs)="f(xs)"
  return(unname(fxs))
}

