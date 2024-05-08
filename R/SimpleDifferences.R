#' @title Simple Finite Differences
#' @description
#' This function calculates the simple finite differences of a specified order for the function values represented in matrix \code{A}.
#' The order of differences determines how many times the differences are taken consecutively.
#' The result is returned as the value of the highest order difference.
#'
#' @param A Matrix containing the function values. The first column should represent the x-values and the second column the corresponding function values.
#' @param order Order of the finite differences to be computed.
#' @param showTable Logical indicating whether to display the differences table (default: FALSE).
#'
#' @return The computed finite difference value of the specified order.
#'
#' @export simpleDiff
#'
#' @examples
#' # Example: Compute simple finite differences
#' x <- c(1, 2, 3, 4, 5, 6, 7)
#' fx <- c(144, 56, 35, 22, 78, 3, 17)
#' A <- cbind(x, fx)
#' simpleDiff(A, 6, showTable = TRUE)

simpleDiff = function(A, order, showTable = F){
  if (nrow(A) <= order)
    return("The order of differences cannot be greater than the number of rows in the matrix.")

  if (order < 1)
    return("Order should be greater than 0.")

  separatorLine = paste(rep("-", 10), collapse = "")
  simpleDiff = A

  for (i in 1:order){
    result = rep(0, nrow(A))

    for (j in 1:(nrow(A) - i)){
      currentDiff = simpleDiff[j + 1, ncol(simpleDiff)]
      previousDiff = simpleDiff[j, ncol(simpleDiff)]

      result[j] = currentDiff - previousDiff
    }

    simpleDiff = cbind(simpleDiff, result)

    colnames(simpleDiff)[i+2] = paste("D", i)
  }

  if (showTable){
    print("Differences Table: ")
    print(simpleDiff)
    print(separatorLine)
  }
  return(unname(simpleDiff[1, ncol(simpleDiff)]))
}
