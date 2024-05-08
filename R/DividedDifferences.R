#' @title Calculate divided differences of a given dataset up to a specified order.
#' @description This function calculates the divided differences of a given dataset up to a specified order using Newton's divided differences method.
#'
#' @param A A matrix containing the dataset, where the first column represents the x-values and the second column represents the corresponding function values (f(x)).
#' @param order The order of divided differences to calculate.
#' @param showTable Logical indicating whether to display the divided differences table (default: FALSE).
#'
#' @return If `showTable` is `TRUE`, prints the divided differences table and returns the highest order divided difference (which corresponds to the Newton interpolating polynomial coefficient).
#'         If `showTable` is `FALSE`, returns the highest order divided difference directly.
#'
#' @export divDiff
#'
#' @examples
#' # Calculate and display the divided differences of a dataset up to order 6
#' x <- c(1, 2, 3, 4, 5, 6, 7)
#' fx <- c(144, 56, 35, 22, 78, 3, 17)
#' divDiff(cbind(x, fx), 6, showTable = TRUE)
divDiff = function(A, order, showTable = F){
  if (nrow(A) <= order)
    return("The order of differences cannot be greater than the number of rows in the matrix.")

  if (order < 1)
    return("Order should be greater than 0.")

  separatorLine = paste(rep("-", 10), collapse = "")
  divDiff = A

  for (i in 1:order){
    result = rep(0, nrow(A))

    for (j in 1:(nrow(A) - i)){
      currentDiff = divDiff[j + 1, ncol(divDiff)]
      previousDiff = divDiff[j, ncol(divDiff)]

      currentX = divDiff[j + i, 1]
      prevX = divDiff[j, 1]

      result[j] = (currentDiff - previousDiff) / (currentX - prevX)
    }

    divDiff = cbind(divDiff, result)

    colnames(divDiff)[i+2] = paste("DD", i)
  }

  if (showTable){
    print("Divided differences Table: ")
    print(divDiff)
    print(separatorLine)
  }

  return(unname(divDiff[1, ncol(divDiff)]))
}
