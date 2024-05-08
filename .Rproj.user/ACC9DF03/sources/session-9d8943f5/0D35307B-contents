#' @title Hermite Polynomial Interpolation
#' @description
#' This function computes the Hermite polynomial interpolation based on the data matrix \code{A}, which contains data points (x, fx) along with their derivatives (dx, dxx).
#' It evaluates the Hermite polynomial at the specified values \code{xs}.
#'
#' @param A Matrix containing the data points (x, fx) and their derivatives (dx, dxx).
#' @param xs Values at which to evaluate the Hermite polynomial.
#' @param showTable Logical indicating whether to display the divided differences table (default: FALSE).
#' @param showExpr Logical indicating whether to display the interpolated Hermite polynomial expression (default: FALSE).
#' @param showMatrix Logical indicating whether to display the coefficient matrix (default: FALSE).
#' @param showGraph Logical indicating whether to plot the interpolated Hermite polynomial (default: FALSE).
#'
#' @return A vector containing the interpolated function values at specified \code{xs}.
#'
#' @export interpHermite
#'
#' @examples
#' # Example 1: Hermite interpolation with given derivatives
#' x <- c(3, 2, 5, 6, 9)
#' fx <- c(27, 6, 135, 234, 783)
#' dx <- c(33, NA, 85, NA, NA)
#' dxx <- c(20, NA, NA, 15, NA)
#' A <- cbind(x, fx, dx, dxx)
#' interpHermite(A, 4.5, showTable = TRUE, showExpr = TRUE, showMatrix = TRUE, showGraph = TRUE)

#' # Example 2: Hermite interpolation with missing derivatives
#' x <- c(-1, 1, 2)
#' fx <- c(3, 3, 6)
#' dx <- c(1, NA, -1)
#' dxx <- c(NA, NA, 10)
#' A <- cbind(x, fx, dx, dxx)
#' interpHermite(A, 1.5, showTable = TRUE, showExpr = TRUE, showMatrix = TRUE, showGraph = TRUE)

interpHermite = function(A, xs, showTable = F, showExpr = F, showMatrix = F, showGraph = F){
  A = A[order(A[, 1]), ]

  expandedTable = A
  separatorLine = paste(rep("-", 10), collapse = "")

  pos = 0
  for (i in 1:nrow(A)){
    for (j in ncol(A):3){
      if (!is.na(A[i,j])){
        if ((sum(is.na(A[i,3:j])) == 0)){
          newRow = A[i,]
          newRow[(j):length(newRow)] = NA
          isLastRow = (i + 1 + pos) > nrow(expandedTable)
          if (isLastRow){
            expandedTable = rbind(expandedTable, newRow)

          } else{
            expandedTable = rbind(expandedTable[1:(i + pos), ], newRow,
                                  expandedTable[(i + 1 + pos):nrow(expandedTable), ])
          }
          pos = pos + 1
        } else{
          expandedTable[(i + pos), j] = NA
        }
      }
    }
  }
  row.names(expandedTable) = NULL


  for (i in 1:nrow(expandedTable)) {
    for (j in 3:ncol(expandedTable)) {
      if (!is.na(expandedTable[i, j])) {
        expandedTable[i, j] <- expandedTable[i, j] / (factorial(j - 2))
      }
    }
  }


  degree = nrow(expandedTable) - 1
  x = xs
  divDiff = expandedTable

  for (i in 1:degree){
    result = rep(0, nrow(expandedTable))

    for (j in 1:(nrow(expandedTable) - i)){
      currentCol = i + 1

      if ((currentCol + 1) > ncol(divDiff) || is.na(divDiff[j, (currentCol + 1)])){
        currentDiff = divDiff[(j + 1), currentCol]
        previousDiff = divDiff[j, currentCol]
        currentX = divDiff[(j + i), 1]
        prevX = divDiff[j, 1]

        result[j] = (currentDiff - previousDiff) / (currentX - prevX)
      } else{
        result[j] = divDiff[j, (currentCol + 1)]
      }
    }

    if (currentCol < ncol(expandedTable)){
      divDiff[, (currentCol + 1)] = result
    } else{
      divDiff = cbind(divDiff, result)
    }

    colnames(divDiff)[i+2] = paste("DD", i)
  }

  if (showTable){
    print("Divided differences: ")
    print(divDiff)
    print(separatorLine)
  }

  expr = as.character(divDiff[1, 2])
  matrix = matrix(c(0, divDiff[1, 2], divDiff[1, 2], divDiff[1, 2]),
                  nrow = nrow(expandedTable), ncol = 4, byrow = TRUE)
  colnames(matrix) = list("Order", "Coef", "Term value", "Aprox")

  for (i in 1:degree) {
    term = "1"

    for (j in 1:i) {
      term = paste(term, "*", paste("(x - ", divDiff[j, 1], ")", sep = ""))
    }

    expr = paste(expr, "+", as.character(divDiff[1, i + 2]), "*", term)
    aprox = eval(parse(text = expr))
    coefValue = eval(parse(text = paste(as.character(divDiff[1, i + 2]), "*", term)))
    matrix[i+1,] = c(i, divDiff[1, (i + 2)], coefValue, aprox)
  }
  expr = gsub(" 1 \\*", "", expr)

  if (showExpr){
    cat("P(x) = ", expr, "\n")
    print(separatorLine)
  }

  # Graphs
  polyFunc = function(x) {
    eval(parse(text = expr))
  }

  if (showGraph){
    curve_results = curve(polyFunc, from = min(expandedTable[,1]),
                          to = max(expandedTable[,1]))

    plot(expandedTable[,1], expandedTable[,2], pch = 21, cex = 1.5, lwd = 0.5,
         bg = "red", main = "Hermite Polynomial", xlab = "X", ylab = "P(x)",
         ylim = c(min(curve_results$y), max(curve_results$y)))

    lines(curve_results, col = "black", lwd = 2)

    points(xs, polyFunc(xs), pch = 21, cex = 1.8, lwd = 0.9, bg = "yellow")
  }


  if (showMatrix){
    print("Matrix: ")
    print(matrix)
    print(separatorLine)
  }

  cat("Aprox: ")

  return(polyFunc(xs))
}
