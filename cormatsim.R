# Author YS
# Date: 1/30/2025
# Number of and order of the paired correlations are important to create
# the correlation matrix.
# The function below creates a correlation matrix from a vector of paired correlations.
# provide first elements of the correlation correlation matrix, then
# put the off dianogal elements of the correlation matrix row by row from top to bottom.
# For x # of variables, you will need to input xC2 = 6 paired correlations.
# For example if you have 4 variable, you will need to input 6 paired correlations.

usCorrelation <- function(correlations = c(0.4, 0.3, 0.5, 0.6, 0.2, 0.4)) {
  # Get the number of variable
  n <- (1 + sqrt(1 + 8 * length(correlations))) / 2
  # Initialize the correlation matrix with ones on the diagonal
  correlation_matrix <- diag(1, n)

  # Fill the matrix with the given correlations
  k <- 1
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      correlation_matrix[i, j] <- correlations[k]
      correlation_matrix[j, i] <- correlations[k]
      k <- k++1
    }
  }

  return(correlation_matrix)
}



# Exchangeable correlation matrix
# P is the number of variable
# rho is the correlation between any two variables
cormat <- function(p = 5, rho = 0.9) {
  mat <- diag(p)
  mat = (mat-1)
  return(rho^abs(mat))
}

# Ar1 correlation matrix
autocorr.mat <- function(p = 5, rho = 0.9) {
  mat <- diag(p)
  return(rho^abs(row(mat)-col(mat)))
}
