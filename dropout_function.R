# Author YS
# Date: 01/28/2025
# function for dorpout handling,
# M = vector of visit time points
# Missing at random with monotone missing pattern for longitudinal data
# Data needs to contain the following columns: id, time, response(outVariable)
# assign probability weight for each time point to be missing

introduce_missing <- function(df, outVariable = "response",
                              missing_percentage = 0.1,
                              prob = c(0.05, 0.1, 0.2, 0.3, 0.4)) {
  set.seed(123)  # For reproducibility
  unique_subjects <- unique(df$id)
  n <- length(unique_subjects)
  missing_subjects <- sample(unique_subjects, size = round(missing_percentage * n))

  for (id in missing_subjects) {
    subject_rows <- which(df$id == id)
    missing_start <- sample(subject_rows, 1, prob = prob)
    df[missing_start:max(subject_rows), variable] <- NA
  }
  return(df)
}

# Example usage
# df <- data.frame(
#   subject = rep(1:100, each = 5),
#   time = rep(1:5, times = 10),
#   response = rnorm(50)
# )
# df_with_missing_monotone <- introduce_missing_monotone_longitudinal(data, "chg")
# print(df_with_missing_monotone)


