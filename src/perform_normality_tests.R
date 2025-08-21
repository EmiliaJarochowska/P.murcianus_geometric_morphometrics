# Function to perform Shapiro-Wilk normality test of length by groups 
perform_normality_tests <- function(data, group_var) {
  groups <- unique(data[[group_var]])
  
  results <- data.frame(Group = character(0), W = numeric(0), p_value = numeric(0))
  
  for (g in groups) {
    subset_data <- data$Length[data[[group_var]] == g]
    if(length(subset_data) >= 3) {  # Shapiro needs at least 3 observations
      test <- shapiro.test(subset_data)
      results <- rbind(results, data.frame(Group = g, W = test$statistic, p_value = test$p.value))
    } else {
      results <- rbind(results, data.frame(Group = g, W = NA, p_value = NA))
    }
  }
  print(results)
}