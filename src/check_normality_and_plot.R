check_normality_and_plot <- function(data, group_var) {
  cat("\n=== Normality Test for PC1 by", group_var, "===\n")
  groups <- unique(data[[group_var]])
  
  for (grp in groups) {
    grp_data <- data[data[[group_var]] == grp, "PC1", drop = TRUE]
    grp_data <- na.omit(grp_data)
    
    if (length(grp_data) >= 3) {
      test <- shapiro.test(grp_data)
      cat(sprintf("\nGroup: %s\n", grp))
      print(test)
      
      p <- ggplot(data.frame(PC1 = grp_data), aes(x = PC1)) +
        geom_histogram(aes(y = after_stat(density)), bins = 15, fill = "skyblue", color = "black", alpha = 0.7) +
        stat_function(fun = dnorm, args = list(mean = mean(grp_data), sd = sd(grp_data)), color = "red", size = 1.2) +
        labs(title = paste("PC1 Distribution -", grp), x = "PC1", y = "Density") +
        theme_minimal()
      
      print(p)
    } else {
      cat(sprintf("\nGroup: %s — Not enough data for Shapiro-Wilk test\n", grp))
    }
  }
}