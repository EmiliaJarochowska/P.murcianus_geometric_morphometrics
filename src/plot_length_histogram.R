plot_length_histogram <- function(data, group_var) {
  ggplot(data, aes_string(x = "Length")) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
    stat_function(fun = dnorm, 
                  args = list(mean = mean(data$Length, na.rm = TRUE), 
                              sd = sd(data$Length, na.rm = TRUE)), 
                  color = "red", size = 1) +
    facet_wrap(as.formula(paste("~", group_var))) +
    theme_minimal() +
    labs(x = "Length [μm]", y = "Density") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}