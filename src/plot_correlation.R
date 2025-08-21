# Function to create correlation plots
plot_correlation <- function(data, x_var, y_var, group_var, colors, title) {
  # Calculate correlation
  cor_test <- stats::cor.test(data[[x_var]], data[[y_var]])
  
  p <- ggplot2::ggplot(data, ggplot2::aes_string(x = x_var, y = y_var)) +
    ggplot2::geom_point(ggplot2::aes_string(color = group_var), size = 3, alpha = 0.8) +
    ggplot2::geom_smooth(ggplot2::aes_string(color = group_var), method = "lm", se = FALSE, linetype = "solid") +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_x_continuous(
      name = paste(x_var, "Score"),
      labels = scales::number_format(accuracy = 0.1),
      breaks = scales::pretty_breaks()
    ) +
    ggplot2::scale_y_continuous(
      name = "Length (μm)",
      labels = scales::number_format(accuracy = 0.1),
      breaks = scales::pretty_breaks()
    ) +
    ggplot2::labs(
      title = title,
      subtitle = paste("R =", round(cor_test$estimate, 3), 
                       "| p-value =", signif(cor_test$p.value, 3)),
      color = group_var
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = "right",
      panel.grid = ggplot2::element_line(color = "grey90")
    )
  
  return(p)
}