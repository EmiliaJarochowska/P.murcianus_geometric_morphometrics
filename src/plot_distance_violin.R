plot_distance_violin <- function(data, group_var, colors, title = NULL) {
  ggplot(data, aes_string(x = group_var, y = "Length", fill = group_var)) +
    geom_violin(alpha = 0.7, trim = FALSE) +
    geom_boxplot(
      width = 0.15,
      alpha = 0.8,
      outlier.shape = 16,
      outlier.size = 1.5
    ) +
    theme_minimal() +
    scale_fill_manual(values = colors) +
    scale_y_continuous(
      limits = c(100, 950),
      breaks = c(200, 400, 600, 800)
    ) +
    labs(
      x = group_var,
      y = "Length [μm]",
      title = title
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}