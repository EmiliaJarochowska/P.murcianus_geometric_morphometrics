plot_distance_boxplot <- function(data, group_var, colors, title) {
  ggplot(data, aes_string(x = group_var, y = "Length", fill = group_var)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
    theme_minimal() +
    scale_fill_manual(values = colors) +
    labs(x = group_var, y = "Length [um]") +
    theme(legend.position = "none") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}