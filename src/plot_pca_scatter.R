plot_pca_scatter <- function(data, group_var, colors) {
  ggplot(data, aes_string(x = "PC1", y = "PC2", fill = group_var, shape = group_var)) +
    geom_point(size = 2, alpha = 0.7, color = "black") +
    scale_fill_manual(values = colors) +
    scale_shape_manual(values = seq(from = 21, length.out = length(unique(data[[group_var]])))) +
    labs(x = "PC1", y = "PC2", fill = group_var, shape = group_var) +
    theme_minimal(base_size = 9) +
    theme(
      legend.position = "bottom",
      panel.grid = element_line(color = "grey90")
    )
}
