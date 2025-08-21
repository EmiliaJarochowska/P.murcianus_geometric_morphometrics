plot_pca_scatter <- function(data, group_var, colors, title) {
  ggplot(data, aes_string(x = "PC1", y = "PC2", color = group_var, shape = group_var)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = seq(from=15, length.out = length(unique(data[[group_var]])))) +
    labs(title = title, x = "PC1", y = "PC2", color = group_var, shape = group_var) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      panel.grid = element_line(color = "grey90")
    )
}
