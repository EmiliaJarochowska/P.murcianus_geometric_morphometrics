plot_pca_scatter <- function(data, group_var, colors, title) {
  ggplot2::ggplot(data, ggplot2::aes_string(x = "PC1", y = "PC2", color = group_var)) +
    ggplot2::geom_point(size = 3, alpha = 0.8) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(title = title, x = "PC1", y = "PC2", color = group_var) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = "right",
      panel.grid = ggplot2::element_line(color = "grey90")
    )
}