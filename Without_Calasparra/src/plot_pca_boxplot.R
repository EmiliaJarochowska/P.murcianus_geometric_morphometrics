# Function to create PCA boxplots
plot_pca_boxplot <- function(data, group_var, colors, title) {
  pca_long <- data %>%
    dplyr::select(PC1, PC2, sym(group_var)) %>%
    tidyr::pivot_longer(cols = c("PC1", "PC2"), names_to = "PC", values_to = "Score")
  
  ggplot2::ggplot(pca_long, ggplot2::aes_string(x = group_var, y = "Score", fill = group_var)) +
    ggplot2::geom_boxplot(alpha = 0.7) +
    ggplot2::facet_wrap(~ PC, scales = "free_y") +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = group_var, y = "PC Score", fill = group_var) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}