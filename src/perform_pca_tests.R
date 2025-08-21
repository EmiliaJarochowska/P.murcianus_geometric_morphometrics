# Function to perform statistical tests on PCA scores
perform_pca_tests <- function(data, group_var) {
  pca_long <- data %>%
    dplyr::select(PC1, PC2, sym(group_var)) %>%
    tidyr::pivot_longer(cols = c("PC1", "PC2"), names_to = "PC", values_to = "Score")
  
  for (pc in c("PC1", "PC2")) {
    pc_data <- dplyr::filter(pca_long, PC == pc)
    kw <- stats::kruskal.test(Score ~ get(group_var), data = pc_data)
    print(kw)
    if (kw$p.value < 0.05) {
      dunn <- dunn.test::dunn.test(pc_data$Score, pc_data[[group_var]], method = "bonferroni")
      print(dunn)
    } else {
      cat("No significant differences for", pc, "\n")
    }
  }
}