#### Dependencies ####
library(vegan)
library(geomorph)
library(dunn.test)
library(egg)
library(ggpubr)
library(dplyr)

#### Statistical Analysis and Visualization Script ####

# Load and clean up data, set global variables
source("src/import_data.R")

# Or load from saved file if import_data.R had run previously
# load("data/processed_data.RData")

#### Load helper functions ####

# Function to create PCA scatter plots
source("src/plot_pca_scatter.R")
source("src/plot_pca_boxplot.R")
source("src/perform_pca_tests.R")
source("src/plot_correlation.R")

#### Statistical Analysis of Chirality ####

# PCA on chirality-filtered coordinates
pca_chirality <- geomorph::gm.prcomp(landmarks.gpa$coords)

# Test if PC scores differ between Left and Right
morphol.disparity(coords ~ 1, groups = specimen_info_matched$Chirality, 
                  data = landmarks.gpa,  
                  print.progress = TRUE)

# Plot PCA colored by chirality
par(mfrow = c(1, 1))
cols <- c("Left" = "blue", "Right" = "red")
plot(pca_chirality$x[, 1], pca_chirality$x[, 2], col = cols[specimen_info_matched$Chirality], pch = 19,
     xlab = "PC1", ylab = "PC2", main = "PCA of Shape by Chirality")
legend("topright", legend = levels(specimen_info_matched$Chirality), col = cols, pch = 19)

# Calculate mean shapes for each group
mean_shape_left <- geomorph::mshape(landmarks.gpa$coords[, , specimen_info_matched$Chirality == "Left"])
mean_shape_right <- geomorph::mshape(landmarks.gpa$coords[, , specimen_info_matched$Chirality == "Right"])

# Plot mean shapes with deformation grids
par(mfrow = c(1, 2))
geomorph::plotRefToTarget(mean_shape_left, mean_shape_right, method = "points",
                          main = "Mean Shape: Left -> Right")
geomorph::plotRefToTarget(mean_shape_right, mean_shape_left, method = "points",
                          main = "Mean Shape: Right -> Left")
par(mfrow = c(1, 1))

#### Analysis of Size Differences ####

combined_data_chirality <- data_combined %>% 
  dplyr::mutate(Chirality = ifelse(data_combined$Chirality == "R", "Right", "Left"))

print(table(combined_data_chirality$Chirality))
mean(combined_data_chirality$Length[combined_data_chirality$Chirality == "Right"])
mean(combined_data_chirality$Length[combined_data_chirality$Chirality == "Left"])

# Test for size differences
stats::kruskal.test(Length ~ Chirality, data = combined_data_chirality)

# Count specimens by Country and Group (Left/Right)
counts <- data_combined %>%
  dplyr::group_by(Country, Chirality) %>%
  dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
  dplyr::arrange(Country, Chirality)

print(counts)

##### Normality of the lengths #####

# Perform Shapiro-Wilk tests by Country, Section, Region on Length
perform_normality_tests(data_combined, "Country")
perform_normality_tests(data_combined, "Section")
perform_normality_tests(data_combined, "Region")

plot_length_histogram(data_combined, "Country", "Length by Country with Normal Distribution")
plot_length_histogram(data_combined, "Section", "Length by Section with Normal Distribution")
plot_length_histogram(data_combined, "Region", "Length by Region with Normal Distribution")

##### Length plots #####

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

FZ_length <- plot_distance_boxplot(data_combined, "FaciesZone", colors_list$FaciesZone,)

plot_distance_boxplot(data_combined, "Country", colors_list$Country)

Region_length <- plot_distance_boxplot(data_combined, "Region", colors_list$Region)

Section_length <- plot_distance_boxplot(data_combined, "Section", colors_list$Section)


ggarrange(Section_length, Region_length, FZ_length, 
          ncol=3, widths = c(3,2,2), 
          labels = c("A", "B", "C")) %>%
  ggsave(filename = "figs/lengths.pdf", width = 170, units = "mm")

data_combined %>%
  group_by(Section) %>%
  summarise(
    Median = median(Length),
    Mean = mean(Length)
  )

##### Kruskal-Wallis test #####

kruskal.test(Length ~ Section, data = data_combined)
dunn.test(data_combined$Length, g=data_combined$Section)

### Normality ####

check_normality_and_plot(data_combined, "Country")
check_normality_and_plot(data_combined, "Section")
check_normality_and_plot(data_combined, "Region")
check_normality_and_plot(data_combined, "Part")

#### PCA Analysis by Grouping Variables ####

grouping_vars <- c("FaciesZone", "Country", "Region", "Section", "Part")

for (var in grouping_vars) {
  if (var %in% names(data_combined) && var %in% names(colors_list)) {
    cat("\n--- Analysis for", var, "---\n")
    
    # Create plots
    p1 <- plot_pca_scatter(data_combined, var, colors_list[[var]], paste("PCA by", var))
    p2 <- plot_pca_boxplot(data_combined, var, colors_list[[var]], paste("PC Scores by", var))
    
    # Display plots
    print(p1)
    print(p2)
    
    # Perform statistical tests
    perform_pca_tests(data_combined, var)
  }
}

#### Correlation Analysis ####

# Overall correlation
cor_test_overall <- stats::cor.test(data_combined$PC1, data_combined$Length)
cat("Pearson correlation for entire dataset:\n")
cat("  R =", round(cor_test_overall$estimate, 3), 
    "| p-value =", signif(cor_test_overall$p.value, 3), "\n\n")

# Correlation plots by Region
if ("Region" %in% names(data_combined)) {
  p_region <- plot_correlation(data_combined, "PC1", "Length", "Region", 
                               colors_list$Region, "Relationship Between PC1 and Length by Region")
  print(p_region)
  
  # Calculate slopes by region
  calculate_slopes(data_combined, "PC1", "Length", "Region")
  
  # Kruskal-Wallis test by Region
  kruskal_result <- stats::kruskal.test(PC1 ~ Region, data = data_combined)
  print(kruskal_result)
  
  if (kruskal_result$p.value < 0.05) {
    dunn_result <- dunn.test::dunn.test(data_combined$PC1, data_combined$Region, method = "bonferroni")
    print(dunn_result)
  }
}

#### Shape Visualization ####

# Variance explained
variance_explained <- PCA$sdev^2 / sum(PCA$sdev^2)
cat(sprintf("PC1: %.2f%%\nPC2: %.2f%%\n", 
            variance_explained[1] * 100, variance_explained[2] * 100))

# Plot shape changes along PC axes
par(mfrow = c(2, 2))

# PC1 shape changes
geomorph::plotRefToTarget(PCA$shapes$shapes.comp1$min, msho, method = "vector", 
                          main = "PC1 Minimum")
geomorph::plotRefToTarget(PCA$shapes$shapes.comp1$max, msho, method = "vector",
                          main = "PC1 Maximum")

# PC2 shape changes  
geomorph::plotRefToTarget(PCA$shapes$shapes.comp2$min, msho, method = "vector",
                          main = "PC2 Minimum")
geomorph::plotRefToTarget(PCA$shapes$shapes.comp2$max, msho, method = "vector",
                          main = "PC2 Maximum")

par(mfrow = c(1, 1))

#### Summary ####

cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Total specimens analyzed:", nrow(data_combined), "\n")
cat("Specimens with valid chirality:", sum(valid_rows), "\n")
cat("PC1 variance explained:", round(variance_explained[1] * 100, 2), "%\n")
cat("PC2 variance explained:", round(variance_explained[2] * 100, 2), "%\n")
cat("Overall PC1-Length correlation:", round(cor_test_overall$estimate, 3), "\n")

cat("\nAnalysis complete!\n")