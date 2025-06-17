### STATISTICAL DIFFERENCE IN SHAPE BETWEEN LEFT AND RIGHT ELEMENTS ###
# Load required libraries
library(geomorph)
library(readxl)

# 1. Load TPS file and metadata Excel file
tps_file <- "All_sections.TPS"
specimen_info_path <- "Specimens_info.xlsx"

landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)
specimen_info <- read_excel(specimen_info_path)

# 2. Clean specimen IDs: trim whitespace and unify case
specimen_info$ID <- trimws(specimen_info$ID)
specimen_info$ID <- toupper(specimen_info$ID)

tps_ids <- dimnames(landmarks)[[3]]
tps_ids <- trimws(tps_ids)
tps_ids <- toupper(tps_ids)
dimnames(landmarks)[[3]] <- tps_ids  # overwrite with cleaned IDs

# 3. Check for unmatched specimen IDs between TPS and Excel
missing_in_excel <- tps_ids[!tps_ids %in% specimen_info$ID]
missing_in_tps <- specimen_info$ID[!specimen_info$ID %in% tps_ids]

if (length(missing_in_excel) > 0) {
  cat("The following TPS specimen IDs are NOT found in specimen_info:\n")
  print(missing_in_excel)
}

if (length(missing_in_tps) > 0) {
  cat("The following specimen_info IDs are NOT found in TPS:\n")
  print(missing_in_tps)
}

# Stop if any missing IDs to avoid misalignment
if (length(missing_in_excel) > 0 || length(missing_in_tps) > 0) {
  stop("Mismatch between TPS IDs and specimen_info IDs. Please fix and rerun.")
}

# 4. Match specimen_info order to TPS specimen IDs
match_rows <- match(tps_ids, specimen_info$ID)
specimen_info_matched <- specimen_info[match_rows, ]

# 5. Define sliding landmarks (adjust indices to your dataset)
sliders <- define.sliders(3:138)

# 6. Perform Generalized Procrustes Analysis (GPA)
gpa <- gpagen(landmarks, curves = sliders)

# 7. Filter specimens with Chirality = "L" or "R"
valid_rows <- specimen_info_matched$Chirality %in% c("L", "R")
specimen_info_filtered <- specimen_info_matched[valid_rows, ]
gpa_coords_filtered <- gpa$coords[,,valid_rows]

# 8. Create grouping factor
group <- factor(ifelse(specimen_info_filtered$Chirality == "R", "Right", "Left"))

# 9. PCA on GPA-aligned coordinates
pca <- gm.prcomp(gpa_coords_filtered)

# 10. Plot PCA colored by chirality
cols <- c("Left" = "blue", "Right" = "red")
plot(pca$x[,1], pca$x[,2], col = cols[group], pch = 19,
     xlab = "PC1", ylab = "PC2", main = "PCA of Shape by Chirality")
legend("topright", legend = levels(group), col = cols, pch = 19)

# 11. Calculate mean shapes for each group
mean_shape_left <- mshape(gpa_coords_filtered[,,group == "Left"])
mean_shape_right <- mshape(gpa_coords_filtered[,,group == "Right"])

# 12. Plot mean shapes with deformation grids
par(mfrow = c(1, 2))
plotRefToTarget(mean_shape_left, mean_shape_right, method = "points",
                main = "Mean Shape: Left -> Right")
plotRefToTarget(mean_shape_right, mean_shape_left, method = "points",
                main = "Mean Shape: Right -> Left")
par(mfrow = c(1, 1))

# 13. Procrustes ANOVA test for shape difference by Chirality
gdf <- geomorph.data.frame(coords = gpa_coords_filtered, group = group)
proc_test <- procD.lm(coords ~ group, data = gdf, iter = 1000)
summary(proc_test)





### STATISTICAL DIFFERENCE IN LENGTH BETWEENLEFT AND RIGHT ELEMENTS ###

library(dplyr)
library(readxl)

# Function to calculate Euclidean distance
calculate_distance <- function(x1, y1, x2, y2) {
  return(sqrt((x2 - x1)^2 + (y2 - y1)^2))
}

# Function to process TPS file and compute scaled distances
process_data <- function(file_path) {
  lines <- readLines(file_path)
  lm_indices <- which(grepl("LM=", lines))
  scale_indices <- which(grepl("SCALE=", lines))
  id_indices <- which(grepl("ID=", lines))
  
  if(length(lm_indices) == 0 || length(scale_indices) == 0 || length(id_indices) == 0) {
    stop("Missing LM=, SCALE=, or ID= fields in TPS file.")
  }
  
  results <- data.frame(ID = character(), Mean_Distance = numeric(), stringsAsFactors = FALSE)
  
  for (i in seq_along(lm_indices)) {
    num_landmarks <- as.numeric(gsub("LM=", "", lines[lm_indices[i]]))
    scale <- as.numeric(gsub("SCALE=", "", lines[scale_indices[i]]))
    specimen_id <- gsub("ID=", "", lines[id_indices[i]])
    specimen_id <- gsub("\\s+", "", specimen_id)  # clean spaces
    
    landmark_lines <- lines[(lm_indices[i] + 1):(lm_indices[i] + num_landmarks)]
    landmarks <- do.call(rbind, strsplit(landmark_lines, "\\s+"))
    landmarks <- as.data.frame(landmarks, stringsAsFactors = FALSE)
    landmarks <- mutate_all(landmarks, as.numeric)
    
    if(nrow(landmarks) < 2 || is.na(scale)) {
      next 
    }
    
    # Calculate all pairwise distances
    distances <- numeric()
    for (j in 1:(nrow(landmarks) - 1)) {
      for (k in (j + 1):nrow(landmarks)) {
        dist <- calculate_distance(landmarks[j, 1], landmarks[j, 2], landmarks[k, 1], landmarks[k, 2])
        scaled_distance <- dist * scale
        distances <- c(distances, scaled_distance)
      }
    }
    
    mean_dist <- mean(distances)
    results <- rbind(results, data.frame(ID = specimen_id, Mean_Distance = mean_dist, stringsAsFactors = FALSE))
  }
  
  return(results)
}

# Load specimen metadata from Excel
specimen_info <- read_excel("Specimens_info.xlsx")
specimen_info$ID <- gsub("\\s+", "", specimen_info$ID)
specimen_info$Chirality <- toupper(specimen_info$Chirality)

# Process TPS file
tps_results <- process_data("All_sections.TPS")

# Merge distances with metadata
combined_data <- merge(tps_results, specimen_info[, c("ID", "Chirality")], by = "ID")

# Filter for valid chirality values
combined_data <- combined_data %>% filter(Chirality %in% c("L", "R"))
combined_data$Group <- ifelse(combined_data$Chirality == "R", "Right", "Left")

# Summary
print(table(combined_data$Group))
cat("Mean Distance for Right group:", mean(combined_data$Mean_Distance[combined_data$Group == "Right"]), "\n")
cat("Mean Distance for Left group:", mean(combined_data$Mean_Distance[combined_data$Group == "Left"]), "\n")

# Kruskal-Wallis test
if (length(unique(combined_data$Group)) >= 2) {
  kruskal_test <- kruskal.test(Mean_Distance ~ Group, data = combined_data)
  print(kruskal_test)
} else {
  cat("Not enough groups for Kruskal-Wallis test.\n")
}



### number of lft and right specimens for each Country ###
# Load dplyr for convenient data manipulation
library(dplyr)

# Assuming specimen_info_matched has columns: Country, Chirality
# Filter for valid chirality values
filtered_data <- specimen_info_matched %>% 
  filter(Chirality %in% c("L", "R")) %>%
  mutate(Group = ifelse(Chirality == "R", "Right", "Left"))

# Count specimens by Country and Group (Left/Right)
counts <- filtered_data %>%
  group_by(Country, Group) %>%
  summarise(Count = n()) %>%
  arrange(Country, Group)

print(counts)




#### SHAPE ####

# --- Load Libraries ---
library(geomorph)
library(readxl)
library(dplyr)
library(ggplot2)
library(FSA)
library(dunn.test)
library(tidyr)
library(ggpubr)
# --- Distance Function ---
calculate_distance <- function(x1, y1, x2, y2) {
  sqrt((x2 - x1)^2 + (y2 - y1)^2)
}

# --- Process TPS File to Calculate Scaled Mean Distances ---
process_data <- function(file_path) {
  lines <- readLines(file_path)
  lm_indices <- which(grepl("LM=", lines))
  scale_indices <- which(grepl("SCALE=", lines))
  id_indices <- which(grepl("ID=", lines))
  
  results <- data.frame(ID = character(), Mean_Distance = numeric(), stringsAsFactors = FALSE)
  
  for (i in seq_along(lm_indices)) {
    num_landmarks <- as.numeric(gsub("LM=", "", lines[lm_indices[i]]))
    scale <- as.numeric(gsub("SCALE=", "", lines[scale_indices[i]]))
    specimen_id <- gsub("ID=", "", lines[id_indices[i]])
    specimen_id <- gsub("\\s+", "", specimen_id)
    
    landmark_lines <- lines[(lm_indices[i] + 1):(lm_indices[i] + num_landmarks)]
    landmarks <- do.call(rbind, strsplit(landmark_lines, "\\s+"))
    landmarks <- as.data.frame(landmarks, stringsAsFactors = FALSE)
    landmarks <- mutate_all(landmarks, as.numeric)
    
    if(nrow(landmarks) < 2 || is.na(scale)) next
    
    distances <- numeric()
    for (j in 1:(nrow(landmarks) - 1)) {
      for (k in (j + 1):nrow(landmarks)) {
        dist <- calculate_distance(landmarks[j, 1], landmarks[j, 2], landmarks[k, 1], landmarks[k, 2])
        distances <- c(distances, dist * scale)
      }
    }
    
    mean_dist <- mean(distances)
    results <- rbind(results, data.frame(ID = specimen_id, Mean_Distance = mean_dist, stringsAsFactors = FALSE))
  }
  return(results)
}
# --- File Paths ---
tps_file <- "All_sections.TPS"
metadata_path <- "Specimens_info.xlsx"


landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)
specimen_info <- read_excel(metadata_path)

# --- Clean and Sync IDs ---
specimen_info$ID <- toupper(trimws(specimen_info$ID))
tps_ids <- toupper(trimws(dimnames(landmarks)[[3]]))
dimnames(landmarks)[[3]] <- tps_ids
match_idx <- match(tps_ids, specimen_info$ID)
specimen_info <- specimen_info[match_idx, ]

# --- GPA and PCA ---
sliders <- define.sliders(3:138)
landmarks.gpa <- gpagen(landmarks, curves = sliders)
PCA <- gm.prcomp(landmarks.gpa$coords)

# --- Calculate Mean Distances ---
mean_distances <- process_data(tps_file)
mean_distances$ID <- toupper(trimws(mean_distances$ID))

# --- Combine All Data ---
data_combined <- specimen_info %>%
  left_join(mean_distances, by = "ID")

# --- Fill NA Mean_Distance with Median ---
data_combined$Mean_Distance[is.na(data_combined$Mean_Distance)] <- median(data_combined$Mean_Distance, na.rm = TRUE)

# --- Add PCA Scores ---
data_combined$PC1 <- PCA$x[, 1]
data_combined$PC2 <- PCA$x[, 2]

# --- Recode Region from Country ---
data_combined$Region <- dplyr::recode_factor(
  data_combined$Country,
  "Slovenia" = "North-Eastern part of Sephardic Province",
  "Bosnia and Herzegovina" = "North-Eastern part of Sephardic Province",
  "Spain" = "Western part of Sephardic Province",
  .default = NA_character_
)

# --- Define Section colors and custom order ---
section_order <- c("Calasparra", "Henarejos", "Libros", "Bugarra", "Prikrnica", "Drežnica")
section_colors <- c(
  "Calasparra" = "blue",
  "Henarejos" = "#D4A017",
  "Libros" = "darkgreen",
  "Bugarra" = "#17BCC4",
  "Prikrnica" = "#D91C93",
  "Drežnica" = "lightgreen"
)

# Apply factor order
data_combined$Section <- factor(data_combined$Section, levels = section_order)

# --- Color Palettes ---
colors_list <- list(
  FaciesZone = c("FZ3" = "grey", "FZ7" = "#0033A0", "FZ8" = "#D81B60"),
  Country = c("Slovenia" = "grey", "Bosnia and Herzegovina" = "#0033A0", "Spain" = "#D81B60"),
  Region = c("North-Eastern part of Sephardic Province" = "grey", "Western part of Sephardic Province" = "#0033A0"),
  Section = section_colors
)

# --- Define Par variable as subregion of Sephardic Province ---
data_combined$Par <- NA_character_

data_combined$Par[data_combined$Section %in% c("Calasparra", "Henarejos", "Libros", "Bugarra")] <- "Western Subprovince"
data_combined$Par[data_combined$Section %in% c("Prikrnica", "Drežnica")] <- "North-Eastern Subprovince"

data_combined$Par <- factor(data_combined$Par)

# --- Normality Test and Distribution Plots for PC1 by Grouping Variables ---
check_normality_and_plot <- function(data, group_var) {
  cat("\n=== Normality Test for PC1 by", group_var, "===\n")
  groups <- unique(data[[group_var]])
  
  for (grp in groups) {
    grp_data <- data[data[[group_var]] == grp, "PC1", drop = TRUE]
    grp_data <- na.omit(grp_data)
    
    if (length(grp_data) >= 3) {
      test <- shapiro.test(grp_data)
      cat(sprintf("\nGroup: %s\n", grp))
      print(test)
      
      p <- ggplot(data.frame(PC1 = grp_data), aes(x = PC1)) +
        geom_histogram(aes(y = after_stat(density)), bins = 15, fill = "skyblue", color = "black", alpha = 0.7) +
        stat_function(fun = dnorm, args = list(mean = mean(grp_data), sd = sd(grp_data)), color = "red", size = 1.2) +
        labs(title = paste("PC1 Distribution -", grp), x = "PC1", y = "Density") +
        theme_minimal()
      
      print(p)
    } else {
      cat(sprintf("\nGroup: %s — Not enough data for Shapiro-Wilk test\n", grp))
    }
  }
}

check_normality_and_plot(data_combined, "Country")
check_normality_and_plot(data_combined, "Section")
check_normality_and_plot(data_combined, "Region")
check_normality_and_plot(data_combined, "Par")

# --- Plotting Functions ---
plot_pca_scatter <- function(data, group_var, colors, title) {
  ggplot(data, aes_string(x = "PC1", y = "PC2", color = group_var, size = "Mean_Distance")) +
    geom_point(alpha = 0.8) +
    stat_ellipse(aes_string(group = group_var), linewidth = 1, show.legend = FALSE) +
    scale_color_manual(values = colors) +
    scale_size_continuous(range = c(2, 7)) +
    theme_minimal() +
    labs(title = title, x = "PC1", y = "PC2", color = group_var, size = "Mean Scaled Distance") +
    guides(color = guide_legend(override.aes = list(size = 4)), size = guide_legend())
}

plot_pca_boxplot <- function(data, group_var, colors, title) {
  pca_long <- data %>%
    select(PC1, PC2, !!sym(group_var)) %>%
    pivot_longer(cols = c("PC1", "PC2"), names_to = "PC", values_to = "Score")
  
  ggplot(pca_long, aes_string(x = group_var, y = "Score", fill = group_var)) +
    geom_boxplot(alpha = 0.7) +
    facet_wrap(~ PC, scales = "free_y") +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    labs(title = title, x = group_var, y = "PC Score", fill = group_var) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

perform_pca_tests <- function(data, group_var) {
  pca_long <- data %>%
    select(PC1, PC2, !!sym(group_var)) %>%
    pivot_longer(cols = c("PC1", "PC2"), names_to = "PC", values_to = "Score")
  
  for (pc in c("PC1", "PC2")) {
    cat("\nTesting", pc, "by", group_var, "\n")
    pc_data <- filter(pca_long, PC == pc)
    kw <- kruskal.test(Score ~ get(group_var), data = pc_data)
    print(kw)
    if (kw$p.value < 0.05) {
      dunn <- dunn.test(pc_data$Score, pc_data[[group_var]], method = "bonferroni")
      print(dunn)
    } else {
      cat("No significant differences for", pc, "\n")
    }
  }
}

# --- Plot and Test for Each Grouping Variable ---
plot_pca_scatter(data_combined, "FaciesZone", colors_list$FaciesZone, "PCA by Facies Zone")
plot_pca_boxplot(data_combined, "FaciesZone", colors_list$FaciesZone, "PC Scores by Facies Zone")
perform_pca_tests(data_combined, "FaciesZone")

plot_pca_scatter(data_combined, "Country", colors_list$Country, "PCA by Country")
plot_pca_boxplot(data_combined, "Country", colors_list$Country, "PC Scores by Country")
perform_pca_tests(data_combined, "Country")

plot_pca_scatter(data_combined, "Region", colors_list$Region, "PCA by Region")
plot_pca_boxplot(data_combined, "Region", colors_list$Region, "PC Scores by Region")
perform_pca_tests(data_combined, "Region")

plot_pca_scatter(data_combined, "Section", colors_list$Section, "PCA by Section")
plot_pca_boxplot(data_combined, "Section", colors_list$Section, "PC Scores by Section")
perform_pca_tests(data_combined, "Section")

plot_pca_scatter(data_combined, "Par", c("Western Subprovince" = "#0033A0", "North-Eastern Subprovince" = "grey"), "PCA by Par Subprovince")
plot_pca_boxplot(data_combined, "Par", c("Western Subprovince" = "#0033A0", "North-Eastern Subprovince" = "grey"), "PC Scores by Par Subprovince")
perform_pca_tests(data_combined, "Par")

# --- PCA Shape Visualization ---
msho <- mshape(landmarks.gpa$coords)
summary(msho)

plotRefToTarget(PCA$shapes$shapes.comp1$min, msho, method = "vector")
plotRefToTarget(PCA$shapes$shapes.comp1$max, msho, method = "vector")
plotRefToTarget(PCA$shapes$shapes.comp2$min, msho, method = "vector")
plotRefToTarget(PCA$shapes$shapes.comp2$max, msho, method = "vector")

# --- Variance Explained ---
variance_explained <- PCA$sdev^2 / sum(PCA$sdev^2)
cat(sprintf("\nPC1: %.2f%%\nPC2: %.2f%%\n", variance_explained[1]*100, variance_explained[2]*100))

# --- Working on residuals PC1 ---
# Calculate PC1 residuals after removing effect of Mean_Distance
lm_pc1 <- lm(PC1 ~ Mean_Distance, data = data_combined)
data_combined$PC1_residual <- residuals(lm_pc1)

kruskal_dunn_plot <- function(data, group_var, colors = NULL) {
  cat(sprintf("\n--- Kruskal-Wallis test on PC1 residuals by %s ---\n", group_var))
  kruskal_res <- kruskal.test(PC1_residual ~ get(group_var), data = data)
  print(kruskal_res)
  
  if (kruskal_res$p.value < 0.05) {
    cat("Significant differences found, performing Dunn's post hoc test:\n")
    dunn_res <- dunn.test(data$PC1_residual, data[[group_var]], method = "bonferroni")
    print(dunn_res)
  } else {
    cat("No significant differences found in Kruskal-Wallis test.\n")
  }
  
  p <- ggplot(data, aes_string(x = group_var, y = "PC1_residual", fill = group_var)) +
    geom_boxplot(alpha = 0.7) +
    labs(title = paste("PC1 Residuals by", group_var), x = group_var, y = "PC1 Residuals") +
    theme_minimal()
  
  if (!is.null(colors)) {
    p <- p + scale_fill_manual(values = colors)
  }
  
  print(p)
}

kruskal_dunn_plot(data_combined, "Section", colors_list$Section)
kruskal_dunn_plot(data_combined, "FaciesZone", colors_list$FaciesZone)
kruskal_dunn_plot(data_combined, "Country", colors_list$Country)
kruskal_dunn_plot(data_combined, "Region", colors_list$Region)
kruskal_dunn_plot(data_combined, "Par", c("Western Subprovince" = "#0033A0", "North-Eastern Subprovince" = "grey"))







#### Length of elements ####


library(ggplot2)
library(dunn.test)
# Function to plot histogram + normal distribution overlay for Mean_Distance (Length proxy)
plot_length_histogram <- function(data, group_var, title) {
  ggplot(data, aes_string(x = "Mean_Distance")) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
    stat_function(fun = dnorm, 
                  args = list(mean = mean(data$Mean_Distance, na.rm = TRUE), 
                              sd = sd(data$Mean_Distance, na.rm = TRUE)), 
                  color = "red", size = 1) +
    facet_wrap(as.formula(paste("~", group_var))) +
    theme_minimal() +
    labs(title = title, x = "Mean Scaled Distance (Length)", y = "Density") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Function to perform Shapiro-Wilk normality test by groups on Mean_Distance
perform_normality_tests <- function(data, group_var) {
  cat("\n--- Shapiro-Wilk Normality Test for Mean_Distance by", group_var, "---\n")
  data_filtered <- data[!is.na(data[[group_var]]) & !is.na(data$Mean_Distance), ]
  
  groups <- unique(data_filtered[[group_var]])
  
  results <- data.frame(Group = character(0), W = numeric(0), p_value = numeric(0))
  
  for (g in groups) {
    subset_data <- data_filtered$Mean_Distance[data_filtered[[group_var]] == g]
    if(length(subset_data) >= 3) {  # Shapiro needs at least 3 observations
      test <- shapiro.test(subset_data)
      results <- rbind(results, data.frame(Group = g, W = test$statistic, p_value = test$p.value))
    } else {
      results <- rbind(results, data.frame(Group = g, W = NA, p_value = NA))
    }
  }
  print(results)
}

# Plot histograms + normal distribution by Country, Section, Region using corrected column
plot_length_histogram(data_combined, "Country", "Length (Mean Scaled Distance) by Country with Normal Distribution")
plot_length_histogram(data_combined, "Section", "Length (Mean Scaled Distance) by Section with Normal Distribution")
plot_length_histogram(data_combined, "Region", "Length (Mean Scaled Distance) by Region with Normal Distribution")

# Perform Shapiro-Wilk tests by Country, Section, Region on Length
perform_normality_tests(data_combined, "Country")
perform_normality_tests(data_combined, "Section")
perform_normality_tests(data_combined, "Region")


plot_distance_boxplot <- function(data, group_var, colors, title) {
  ggplot(data, aes_string(x = group_var, y = "Mean_Distance", fill = group_var)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    labs(title = title, x = group_var, y = "Mean Scaled Distance", fill = group_var) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}


# --- Plot Scaled Mean Distance by Facies Zone ---
plot_distance_boxplot(data_combined, "FaciesZone", colors_list$FaciesZone, "Mean Scaled Distance by Facies Zone")

# --- Plot Scaled Mean Distance by Country ---
plot_distance_boxplot(data_combined, "Country", colors_list$Country, "Mean Scaled Distance by Country")

# --- Plot Scaled Mean Distance by Sephardic Province Region ---
plot_distance_boxplot(data_combined, "Region", colors_list$Region, "Mean Scaled Distance by Part of Sephardic Province")

# --- Plot Scaled Mean Distance by Section ---
plot_distance_boxplot(data_combined, "Section", colors_list$Section, "Mean Scaled Distance by Section")


perform_distance_tests <- function(data, group_var) {
  cat("\n--- Kruskal-Wallis and Dunn Test for", group_var, "---\n")
  
  # Subset non-missing data
  data_filtered <- data[!is.na(data[[group_var]]), ]
  
  # Kruskal-Wallis test
  kruskal <- kruskal.test(data_filtered$Mean_Distance ~ data_filtered[[group_var]])
  print(kruskal)
  
  # Dunn test if significant
  if (kruskal$p.value < 0.05) {
    dunn <- dunn.test(data_filtered$Mean_Distance, data_filtered[[group_var]], method = "bonferroni")
    print(dunn)
  } else {
    cat("No significant differences found (p =", round(kruskal$p.value, 4), ")\n")
  }
}


perform_distance_tests(data_combined, "FaciesZone")
perform_distance_tests(data_combined, "Country")
perform_distance_tests(data_combined, "Region")
perform_distance_tests(data_combined, "Section")






### relationship between Length and PC1 scores ###

# For simplicity, let's rename Mean_Distance as Length (or create Length if needed)
data_combined$Length <- data_combined$Mean_Distance

# Optional: Create a grouping variable similar to your example (adjust country names accordingly)
data_combined$group <- ifelse(data_combined$Country %in% c("Slovenia", "Bosnia and Herzegovina"), 
                              "North-Eastern part", "Western part")

# Plot Length vs PC1 with colors by group and shapes by Section
ggplot(data_combined, aes(x = PC1, y = Length, color = group, shape = Section)) +
  geom_point(aes(size = Length), alpha = 0.8) +        # Points sized by Length
  geom_smooth(aes(group = group), method = "lm", se = FALSE) +  # Linear trend lines per group
  scale_color_manual(
    name = "Sephardic Province Region",
    values = c("North-Eastern part" = "black", "Western part" = "lightblue")
  ) +
  scale_shape_manual(
    name = "Section",
    values = c(16, 17, 18, 19, 15, 13),  # Use different point shapes
    breaks = levels(data_combined$Section)
  ) +
  scale_size_continuous(range = c(2, 6)) +
  labs(
    title = "Relationship between Length and PC1",
    x = "PC1 Score",
    y = "Length (Mean Distance)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right"
  )

# Optional: Calculate correlation and linear model slopes per group
library(dplyr)

# Correlation test for entire dataset
cor_test <- cor.test(data_combined$PC1, data_combined$Length)
print(paste("Pearson correlation: ", round(cor_test$estimate, 3), "p-value:", signif(cor_test$p.value, 3)))

# Fit linear models per group
lm_north <- lm(Length ~ PC1, data = filter(data_combined, group == "North-Eastern part"))
lm_west <- lm(Length ~ PC1, data = filter(data_combined, group == "Western part"))

cat("Slope for North-Eastern part:", coef(lm_north)["PC1"], "\n")
cat("Slope for Western part:", coef(lm_west)["PC1"], "\n")


## residuals PC1 VS Length ##

# Fit linear model to remove effect of size on PC1
lm_pc1_size <- lm(PC1 ~ Mean_Distance, data = data_combined)

# Add residuals as a new column in your dataframe
data_combined$PC1_residuals <- residuals(lm_pc1_size)
ggplot(data_combined, aes(x = PC1_residuals, y = Length, color = group, shape = Section)) +
  geom_point(aes(size = Length), alpha = 0.8) +
  geom_smooth(aes(group = group), method = "lm", se = FALSE) +  # explicit grouping, no warning
  scale_color_manual(
    name = "Sephardic Province Region",
    values = c("North-Eastern part" = "black", "Western part" = "lightblue")
  ) +
  scale_shape_manual(
    name = "Section",
    values = c(16, 17, 18, 19, 15, 13),
    breaks = levels(data_combined$Section)
  ) +
  scale_size_continuous(range = c(2, 6)) +
  labs(
    title = "Relationship between Length and PC1 Residuals (Size-Corrected)",
    x = "PC1 Residuals",
    y = "Length (Mean Distance)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right"
  )
