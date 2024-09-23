### IN ANAYSE ARE INCORPORATE ELEMENTS FROM SLOVENIA, SPAIN AND BOSNIA & HERZEGOVINA ###

library(geomorph)
tps_file <- "skupno.TPS"
landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)
sliders <- define.sliders(3:138) 

### Mean shape, PC1 min and max, PC2 min and max ###

#procrusteranalyse
landmarks.gpa<-gpagen(landmarks, curves = sliders)
plot(landmarks.gpa)
PCA <- gm.prcomp(landmarks.gpa$coords) 
plot(PCA)
PCA
msho <- mshape(landmarks.gpa$coords) 
plot(msho)
summary(msho)
plotRefToTarget(PCA$shapes$shapes.comp1$min, msho,                 
                method = "vector")
plotRefToTarget(PCA$shapes$shapes.comp1$max, msho,                 
                method = "vector")
plotRefToTarget(PCA$shapes$shapes.comp2$min, msho,                 
                method = "vector")
plotRefToTarget(PCA$shapes$shapes.comp2$max, msho,                 
                method = "vector")
require(geomorph)


### Principal Component Analysis ####

# Install and load necessary packages
install.packages("ggplot2")
library(ggplot2)
library(geomorph)

tps_file <- "skupno.TPS"
landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)
sliders <- define.sliders(3:138) 

# Assuming `landmarks` is already defined and loaded
samples <- data.frame(name = rep(NA, dim(landmarks)[3]),
                      country = rep(NA, dim(landmarks)[3]),
                      variable1 = rep(NA, dim(landmarks)[3]),
                      variable2 = rep(NA, dim(landmarks)[3]),
                      variable3 = rep(NA, dim(landmarks)[3]))

for (i in 1:dim(landmarks)[3]) {
  samples$name[i] <- dimnames(landmarks)[[3]][i]
  samples$country[i] <- unlist(strsplit(samples$name[i], "_"))[1]
  samples$variable1[i] <- unlist(strsplit(samples$name[i], "_"))[2]
  samples$variable2[i] <- unlist(strsplit(samples$name[i], "_"))[3]
  samples$variable3[i] <- unlist(strsplit(samples$name[i], "_"))[4]
}

samples$country <- as.factor(samples$country)
samples$variable1 <- as.factor(samples$variable1)

# Perform Procrustes analysis
sliders <- define.sliders(3:138)
landmarks.gpa <- gpagen(landmarks, curves = sliders)
plot(landmarks.gpa)
plotAllSpecimens(A = landmarks.gpa$coords, mean = TRUE, label = FALSE)

# Perform PCA
PCA <- gm.prcomp(landmarks.gpa$coords)

# Extract PCA scores and create a data frame
pca_data <- data.frame(PC1 = PCA$x[, 1], PC2 = PCA$x[, 2], Country = samples$country)
# Perform Kruskal-Wallis test on PC1 scores by country
kruskal_test <- kruskal.test(PC1 ~ Country, data = pca_data)
print(kruskal_test)
# Perform Dunn's test for pairwise comparison between countries
install.packages("dunn.test")
library(dunn.test)

# Perform Dunn's test on PC1 scores
dunn_test <- dunn.test(pca_data$PC1, pca_data$Country, method = "bonferroni")
print(dunn_test)



### Print % of explained PC scores ###

# Assuming PCA and samples are already defined and PCA$x contains the PCA scores

# Create a data frame for PCA scores
pca_scores_df <- data.frame(
  PC1 = PCA$x[, 1],
  PC2 = PCA$x[, 2],
  country = samples$country  # Assuming you have the samples dataframe ready
)

# Extract variance explained
variance_explained <- PCA$sdev^2 / sum(PCA$sdev^2)  # Eigenvalues / Total eigenvalues
variance_explained_percent <- variance_explained * 100

# Print the percentage of variance explained
cat("Percentage variance explained by each principal component:\n")
for (i in 1:length(variance_explained_percent)) {
  cat(sprintf("PC%d: %.2f%%\n", i, variance_explained_percent[i]))
}


### Principal Component Analysis with confidence ellipse ####

# Define country full names
country_full_names <- c(
  "SL" = "Slovenia",
  "SP" = "Spain",
  "BAH" = "Bosnia and Herzegovina"
)

# Convert country codes in samples$country to full names
samples$country <- country_full_names[samples$country]

# Ensure factor levels in pca_scores_df match the names in scale_color_manual
levels(pca_scores_df$country) <- country_full_names

# Create the plot
ggplot(pca_scores_df, aes(x = PC1, y = PC2, color = country, shape = country)) +
  geom_point(size = 3) +  # Adjust size if needed
  stat_ellipse(level = 0.95) +  # 95% confidence ellipse
  labs(x = "PC 1", y = "PC 2", title = "Morphospace of the aboral side of P. murcianus") +
  scale_color_manual(values = c(
    "Bosnia and Herzegovina" = "#D81B60",  # Dark yellow
    "Slovenia" = "grey",                        # Black
    "Spain" = "#0033A0"                            # Blue
  )) +
  scale_shape_manual(values = 15 + 0:(length(country_full_names) - 1)) +
  theme_minimal() +
  theme(legend.position = "bottom",  # Place legend at the bottom
        legend.box = "horizontal",  # Arrange legend horizontally
        legend.background = element_rect(fill = "transparent")) +  # Transparent background for legend
  annotate("text", x = Inf, y = -Inf, label = sprintf("PC1: %.2f%%", variance_explained_percent[1]), 
           hjust = 1.1, vjust = -0.5, size = 4, fontface = "italic") +
  annotate("text", x = -Inf, y = Inf, label = sprintf("PC2: %.2f%%", variance_explained_percent[2]), 
           hjust = -0.1, vjust = 1.5, size = 4, fontface = "italic")




### Print the ID od specific specimens who are closer to the PC1 and PC2 min and max ####

# Create a sample data frame assuming you have 'ID' or 'IMAGE' in your dataset
samples <- data.frame(
  ID = dimnames(landmarks)[[3]],  # Replace this with appropriate identifier if different
  PC1 = PCA$x[, 1],
  PC2 = PCA$x[, 2]
)

# Print PCA results
print(PCA)
# Find the IDs with minimum and maximum values for PC1
min_PC1_idx <- which.min(PCA$x[, 1])
max_PC1_idx <- which.max(PCA$x[, 1])

# Find the IDs with minimum and maximum values for PC2
min_PC2_idx <- which.min(PCA$x[, 2])
max_PC2_idx <- which.max(PCA$x[, 2])

# Get corresponding IDs
min_PC1_ID <- samples$ID[min_PC1_idx]
max_PC1_ID <- samples$ID[max_PC1_idx]

min_PC2_ID <- samples$ID[min_PC2_idx]
max_PC2_ID <- samples$ID[max_PC2_idx]

# Calculate the mean of PC1 scores
mean_PC1 <- mean(PCA$x[, 1])

# Find the ID closest to the mean PC1 score
mean_PC1_idx <- which.min(abs(PCA$x[, 1] - mean_PC1))
mean_PC1_ID <- samples$ID[mean_PC1_idx]

# Print the min, max, and mean values and corresponding IDs for PC1
cat("PC1 Minimum, Maximum, and Mean values:\n")
cat("Min PC1: ", min(PCA$x[, 1]), " (ID: ", min_PC1_ID, ")\n")
cat("Max PC1: ", max(PCA$x[, 1]), " (ID: ", max_PC1_ID, ")\n")
cat("Mean PC1: ", mean_PC1, " (ID: ", mean_PC1_ID, ")\n\n")

# Print the min and max values and corresponding IDs for PC2
cat("PC2 Minimum and Maximum values:\n")
cat("Min PC2: ", min(PCA$x[, 2]), " (ID: ", min_PC2_ID, ")\n")
cat("Max PC2: ", max(PCA$x[, 2]), " (ID: ", max_PC2_ID, ")\n\n")




### IN ANAYSE ARE INCORPORATE ELEMENTS FROM SLOVENIA AND SPAIN ###

### Principal Component Analysis with confidence ellipse ####

# Install and load necessary packages
install.packages("ggplot2")
library(ggplot2)
library(geomorph)

tps_file <- "skupno - Copy.TPS"
landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)
sliders <- define.sliders(3:138) 

# Assuming `landmarks` is already defined and loaded
samples <- data.frame(name = rep(NA, dim(landmarks)[3]),
                      country = rep(NA, dim(landmarks)[3]),
                      variable1 = rep(NA, dim(landmarks)[3]),
                      variable2 = rep(NA, dim(landmarks)[3]),
                      variable3 = rep(NA, dim(landmarks)[3]))

for (i in 1:dim(landmarks)[3]) {
  samples$name[i] <- dimnames(landmarks)[[3]][i]
  samples$country[i] <- unlist(strsplit(samples$name[i], "_"))[1]
  samples$variable1[i] <- unlist(strsplit(samples$name[i], "_"))[2]
  samples$variable2[i] <- unlist(strsplit(samples$name[i], "_"))[3]
  samples$variable3[i] <- unlist(strsplit(samples$name[i], "_"))[4]
}

samples$country <- as.factor(samples$country)
samples$variable1 <- as.factor(samples$variable1)

# Perform Procrustes analysis
sliders <- define.sliders(3:138)
landmarks.gpa <- gpagen(landmarks, curves = sliders)
plot(landmarks.gpa)
plotAllSpecimens(A = landmarks.gpa$coords, mean = TRUE, label = FALSE)

# Perform PCA
PCA <- gm.prcomp(landmarks.gpa$coords)

# Calculate the percentage of variance for each PC
variance_explained <- PCA$sdev^2 / sum(PCA$sdev^2) * 100
pc1_var <- round(variance_explained[1], 2)
pc2_var <- round(variance_explained[2], 2)

# Generate a distinct color and symbol for each country
unique_countries <- levels(samples$country)
num_countries <- length(unique_countries)
colors <- rainbow(num_countries)
symbols <- 15 + 0:(num_countries - 1)

# Map colors and symbols to the data
pca_data$Color <- colors[as.numeric(pca_data$Country)]
pca_data$Symbol <- symbols[as.numeric(pca_data$Country)]

# Create PCA data for plotting
pca_data <- data.frame(
  PC1 = PCA$x[, 1],
  PC2 = PCA$x[, 2],
  Country = samples$country
)

# Map colors and symbols based on the country
pca_data$Color <- colors[as.numeric(pca_data$Country)]
pca_data$Symbol <- symbols[as.numeric(pca_data$Country)]

# Define country full names
country_full_names <- c(
  "SL" = "Slovenia",
  "SP" = "Spain"
)

# Plot with ggplot2
ggplot(pca_data, aes(x = PC1, y = PC2, color = Country, shape = Country)) +
  geom_point(size = 2) +
  scale_shape_manual(values = symbols, labels = country_full_names) +  # Update with full names
  scale_color_manual(values = c(
    "SL" = "grey",  # Grey for Slovenia
    "SP" = "#0033A0" # Blue for Spain
  ), labels = country_full_names) +  # Update with full names
  stat_ellipse(type = "t", level = 0.95, segments = 51, na.rm = FALSE) +
  labs(title = "Morphospace of the aboral side of P. murcianus",
       x = paste0("PC 1 (", pc1_var, "% variance)"),
       y = paste0("PC 2 (", pc2_var, "% variance)")) +
  theme_minimal() +
  theme(
    legend.position = "bottom",  # Place legend at the bottom
    legend.box = "horizontal",   # Arrange legend horizontally
    plot.title = element_text(size = 20, face = "bold"),  # Bigger plot title
    axis.title = element_text(size = 16),  # Bigger axis titles
    axis.text = element_text(size = 14),   # Bigger axis text
    legend.text = element_text(size = 14), # Bigger legend text
    legend.title = element_text(size = 16) # Bigger legend title
  )


### Print the ID od specific specimens who are closer to the PC1 and PC2 min and max ####

# Create a sample data frame assuming you have 'ID' or 'IMAGE' in your dataset
samples <- data.frame(
  ID = dimnames(landmarks)[[3]],  # Replace this with appropriate identifier if different
  PC1 = PCA$x[, 1],
  PC2 = PCA$x[, 2]
)

# Print PCA results
print(PCA)
# Find the IDs with minimum and maximum values for PC1
min_PC1_idx <- which.min(PCA$x[, 1])
max_PC1_idx <- which.max(PCA$x[, 1])

# Find the IDs with minimum and maximum values for PC2
min_PC2_idx <- which.min(PCA$x[, 2])
max_PC2_idx <- which.max(PCA$x[, 2])

# Get corresponding IDs
min_PC1_ID <- samples$ID[min_PC1_idx]
max_PC1_ID <- samples$ID[max_PC1_idx]

min_PC2_ID <- samples$ID[min_PC2_idx]
max_PC2_ID <- samples$ID[max_PC2_idx]

# Calculate the mean of PC1 scores
mean_PC1 <- mean(PCA$x[, 1])

# Find the ID closest to the mean PC1 score
mean_PC1_idx <- which.min(abs(PCA$x[, 1] - mean_PC1))
mean_PC1_ID <- samples$ID[mean_PC1_idx]

# Print the min, max, and mean values and corresponding IDs for PC1
cat("PC1 Minimum, Maximum, and Mean values:\n")
cat("Min PC1: ", min(PCA$x[, 1]), " (ID: ", min_PC1_ID, ")\n")
cat("Max PC1: ", max(PCA$x[, 1]), " (ID: ", max_PC1_ID, ")\n")
cat("Mean PC1: ", mean_PC1, " (ID: ", mean_PC1_ID, ")\n\n")

# Print the min and max values and corresponding IDs for PC2
cat("PC2 Minimum and Maximum values:\n")
cat("Min PC2: ", min(PCA$x[, 2]), " (ID: ", min_PC2_ID, ")\n")
cat("Max PC2: ", max(PCA$x[, 2]), " (ID: ", max_PC2_ID, ")\n\n")











### bigger the specimes bigger the dot ###







### GROUPING DATA IN TO WESTERN (SPAIN) AND EASTERN TETHYS (SLOVENIA AND BOSNIA AND HERZEGOVINA) ###
  

### Length distribution ###

### Eastern Tetys ###

require(geomorph)
install.packages("ggplot2")
install.packages("dplyr")
install.packages("tidyverse")
library(ggplot2)
library(geomorph)
library(dplyr) # dplyr don't work now, and I don't know why


tps_file <- "eastern.TPS"
landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)
sliders <- define.sliders(3:138)


# Function to calculate Euclidean distance
calculate_distance <- function(x1, y1, x2, y2) {
  return(sqrt((x2 - x1)^2 + (y2 - y1)^2))
}

# Function to process the data and calculate scaled distance
process_data <- function(file_path) {
  lines <- readLines(file_path)
  lm_indices <- which(grepl("LM=", lines))
  scale_indices <- which(grepl("SCALE=", lines))
  
  if(length(lm_indices) == 0 || length(scale_indices) == 0) {
    stop("LM= or SCALE= not found in the data file.")
  }
  
  results <- list()
  
  for (i in seq_along(lm_indices)) {
    num_landmarks <- as.numeric(gsub("LM=", "", lines[lm_indices[i]]))
    scale <- as.numeric(gsub("SCALE=", "", lines[scale_indices[i]]))
    
    landmark_lines <- lines[(lm_indices[i] + 1):(lm_indices[i] + num_landmarks)]
    landmarks <- do.call(rbind, strsplit(landmark_lines, "\\s+"))
    landmarks <- as.data.frame(landmarks, stringsAsFactors = FALSE)
    landmarks <- mutate_all(landmarks, as.numeric)
    
    if(nrow(landmarks) < 2) {
      next 
    }
    
    # Calculate distances between all pairs of landmarks
    num_landmarks <- nrow(landmarks)
    distances <- numeric()
    for (j in 1:(num_landmarks - 1)) {
      for (k in (j + 1):num_landmarks) {
        dist <- calculate_distance(landmarks[j, 1], landmarks[j, 2], landmarks[k, 1], landmarks[k, 2])
        scaled_distance <- dist * scale
        distances <- c(distances, scaled_distance)
      }
    }
    
    results[[length(results) + 1]] <- list(
      block_index = i,
      scale = scale,
      distances = distances
    )
  }
  
  return(results)
}

# Load the TPS file using geomorph's readland.tps function
tps_file <- "eastern.TPS"
landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)
sliders <- define.sliders(3:138)

# Process the TPS file to get scaled distances
scaled_distances_results <- tryCatch({
  process_data(tps_file)
}, error = function(e) {
  cat("Error:", e$message, "\n")
  return(NULL)
})

if (!is.null(scaled_distances_results)) {
  
  # Convert scaled distances results into a data frame
  distances_data <- do.call(rbind, lapply(seq_along(scaled_distances_results), function(i) {
    result <- scaled_distances_results[[i]]
    data.frame(
      block_index = i,
      scale = result$scale,
      distance = unlist(result$distances)
    )
  }))
  
  # Add the data to the samples data frame
  samples <- data.frame(
    name = dimnames(landmarks)[[3]],
    region = "Eastern Tethys",  # All samples are from Eastern Tethys
    variable1 = NA,
    variable2 = NA,
    variable3 = NA
  )
  
  for (i in 1:nrow(samples)) {
    split_name <- unlist(strsplit(samples$name[i], "_"))
    samples$variable1[i] <- split_name[2]
    samples$variable2[i] <- split_name[3]
    samples$variable3[i] <- split_name[4]
  }
  
  samples$variable1 <- factor(samples$variable1)
  
  # Merge using valid index, assuming 'block_index' corresponds to sample rows
  samples_with_distances <- merge(samples, distances_data, by.x = "row.names", by.y = "block_index", all.x = TRUE)
  colnames(samples_with_distances)[1] <- "specimen_id"  # Rename for clarity
  
  # Function to compute Shapiro-Wilk test results for combined data
  compute_shapiro_results <- function(data) {
    shapiro_test <- shapiro.test(data$distance)
    return(paste0("W = ", round(shapiro_test$statistic, 5), 
                  "\nP-value = ", format(shapiro_test$p.value, digits = 3, scientific = TRUE)))
  }
  
  # Compute Shapiro-Wilk results for the combined data
  shapiro_results <- compute_shapiro_results(samples_with_distances)
  
  # Create the combined histogram plot for Eastern Tethys
  plot <- ggplot(samples_with_distances, aes(x = distance)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    stat_function(
      fun = dnorm, 
      args = list(mean = mean(samples_with_distances$distance, na.rm = TRUE), 
                  sd = sd(samples_with_distances$distance, na.rm = TRUE)),
      color = "red", 
      linewidth = 1
    ) +
    labs(
      title = "Histogram of Length for Eastern Tethys with Gaussian Curve",
      x = "Length", 
      y = "Density"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 22, face = "bold"),   # Increase title size
      axis.title.x = element_text(size = 18),                # Increase X-axis label size
      axis.title.y = element_text(size = 18),                # Increase Y-axis label size
      axis.text.x = element_text(size = 16),                 # Increase X-axis tick size
      axis.text.y = element_text(size = 16),                 # Increase Y-axis tick size
      legend.title = element_text(size = 16),                # Increase legend title size (if applicable)
      legend.text = element_text(size = 14)                  # Increase legend text size (if applicable)
    ) +
    annotate(
      "text", 
      x = Inf, 
      y = Inf, 
      label = shapiro_results,
      hjust = 1.1, 
      vjust = 1.5, 
      size = 6,       # Increase size of the Shapiro-Wilk annotation text
      color = "blue"
    )
  
  # Print the plot
  print(plot)
  
} else {
  cat("No scaled distances were computed.\n")
}




### by sections ###
# Load necessary libraries
library(geomorph)
library(ggplot2)
library(dplyr)

# Function to calculate Euclidean distance
calculate_distance <- function(x1, y1, x2, y2) {
  return(sqrt((x2 - x1)^2 + (y2 - y1)^2))
}

# Function to process the data and calculate scaled distance
process_data <- function(file_path) {
  lines <- readLines(file_path)
  lm_indices <- which(grepl("LM=", lines))
  scale_indices <- which(grepl("SCALE=", lines))
  
  if(length(lm_indices) == 0 || length(scale_indices) == 0) {
    stop("LM= or SCALE= not found in the data file.")
  }
  
  results <- list()
  
  for (i in seq_along(lm_indices)) {
    num_landmarks <- as.numeric(gsub("LM=", "", lines[lm_indices[i]]))
    scale <- as.numeric(gsub("SCALE=", "", lines[scale_indices[i]]))
    
    landmark_lines <- lines[(lm_indices[i] + 1):(lm_indices[i] + num_landmarks)]
    landmarks <- do.call(rbind, strsplit(landmark_lines, "\\s+"))
    landmarks <- as.data.frame(landmarks, stringsAsFactors = FALSE)
    landmarks <- mutate_all(landmarks, as.numeric)
    
    if(nrow(landmarks) < 2) {
      next 
    }
    
    # Calculate distances between all pairs of landmarks
    num_landmarks <- nrow(landmarks)
    distances <- numeric()
    for (j in 1:(num_landmarks - 1)) {
      for (k in (j + 1):num_landmarks) {
        dist <- calculate_distance(landmarks[j, 1], landmarks[j, 2], landmarks[k, 1], landmarks[k, 2])
        scaled_distance <- dist * scale
        distances <- c(distances, scaled_distance)
      }
    }
    
    results[[length(results) + 1]] <- list(
      block_index = i,
      scale = scale,
      distances = distances
    )
  }
  
  return(results)
}

# Load the TPS file
tps_file <- "eastern.TPS"
landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)
sliders <- define.sliders(3:138)

# Process the TPS file to get scaled distances
scaled_distances_results <- tryCatch({
  process_data(tps_file)
}, error = function(e) {
  cat("Error:", e$message, "\n")
  return(NULL)
})

if (!is.null(scaled_distances_results)) {
  
  # Convert scaled distances results into a data frame
  distances_data <- do.call(rbind, lapply(scaled_distances_results, function(result) {
    data.frame(
      block_index = result$block_index,
      scale = result$scale,
      distance = unlist(result$distances)
    )
  }))
  
  # Add the data to the samples data frame
  samples <- data.frame(
    name = dimnames(landmarks)[[3]],
    country = NA,
    variable1 = NA,
    variable2 = NA,
    variable3 = NA
  )
  
  for (i in 1:nrow(samples)) {
    split_name <- unlist(strsplit(samples$name[i], "_"))
    samples$country[i] <- split_name[1]
    samples$variable1[i] <- split_name[2]
    samples$variable2[i] <- split_name[3]
    samples$variable3[i] <- split_name[4]
  }
  
  samples$country <- as.factor(samples$country)
  samples$variable1 <- factor(samples$variable1)
  
  # Add distances to samples data frame
  samples_with_distances <- merge(samples, distances_data, by.x = "row.names", by.y = "block_index")
  colnames(samples_with_distances)[1] <- "specimen_id"
  
  # Function to compute Shapiro-Wilk test results for each level of variable1
  compute_shapiro_results <- function(data) {
    shapiro_test <- shapiro.test(data$distance)
    return(paste0("W = ", round(shapiro_test$statistic, 5), 
                  "\nP-value = ", format(shapiro_test$p.value, digits = 3, scientific = TRUE)))
  }
  
  # Create a list to store the plots
  plot_list <- list()
  
  # Loop through each level of variable1 and create a plot
  for (level in levels(samples_with_distances$variable1)) {
    # Filter data for the current level
    subset_data <- samples_with_distances %>% filter(variable1 == level)
    
    # Compute Shapiro-Wilk results for the subset
    shapiro_results <- compute_shapiro_results(subset_data)
    
    # Create the histogram plot
    plot <- ggplot(subset_data, aes(x = distance)) +
      geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
      stat_function(
        fun = dnorm, 
        args = list(mean = mean(subset_data$distance), sd = sd(subset_data$distance)),
        color = "black", 
        linewidth = 1
      ) +
      labs(
        title = paste("Histogram of Length for", level, "with Gaussian Curve"),
        x = "Length", 
        y = "Density"
      ) +
      theme_minimal() +
      annotate(
        "text", 
        x = Inf, 
        y = Inf, 
        label = shapiro_results,
        hjust = 1.1, 
        vjust = 1.5, 
        size = 5, 
        color = "blue"
      )
    
    # Add the plot to the list
    plot_list[[level]] <- plot
  }
  
  # Print all plots
  for (p in plot_list) {
    print(p)
  }
  
} else {
  cat("No scaled distances were computed.\n")
}



# Recode 'variable1' to replace abbreviations with full names
samples_with_distances$variable1 <- recode(samples_with_distances$variable1,
                                           "Clp" = "Calaspara",
                                           "Li" = "Libros",
                                           "Bu" = "Bugarra",
                                           "Pr" = "Prikrnica",
                                           "He" = "Henarejos")

# Loop through each level of variable1 and create a plot with larger labels and text
for (level in levels(samples_with_distances$variable1)) {
  # Filter data for the current level
  subset_data <- samples_with_distances %>% filter(variable1 == level)
  
  # Compute Shapiro-Wilk results for the subset
  shapiro_results <- compute_shapiro_results(subset_data)
  
  # Create the histogram plot with larger labels and text
  plot <- ggplot(subset_data, aes(x = distance)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    stat_function(
      fun = dnorm, 
      args = list(mean = mean(subset_data$distance), sd = sd(subset_data$distance)),
      color = "black", 
      linewidth = 1
    ) +
    labs(
      title = paste(level),
      x = "Length", 
      y = "Density"
    ) +
    theme_minimal() +
    annotate(
      "text", 
      x = Inf, 
      y = Inf, 
      label = shapiro_results,
      hjust = 1.1, 
      vjust = 1.5, 
      size = 5, 
      color = "blue"
    ) +
    theme(
      plot.title = element_text(size = 20, face = "bold"),  # Title text size
      axis.title.x = element_text(size = 16),               # X-axis label size
      axis.title.y = element_text(size = 16),               # Y-axis label size
      axis.text.x = element_text(size = 14),                # X-axis tick size
      axis.text.y = element_text(size = 14),                # Y-axis tick size
      legend.title = element_text(size = 16),               # Legend title size
      legend.text = element_text(size = 14)                 # Legend text size
    )
  
  # Add the plot to the list
  plot_list[[level]] <- plot
}

# Print all plots
for (p in plot_list) {
  print(p)
}


### Length by sections (nonparametric test) ###
# Install and load required packages
if (!require("ca")) {
  install.packages("ca")
}
if (!require("ggplot2")) {
  install.packages("ggplot2")
}
if (!require("dplyr")) {
  install.packages("dplyr")
}
library(ca)
library(ggplot2)
library(dplyr)

tps_file <- "skupno - Copy.TPS"
landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)
sliders <- define.sliders(3:138)

# Placeholder function to calculate length from landmarks
calculate_length <- function(landmark_data) {
  coords <- as.matrix(landmark_data)
  sqrt(sum((coords[1, ] - coords[nrow(coords), ])^2))
}

# Calculate length for each sample
lengths <- sapply(1:dim(landmarks)[3], function(i) calculate_length(landmarks[,,i]))

# Generate a sample properties data frame
samples <- data.frame(name = dimnames(landmarks)[[3]],
                      country = as.factor(sapply(dimnames(landmarks)[[3]], function(x) unlist(strsplit(x, "_"))[1])),
                      variable1 = as.factor(sapply(dimnames(landmarks)[[3]], function(x) unlist(strsplit(x, "_"))[2])),
                      variable2 = as.factor(sapply(dimnames(landmarks)[[3]], function(x) unlist(strsplit(x, "_"))[3])),
                      variable3 = as.factor(sapply(dimnames(landmarks)[[3]], function(x) unlist(strsplit(x, "_"))[4])))

# Add length variable to the samples data frame
samples$length <- lengths

# Perform the Kruskal-Wallis test
kruskal_test <- kruskal.test(length ~ variable1, data = samples)
print(kruskal_test)

# Visualize the results with boxplots
ggplot(samples, aes(x = variable1, y = length, fill = variable1)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) +
  labs(title = "Length Distribution by Section",
       x = "Section",
       y = "Length") +
  scale_fill_manual(values = c("Pr" = "red", "Li" = "blue", "Bu" = "green", "He" = "orange"),
                    labels = section_names) +
  theme_minimal() +
  theme(legend.position = "none")  # Remove legend if desired


### Length by regions ###

# Install and load required packages
if (!require("ca")) {
  install.packages("ca")
}
if (!require("ggplot2")) {
  install.packages("ggplot2")
}
library(ca)
library(ggplot2)

# Reading landmarks data from TPS file
tps_file <- "skupno - Copy.TPS"
landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)
sliders <- define.sliders(3:138)

# Placeholder function to calculate length from landmarks
# Replace this with your actual method to compute lengths
calculate_length <- function(landmark_data) {
  # Example: Calculate Euclidean distance between first and last landmark points
  coords <- as.matrix(landmark_data)  # Convert landmark data to matrix
  sqrt(sum((coords[1, ] - coords[nrow(coords), ])^2))  # Euclidean distance
}

# Calculate length for each sample
lengths <- sapply(1:dim(landmarks)[3], function(i) calculate_length(landmarks[,,i]))

# Generate a sample properties data frame
samples <- data.frame(name = dimnames(landmarks)[[3]],
                      country = as.factor(sapply(dimnames(landmarks)[[3]], function(x) unlist(strsplit(x, "_"))[1])),
                      variable1 = as.factor(sapply(dimnames(landmarks)[[3]], function(x) unlist(strsplit(x, "_"))[2])),
                      variable2 = as.factor(sapply(dimnames(landmarks)[[3]], function(x) unlist(strsplit(x, "_"))[3])),
                      variable3 = as.factor(sapply(dimnames(landmarks)[[3]], function(x) unlist(strsplit(x, "_"))[4])))

# Add length variable to the samples data frame
samples$length <- lengths

# Group variable1 into Western and Eastern Tethys
samples$tethys_group <- ifelse(samples$variable1 %in% c("Bu", "Li", "He", "Clp"), "Western Tethys", "Eastern Tethys")

# Reorder levels of tethys_group so that Eastern Tethys comes first
samples$tethys_group <- factor(samples$tethys_group, levels = c("Western Tethys", "Eastern Tethys"))

# Fit a linear model using length as the response variable and tethys_group as the only predictor
model <- lm(length ~ tethys_group, data = samples)

# Perform ANOVA with 1 degree of freedom for the Tethys group comparison
anova_results <- anova(model)
print(anova_results)

# Visualize the results with a box-and-whiskers plot for Western and Eastern Tethys
ggplot(samples, aes(x = tethys_group, y = length, fill = tethys_group)) +
  geom_boxplot() +
  labs(title = "Length by Tethys Group",
       x = "Sephardic province",
       y = "Length") +
  scale_fill_manual(values = c("Western Tethys" = "lightblue", "Eastern Tethys" = "gray")) + 
  theme_minimal() +
  theme(legend.position = "none")






### PC1 scores ###

### Weastern Tetys (when analysing the Eastern Tethys just change the file to eastern) ###

# Load necessary libraries
library(geomorph)
library(ggplot2)
library(dplyr) ##same proble it doesn't work

# Load the TPS file for Eastern Tethys
tps_file <- "C:/Users/katja.oselj/Desktop/DOKTORSKA DISERTACIJA/MORFOMETRIČNE ANALIZE/REGIONS/Western T.TPS"
landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)
sliders <- define.sliders(3:138)

# Check if 'landmarks' is a 3D array and if it has dimnames
if (is.array(landmarks) && !is.null(dimnames(landmarks)[[3]])) {
  
  # Initialize the data frame with the correct number of rows
  samples <- data.frame(
    name = rep(NA, dim(landmarks)[3]),
    country = rep(NA, dim(landmarks)[3]),
    variable1 = rep(NA, dim(landmarks)[3]),
    variable2 = rep(NA, dim(landmarks)[3]),
    variable3 = rep(NA, dim(landmarks)[3])
  )
  
  # Populate the samples data frame
  for (i in 1:dim(landmarks)[3]) {
    samples$name[i] <- dimnames(landmarks)[[3]][i]
    
    # Split the name based on underscores
    split_name <- unlist(strsplit(samples$name[i], "_"))
    
    # Assign values based on the split
    samples$country[i] <- split_name[1]
    samples$variable1[i] <- split_name[2]
    samples$variable2[i] <- split_name[3]
    samples$variable3[i] <- split_name[4]
  }
  
  # Convert to factors
  samples$country <- as.factor(samples$country)
  samples$variable1 <- factor(samples$variable1)
  
  # Perform Generalized Procrustes Analysis (GPA) to align shapes
  gpa_result <- gpagen(landmarks, curves = sliders, PrinAxes = TRUE)
  
  # Perform PCA on the aligned shape data
  pca_result <- gm.prcomp(gpa_result$coords)
  
  # Extract PC1 scores
  samples$PC1 <- pca_result$x[,1]
  
  # Compute Shapiro-Wilk test for the overall PC1 distribution
  shapiro_test <- shapiro.test(samples$PC1)
  shapiro_results <- paste0("W = ", round(shapiro_test$statistic, 5), 
                            "\nP-value = ", format(shapiro_test$p.value, digits = 3, scientific = TRUE))
  
  # Create a combined histogram plot for Eastern Tethys with Gaussian curve
  combined_plot <- ggplot(samples, aes(x = PC1)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    stat_function(
      fun = dnorm, 
      args = list(mean = mean(samples$PC1), sd = sd(samples$PC1)),
      color = "black", 
      linewidth = 1
    ) +
    labs(
      title = "Western part",
      x = "PC1 Scores", 
      y = "Density"
    ) +
    theme_minimal() +
    annotate(
      "text", 
      x = Inf, 
      y = Inf, 
      label = shapiro_results,
      hjust = 1.1, 
      vjust = 1.5, 
      size = 3, 
      color = "blue",
      fontface = "italic"
    )
  
  # Print the combined plot
  print(combined_plot)
  
} else {
  stop("The 'landmarks' object is either not a 3D array or does not have dimnames.")
}




# Create a combined histogram plot for Eastern Tethys with Gaussian curve
combined_plot <- ggplot(samples, aes(x = PC1)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
  stat_function(
    fun = dnorm, 
    args = list(mean = mean(samples$PC1), sd = sd(samples$PC1)),
    color = "black", 
    linewidth = 1
  ) +
  labs(
    title = "North-Eastern part",
    x = "PC1 Scores", 
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  ) +
  annotate(
    "text", 
    x = Inf, 
    y = Inf, 
    label = shapiro_results,
    hjust = 1.1, 
    vjust = 1.5, 
    size = 3, 
    color = "blue",
    fontface = "italic"
  )

# Print the combined plot
print(combined_plot)








### Sections ###
  
  # Function to compute Shapiro-Wilk test results for each level of variable1
  compute_shapiro_results <- function(data) {
    shapiro_test <- shapiro.test(data$PC1)
    return(paste0("W = ", round(shapiro_test$statistic, 5), 
                  "\nP-value = ", format(shapiro_test$p.value, digits = 3, scientific = TRUE)))
  }
  
  # Create a list to store the plots
  plot_list <- list()
  
  # Loop through each level of variable1 and create a plot
  for (level in levels(samples$variable1)) {
    # Filter data for the current level
    subset_data <- samples %>% filter(variable1 == level)
    
    # Compute Shapiro-Wilk results for the subset
    shapiro_results <- compute_shapiro_results(subset_data)
    
    # Create the histogram plot
    plot <- ggplot(subset_data, aes(x = PC1)) +
      geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
      stat_function(
        fun = dnorm, 
        args = list(mean = mean(subset_data$PC1), sd = sd(subset_data$PC1)),
        color = "red", 
        linewidth = 1
      ) +
      labs(
        title = paste(level,),
        x = "PC1 Scores", 
        y = "Density"
      ) +
      theme_minimal() +
      annotate(
        "text", 
        x = Inf, 
        y = Inf, 
        label = shapiro_results,
        hjust = 1.1, 
        vjust = 1.5, 
        size = 5, 
        color = "blue",
        fontface = "italic"
      )
    
    # Add the plot to the list
    plot_list[[level]] <- plot
  }
  
  # Print all plots
  for (p in plot_list) {
    print(p)
  }
  
} else {
  stop("The 'landmarks' object is either not a 3D array or does not have dimnames.")
}



# Define level mappings
level_names <- c(
  "Pr" = "Prikrnica",
  "Dr" = "Drežnica",
  "Bu" = "Bugarra",
  "Li" = "Libros",
  "He" = "Henarejos",
  "Clp" = "Calasparra"
)

# Function to compute Shapiro-Wilk test results for each level of variable1
compute_shapiro_results <- function(data) {
  shapiro_test <- shapiro.test(data$PC1)
  return(paste0("W = ", round(shapiro_test$statistic, 5), 
                "\nP-value = ", format(shapiro_test$p.value, digits = 3, scientific = TRUE)))
}

# Create a list to store the plots
plot_list <- list()

# Loop through each level of variable1 and create a plot
for (level in levels(samples$variable1)) {
  # Filter data for the current level
  subset_data <- samples %>% filter(variable1 == level)
  
  # Compute Shapiro-Wilk results for the subset
  shapiro_results <- compute_shapiro_results(subset_data)
  
  # Get the descriptive name for the level
  level_name <- level_names[level]
  
  # Create the histogram plot
  plot <- ggplot(subset_data, aes(x = PC1)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    stat_function(
      fun = dnorm, 
      args = list(mean = mean(subset_data$PC1), sd = sd(subset_data$PC1)),
      color = "black", 
      linewidth = 1
    ) +
    labs(
      title = level_name,
      x = "PC1 Scores", 
      y = "Density"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14)
    ) +
    annotate(
      "text", 
      x = Inf, 
      y = Inf, 
      label = shapiro_results,
      hjust = 1.1, 
      vjust = 1.5, 
      size = 4, 
      color = "blue",
      fontface = "italic"
    )
  
  # Add the plot to the list
  plot_list[[level]] <- plot
}

# Print all plots
for (p in plot_list) {
  print(p)
}







### PC1 by coutry ###

### ANOVA and Tukey's analyse for coutry ###
tps_file <- "skupno - Copy.TPS"
landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)
sliders <- define.sliders(3:138) 


# Check if 'pca_samples' contains 'country' and 'PC1' columns
if ("country" %in% colnames(pca_samples) && "PC1" %in% colnames(pca_samples)) {
  
  # Perform ANOVA for PC1 by country
  anova_results <- aov(PC1 ~ country, data = pca_samples)
  
  # Extract summary of ANOVA
  anova_summary <- summary(anova_results)
  
  # Extract p-value from the ANOVA summary
  p_value <- anova_summary[[1]]$`Pr(>F)`[1]
  
  # Print p-value
  cat("ANOVA p-value for PC1 across countries:", p_value, "\n")
  
} else {
  stop("Error: `pca_samples` should contain both `country` and `PC1` columns.")
}

# Perform Tukey's HSD post-hoc test
tukey_results <- TukeyHSD(anova_results)

# Print the Tukey's HSD results
print(tukey_results)


# Load ggplot2 for plotting
library(ggplot2)

# Create a boxplot of PC1 by country
ggplot(pca_samples, aes(x = country, y = PC1, fill = country)) +
  geom_boxplot() +
  labs(x = "Country", y = "PC1", title = "Distribution of PC1 by Country") +
  scale_color_manual(values = c(
    "Bosnia and Herzegovina" = "#D81B60",  # Dark pink
    "Slovenia" = "grey",                        # Grey
    "Spain" = "#0033A0"                         # Blue
  )) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability






### ANOVA and Tukey's analyse for section ###
# Load necessary libraries
library(ggplot2)
library(geomorph)

# Load the TPS file
tps_file <- "skupno - Copy.TPS"
landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)
sliders <- define.sliders(3:138)

# Check if 'landmarks' is a 3D array and if it has dimnames
if (is.array(landmarks) && !is.null(dimnames(landmarks)[[3]])) {
  
  # Initialize the data frame with the correct number of rows
  samples <- data.frame(
    name = rep(NA, dim(landmarks)[3]),
    country = rep(NA, dim(landmarks)[3]),
    variable1 = rep(NA, dim(landmarks)[3]),
    variable2 = rep(NA, dim(landmarks)[3]),
    variable3 = rep(NA, dim(landmarks)[3])
  )
  
  # Populate the samples data frame
  for (i in 1:dim(landmarks)[3]) {
    samples$name[i] <- dimnames(landmarks)[[3]][i]
    
    # Split the name based on underscores
    split_name <- unlist(strsplit(samples$name[i], "_"))
    
    # Assign values based on the split
    samples$country[i] <- split_name[1]
    samples$variable1[i] <- split_name[2]
    samples$variable2[i] <- split_name[3]
    samples$variable3[i] <- split_name[4]
  }
  
  # Convert to factors
  samples$country <- as.factor(samples$country)
  
  # Manually set the levels for variable1 to match the desired section order
  samples$variable1 <- factor(samples$variable1, levels = c("Clp", "He", "Li", "Bu", "Dr", "Pr"))
  
  # Perform Generalized Procrustes Analysis (GPA) to align shapes
  gpa_result <- gpagen(landmarks, curves = sliders, PrinAxes = TRUE)
  
  # Perform PCA on the aligned shape data
  pca_result <- gm.prcomp(gpa_result$coords)
  
  # Extract PC1 scores
  samples$PC1 <- pca_result$x[,1]
  
  # Plot PC1 scores by variable1 (sections) in the specified order
  ggplot(samples, aes(x = variable1, y = PC1, fill = variable1)) +
    geom_boxplot() +
    labs(x = "Section (Variable1)", y = "PC1 Scores", title = "PC1 Scores by Section (Ordered)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability
  
  # Fit an ANOVA model using PC1 as the dependent variable
  anova_model <- aov(PC1 ~ variable1, data = samples)
  
  # Check the ANOVA summary
  summary(anova_model)
  
  # Perform Tukey's HSD test
  tukey_results <- TukeyHSD(anova_model)
  
  # View Tukey's HSD results
  print(tukey_results)
  
  # Plot the Tukey's HSD results
  plot(tukey_results)
  
} else {
  stop("The 'landmarks' object is either not a 3D array or does not have dimnames.")
}


  
  
  # Load necessary libraries
  library(ggplot2)

# Verify and align PCA results and country labels
# Assuming PCA and samples are already defined

# Check dimensions
cat("Dimensions of PCA$x:", dim(PCA$x), "\n")
cat("Length of samples$country:", length(samples$country), "\n")

# Ensure samples$country matches the number of rows in PCA$x
if (nrow(PCA$x) != length(samples$country)) {
  stop("Error: Number of rows in PCA$x does not match length of samples$country")
}

# Create a data frame for PCA scores and country information
pca_scores_df <- data.frame(
  PC1 = PCA$x[, 1],  # Ensure PCA$x has the correct dimensions
  country = samples$country  # Ensure samples$country is aligned
)

# Define country full names
country_full_names <- c(
  "BAH" = "Bosnia and Herzegovina",
  "SL" = "Slovenia",
  "SP" = "Spain"
)

# Convert country codes in pca_scores_df$country to full names
pca_scores_df$country <- country_full_names[pca_scores_df$country]
# Reorder countries with Spain on the left
pca_scores_df$country <- factor(pca_scores_df$country, levels = c("Spain", "Bosnia and Herzegovina", "Slovenia"))
# Create the plot
ggplot(pca_scores_df, aes(x = country, y = PC1, fill = country)) +
  geom_boxplot() +
  labs(x = "Country", y = "PC1", title = "Distribution of PC1 by Country") +
  scale_fill_manual(values = c(
    "Bosnia and Herzegovina" = "lightgreen",  # Dark pink
    "Slovenia" = "lightblue",                 # Grey
    "Spain" = "grey"                  # Blue
  )) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Increase size of x-axis labels
    axis.text.y = element_text(size = 14),  # Increase size of y-axis labels
    axis.title.x = element_text(size = 16),  # Increase size of x-axis title
    axis.title.y = element_text(size = 16),  # Increase size of y-axis title
    plot.title = element_text(size = 18, face = "bold"),  # Increase size of plot title
    legend.title = element_blank(),  # Remove legend title
    legend.text = element_text(size = 14)  # Increase size of legend text
  )
  



### N-E and W part ###


# Load necessary libraries
library(ggplot2)
library(geomorph)

# Load the TPS file
tps_file <- "C:/Users/katja.oselj/Desktop/DOKTORSKA DISERTACIJA/TPS files - prba/skupno - Copy.TPS"
landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)
sliders <- define.sliders(3:138)

# Check if 'landmarks' is a 3D array and if it has dimnames
if (is.array(landmarks) && !is.null(dimnames(landmarks)[[3]])) {
  
  # Initialize the data frame with the correct number of rows
  samples <- data.frame(
    name = rep(NA, dim(landmarks)[3]),
    country = rep(NA, dim(landmarks)[3]),
    variable1 = rep(NA, dim(landmarks)[3]),
    variable2 = rep(NA, dim(landmarks)[3]),
    variable3 = rep(NA, dim(landmarks)[3])
  )
  
  # Populate the samples data frame
  for (i in 1:dim(landmarks)[3]) {
    samples$name[i] <- dimnames(landmarks)[[3]][i]
    
    # Split the name based on underscores
    split_name <- unlist(strsplit(samples$name[i], "_"))
    
    # Assign values based on the split
    samples$country[i] <- split_name[1]
    samples$variable1[i] <- split_name[2]
    samples$variable2[i] <- split_name[3]
    samples$variable3[i] <- split_name[4]
  }
  
  # Convert to factors
  samples$country <- as.factor(samples$country)
  
  # Define Western and North-Eastern sections
  samples$section <- ifelse(samples$variable1 %in% c("DR", "Pr"), "North-Eastern", "Western")
  
  # Convert 'section' to a factor with Western on the left and North-Eastern on the right
  samples$section <- factor(samples$section, levels = c("Western", "North-Eastern"))
  
  # Perform Generalized Procrustes Analysis (GPA) to align shapes
  gpa_result <- gpagen(landmarks, curves = sliders, PrinAxes = TRUE)
  
  # Perform PCA on the aligned shape data
  pca_result <- gm.prcomp(gpa_result$coords)
  
  # Extract PC1 scores
  samples$PC1 <- pca_result$x[,1]
  
  # Plot PC1 scores by section
  ggplot(samples, aes(x = section, y = PC1, fill = section)) +
    geom_boxplot() +
    labs(x = "Section", y = "PC1 Scores", title = "PC1 Scores by Section (Western vs. North-Eastern)") +
    scale_fill_manual(values = c("North-Eastern part" = "lightblue", "Western part" = "gray")) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Increase size of x-axis labels
      axis.text.y = element_text(size = 14),  # Increase size of y-axis labels
      axis.title.x = element_text(size = 16),  # Increase size of x-axis title
      axis.title.y = element_text(size = 16),  # Increase size of y-axis title
      plot.title = element_text(size = 18, face = "bold"),  # Increase size of plot title
      legend.title = element_blank(),  # Remove legend title
      legend.text = element_text(size = 14)  # Increase size of legend text
    )
  
  # Fit an ANOVA model using PC1 as the dependent variable
  anova_model <- aov(PC1 ~ section, data = samples)
  
  # Check the ANOVA summary
  summary(anova_model)
  
  # Perform Tukey's HSD test
  tukey_results <- TukeyHSD(anova_model)
  
  # View Tukey's HSD results
  print(tukey_results)
  
  # Plot the Tukey's HSD results
  plot(tukey_results)
  
} else {
  stop("The 'landmarks' object is either not a 3D array or does not have dimnames.")
}

# Load necessary libraries
library(geomorph)

# Load the TPS file
tps_file <- "C:/Users/katja.oselj/Desktop/DOKTORSKA DISERTACIJA/MORFOMETRIČNE ANALIZE/REGIONS/Western T.TPS"
landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)

# Perform Procrustes alignment (this step is often required before calculating the mean shape)
gpa <- gpagen(landmarks)

# Calculate the mean shape
msho <- mshape(gpa$coords)

# Plot the mean shape
plot(msho, main = "Mean Shape")





### Allometric relationship between size and shape ###

### Eastern Tethys ###
# Load the necessary libraries
library(geomorph)
library(dplyr)
library(ggplot2)

# Path to your data file
file_path <- "C:/Users/katja.oselj/Desktop/DOKTORSKA DISERTACIJA/TPS files - prba/skupno - Copy.TPS"
tps_file <- file_path

# Load landmarks and define sliders
landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)
sliders <- define.sliders(3:138)

# Function to calculate Euclidean distance
calculate_distance <- function(x1, y1, x2, y2) {
  return(sqrt((x2 - x1)^2 + (y2 - y1)^2))
}

# Function to process the data and calculate scaled distance
process_data <- function(file_path) {
  lines <- readLines(file_path)
  lm_indices <- which(grepl("LM=", lines))
  scale_indices <- which(grepl("SCALE=", lines))
  
  if(length(lm_indices) == 0 || length(scale_indices) == 0) {
    stop("LM= or SCALE= not found in the data file.")
  }
  
  results <- list()
  
  for (i in seq_along(lm_indices)) {
    num_landmarks <- as.numeric(gsub("LM=", "", lines[lm_indices[i]]))
    scale <- as.numeric(gsub("SCALE=", "", lines[scale_indices[i]]))
    
    landmark_lines <- lines[(lm_indices[i] + 1):(lm_indices[i] + num_landmarks)]
    landmarks <- do.call(rbind, strsplit(landmark_lines, "\\s+"))
    landmarks <- as.data.frame(landmarks, stringsAsFactors = FALSE)
    landmarks <- mutate_all(landmarks, as.numeric)
    
    if(nrow(landmarks) < 2) {
      next 
    }
    
    distance <- calculate_distance(landmarks[1, 1], landmarks[1, 2], landmarks[2, 1], landmarks[2, 2])
    scaled_distance <- distance * scale
    
    results[[length(results) + 1]] <- list(
      block_index = i,
      scale = scale,
      scaled_distance = scaled_distance
    )
  }
  
  return(results)
}

# Calculate the scaled distances
tryCatch({
  scaled_distances <- process_data(file_path)
  
  if (length(scaled_distances) > 0) {
    total_distance <- sum(sapply(scaled_distances, function(x) x$scaled_distance))
    average_distance <- total_distance / length(scaled_distances)
    
    for (result in scaled_distances) {
      cat("Block Index:", result$block_index, "\n")
      cat("Scale:", result$scale, "\n")
      cat("Scaled Distance:", result$scaled_distance, "\n\n")
    }
    
    cat("Average Scaled Distance:", average_distance, "\n")
  } else {
    cat("No valid scaled distances were computed.\n")
  }
}, error = function(e) {
  cat("Error:", e$message, "\n")
})

# Perform GPA and PCA
gpa_results <- gpagen(landmarks, curves = sliders, ProcD = FALSE)
pca_results <- gm.prcomp(gpa_results$coords)
pc1_scores <- pca_results$x[,1]

# Create and populate the samples data frame
samples <- data.frame(name = rep(NA, dim(landmarks)[3]),
                      country = rep(NA, dim(landmarks)[3]),
                      variable1 = rep(NA, dim(landmarks)[3]),
                      Section = rep(NA, dim(landmarks)[3]),
                      variable3 = rep(NA, dim(landmarks)[3]))

for (i in 1:dim(landmarks)[3]) {
  samples$name[i] <- dimnames(landmarks)[[3]][i]
  samples$country[i] <- unlist(strsplit(samples$name[i], "_"))[1]
  samples$variable1[i] <- unlist(strsplit(samples$name[i], "_"))[2]
  samples$Section[i] <- unlist(strsplit(samples$name[i], "_"))[3]
  samples$variable3[i] <- unlist(strsplit(samples$name[i], "_"))[4]
}

# Convert columns to factors
samples$country <- as.factor(samples$country)
samples$variable1 <- as.factor(samples$variable1)
samples$Section <- factor(samples$Section)

# Add PC1 scores to the samples data frame
samples$PC1 <- pc1_scores

# Ensure that scaled_distances has the correct number of entries
if (length(scaled_distances) != dim(landmarks)[3]) {
  stop("Mismatch between the number of scaled distances and samples.")
}

# Assign scaled distances to samples
samples$Length <- sapply(1:dim(landmarks)[3], function(i) {
  if (i <= length(scaled_distances)) {
    return(scaled_distances[[i]]$scaled_distance)
  } else {
    return(NA)
  }
})

# Calculate the correlation and p-value
correlation_test <- cor.test(samples$PC1, samples$Length, method = "pearson")
correlation_value <- correlation_test$estimate
p_value <- correlation_test$p.value
library(ggplot2)
library(ggplot2)



samples$variable_grouped <- ifelse(samples$variable1 %in% c("He", "Li", "BU"), "Combined", samples$variable1)

ggplot(samples, aes(x = PC1, y = Length, color = variable_grouped)) +
  geom_point() +  # Points colored by variable1
  geom_smooth(method = "lm", se = FALSE, aes(group = variable1), color = "black") +  # Add separate regression lines for each section
  theme_minimal() +
  labs(title = "PC1 vs Length",
       x = "PC1 Score",
       y = "Length") +
  scale_color_discrete(name = "Section",
                       labels = c("Bu" = "Bugarra", "He" = "Henarejos", "Li" = "Libros", "Pr" = "Prikrnica")) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )








library(ggplot2)

# Plot with varying dot sizes
ggplot(samples, aes(x = PC1, y = Length)) +
  geom_point(aes(color = variable1, size = Length)) +  # Points colored by variable1 and sized by Length
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add linear regression line
  theme_minimal() +
  labs(title = "PC1 vs Length",
       x = "PC1 Score",
       y = "Length") +
  scale_color_discrete(name = "Section",
                       labels = c("Dr" = "Drežnica", "Pr" = "Prikrnica")) +
  scale_size_continuous(name = "Length", range = c(1, 5)) +  # Adjust the size range as needed
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )


# Perform linear regression
lm_model <- lm(Length ~ PC1, data = samples)

# Print the summary of the regression model
summary(lm_model)

# Extract and print the slope (inclination) of the regression line
slope <- coef(lm_model)[2]
cat("Inclination (Slope) of the regression line:", slope, "\n")





### Eastern Tethys ###
# Load the necessary libraries
library(geomorph)
library(dplyr)
library(ggplot2)

# Path to your data file
file_path <- "C:/Users/katja.oselj/Desktop/DOKTORSKA DISERTACIJA/MORFOMETRIČNE ANALIZE/REGIONS/eastern.TPS"
tps_file <- file_path

# Load landmarks and define sliders
landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)
sliders <- define.sliders(3:138)

# Function to calculate Euclidean distance
calculate_distance <- function(x1, y1, x2, y2) {
  return(sqrt((x2 - x1)^2 + (y2 - y1)^2))
}

# Function to process the data and calculate scaled distance
process_data <- function(file_path) {
  lines <- readLines(file_path)
  lm_indices <- which(grepl("LM=", lines))
  scale_indices <- which(grepl("SCALE=", lines))
  
  if(length(lm_indices) == 0 || length(scale_indices) == 0) {
    stop("LM= or SCALE= not found in the data file.")
  }
  
  results <- list()
  
  for (i in seq_along(lm_indices)) {
    num_landmarks <- as.numeric(gsub("LM=", "", lines[lm_indices[i]]))
    scale <- as.numeric(gsub("SCALE=", "", lines[scale_indices[i]]))
    
    landmark_lines <- lines[(lm_indices[i] + 1):(lm_indices[i] + num_landmarks)]
    landmarks <- do.call(rbind, strsplit(landmark_lines, "\\s+"))
    landmarks <- as.data.frame(landmarks, stringsAsFactors = FALSE)
    landmarks <- mutate_all(landmarks, as.numeric)
    
    if(nrow(landmarks) < 2) {
      next 
    }
    
    distance <- calculate_distance(landmarks[1, 1], landmarks[1, 2], landmarks[2, 1], landmarks[2, 2])
    scaled_distance <- distance * scale
    
    results[[length(results) + 1]] <- list(
      block_index = i,
      scale = scale,
      scaled_distance = scaled_distance
    )
  }
  
  return(results)
}

# Calculate the scaled distances
tryCatch({
  scaled_distances <- process_data(file_path)
  
  if (length(scaled_distances) > 0) {
    total_distance <- sum(sapply(scaled_distances, function(x) x$scaled_distance))
    average_distance <- total_distance / length(scaled_distances)
    
    for (result in scaled_distances) {
      cat("Block Index:", result$block_index, "\n")
      cat("Scale:", result$scale, "\n")
      cat("Scaled Distance:", result$scaled_distance, "\n\n")
    }
    
    cat("Average Scaled Distance:", average_distance, "\n")
  } else {
    cat("No valid scaled distances were computed.\n")
  }
}, error = function(e) {
  cat("Error:", e$message, "\n")
})

# Perform GPA and PCA
gpa_results <- gpagen(landmarks, curves = sliders, ProcD = FALSE)
pca_results <- gm.prcomp(gpa_results$coords)
pc1_scores <- pca_results$x[,1]

# Create and populate the samples data frame
samples <- data.frame(name = rep(NA, dim(landmarks)[3]),
                      country = rep(NA, dim(landmarks)[3]),
                      variable1 = rep(NA, dim(landmarks)[3]),
                      Section = rep(NA, dim(landmarks)[3]),
                      variable3 = rep(NA, dim(landmarks)[3]))

for (i in 1:dim(landmarks)[3]) {
  samples$name[i] <- dimnames(landmarks)[[3]][i]
  samples$country[i] <- unlist(strsplit(samples$name[i], "_"))[1]
  samples$variable1[i] <- unlist(strsplit(samples$name[i], "_"))[2]
  samples$Section[i] <- unlist(strsplit(samples$name[i], "_"))[3]
  samples$variable3[i] <- unlist(strsplit(samples$name[i], "_"))[4]
}

# Convert columns to factors
samples$country <- as.factor(samples$country)
samples$variable1 <- as.factor(samples$variable1)
samples$Section <- factor(samples$Section)

# Add PC1 scores to the samples data frame
samples$PC1 <- pc1_scores

# Ensure that scaled_distances has the correct number of entries
if (length(scaled_distances) != dim(landmarks)[3]) {
  stop("Mismatch between the number of scaled distances and samples.")
}

# Assign scaled distances to samples
samples$Length <- sapply(1:dim(landmarks)[3], function(i) {
  if (i <= length(scaled_distances)) {
    return(scaled_distances[[i]]$scaled_distance)
  } else {
    return(NA)
  }
})

# Calculate the correlation and p-value
correlation_test <- cor.test(samples$PC1, samples$Length, method = "pearson")
correlation_value <- correlation_test$estimate
p_value <- correlation_test$p.value



ggplot(samples, aes(x = PC1, y = Length)) +
  geom_point(aes(color = variable1)) +  # Points colored by variable1
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add linear regression line
  theme_minimal() +
  labs(title = "PC1 vs Length",
       x = "PC1 Score",
       y = "Length") +
  scale_color_discrete(name = "Section",
                       labels = c("He" = "Henarejos", "Li" = "Libros", "Bu" = "Bugarra", "Clp" = "Calasparra")) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )


ggplot(samples, aes(x = PC1, y = Length)) +
  geom_point(aes(color = variable1, size = Length)) +  # Points colored by variable1 and sized by Length
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add linear regression line
  theme_minimal() +
  labs(title = "PC1 vs Length",
       x = "PC1 Score",
       y = "Length") +
  scale_color_discrete(name = "Section",
                       labels = c("He" = "Henarejos", "Li" = "Libros", "Bu" = "Bugarra", "Clp" = "Calasparra")) +
  scale_size_continuous(name = "Length", range = c(1, 5)) +  # Adjust the size range as needed
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )


# Perform linear regression
lm_model <- lm(Length ~ PC1, data = samples)

# Print the summary of the regression model
summary(lm_model)

# Extract and print the slope (inclination) of the regression line
slope <- coef(lm_model)[2]
cat("Inclination (Slope) of the regression line:", slope, "\n")






# Add PC2 scores to the samples data frame
pc2_scores <- pca_results$x[,2]  # Extract PC2 scores
samples$PC2 <- pc2_scores         # Add PC2 scores to the samples data frame

# Plot PC2 vs Length with updated labels and legend
ggplot(samples, aes(x = PC2, y = Length)) +
  geom_point(aes(color = variable1)) +  # Points colored by variable1
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add linear regression line
  theme_minimal() +
  labs(title = "PC2 vs Length",     # Updated title to PC2
       x = "PC2 Score",             # X-axis is now PC2
       y = "Length") +
  scale_color_discrete(name = "Section",
                       labels = c("Dr" = "Drežnica", "Pr" = "Prikrnica")) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )















### mshape ###
# Load the required library
library(geomorph)

# Load TPS file and read landmarks
tps_file <- "C:/Users/katja.oselj/Desktop/DOKTORSKA DISERTACIJA/TPS files - prba/skupno - Copy.TPS"
landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)

# Define sliders for the landmarks
sliders <- define.sliders(3:138)

# Perform Procrustes analysis (Generalized Procrustes Analysis - GPA)
landmarks.gpa <- gpagen(landmarks, curves = sliders)

# Perform PCA on the aligned coordinates
PCA <- gm.prcomp(landmarks.gpa$coords)

# Define variable1 based on the sample names (as in your previous code)
samples <- data.frame(name = dimnames(landmarks)[[3]],
                      country = rep(NA, dim(landmarks)[3]),
                      variable1 = rep(NA, dim(landmarks)[3]),
                      variable2 = rep(NA, dim(landmarks)[3]),
                      variable3 = rep(NA, dim(landmarks)[3]))

for (i in 1:dim(landmarks)[3]) {
  samples$name[i] <- dimnames(landmarks)[[3]][i]
  
  # Split the name based on underscores
  split_name <- unlist(strsplit(samples$name[i], "_"))
  
  # Assign values based on the split
  samples$country[i] <- split_name[1]
  samples$variable1[i] <- split_name[2]
  samples$variable2[i] <- split_name[3]
  samples$variable3[i] <- split_name[4]
}

samples$variable1 <- factor(samples$variable1, levels = c("He", "Li", "Bu", "Clp", "Dr", "Pr"))

# Initialize a list to store mean shapes for each variable1 level
mean_shapes <- list()

# Loop over each level of variable1
for (level in levels(samples$variable1)) {
  # Subset the indices for the current level of variable1
  subset_indices <- which(samples$variable1 == level)
  
  # Subset the GPA coordinates and PCA results for the current level
  subset_coords <- landmarks.gpa$coords[, , subset_indices]
  subset_PCA <- gm.prcomp(subset_coords)
  
  # Compute the mean shape for this subset based on PC1
  mean_shape <- mshape(subset_coords)
  
  # Store the mean shape in the list
  mean_shapes[[level]] <- mean_shape
  
  # Plot the mean shape
  plot(mean_shape, main = paste("Mean Shape for", level, "based on PC1"), col = "blue", pch = 16, 
       xlab = "X Coordinate", ylab = "Y Coordinate")
  
  # Add text annotation to the plot
  text(x = mean_shape[1, 1], y = mean_shape[1, 2], labels = level, pos = 4, cex = 1.2, col = "red")
}

# Optionally, you can save these plots to files
for (level in names(mean_shapes)) {
  png_filename <- paste0("mean_shape_", level, "_PC1.png")
  png(filename = png_filename, width = 800, height = 600)
  
  # Plot the mean shape and save
  plot(mean_shapes[[level]], main = paste("Mean Shape for", level, "based on PC1"), col = "blue", pch = 16, 
       xlab = "X Coordinate", ylab = "Y Coordinate")
  
  # Add text annotation to the plot
  text(x = mean_shapes[[level]][1, 1], y = mean_shapes[[level]][1, 2], labels = level, pos = 4, cex = 1.2, col = "red")
  
  dev.off()
}




# Check if 'landmarks' is a 3D array and if it has dimnames
if (is.array(landmarks) && !is.null(dimnames(landmarks)[[3]])) {
  
  # Initialize the data frame with the correct number of rows
  samples <- data.frame(name = rep(NA, dim(landmarks)[3]),
                        country = rep(NA, dim(landmarks)[3]),
                        variable1 = rep(NA, dim(landmarks)[3]),
                        variable2 = rep(NA, dim(landmarks)[3]),
                        variable3 = rep(NA, dim(landmarks)[3]))
  
  for (i in 1:dim(landmarks)[3]) {
    samples$name[i] <- dimnames(landmarks)[[3]][i]
    
    # Split the name based on underscores
    split_name <- unlist(strsplit(samples$name[i], "_"))
    
    # Assign values based on the split
    samples$country[i] <- split_name[1]
    samples$variable1[i] <- split_name[2]
    samples$variable2[i] <- split_name[3]
    samples$variable3[i] <- split_name[4]
  }
  
  # Convert to factors and reorder variable1
  samples$country <- as.factor(samples$country)
  samples$variable1 <- factor(samples$variable1, levels = c("He", "Li", "Bu", "Pr", "Dr"))
  
} else {
  stop("Error: 'landmarks' should be a 3D array with appropriate dimnames.")
}

# Ensure pca_data contains 'name' column for merging
if (!"name" %in% colnames(pca_data)) {
  stop("Error: `pca_data` should contain a 'name' column for merging.")
}

# Merge `pca_data` with `samples` data frame
pca_samples <- merge(samples, pca_data, by = "name")

# Load ggplot2 for plotting
library(ggplot2)

# Get the list of unique countries
countries <- unique(pca_samples$country)

# Initialize an empty list to store plots
plot_list <- list()

# Generate a plot for each country
for (country in countries) {
  country_data <- subset(pca_samples, country == country)
  
  # Fit a linear model
  lm_model <- lm(PC1 ~ variable1, data = country_data)
  
  # Print the summary of the linear model
  cat("\nSummary of regression model for", country, ":\n")
  print(summary(lm_model))
  
  # Create plot for the current country
  plot <- ggplot(country_data, aes(x = variable1, y = PC1, fill = variable1, group = variable1)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, color = "black", alpha = 0.5) +
    geom_smooth(method = "lm", aes(group = 1), color = "blue", se = FALSE) +  # Add regression line
    labs(x = "Section", y = "PC1", title = paste("PC1 Distribution in sections", country)) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate variable names for better readability
  
  # Add plot to the list
  plot_list[[country]] <- plot
}

# Print each plot
for (country in countries) {
  print(plot_list[[country]])
}




