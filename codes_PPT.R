### STATISTICAL DIFFERENCE IN SHAPE BETWEEN LEFT AND RIGHT ELEMENTS ###
library(geomorph)

# Load TPS file
tps_file <- "group_analise.TPS"
landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)

# Define sliding landmarks
sliders <- define.sliders(3:138)

# Perform Generalized Procrustes Analysis (GPA)
landmarks.gpa <- gpagen(landmarks, curves = sliders)
PCA <- gm.prcomp(landmarks.gpa$coords) 

# Compute Mean Shape for each specimen
mean_shapes <- apply(landmarks.gpa$coords, 3, function(x) mean(x)) 

# Extract specimen IDs
specimen_IDs <- dimnames(landmarks.gpa$coords)[[3]]

# Create a grouping factor (assumes IDs containing 'R' belong to one group)
group <- ifelse(grepl("R", specimen_IDs), "R", "No_R")

# Perform Kruskal-Wallis test
kruskal_test <- kruskal.test(mean_shapes ~ group)

# Print results
print(kruskal_test)


### STATISTICAL DIFFERENCE IN LENGTH BETWEENLEFT AND RIGHT ELEMENTS ###

library(dplyr)

# Function to calculate Euclidean distance
calculate_distance <- function(x1, y1, x2, y2) {
  return(sqrt((x2 - x1)^2 + (y2 - y1)^2))
}

# Function to process TPS file and compute mean Euclidean distances
process_data <- function(file_path, group_label) {
  lines <- readLines(file_path)
  lm_indices <- which(grepl("LM=", lines))
  scale_indices <- which(grepl("SCALE=", lines))
  id_indices <- which(grepl("ID=", lines))  # Extract specimen IDs
  
  if (length(lm_indices) == 0 || length(scale_indices) == 0 || length(id_indices) == 0) {
    stop("LM=, SCALE=, or ID= not found in the data file.")
  }
  
  results <- data.frame(ID = character(), Mean_Distance = numeric(), Group = character(), stringsAsFactors = FALSE)
  
  for (i in seq_along(lm_indices)) {
    num_landmarks <- as.numeric(gsub("LM=", "", lines[lm_indices[i]]))
    scale <- as.numeric(gsub("SCALE=", "", lines[scale_indices[i]]))
    
    # Extract specimen ID
    specimen_ID <- gsub("ID=", "", lines[id_indices[i]])
    
    landmark_lines <- lines[(lm_indices[i] + 1):(lm_indices[i] + num_landmarks)]
    landmarks <- do.call(rbind, strsplit(landmark_lines, "\\s+"))
    landmarks <- as.data.frame(landmarks, stringsAsFactors = FALSE)
    landmarks <- mutate_all(landmarks, as.numeric)
    
    if (nrow(landmarks) < 2) {
      next 
    }
    
    # Compute Euclidean distances and scale them
    num_landmarks <- nrow(landmarks)
    distances <- numeric()
    for (j in 1:(num_landmarks - 1)) {
      for (k in (j + 1):num_landmarks) {
        dist <- calculate_distance(landmarks[j, 1], landmarks[j, 2], landmarks[k, 1], landmarks[k, 2])
        scaled_distance <- dist * scale  # Use scaled distance
        distances <- c(distances, scaled_distance)
      }
    }
    
    # Compute mean scaled Euclidean distance for this specimen
    mean_dist <- mean(distances)
    
    # Store results with the specified group
    results <- rbind(results, data.frame(ID = specimen_ID, Mean_Distance = mean_dist, Group = group_label, stringsAsFactors = FALSE))
  }
  
  return(results)
}

# File paths for right and left TPS files
right_file_path <- "right.TPS"
left_file_path <- "Left.TPS"

# Process the right and left files
right_data <- process_data(right_file_path, "Right")
left_data <- process_data(left_file_path, "Left")

# Combine the data from both files into one data frame
combined_data <- rbind(right_data, left_data)

# Check the grouping
print(table(combined_data$Group))
# Compute the mean distance for each file (group)
mean_right <- mean(right_data$Mean_Distance)
mean_left <- mean(left_data$Mean_Distance)

# Print the mean distances for right and left files
cat("Mean Distance for Right file: ", mean_right, "\n")
cat("Mean Distance for Left file: ", mean_left, "\n")

# Ensure we have at least two groups before running Kruskal-Wallis test
if (length(unique(combined_data$Group)) < 2) {
  stop("Kruskal-Wallis test requires at least two groups.")
}

# Perform Kruskal-Wallis test using scaled mean distances
kruskal_test <- kruskal.test(Mean_Distance ~ Group, data = combined_data)

# Print test result
print(kruskal_test)


### IN ANAYSE ARE INCORPORATE ELEMENTS FROM SLOVENIA, SPAIN AND BOSNIA & HERZEGOVINA ###

library(geomorph)
tps_file <- "All_sections.TPS"
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

tps_file <- "All_sections.TPS"
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

# Install and load necessary packages
install.packages("ggplot2")
library(ggplot2)
library(geomorph)

tps_file <- "All_sections.TPS"
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
# Map abbreviations to full country names
samples$country <- factor(samples$country,
                          levels = c("SL", "SP", "BAH"),
                          labels = c("Slovenia", "Spain", "Bosnia and Herzegovina"))

# Perform Procrustes analysis (as already defined)
landmarks.gpa <- gpagen(landmarks, curves = sliders)
plot(landmarks.gpa)
plotAllSpecimens(A = landmarks.gpa$coords, mean = TRUE, label = FALSE)

# Perform PCA
PCA <- gm.prcomp(landmarks.gpa$coords)

# Extract PCA scores and create a data frame
pca_data <- data.frame(PC1 = PCA$x[, 1], PC2 = PCA$x[, 2], Country = samples$country)

# Define specific colors for each country
colors <- c("Slovenia" = "grey", "Spain" = "#0033A0", "Bosnia and Herzegovina" = "#D81B60")

# Define symbols for each country
symbols <- c("Slovenia" = 16, "Spain" = 17, "Bosnia and Herzegovina" = 18)

# Calculate the explained variance
eigenvalues <- PCA$sdev^2  # Variance is the square of standard deviation
total_variance <- sum(eigenvalues)  # Total variance
explained_variance <- (eigenvalues / total_variance) * 100  # Percentage of variance explained

# Extract percentage for PC1 and PC2
pc1_variance <- round(explained_variance[1], 2)  # PC1 variance
pc2_variance <- round(explained_variance[2], 2)  # PC2 variance

# Create axis labels with explained variance
x_label <- paste("PC 1 (", pc1_variance, "%)", sep = "")
y_label <- paste("PC 2 (", pc2_variance, "%)", sep = "")

# Plot with ggplot2
ggplot(pca_data, aes(x = PC1, y = PC2, color = Country, shape = Country)) +
  geom_point(size = 2) +
  scale_shape_manual(values = symbols) +
  scale_color_manual(values = colors) +
  stat_ellipse(type = "t", level = 0.95, segments = 51, na.rm = FALSE) +
  labs(title = "Morphospace of the aboral side of P. murcianus",
       x = x_label,
       y = y_label) +
  theme_minimal() +
  theme(legend.position = c(0.1, 0.1))

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








### PCA: GROUPING DATA IN TO WESTERN (SPAIN) AND Norteastern_sections TETHYS (SLOVENIA AND BOSNIA AND HERZEGOVINA) ###
  

library(ggplot2)
library(geomorph)

tps_file <- "All_sections.TPS"
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
# Map abbreviations to full country names
samples$country <- factor(samples$country,
                          levels = c("SL", "SP", "BAH"),
                          labels = c("Slovenia", "Spain", "Bosnia and Herzegovina"))

# Perform Procrustes analysis (as already defined)
landmarks.gpa <- gpagen(landmarks, curves = sliders)
plot(landmarks.gpa)
plotAllSpecimens(A = landmarks.gpa$coords, mean = TRUE, label = FALSE)

# Perform PCA
PCA <- gm.prcomp(landmarks.gpa$coords)

# Extract PCA scores and create a data frame
pca_data <- data.frame(PC1 = PCA$x[, 1], PC2 = PCA$x[, 2], Country = samples$country)
# Modify the samples data frame to include region
samples$region <- ifelse(samples$country %in% c("Slovenia", "Bosnia and Herzegovina"), 
                         "North-Norteastern_sections", "Western")
samples$region <- as.factor(samples$region)

# Create PCA data frame including the region
pca_data <- data.frame(PC1 = PCA$x[, 1], PC2 = PCA$x[, 2], 
                       Country = samples$country, Region = samples$region)

# Define colors for each region
region_colors <- c("North-Norteastern_sections" = "grey", "Western" = "#0033A0")


# Plot with ggplot2 using region as color and removing country from the legend
ggplot(pca_data, aes(x = PC1, y = PC2, color = Region)) +
  geom_point(aes(shape = Country), size = 2, show.legend = FALSE) +  # Hide Country from legend
  scale_shape_manual(values = symbols) +
  scale_color_manual(values = region_colors) +
  # Create a single ellipse for the North-Eastern part
  stat_ellipse(data = pca_data[pca_data$Region == "North-Norteastern_sections", ], 
               aes(color = Region), 
               type = "t", level = 0.95, segments = 51, na.rm = FALSE) +
  # Create individual ellipses for the Western region
  stat_ellipse(data = pca_data[pca_data$Region == "Western", ], 
               aes(color = Region), 
               type = "t", level = 0.95, segments = 51, na.rm = FALSE) +
  labs(title = "Morphospace of the aboral side of P. murcianus",
       x = x_label,
       y = y_label) +
  theme_minimal() +
  theme(
    legend.position = "right",         # Move legend to the right
    legend.margin = margin(0, 20, 0, 0)  # Add space between the plot and the legend
  ) +
  guides(shape = guide_legend(override.aes = list(size = 4)))  # Adjust shape legend




### longer the specimen bigger the symbol ###
# Install required packages if they are not installed
if (!require("geomorph")) install.packages("geomorph", dependencies = TRUE)
if (!require("ggplot2")) install.packages("ggplot2", dependencies = TRUE)

# Load necessary libraries
library(geomorph)
library(ggplot2)

# Function to calculate Euclidean distance
calculate_distance <- function(x1, y1, x2, y2) {
  return(sqrt((x2 - x1)^2 + (y2 - y1)^2))
}

# Function to process the data and calculate scaled distance (using base R)
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
    
    # Use apply to convert all columns to numeric
    landmarks <- as.data.frame(apply(landmarks, 2, as.numeric))
    
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

# Process the data and calculate scaled distances
tps_file <- "All_sections.TPS"  
landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)
sliders <- define.sliders(3:138) 

# Assuming the process_data function calculates scaled distances for each specimen
scaled_distances <- process_data(tps_file)

# Extract the mean scaled distance for each specimen and add to the samples data frame
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

# Map abbreviations to full country names
samples$country <- factor(samples$country,
                          levels = c("SL", "SP", "BAH"),
                          labels = c("Slovenia", "Spain", "Bosnia and Herzegovina"))

# Add scaled distance to the samples data frame
samples$scaled_distance <- sapply(scaled_distances, function(x) mean(x$distances))

# Add region data
samples$region <- ifelse(samples$country %in% c("Slovenia", "Bosnia and Herzegovina"), 
                         "North-Eastern", "Western")
samples$region <- as.factor(samples$region)

# Perform Procrustes analysis (as already defined)
landmarks.gpa <- gpagen(landmarks, curves = sliders)

# Perform PCA
PCA <- gm.prcomp(landmarks.gpa$coords)

# Extract PCA scores and create a data frame
pca_data <- data.frame(PC1 = PCA$x[, 1], PC2 = PCA$x[, 2], 
                       Country = samples$country, Region = samples$region,
                       Scaled_Distance = samples$scaled_distance)  # Add scaled distance

# Define colors for each region
region_colors <- c("North-Eastern" = "grey", "Western" = "#0033A0") 

# Plot with ggplot2 using region as color and removing country from the legend
ggplot(pca_data, aes(x = PC1, y = PC2, color = Region, size = Scaled_Distance)) +  # Size by scaled distance
  geom_point(aes(shape = Country), show.legend = FALSE) +  # Hide Country from legend
  scale_shape_manual(values = c(16, 17, 18)) +  # Replace with actual symbols if needed
  scale_color_manual(values = region_colors) +
  # Create a single ellipse for the North-Eastern part
  stat_ellipse(data = pca_data[pca_data$Region == "North-Eastern", ], 
               aes(color = Region), 
               type = "t", level = 0.95, segments = 51, na.rm = FALSE) +
  # Create individual ellipses for the Western region
  stat_ellipse(data = pca_data[pca_data$Region == "Western", "North-Eastern"], 
               aes(color = Region), 
               type = "t", level = 0.95, segments = 51, na.rm = FALSE) +
  labs(title = "Morphospace of the aboral side of P. murcianus",
       x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(
    legend.position = "right",         # Move legend to the right
    legend.margin = margin(0, 20, 0, 0)  # Add space between the plot and the legend
  ) +
  guides(shape = guide_legend(override.aes = list(size = 4)), 
         size = guide_legend(title = "Length (µm)"))  # Change size legend title to Length






### normal distribution of PC1 scores ###

### Weastern_sections part ###

# Load necessary libraries
library(geomorph)
library(ggplot2)
library(dplyr) ##same proble it doesn't work

# Load the TPS file for Norteastern_sections Tethys
tps_file <- "Western_sections.TPS"
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
  
  # Create a combined histogram plot for Norteastern_sections Tethys with Gaussian curve
  combined_plot <- ggplot(samples, aes(x = PC1)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    stat_function(
      fun = dnorm, 
      args = list(mean = mean(samples$PC1), sd = sd(samples$PC1)),
      color = "black", 
      linewidth = 1
    ) +
    labs(
      title = "Western part of Sephardic province",
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
      color = "black",
      fontface = "italic"
    )
  
  # Print the combined plot
  print(combined_plot)
  
} else {
  stop("The 'landmarks' object is either not a 3D array or does not have dimnames.")
}


# Create a combined histogram plot for Norteastern_sections Tethys with Gaussian curve
combined_plot <- ggplot(samples, aes(x = PC1)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
  stat_function(
    fun = dnorm, 
    args = list(mean = mean(samples$PC1), sd = sd(samples$PC1)),
    color = "black", 
    linewidth = 1
  ) +
  labs(
    title = "Western part of Sephardic province",
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
    color = "black",
    fontface = "italic"
  )

# Print the combined plot
print(combined_plot)

### Sections in western part ###

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
      color = "black", 
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
      color = "black",
      fontface = "italic"
    )
  
  # Add the plot to the list
  plot_list[[level]] <- plot
}

# Print all plots
for (p in plot_list) {
  print(p)
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
      color = "black",
      fontface = "italic"
    )
  
  # Add the plot to the list
  plot_list[[level]] <- plot
}

# Print all plots
for (p in plot_list) {
  print(p)
}





### North-Eastern part ###

# Load necessary libraries
library(geomorph)
library(ggplot2)
library(dplyr) ##same proble it doesn't work

# Load the TPS file for Norteastern_sections Tethys
tps_file <- "Norteastern_sections.TPS"
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
  
  # Create a combined histogram plot for Norteastern_sections Tethys with Gaussian curve
  combined_plot <- ggplot(samples, aes(x = PC1)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    stat_function(
      fun = dnorm, 
      args = list(mean = mean(samples$PC1), sd = sd(samples$PC1)),
      color = "black", 
      linewidth = 1
    ) +
    labs(
      title = "North-Eastern part of Sephardic province",
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
      color = "black",
      fontface = "italic"
    )
  
  # Print the combined plot
  print(combined_plot)
  
} else {
  stop("The 'landmarks' object is either not a 3D array or does not have dimnames.")
}

# Create a combined histogram plot for Norteastern_sections Tethys with Gaussian curve
combined_plot <- ggplot(samples, aes(x = PC1)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
  stat_function(
    fun = dnorm, 
    args = list(mean = mean(samples$PC1), sd = sd(samples$PC1)),
    color = "black", 
    linewidth = 1
  ) +
  labs(
    title = "North-Eastern part of Sephardic province",
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
    color = "black",
    fontface = "italic"
  )

# Print the combined plot
print(combined_plot)



### Sections in North-Eastern part ###
  
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
        color = "black", 
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
        color = "black",
        fontface = "italic"
      )
    
    # Add the plot to the list
    plot_list[[level]] <- plot
  }
  
  # Print all plots
  for (p in plot_list) {
    print(p)
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
      color = "black",
      fontface = "italic"
    )
  
  # Add the plot to the list
  plot_list[[level]] <- plot
}

# Print all plots
for (p in plot_list) {
  print(p)
}



### nonparametric test for section ###

### PC sections ###
# Load necessary libraries
library(geomorph)  # For GPA, PCA, and shape analysis
library(ggplot2)   # For plotting
library(MASS)      # For Kernel Density Estimation
library(viridis)   # For colorblind-friendly palettes
library(dplyr)     # For data manipulation

# 1. Read the TPS file and define the sliders
tps_file <- "All_sections.TPS"
landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)
sliders <- define.sliders(3:138) 

# 2. Create the samples data frame
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
# Load the necessary library for data manipulation
library(dplyr)

# 3. Perform GPA (Generalized Procrustes Analysis) to align the shape data
gpa <- gpagen(landmarks, curves = sliders, print.progress = FALSE) # GPA alignment

# 4. Perform Principal Component Analysis (PCA)
pca_result <- gm.prcomp(gpa$coords) # PCA on GPA aligned data

# 5. Extract the first Principal Component (PC1) and add it to the samples data frame
samples$PC1 <- pca_result$x[, 1] # First Principal Component

# 6. Reorder the levels of variable1 to ensure the desired order in the boxplot with full names
samples$variable1 <- factor(samples$variable1, 
                            levels = c("Clp", "He", "Li", "Bu", "Dr", "Pr"),
                            labels = c("Calasparra", "Henarejos", "Libros", "Bugarra", "Drežnica", "Prikrnica"))

# 7. Filter the samples data frame to include only the desired sections
filtered_samples <- samples %>% 
  filter(variable1 %in% c("Henarejos", "Libros", "Bugarra", "Prikrnica"))

# Define custom colors for sections
section_colors <- c(
  "Henarejos" = "#D4A017",   # Yellow for Henarejos
  "Libros" = "#17A909",      # Green for Libros
  "Bugarra" = "#17BCC4",     # Cyan/Teal for Bugarra
  "Prikrnica" = "#D91C93"    # Pink for Prikrnica
)

# 8. Create a boxplot for PC1 by sections (variable1) with ordered levels
ggplot(filtered_samples, aes(x = variable1, y = PC1, fill = variable1)) +
  geom_boxplot(alpha = 0.6) +         # Plot the boxplot with some transparency
  scale_fill_manual(values = section_colors, drop = FALSE) +  # Use drop = FALSE to keep unused factor levels
  labs(title = "Boxplot of PC1 by Sections", 
       x = "Section", 
       y = "PC1 scores") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),  # Remove legend title
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
    axis.text.y = element_text(size = 12),  # Increase y-axis text size
    axis.title.x = element_text(size = 14),  # Increase x-axis title size
    axis.title.y = element_text(size = 14),  # Increase y-axis title size
    plot.title = element_text(size = 16, face = "bold")  # Increase plot title size and make it bold
  ) 


# Load the FSA package for Dunn test
install.packages("FSA")
library(FSA)

# Perform Kruskal-Wallis test to compare PC1 values between sections (variable1)
kruskal_test <- kruskal.test(PC1 ~ variable1, data = samples)
print(kruskal_test)

# If Kruskal-Wallis is significant, perform pairwise comparisons using Dunn test
if (kruskal_test$p.value < 0.05) {
  dunn_test <- dunnTest(PC1 ~ variable1, data = samples, method = "bonferroni")
  print(dunn_test)
}



### nonparametric test for regions in Sephardic province ###

# Load necessary libraries
library(ggplot2)
library(geomorph)
library(dplyr)
library(FSA)  # For Dunn's test

# Read landmarks from TPS file
tps_file <- "All_sections.TPS"
landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)
sliders <- define.sliders(3:138)

# Initialize samples data frame
samples <- data.frame(name = rep(NA, dim(landmarks)[3]),
                      country = rep(NA, dim(landmarks)[3]),
                      variable1 = rep(NA, dim(landmarks)[3]),
                      variable2 = rep(NA, dim(landmarks)[3]),
                      variable3 = rep(NA, dim(landmarks)[3]))

# Populate samples data frame
for (i in 1:dim(landmarks)[3]) {
  samples$name[i] <- dimnames(landmarks)[[3]][i]
  samples$country[i] <- unlist(strsplit(samples$name[i], "_"))[1]
  samples$variable1[i] <- unlist(strsplit(samples$name[i], "_"))[2]
  samples$variable2[i] <- unlist(strsplit(samples$name[i], "_"))[3]
  samples$variable3[i] <- unlist(strsplit(samples$name[i], "_"))[4]
}

# Factorize country and variable1
samples$country <- as.factor(samples$country)
samples$variable1 <- as.factor(samples$variable1)

# Map country abbreviations to full names
samples$country <- factor(samples$country,
                          levels = c("SL", "SP", "BAH"),
                          labels = c("Slovenia", "Spain", "Bosnia and Herzegovina"))

# Perform Procrustes analysis
landmarks.gpa <- gpagen(landmarks, curves = sliders)

# Perform PCA
PCA <- gm.prcomp(landmarks.gpa$coords)

# Modify the samples data frame to include region
samples$region <- ifelse(samples$country %in% c("Slovenia", "Bosnia and Herzegovina"), 
                         "North-Eastern", "Western")
samples$region <- as.factor(samples$region)

# Create PCA data frame including the region
pca_data <- data.frame(PC1 = PCA$x[, 1], PC2 = PCA$x[, 2], 
                       Country = samples$country, Region = samples$region)

# 1. Conduct Kruskal-Wallis Test on PC1
kruskal_test <- kruskal.test(PC1 ~ Region, data = pca_data)
print(kruskal_test)

# 2. If significant, conduct Dunn's Test for pairwise comparisons
if (kruskal_test$p.value < 0.05) {
  dunn_test <- dunnTest(PC1 ~ Region, data = pca_data, method = "bonferroni")
  print(dunn_test)
} else {
  cat("No significant differences between regions.\n")
}


# Load necessary libraries
library(ggplot2)
library(geomorph)
library(dplyr)
library(FSA)  # For Dunn's test

# Read landmarks from TPS file
tps_file <- "All_sections.TPS"
landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)
sliders <- define.sliders(3:138)

# Initialize samples data frame
samples <- data.frame(name = rep(NA, dim(landmarks)[3]),
                      country = rep(NA, dim(landmarks)[3]),
                      variable1 = rep(NA, dim(landmarks)[3]),
                      variable2 = rep(NA, dim(landmarks)[3]),
                      variable3 = rep(NA, dim(landmarks)[3]))

# Populate samples data frame
for (i in 1:dim(landmarks)[3]) {
  samples$name[i] <- dimnames(landmarks)[[3]][i]
  samples$country[i] <- unlist(strsplit(samples$name[i], "_"))[1]
  samples$variable1[i] <- unlist(strsplit(samples$name[i], "_"))[2]
  samples$variable2[i] <- unlist(strsplit(samples$name[i], "_"))[3]
  samples$variable3[i] <- unlist(strsplit(samples$name[i], "_"))[4]
}

# Factorize country and variable1
samples$country <- as.factor(samples$country)
samples$variable1 <- as.factor(samples$variable1)

# Map country abbreviations to full names
samples$country <- factor(samples$country,
                          levels = c("SL", "SP", "BAH"),
                          labels = c("Slovenia", "Spain", "Bosnia and Herzegovina"))

# Perform Procrustes analysis
landmarks.gpa <- gpagen(landmarks, curves = sliders)

# Perform PCA
PCA <- gm.prcomp(landmarks.gpa$coords)

# Modify the samples data frame to include region
samples$region <- ifelse(samples$country %in% c("Slovenia", "Bosnia and Herzegovina"), 
                         "North-Eastern", "Western")
samples$region <- as.factor(samples$region)

# Set the order of regions: Western on the left, North-Eastern on the right
samples$region <- factor(samples$region, levels = c("Western", "North-Eastern"))
# Create PCA data frame including the region
pca_data <- data.frame(PC1 = PCA$x[, 1], PC2 = PCA$x[, 2], 
                       Country = samples$country, Region = samples$region)

# 1. Conduct Kruskal-Wallis Test on PC1
kruskal_test <- kruskal.test(PC1 ~ Region, data = pca_data)
print(kruskal_test)

# 2. If significant, conduct Dunn's Test for pairwise comparisons
if (kruskal_test$p.value < 0.05) {
  dunn_test <- dunnTest(PC1 ~ Region, data = pca_data, method = "bonferroni")
  print(dunn_test)
} else {
  cat("No significant differences between regions.\n")
}

# 3. Create a boxplot for PC1 by region
ggplot(pca_data, aes(x = Region, y = PC1, fill = Region)) +
  geom_boxplot(alpha = 0.6) +  # Add boxplot with some transparency
  scale_fill_manual(values = c("Western" = "lightblue", "North-Eastern" = "gray")) + 
  labs(title = "Boxplot of PC1 by part of Sephardic province", 
       x = "Part of Sephardic province", 
       y = "PC1") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),  # Remove legend title
    axis.text.x = element_text(size = 12),  # Increase x-axis text size
    axis.text.y = element_text(size = 12),  # Increase y-axis text size
    axis.title.x = element_text(size = 14),  # Increase x-axis title size
    axis.title.y = element_text(size = 14),  # Increase y-axis title size
    plot.title = element_text(size = 16, face = "bold")  # Increase plot title size and make it bold
  ) 






### plots for mean shape for western and Norteastern_sections part ###
# Load the TPS file
tps_file <- "Western_sections.TPS"
landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)

# Perform Procrustes alignment (this step is often required before calculating the mean shape)
gpa <- gpagen(landmarks)

# Calculate the mean shape
msho <- mshape(gpa$coords)

# Plot the mean shape
plot(msho, main = "Mean Shape")



# Load the TPS file
tps_file <- "Norteastern_sections.TPS"
landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)

# Perform Procrustes alignment (this step is often required before calculating the mean shape)
gpa <- gpagen(landmarks)

# Calculate the mean shape
msho <- mshape(gpa$coords)

# Plot the mean shape
plot(msho, main = "Mean Shape")





### Length distribution ###

### Norteastern_sections Tetys ###

require(geomorph)
install.packages("ggplot2")
install.packages("dplyr")
install.packages("tidyverse")
library(ggplot2)
library(geomorph)
library(dplyr) # dplyr don't work now, and I don't know why


tps_file <- "Norteastern_sections.TPS"
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
tps_file <- "Norteastern_sections.TPS"
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
    region = "North-Eastern part of Sephardic province",  # All samples are from Norteastern_sections Tethys
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
  
  # Create the combined histogram plot for Norteastern_sections Tethys
  plot <- ggplot(samples_with_distances, aes(x = distance)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    stat_function(
      fun = dnorm, 
      args = list(mean = mean(samples_with_distances$distance, na.rm = TRUE), 
                  sd = sd(samples_with_distances$distance, na.rm = TRUE)),
      color = "black", 
      linewidth = 1
    ) +
    labs(
      title = "North-Eastern part of Sephardic province",
      x = "Length (µm)", 
      y = "Density"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 18, face = "bold"),   # Increase title size
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
      size = 5,       # Increase size of the Shapiro-Wilk annotation text
      color = "black"
    )
  
  # Print the plot
  print(plot)
  
} else {
  cat("No scaled distances were computed.\n")
}






### western part ###


tps_file <- "Western_sections.TPS"
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
tps_file <- "Western_sections.TPS"
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
    region = "Western part of Sephardic Province",  # All samples are from Norteastern_sections Tethys
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
  
  # Create the combined histogram plot for Norteastern_sections Tethys
  plot <- ggplot(samples_with_distances, aes(x = distance)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    stat_function(
      fun = dnorm, 
      args = list(mean = mean(samples_with_distances$distance, na.rm = TRUE), 
                  sd = sd(samples_with_distances$distance, na.rm = TRUE)),
      color = "black", 
      linewidth = 1
    ) +
    labs(
      title = "Western part of Sephardic province",
      x = "Length (µm)", 
      y = "Density"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 18, face = "bold"),   # Increase title size
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
      size = 5,       # Increase size of the Shapiro-Wilk annotation text
      color = "black"
    )
  
  # Print the plot
  print(plot)
  
} else {
  cat("No scaled distances were computed.\n")
}



### Normal distribution of length for sections in NE part of sephardic province ###

# Load necessary libraries
library(geomorph)
library(ggplot2)
library(dplyr)
# Process the data and calculate scaled distances
tps_file <- "All_sections.TPS"  # Replace with the actual file path of your TPS file
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

# Load the TPS file
tps_file <- "Norteastern_sections.TPS"
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
        x = "Length (µm)", 
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
        color = "black"
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
                                           "Pr" = "Prikrnica", "Dr" = "Drežnica")

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
      x = "Length (µm)", 
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
      color = "black"
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




### Normal distribution of length for sections in W part of sephardic province ###

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
tps_file <- "Western_sections.TPS"
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
        x = "Length (µm)", 
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
        color = "black"
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
                                           "Bu" = "Bugarra",
                                           "Li" = "Libros",
                                           "Clp" = "Calasparra",
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
      x = "Length (µm)", 
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
      color = "black"
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




### Length by regions in Sephardis province and statistically significant ###
# Load required packages
install.packages("ggpubr")
library(dplyr)
library(FSA)
library(ggplot2)
library(ggpubr)
library(geomorph)
# Process the data and calculate scaled distances
tps_file <- "All_sections.TPS"  # Replace with the actual file path of your TPS file
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

# Group variable1 into Western and Norteastern_sections part
samples$Sephardic_group <- ifelse(samples$variable1 %in% c("Bu", "Li", "He", "Clp"), "Western part", "North-Eastern part")

# Reorder levels of tethys_group
samples$Sephardic_group <- factor(samples$Sephardic_group, levels = c("Western part", "North-Eastern part"))

# Perform Kruskal-Wallis test
kruskal_test <- kruskal.test(length ~ Sephardic_group, data = samples)
print(kruskal_test)

# Perform Dunn's test
dunn_results <- dunnTest(length ~ Sephardic_group, data = samples, method = "bonferroni")
print(dunn_results)

# Visualize the results with a box-and-whiskers plot
ggplot(samples, aes(x = Sephardic_group, y = length, fill = Sephardic_group)) +
  geom_boxplot() +
  labs(title = "Length by part of Sephardic province",
       x = "Part of Sephardic province",
       y = "Length (µm)") +
  scale_fill_manual(values = c("Western part" = "lightblue", "North-Eastern part" = "gray")) + 
  theme_minimal() +
  theme(
    legend.title = element_blank(),  # Remove legend title
    axis.text.x = element_text(size = 12),  # Increase x-axis text size
    axis.text.y = element_text(size = 12),  # Increase y-axis text size
    axis.title.x = element_text(size = 14),  # Increase x-axis title size
    axis.title.y = element_text(size = 14),  # Increase y-axis title size
    plot.title = element_text(size = 16, face = "bold")  # Increase plot title size and make it bold
  ) 





### length by sections and statistically significant ###

# Load necessary libraries
library(geomorph)
library(dplyr)
library(ggplot2)

# Function to calculate Euclidean distance
calculate_distance <- function(x1, y1, x2, y2) {
  return(sqrt((x2 - x1)^2 + (y2 - y1)^2))
}

# Function to calculate lengths from landmarks
calculate_length <- function(landmark_data) {
  coords <- as.matrix(landmark_data)
  return(sqrt(sum((coords[1, ] - coords[nrow(coords), ])^2)))
}

# Read TPS file for landmarks
tps_file <- "Selected_sections.TPS"  # Replace with your actual file path
landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)
sliders <- define.sliders(3:138)

# Print the dimensions of the landmarks to check the number of samples
cat("Dimensions of landmarks array:\n")
print(dim(landmarks))

# Calculate length for each sample
lengths <- sapply(1:dim(landmarks)[3], function(i) calculate_length(landmarks[,,i]))

# Print the length of the `lengths` vector
cat("Number of calculated lengths: ", length(lengths), "\n")

# Create samples data frame
samples <- data.frame(name = dimnames(landmarks)[[3]],
                      country = sapply(dimnames(landmarks)[[3]], function(x) unlist(strsplit(x, "_"))[1]),
                      variable1 = sapply(dimnames(landmarks)[[3]], function(x) unlist(strsplit(x, "_"))[2]),
                      variable2 = sapply(dimnames(landmarks)[[3]], function(x) unlist(strsplit(x, "_"))[3]),
                      variable3 = sapply(dimnames(landmarks)[[3]], function(x) unlist(strsplit(x, "_"))[4]),
                      stringsAsFactors = FALSE)

# Print the number of rows in the samples data frame
cat("Number of rows in samples: ", nrow(samples), "\n")

# Convert relevant columns to factors
samples$country <- as.factor(samples$country)
samples$variable1 <- as.factor(samples$variable1)

# Update the variable1 names
samples <- samples %>%
  mutate(variable1 = recode(variable1, 
                            "He" = "Henarejos", 
                            "Li" = "Libros", 
                            "Bu" = "Bugarra", 
                            "Pr" = "Prikrnica")) %>%
  filter(variable1 %in% c("Henarejos", "Libros", "Bugarra", "Prikrnica"))

# Ensure lengths match the samples, and filter accordingly
if (length(lengths) == nrow(samples)) {
  samples$length <- lengths
} else {
  warning("Mismatch between number of lengths and samples. Attempting to filter.")
  
  # Assuming dimnames(landmarks)[[3]] correspond to sample names
  valid_samples <- dimnames(landmarks)[[3]][1:length(lengths)]
  
  # Filter samples that have corresponding valid landmarks
  samples <- samples[samples$name %in% valid_samples, ]
  
  # Check if they still mismatch
  cat("Filtered number of rows in samples: ", nrow(samples), "\n")
  cat("Filtered number of lengths: ", length(lengths), "\n")
  
  # Only use the number of lengths that matches the filtered samples
  samples$length <- lengths[1:nrow(samples)]
}

# Set the order of levels for variable1
samples$variable1 <- factor(samples$variable1, 
                            levels = c("Henarejos", "Libros", "Bugarra", "Prikrnica"))

# Define custom colors for sections
section_colors <- c(
  "Henarejos" = "#D4A017",   # Yellow for Henarejos
  "Libros" = "#17A909",      # Green for Libros
  "Bugarra" = "#17BCC4",     # Cyan/Teal for Bugarra
  "Prikrnica" = "#D91C93"    # Pink for Prikrnica
)

# Box and Whisker Plot with custom colors for Length by Section
ggplot(samples, aes(x = variable1, y = length, fill = variable1)) +
  geom_boxplot(alpha = 0.7, outlier.size = 2, outlier.shape = 16) +
  scale_fill_manual(values = section_colors, drop = FALSE) +  # Use drop = FALSE to keep unused factor levels
  labs(title = "Boxplot of Length by Section",
       x = "Section",
       y = "Length (µm)",
       fill = "Sections") +  # Update the legend title to "Sections"
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),   # Title text size
    axis.title.x = element_text(size = 16),                # X-axis label size
    axis.title.y = element_text(size = 16),                # Y-axis label size
    axis.text.x = element_text(size = 14),                 # X-axis tick label size
    axis.text.y = element_text(size = 14),                 # Y-axis tick label size
    legend.title = element_text(size = 16),                # Legend title size
    legend.text = element_text(size = 14)                  # Legend text size
  )


# Load necessary libraries
library(dplyr)
library(FSA)  # For Dunn's test

# 1. Perform the Kruskal-Wallis test
kruskal_test <- kruskal.test(length ~ variable1, data = samples)

# Print the result of the Kruskal-Wallis test
print(kruskal_test)

# 2. Conduct Dunn's Post Hoc Test if the Kruskal-Wallis test is significant
if (kruskal_test$p.value < 0.05) {
  dunn_test <- dunnTest(length ~ variable1, data = samples, method = "bonferroni")
  
  # Print the results of Dunn's test
  print(dunn_test)
} else {
  print("No significant differences found in the Kruskal-Wallis test.")
}






### relationship between Length and PC1 scores ###
### All together ###
# Load the necessary libraries
library(geomorph)
library(dplyr)
library(ggplot2)

# Path to your data file
file_path <- "All_sections - Copy.TPS"
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


# Create the grouping variable with updated names
samples$group <- ifelse(samples$variable1 %in% c("Pr", "Dr"), "North-Eastern part", "Western part")

# Update the plot with specific colors for each group
ggplot(samples, aes(x = PC1, y = Length, color = group, shape = variable1)) +  # Color by group and shape by section
  geom_point(aes(size = Length)) +  # Points sized by Length
  geom_smooth(aes(group = group), method = "lm", se = FALSE) +  # Regression lines grouped by the new 'group' column
  theme_minimal() +
  labs(title = "Relationship between size and shape",
       x = "PC1 Score",
       y = "Length (µm)") +
  scale_color_manual(name = "Sephardic province",  # Legend title updated to 'Sephardic province'
                     values = c("North-Eastern part" = "black", "Western part" = "lightblue"),  # Assigning specific colors
                     labels = c("North-Eastern part", "Western part")) +  # Updated group labels
  scale_shape_manual(name = "Section",
                     values = c(16, 17, 18, 19, 15, 13),  # Assign different shapes to each section
                     labels = c("Pr" = "Prikrnica", "Dr" = "Drežnica", 
                                "Li" = "Libros", "Bu" = "Bugarra", 
                                "He" = "Henarejos", "Clp" = "Calasparra")) +  # Shape labels for sections
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  scale_size_continuous(range = c(1, 5))  # Adjust the range of point sizes as needed

# Subset the data for each group
north_Norteastern_sections <- subset(samples, group == "North-Eastern part")
western <- subset(samples, group == "Western part")

# Fit linear models for each group
lm_north_Norteastern_sections <- lm(Length ~ PC1, data = north_Norteastern_sections)
lm_western <- lm(Length ~ PC1, data = western)

# Extract slopes (coefficients for PC1)
slope_north_Norteastern_sections <- coef(lm_north_Norteastern_sections)["PC1"]
slope_western <- coef(lm_western)["PC1"]

# Print the slopes
print(paste("Slope for North-Eastern part:", slope_north_Norteastern_sections))
print(paste("Slope for Western part:", slope_western))





### working on residuals ###


#### PC1 and PC2 resid. plot ###
# Load necessary libraries
if (!require("geomorph")) install.packages("geomorph", dependencies = TRUE)
if (!require("ggplot2")) install.packages("ggplot2", dependencies = TRUE)

library(geomorph)
library(ggplot2)

# Function to calculate Euclidean distance
calculate_distance <- function(x1, y1, x2, y2) {
  return(sqrt((x2 - x1)^2 + (y1 - y2)^2))
}

# Function to process the data and calculate scaled distance (using base R)
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
    
    # Use apply to convert all columns to numeric
    landmarks <- as.data.frame(apply(landmarks, 2, as.numeric))
    
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

# Path to your data file
tps_file <- "All_sections.TPS"

# Process the data to calculate distances
scaled_distances_results <- process_data(tps_file)

# Conduct Generalized Procrustes Analysis (GPA)
landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)
gpa <- gpagen(landmarks)

# Perform Principal Component Analysis (PCA) on Procrustes-aligned data
pca_results <- gm.prcomp(gpa$coords)

# Extract PC scores (PC1 and PC2)
pc_scores <- pca_results$x[, 1:2]  # Getting PC1 and PC2

# Calculate the residuals
residuals <- array(0, dim = dim(gpa$coords))  # Initialize residual array

for (i in 1:dim(gpa$coords)[3]) {
  residuals[,,i] <- gpa$coords[,,i] - gpa$consensus  # Subtract consensus from each specimen
}

# Sum the squared residuals across landmarks and dimensions (for each specimen)
pc1_residuals <- apply(residuals, 3, function(x) sum(x[,1]^2))  # PC1 residuals
pc2_residuals <- apply(residuals, 3, function(x) sum(x[,2]^2))  # PC2 residuals

# Combine results in a data frame with the specimen names
residuals_data <- data.frame(
  Specimen = dimnames(landmarks)[[3]],
  PC1_Residuals = pc1_residuals,
  PC2_Residuals = pc2_residuals
)

# Create a sample properties data frame
samples <- data.frame(
  Specimen = dimnames(landmarks)[[3]],
  Country = sapply(dimnames(landmarks)[[3]], function(x) strsplit(x, "_")[[1]][1])  # Extract country from specimen name
)

# Map abbreviations to full country names
country_full_names <- c("BAH" = "Bosnia and Herzegovina", 
                        "SP" = "Spain", 
                        "SL" = "Slovenia")

# Create a new column with full country names
samples$FullCountry <- country_full_names[samples$Country]

# Assign groups
samples$Group <- ifelse(samples$Country %in% c("SL", "BAH"), "North-Eastern part", "Western part")

# Merge residuals data with samples data (to include country information)
merged_data <- merge(residuals_data, samples, by = "Specimen")

# Extract the scaled distances from the processed data
scaled_distances <- unlist(lapply(scaled_distances_results, function(x) x$distances))

# Add scaled distances to the merged data frame
merged_data$Scaled_Distance <- scaled_distances[1:nrow(merged_data)]  # Ensure sizes match

# Calculate the variance explained by PC1 and PC2
eigenvalues <- pca_results$sdev^2  # Eigenvalues are the squared singular values
total_variance <- sum(eigenvalues)  # Total variance

# Calculate the proportion of variance explained by PC1 and PC2
variance_explained <- eigenvalues / total_variance
pc1_variance <- variance_explained[1] * 100  # Percentage of variance for PC1
pc2_variance <- variance_explained[2] * 100  # Percentage of variance for PC2

# Calculate explained variance for residuals in PC1 and PC2
residuals_pc1_variance <- var(merged_data$PC1_Residuals)  # Variance of PC1 residuals
residuals_pc2_variance <- var(merged_data$PC2_Residuals)  # Variance of PC2 residuals

# Calculate explained variance percentages for residuals
explained_variance_pc1 <- (residuals_pc1_variance / total_variance) * 100  # Explained variance percentage for PC1
explained_variance_pc2 <- (residuals_pc2_variance / total_variance) * 100  # Explained variance percentage for PC2

# Print variance explained
cat("Variance explained by PC1: ", round(pc1_variance, 2), "%\n")
cat("Variance explained by PC2: ", round(pc2_variance, 2), "%\n")
cat("Explained variance by residuals for PC1: ", round(explained_variance_pc1, 2), "%\n")
cat("Explained variance by residuals for PC2: ", round(explained_variance_pc2, 2), "%\n")

# Plot the residuals for PC1 and PC2, coloring by group and adding confidence ellipses
ggplot(merged_data, aes(x = PC1_Residuals, y = PC2_Residuals, color = Group)) +
  geom_point(aes(size = Scaled_Distance), alpha = 0.7) +  # Use scaled distance to size points
  stat_ellipse(level = 0.95) +  # Add 95% confidence ellipses
  labs(title = "PC1 vs PC2 Residuals with Confidence Ellipses", 
       x = "PC1 Residuals",  # Label for x-axis
       y = "PC2 Residuals") +  # Label for y-axis
  theme_minimal() +  # Use a clean theme
  theme(legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8)) +
  scale_color_manual(name = "Country Group", 
                     values = c("North-Eastern part" = "grey", "Western part" = "#0033A0")) +  # Define colors for groups
  scale_size_continuous(name = "Scaled Distance", range = c(1, 4)) +  # Adjust size range for visibility
  annotate("text", x = max(merged_data$PC1_Residuals) * 0.95, 
           y = max(merged_data$PC2_Residuals) * 0.5, 
           label = paste("Variance explained by PC1:", round(pc1_variance, 2), "%"), 
           size = 4, color = "black", hjust = 0) +  # Annotation for PC1 variance
  annotate("text", x = max(merged_data$PC1_Residuals) * 0.05, 
           y = max(merged_data$PC2_Residuals) * 0.95, 
           label = paste("Variance explained by PC2:", round(pc2_variance, 2), "%"), 
           size = 4, color = "black", vjust = 0) +  # Annotation for PC
  
  print(
    ggplot(merged_data, aes(x = PC1_Residuals, y = PC2_Residuals, color = Group, shape = Group)) +  # Map shape to Group
      geom_point(aes(size = Scaled_Distance)) +  # Fully opaque points
      stat_ellipse(level = 0.95) +
      labs(title = "PC1 vs PC2 Residuals with Confidence Ellipses", 
           x = "PC1 Residuals", 
           y = "PC2 Residuals") +
      theme_minimal() +
      theme(legend.title = element_text(size = 10), 
            legend.text = element_text(size = 8)) +
      scale_color_manual(name = "Country Group", 
                         values = c("North-Eastern part" = "grey", "Western part" = "#0033A0")) +
      scale_shape_manual(name = "Country Group", 
                         values = c("North-Eastern part" = 16, "Western part" = 17)) +  # Circles for North-Eastern, Triangles for Western
      scale_size_continuous(name = "Scaled Distance", range = c(1, 4)) +
      
      # PC1 annotation moved next to the x-axis label
      annotate("text", x = max(merged_data$PC1_Residuals) * 0.95, 
               y = max(merged_data$PC2_Residuals) * -0.03, 
               label = paste0(round(pc1_variance, 0), "%"),  
               size = 4, color = "black", hjust = 0) +
      
      # PC2 annotation moved next to the y-axis label
      annotate("text", x = max(merged_data$PC1_Residuals) * -0.03, 
               y = max(merged_data$PC2_Residuals) * 0.95, 
               label = paste0(round(pc2_variance, 0), "%"),  
               size = 4, color = "black", vjust = 0)
  )



### PC1 resids. by section###
### PC sections ###
# Load necessary libraries
library(geomorph)  # For GPA, PCA, and shape analysis
library(ggplot2)   # For plotting
library(MASS)      # For Kernel Density Estimation
library(viridis)   # For colorblind-friendly palettes
library(dplyr)     # For data manipulation
library(FSA)       # For Dunn test

# 1. Read the TPS file and define the sliders
tps_file <- "All_sections.TPS"
landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)
sliders <- define.sliders(3:138) 

# 2. Create the samples data frame
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

# 3. Perform GPA (Generalized Procrustes Analysis) to align the shape data
gpa <- gpagen(landmarks, curves = sliders, print.progress = FALSE) # GPA alignment

# 4. Perform Principal Component Analysis (PCA)
pca_result <- gm.prcomp(gpa$coords) # PCA on GPA aligned data

# 5. Calculate residuals for PC1
# Assume you have a linear model for PC1, e.g., PC1 ~ variable1
lm_model <- lm(pca_result$x[, 1] ~ samples$variable1) # Linear model for PC1
samples$PC1_residuals <- residuals(lm_model)  # Calculate residuals for PC1

# 6. Reorder the levels of variable1 to ensure the desired order in the boxplot with full names
samples$variable1 <- factor(samples$variable1, 
                            levels = c("Clp", "He", "Li", "Bu", "Dr", "Pr"),
                            labels = c("Calasparra", "Henarejos", "Libros", "Bugarra", "Drežnica", "Prikrnica"))

# 7. Filter the samples data frame to include only the desired sections
filtered_samples <- samples %>% 
  filter(variable1 %in% c("Henarejos", "Libros", "Bugarra", "Prikrnica"))

# Define custom colors for sections
section_colors <- c(
  "Henarejos" = "#D4A017",   # Yellow for Henarejos
  "Libros" = "#17A909",      # Green for Libros
  "Bugarra" = "#17BCC4",     # Cyan/Teal for Bugarra
  "Prikrnica" = "#D91C93"    # Pink for Prikrnica
)

# 8. Create a boxplot for PC1 residuals by sections (variable1) with ordered levels
ggplot(filtered_samples, aes(x = variable1, y = PC1_residuals, fill = variable1)) +
  geom_boxplot(alpha = 0.6) +         # Plot the boxplot with some transparency
  scale_fill_manual(values = section_colors, drop = FALSE) +  # Use drop = FALSE to keep unused factor levels
  labs(title = "Boxplot of PC1 Residuals by Sections", 
       x = "Section", 
       y = "PC1 Residuals") +  # Updated label for y-axis
  theme_minimal() +
  theme(
    legend.title = element_blank(),  # Remove legend title
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
    axis.text.y = element_text(size = 12),  # Increase y-axis text size
    axis.title.x = element_text(size = 14),  # Increase x-axis title size
    axis.title.y = element_text(size = 14),  # Increase y-axis title size
    plot.title = element_text(size = 16, face = "bold")  # Increase plot title size and make it bold
  ) 



### PC1 resids. statist. signi. by sections ###
# Load necessary libraries
library(geomorph)  # For GPA, PCA, and shape analysis
library(ggplot2)   # For plotting
library(MASS)      # For Kernel Density Estimation
library(viridis)   # For colorblind-friendly palettes
library(dplyr)     # For data manipulation
library(FSA)       # For Dunn test

# 1. Read the TPS file and define the sliders
tps_file <- "All_sections.TPS"
landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)
sliders <- define.sliders(3:138) 

# 2. Create the samples data frame
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

# 3. Perform GPA (Generalized Procrustes Analysis) to align the shape data
gpa <- gpagen(landmarks, curves = sliders, print.progress = FALSE) # GPA alignment

# 4. Perform Principal Component Analysis (PCA)
pca_result <- gm.prcomp(gpa$coords) # PCA on GPA aligned data

# 5. Calculate residuals for PC1
# Assume you have a linear model for PC1, e.g., PC1 ~ variable1
lm_model <- lm(pca_result$x[, 1] ~ samples$variable1) # Linear model for PC1
samples$PC1_residuals <- residuals(lm_model)  # Calculate residuals for PC1

# 6. Reorder the levels of variable1 to ensure the desired order in the boxplot with full names
samples$variable1 <- factor(samples$variable1, 
                            levels = c("Clp", "He", "Li", "Bu", "Dr", "Pr"),
                            labels = c("Calasparra", "Henarejos", "Libros", "Bugarra", "Drežnica", "Prikrnica"))


# 8. Perform Kruskal-Wallis test to compare PC1 residuals between all sections (variable1)
kruskal_test <- kruskal.test(PC1_residuals ~ variable1, data = samples)
print(kruskal_test)

# If Kruskal-Wallis is significant, perform pairwise comparisons using Dunn test
if (kruskal_test$p.value < 0.05) {
  dunn_test <- dunnTest(PC1_residuals ~ variable1, data = samples, method = "bonferroni")
  print(dunn_test)
}



### PC1 resids by part of sephard, province ###

# Load necessary libraries
library(geomorph)  # For GPA, PCA, and shape analysis
library(ggplot2)   # For plotting
library(MASS)      # For Kernel Density Estimation
library(viridis)   # For colorblind-friendly palettes
library(dplyr)     # For data manipulation
library(FSA)       # For Dunn test

# 1. Read the TPS file and define the sliders
tps_file <- "All_sections.TPS"
landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)
sliders <- define.sliders(3:138) 

# 2. Create the samples data frame
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

# 3. Perform GPA (Generalized Procrustes Analysis) to align the shape data
gpa <- gpagen(landmarks, curves = sliders, print.progress = FALSE) # GPA alignment

# 4. Perform Principal Component Analysis (PCA)
pca_result <- gm.prcomp(gpa$coords) # PCA on GPA aligned data

# 5. Calculate residuals for PC1
lm_model <- lm(pca_result$x[, 1] ~ samples$variable1) # Linear model for PC1
samples$PC1_residuals <- residuals(lm_model)  # Calculate residuals for PC1

# 6. Reorder the levels of variable1 to ensure the desired order in the boxplot with full names
samples$variable1 <- factor(samples$variable1, 
                            levels = c("Clp", "He", "Li", "Bu", "Dr", "Pr"),
                            labels = c("Calasparra", "Henarejos", "Libros", "Bugarra", "Drežnica", "Prikrnica"))

# 7. Group sections into North-Eastern and Western parts
samples$Group <- ifelse(samples$variable1 %in% c("Prikrnica", "Drežnica"), 
                        "North-Eastern part", 
                        "Western part")


# 8. Create a boxplot for PC1 residuals by groups with reordered levels
samples$Group <- factor(samples$Group, levels = c("Western part", "North-Eastern part"))  # Reorder levels

ggplot(samples, aes(x = Group, y = PC1_residuals, fill = Group)) +
  geom_boxplot(alpha = 0.6) +         # Plot the boxplot with some transparency
  scale_fill_manual(values = c("Western part" = "lightblue", "North-Eastern part" = "grey")) +  # Define custom colors for groups
  labs(title = "Boxplot of PC1 Residuals by Sephardic part", 
       x = "Sephardic part", 
       y = "PC1 Residuals") +  # Updated label for y-axis
  theme_minimal() +
  theme(
    legend.title = element_blank(),  # Remove legend title
    axis.text.x = element_text(size = 12),  # Increase x-axis text size
    axis.text.y = element_text(size = 12),  # Increase y-axis text size
    axis.title.x = element_text(size = 14),  # Increase x-axis title size
    axis.title.y = element_text(size = 14),  # Increase y-axis title size
    plot.title = element_text(size = 16, face = "bold")  # Increase plot title size and make it bold
  ) 


# 9. Perform Kruskal-Wallis test to compare PC1 residuals between groups
kruskal_test <- kruskal.test(PC1_residuals ~ Group, data = samples)
print(kruskal_test)

# If Kruskal-Wallis is significant, perform pairwise comparisons using Dunn test
if (kruskal_test$p.value < 0.05) {
  dunn_test <- dunnTest(PC1_residuals ~ Group, data = samples, method = "bonferroni")
  print(dunn_test)
}




### relationship length and PC1 resid. ###  
### Load necessary libraries
if (!require("geomorph")) install.packages("geomorph", dependencies = TRUE)
if (!require("ggplot2")) install.packages("ggplot2", dependencies = TRUE)
if (!require("dplyr")) install.packages("dplyr", dependencies = TRUE)

library(geomorph)
library(ggplot2)
library(dplyr)

# Function to calculate Euclidean distance
calculate_distance <- function(x1, y1, x2, y2) {
  return(sqrt((x2 - x1)^2 + (y1 - y2)^2))
}

# Function to process the data and calculate scaled distance (using base R)
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
    
    # Use apply to convert all columns to numeric
    landmarks <- as.data.frame(apply(landmarks, 2, as.numeric))
    
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

# Path to your data file
tps_file <- "All_sections.TPS"

# Process the data to calculate distances
scaled_distances_results <- process_data(tps_file)

# Conduct Generalized Procrustes Analysis (GPA)
landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)
gpa <- gpagen(landmarks)

# Perform Principal Component Analysis (PCA) on Procrustes-aligned data
pca_results <- gm.prcomp(gpa$coords)

# Extract PC scores (PC1 and PC2)
pc_scores <- pca_results$x[, 1:2]  # Getting PC1 and PC2

# Calculate the residuals
residuals <- array(0, dim = dim(gpa$coords))  # Initialize residual array

for (i in 1:dim(gpa$coords)[3]) {
  residuals[,,i] <- gpa$coords[,,i] - gpa$consensus  # Subtract consensus from each specimen
}

# Sum the squared residuals across landmarks and dimensions (for each specimen)
pc1_residuals <- apply(residuals, 3, function(x) sum(x[,1]^2))  # PC1 residuals
pc2_residuals <- apply(residuals, 3, function(x) sum(x[,2]^2))  # PC2 residuals

# Combine results in a data frame with the specimen names
residuals_data <- data.frame(
  Specimen = dimnames(landmarks)[[3]],
  PC1_Residuals = pc1_residuals,
  PC2_Residuals = pc2_residuals
)

# Create a sample properties data frame
samples <- data.frame(
  Specimen = dimnames(landmarks)[[3]],
  Country = sapply(dimnames(landmarks)[[3]], function(x) strsplit(x, "_")[[1]][1])  # Extract country from specimen name
)

# Map abbreviations to full country names
country_full_names <- c("BAH" = "Bosnia and Herzegovina", 
                        "SP" = "Spain", 
                        "SL" = "Slovenia")

# Create a new column with full country names
samples$FullCountry <- country_full_names[samples$Country]

# Assign groups
samples$Group <- ifelse(samples$Country %in% c("SL", "BAH"), "North-Eastern part", "Western part")

# Merge residuals data with samples data (to include country information)
merged_data <- merge(residuals_data, samples, by = "Specimen")

# Extract the scaled distances from the processed data
scaled_distances <- unlist(lapply(scaled_distances_results, function(x) x$distances))

# Add scaled distances to the merged data frame
merged_data$Scaled_Distance <- scaled_distances[1:nrow(merged_data)]  # Ensure sizes match

# Calculate regression slopes for each group
slopes <- merged_data %>%
  group_by(Group) %>%
  summarize(slope = coef(lm(Scaled_Distance ~ PC1_Residuals))[2])  # Get slope (coefficient of PC1)

# Print slopes
print(slopes)

# Plot the residuals for PC1 against scaled distances, coloring by group and adding regression lines
ggplot(merged_data, aes(x = PC1_Residuals, y = Scaled_Distance, color = Group, shape = Group)) +
  geom_point(aes(size = Scaled_Distance)) +  # Fully opaque points
  geom_smooth(method = "lm", aes(group = Group), se = FALSE) +  # Add regression line for each group
  scale_color_manual(name = "Country Group", 
                     values = c("North-Eastern part" = "lightblue", "Western part" = "gray")) +  # Define colors for groups
  scale_shape_manual(name = "Country Group", 
                     values = c("North-Eastern part" = 16, "Western part" = 17)) +  # Circle for NE, triangle for Western
  labs(title = "PC1 Residuals vs Length", 
       x = "PC1 Residuals",  # Label for x-axis
       y = "Length (µm)") +  # Updated label for y-axis
  theme_minimal() +  # Use a clean theme
  theme(legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8)) +
  scale_size_continuous(name = "Scaled Distance", range = c(1, 4))  # Adjust size range for visibility


# Subset the data for each group
north_eastern_sections <- subset(merged_data, Group == "North-Eastern part")
western_sections <- subset(merged_data, Group == "Western part")

# Fit linear models for each group
lm_north_eastern <- lm(Scaled_Distance ~ PC1_Residuals, data = north_eastern_sections)
lm_western <- lm(Scaled_Distance ~ PC1_Residuals, data = western_sections)

# Extract slopes (coefficients for PC1)
slope_north_eastern <- coef(lm_north_eastern)["PC1_Residuals"]
slope_western <- coef(lm_western)["PC1_Residuals"]


