#### Import script and data preprocessing ####

# Load required libraries
library(dplyr)
set.seed(42)

#### Import landmarks and specimen information ####

tps_file = "data/All_sections.TPS"


# Import landmarks from TPS file
landmarks <- geomorph::readland.tps(tps_file, specID = "ID", readcurves = TRUE)

# Import specimen information from CSV
specimen_info <- read.csv("data/Specimens_info.csv", header = TRUE)

#### Data cleaning and preparation ####

specimen_info$Region <- trimws(specimen_info$Region)

# Clean specimen IDs: trim whitespace and unify case
specimen_info$ID <- specimen_info$ID %>% 
  trimws() %>% 
  toupper()

tps_ids <- dimnames(landmarks)[[3]] %>% 
  trimws() %>% 
  toupper()

dimnames(landmarks)[[3]] <- tps_ids

# Check for unmatched specimen IDs
missing_in_csv <- tps_ids[!tps_ids %in% specimen_info$ID]
missing_in_tps <- specimen_info$ID[!specimen_info$ID %in% tps_ids]

if (length(missing_in_csv) > 0) {
  cat("TPS specimen IDs NOT found in Specimens_info.csv:\n")
  print(missing_in_csv)
}

if (length(missing_in_tps) > 0) {
  cat("Specimens_info.csv IDs NOT found in TPS:\n")
  print(missing_in_tps)
}

# Stop if any missing IDs to avoid misalignment
if (length(missing_in_csv) > 0 || length(missing_in_tps) > 0) {
  stop("Mismatch between TPS IDs and specimen_info IDs. Please fix and rerun.")
}

# Match specimen_info order to TPS specimen IDs
match_rows <- match(tps_ids, specimen_info$ID)
specimen_info_matched <- specimen_info[match_rows, ]
rm(tps_ids)
rm(missing_in_csv)
rm(missing_in_tps)
rm(match_rows)

#### Geometric morphometric preprocessing ####

# Define sliding landmarks
sliders <- geomorph::define.sliders(3:138)
# Perform Generalized Procrustes Analysis (GPA)
landmarks.gpa <- geomorph::gpagen(landmarks, curves = sliders)


# Perform PCA on GPA-aligned coordinates
PCA <- geomorph::gm.prcomp(landmarks.gpa$coords)
plot(PCA)


# Calculate mean shape
msho <- geomorph::mshape(landmarks.gpa$coords)

#### Calculate mean distances and combine data ####

source("src/calculate_distance.R")
source("src/process_data.R")

lengths <- process_data(tps_file)

# Combine all data
data_combined <- specimen_info_matched %>%
  dplyr::left_join(lengths, by = "ID")

# Add PCA scores
data_combined$PC1 <- PCA$x[, 1]
data_combined$PC2 <- PCA$x[, 2]
data_combined$PC3 <- PCA$x[, 3]
data_combined$PC4 <- PCA$x[, 4]
data_combined$PC5 <- PCA$x[, 5]
data_combined$PC6 <- PCA$x[, 6]
data_combined$PC7 <- PCA$x[, 7]

all_pcs <- PCA$x


# Function to determine the number of informative PCs using the Broken Stick Model 
# After GUENSER et al., 2022
doPcaSignif <- function(eigvals) {
  
  acron <- character()
  criterion <- character()
  nsigncomp <- numeric()
  n <- length(eigvals)
  
  # Percentage variance explained
  pcvars <- 100 * eigvals / sum(eigvals)
  
  # Broken stick model
  bsm <- data.frame(j = seq_len(n), p = 0)
  bsm$p[1] <- 1 / n
  for(i in 2:n){
    bsm$p[i] <- bsm$p[i - 1] + (1 / (n + 1 - i))
  }
  bsm$p <- 100 * bsm$p / n
  bsvars <- rev(bsm$p)
  
  acron <- "BSM"
  criterion <- "Broken stick model"
  nsigncomp <- sum(pcvars >= bsvars)
  
  data.frame(acron, criterion, nsigncomp)
}

# Eigenvalues from your geomorph PCA
pcasign <- doPcaSignif(PCA$d^2)

print(pcasign)


# Plot Eigenvalues of first 20 PCs, to visualise
eigvals_all <- (PCA$d)^2

# Number of PCs to display
npcs <- min(20, length(eigvals_all))

eigvals <- eigvals_all[1:npcs]

# % variance explained
pcvars <- 100 * eigvals / sum(eigvals_all)

# Broken-stick expectations
n <- length(eigvals_all)

bsvars_all <- sapply(1:n, function(i) {
  100 * sum(1/(i:n))/n
})

bsvars <- bsvars_all[1:npcs]

# Plot
plot(pcvars, type = "b", pch = 19,
     xlab = "Principal Component",
     ylab = "Variance explained (%)",
     xlim = c(1, npcs),
     ylim = c(0, max(pcvars, bsvars) * 1.05),
     main = "Scree Plot with Broken Stick Model")

lines(bsvars, type = "b", pch = 16, col = "red", lwd = 2)
legend("topright",
       legend = c("Observed", "Broken stick"),
       fill = c(NA, NA),
       border = c(NA, NA),
       lty = c(1, 1),
       pch = c(16, 16),
       col = c("black", "red"))

informative_pcs <- PCA$x[, 1:pcasign$nsigncomp]

data_combined <- cbind(
  data_combined,
  informative_pcs
)

PC_scores <- data_combined[, c("Comp1", "Comp2", "Comp3", 
                               "Comp4", "Comp5", "Comp6")]


#### Define Region based on Section ####

data_combined$Region <- as.character(data_combined$Region)
data_combined$Region[data_combined$Section %in% c("Henarejos","Libros","Bugarra")] <- "Western"
data_combined$Region[data_combined$Section %in% c("Prikrnica", "Drežnica")] <- "Northeastern"


#### Define factor orders and color schemes ####

# Define Section colors and custom order
section_order <- c("Henarejos", "Libros", "Bugarra", "Prikrnica", "Drežnica")
Region_order <- c("Western", "Northeastern")
FaciesZone_order <- c("FZ8", "FZ3")
data_combined$Section <- factor(data_combined$Section, levels = section_order)
data_combined$Region <- factor(data_combined$Region, levels = Region_order)
data_combined$FaciesZone <- factor(data_combined$FaciesZone, levels = FaciesZone_order)
data_combined$Chirality <- factor(data_combined$Chirality,levels = c("L", "R"), labels = c("Sinistral", "Dextral"))

# Color palettes
colors_list <- list(
  FaciesZone = c("FZ8" = "#fc8d59", "FZ3" = "#91bfdb"),
  Country = c("Slovenia" = "grey", "Bosnia and Herzegovina" = "#0033A0", "Spain" = "#D81B60"),
  Region = c("Western" = "#fc8d59", "Northeastern" = "#91bfdb"),
  Chirality = c("Sinistral" = "#fc8d59", "Dextral" = "#91bfdb" ),
  Section = c(
    "Henarejos" = "#fc8d59", 
    "Libros" = "#fee090",
    "Bugarra" = "#ffffbf",
    "Prikrnica" = "darkblue",
    "Drežnica" = "#4575b4"
  )
)


# Create chirality grouping factor  
specimen_info_matched$Chirality <- factor(ifelse(specimen_info_matched$Chirality == "R", "Dextral", "Sinistral"))

#### Save all processed data and objects ####
save(
  landmarks, specimen_info_matched, landmarks.gpa, PCA, msho,
  data_combined, colors_list, 
  file = "data/processed_data.RData"
)

