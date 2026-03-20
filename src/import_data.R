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
all_pcs <- PCA$x

#### Separate Prikrnica ####
data_combined$Section <- as.character(data_combined$Section)
data_combined$Section[c(21:62, 187:370)] <- "lower part of Prikrnica"
data_combined$Section[c(63:186, 371:384)] <- "upper part of Prikrnica"


#### Define Region based on Section ####

data_combined$Region <- as.character(data_combined$Region)
data_combined$Region[data_combined$Section %in% c("Henarejos","Libros","Bugarra")] <- "Western"
data_combined$Region[data_combined$Section %in% c("lower part of Prikrnica", "upper part of Prikrnica", "Drežnica")] <- "Northeastern"

#### Define factor orders and color schemes ####

# Define Section colors and custom order
section_order <- c("Henarejos", "Libros", "Bugarra", "lower part of Prikrnica", "upper part of Prikrnica", "Drežnica")
Region_order <- c("Western", "Northeastern")
FaciesZone_order <- c("FZ8", "FZ3")
data_combined$Section <- factor(data_combined$Section, levels = section_order)
data_combined$Region <- factor(data_combined$Region, levels = Region_order)
data_combined$FaciesZone <- factor(data_combined$FaciesZone, levels = FaciesZone_order)

# Color palettes
colors_list <- list(
  FaciesZone = c("FZ8" = "#fc8d59", "FZ3" = "#91bfdb"),
  Country = c("Slovenia" = "grey", "Bosnia and Herzegovina" = "#0033A0", "Spain" = "#D81B60"),
  Region = c("Western" = "#fc8d59", "Northeastern" = "#91bfdb"),
  Section = c(
    "Henarejos" = "#fc8d59", 
    "Libros" = "#fee090",
    "Bugarra" = "#ffffbf",
    "lower part of Prikrnica" = "#91bfdb",
    "upper part of Prikrnica" = "darkblue",
    "Drežnica" = "#4575b4"
  )
)

# Create chirality grouping factor  
specimen_info_matched$Chirality <- factor(ifelse(specimen_info_matched$Chirality == "R", "Right", "Left"))

#### Save all processed data and objects ####
save(
  landmarks, specimen_info_matched, landmarks.gpa, PCA, msho,
  data_combined, colors_list, 
  file = "data/processed_data.RData"
)

