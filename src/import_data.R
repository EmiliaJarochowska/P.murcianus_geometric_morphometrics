#### Import script and data preprocessing ####

# Load required libraries
library(dplyr)

#### Import landmarks and specimen information ####

tps_file = "data/All_sections.TPS"

# Import landmarks from TPS file
landmarks <- geomorph::readland.tps(tps_file, specID = "ID", readcurves = TRUE)

# Import specimen information from CSV
specimen_info <- utils::read.csv("data/Specimens_info.csv", header = TRUE)

#### Data cleaning and preparation ####

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

#### Create regional classifications ####

# Recode Region from Country
data_combined$Region <- dplyr::recode_factor(
  data_combined$Country,
  "Slovenia" = "North-Eastern part",
  "Bosnia and Herzegovina" = "North-Eastern part",
  "Spain" = "Western part",
  .default = NA_character_
)

# Define the Part variable as a subregion of the Sephardic Province

data_combined$Part <- dplyr::recode_factor(
  data_combined$Country,
  "Slovenia" = "North-Eastern Subprovince",
  "Bosnia and Herzegovina" = "North-Eastern Subprovince", 
  "Spain" = "Western Subprovince",
  .default = NA_character_
)

#### Define factor orders and color schemes ####

# Define Section colors and custom order
section_order <- c("Calasparra", "Henarejos", "Libros", "Bugarra", "Prikrnica", "Drežnica")
data_combined$Section <- factor(data_combined$Section, levels = section_order)

# Color palettes
colors_list <- list(
  FaciesZone = c("FZ3" = "grey", "FZ7" = "#0033A0", "FZ8" = "#D81B60"),
  Country = c("Slovenia" = "grey", "Bosnia and Herzegovina" = "#0033A0", "Spain" = "#D81B60"),
  Region = c("North-Eastern part" = "grey", 
                   "Western part" = "#0033A0"),
  Section = c(
    "Calasparra" = "#d73027",
    "Henarejos" = "#fc8d59", 
    "Libros" = "#fee090",
    "Bugarra" = "#ffffbf",
    "Prikrnica" = "#91bfdb",
    "Drežnica" = "#4575b4"
  ),
  Part = c("Western Subprovince" = "#0033A0", "North-Eastern Subprovince" = "grey")
)

# Create chirality grouping factor  
specimen_info_matched$Chirality <- factor(ifelse(specimen_info_matched$Chirality == "R", "Right", "Left"))

#### Save all processed data and objects ####
save(
  landmarks, specimen_info_matched, landmarks.gpa, PCA, msho,
  data_combined, colors_list, 
  file = "data/processed_data.RData"
)
