#### Include code for data import ####


require(paleoTS)
library(dplyr)
set.seed(42)

#### Import landmarks and specimen information ####

tps_file = "data/Prikrnica_mode_evolution.TPS"


# Import landmarks from TPS file
landmarks <- geomorph::readland.tps(tps_file, specID = "ID", readcurves = TRUE)

# Import specimen information from CSV
specimen_info <- read.csv("data/Specimens_info_mode_evolution.csv", header = TRUE)

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

samples <- data.frame(name = rep(NA, dim(landmarks)[3]),
                     country = rep(NA, dim(landmarks)[3]),
                     section = rep(NA, dim(landmarks)[3]),
                     rock_sample = rep(NA, dim(landmarks)[3]),
                     element_number = rep(NA, dim(landmarks)[3]))

for (i in 1:dim(landmarks)[3]) {
  samples$name[i] <- dimnames(landmarks)[[3]][i]
  samples$country[i] <- unlist(strsplit(samples$name[i], "_"))[1]
  samples$section[i] <- unlist(strsplit(samples$name[i], "_"))[2]
  samples$rock_sample[i] <- unlist(strsplit(samples$name[i], "_"))[3]
  samples$element_number[i] <- unlist(strsplit(samples$name[i], "_"))[4]
}

samples$country <- as.factor(samples$country)
samples$rock_sample <- as.factor(samples$rock_sample)
samples$country <- factor(samples$country,
                          levels = c("SL", "SP", "BAH"),
                          labels = c("Slovenia", "Spain", "Bosnia and Herzegovina"))


pca_data <- data.frame(PC1 = PCA$x[, 1], 
                       PC2 = PCA$x[, 2], 
                       Country = samples$country,
                       Section = samples$section,
                       Rock_sample = samples$rock_sample,
                       Element_number = samples$element_number)

pca_Pr <- pca_data[pca_data$Section == "PR",]
rm(pca_data)

# Remove A at the end of sample names
# Here we assume that A, B, C etc are subsequent samples taken from the same 
# stratigraphic positions and can be pooled

pca_Pr$Rock_sample <- gsub("\\D$", 
                           x = pca_Pr$Rock_sample,
                           replacement = "")
pca_Pr$Rock_sample <- as.factor(pca_Pr$Rock_sample)

Pr_heights <- read.csv(file = "data/Prikrnica_heights.csv",
                       header = F)
colnames(Pr_heights) <- c("Rock_sample", "Height")

pca_Pr <- pca_Pr[order(pca_Pr$Rock_sample),]
Pr_heights <- Pr_heights[order(Pr_heights$Rock_sample),]

pca_Pr_heights = dplyr::left_join(x=pca_Pr,
                            y=Pr_heights,
                            by="Rock_sample")

PC1_ts <- pca_Pr_heights %>%
  group_by(Height) %>%
  summarize(
    count = n(),
    mean = mean(PC1, na.rm = TRUE),
    var = var(PC1, na.rm = TRUE)
  )
PC1_ts <- PC1_ts %>% relocate("Height", .after = last_col())

# Plot the morphospace to see if there is any relationship between position 
# in the section and in the morphospace
# Point size is proportional to height in the section

plot(pca_Pr$PC1, pca_Pr$PC2,
     cex = (PC1_ts$Height-5)^(1/3),
     xlab = "PC1",
     ylab = "PC2",
     main = "Relationship between position in the section and in the morphospace (Point size is proportional to height in the section)",
     col = "orange",
     pch = 16
)


PC1_ts <- PC1_ts[order(PC1_ts$Height),] 
PC1_ts <- PC1_ts[complete.cases(PC1_ts), ] 

PC1_mode <- paleoTS::as.paleoTS(mm = PC1_ts$mean,
                                vv = PC1_ts$var,
                                nn = PC1_ts$count,
                                tt = PC1_ts$Height,
                                oldest = "first",
                                reset.time = T)  

mode_evolution <- fit9models(PC1_mode, 
                             method = "AD",
                             pool = FALSE)
write.csv(mode_evolution, file="supplementary_material/Tab.S4.csv")

plot(PC1_mode)



# fit a punctuated model from the data
punc <- fitGpunc(PC1_mode, oshare = FALSE, pool = FALSE)

jpeg(file="figs/Fig.4.jpg", width = 2000, height = 1000, res = 300)
plot(PC1_mode, modelFit = punc)
dev.off()

plot(PC1_mode, modelFit = punc)

stasis <- fitSimple(PC1_mode, 
                    model = "Stasis",
                    pool = FALSE)

plot(PC1_mode, modelFit = stasis)
