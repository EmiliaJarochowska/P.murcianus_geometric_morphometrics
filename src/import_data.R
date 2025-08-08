#### Import script and data preprocessing ####

# Imports

library(geomorph)

# Import landmarks

tps_file <- "All_sections.TPS"
landmarks <- readland.tps(tps_file, specID = "ID", readcurves = TRUE)
sliders <- define.sliders(3:138) 

##### Mean shape, PC1 min and max, PC2 min and max #####

landmarks.gpa<-gpagen(landmarks, curves = sliders)
PCA <- gm.prcomp(landmarks.gpa$coords) 
msho <- mshape(landmarks.gpa$coords) 

#### Prepare sample properties dataset ####

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

summary <- summary(PCA)

save(samples, PCA, msho, summary,
     file = "data/imported.RData")
