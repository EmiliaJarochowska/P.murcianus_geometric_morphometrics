#### Include code for data import ####
source("src/import_data.R")
require(dplyr)
require(paleoTS)

##### 
pca_data <- data.frame(PC1 = PCA$x[, 1], 
                       PC2 = PCA$x[, 2], 
                       Country = samples$country,
                       Section = samples$section,
                       Rock_sample = samples$rock_sample,
                       Element_number = samples$element_number)

pca_Pr <- pca_data[pca_data$Section == "Pr",]

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

fit9models(PC1_mode, 
           method = "AD",
           pool = FALSE)

PC1_sstasis <- opt.joint.Stasis(y = PC1_mode,
                                pool = F)
plot(PC1_mode, modelFit = PC1_sstasis)
