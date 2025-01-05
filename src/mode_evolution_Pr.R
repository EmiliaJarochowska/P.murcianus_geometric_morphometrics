#### Include code for data import ####

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

pca_Pr$Rock_sample <- droplevels(pca_Pr$Rock_sample)
Pr_samples <- unique(pca_Pr$Rock_sample)
write.csv(Pr_samples,
          file = "data/Pr_samples.csv",
          sep = ".",
          col.names = T)

Pr_heights <- read.csv(file = "data/Prikrnica_heights.csv",
                       header = F)
colnames(Pr_heights) <- c("Sample", "Height")

pca_Pr <- pca_Pr[order(pca_Pr$Rock_sample),]
Pr_heights <- Pr_heights[order(Pr_heights$Sample),]

pca_Pr_heights <- merge(x = pca_Pr,
                        y = Pr_heights,
                        by.x = "Rock_sample",
                        by.y = "Sample",
                        all.x = all)
