source("src/import_data.R")
require(dplyr)
library("ggplot2")
library(car)


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



PC1_var <- pca_data %>%
  group_by(Rock_sample) %>%
  summarise(
    n = n(),
    PC1_mean = mean(PC1, na.rm = TRUE),
    PC1_var = var(PC1, na.rm = TRUE)
  ) %>%
  ungroup()


pca_data <- merge(pca_data, PC1_var, by = "Rock_sample")


ggplot(pca_data, aes_string(x = "n", y = "PC1_var")) +
  geom_point()

boxplot(PC1 ~ Rock_sample, data = pca_data,
        ylab = "PC1 scores", xlab = "Rock sample",
        cex.axis = 0.7)

# plotted var. for each rock sample
ggplot(pca_data, aes(x = Rock_sample, y = PC1_var)) +
  geom_point() +
  theme_minimal()


# mean and SD
ggplot(pca_data, aes(x = Rock_sample, y = PC1)) +
  geom_point(data = pca_data,aes(x = Rock_sample, y = PC1_mean)) +
  geom_errorbar(data = pca_data, aes(x = Rock_sample, ymin = PC1_mean - sqrt(PC1_var), ymax = PC1_mean + sqrt(PC1_var))
  ) +
  theme_minimal()


var <- ggplot(pca_data, aes(x = Rock_sample, y = PC1)) +
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot() +
  labs(x = "Rock sample", y = "PC1 scores") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45))

ggsave("supplementary_material/Fig.S2.jpg", var, width=170, units="mm", dpi = 300)

shapiro.test(pca_data$PC1)

# data is non-normal distributed
leveneTest(PC1 ~ Rock_sample, data = pca_data)
