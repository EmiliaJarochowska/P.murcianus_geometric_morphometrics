#### Subsampling and Regression Analysis ####
library(ggplot2)
library(dplyr)
library(visreg)
library(tidyr)
library(purrr)

# Load processed data
set.seed(42)
source("src/import_data.R")
data_combined$FaciesZone <- as.factor(data_combined$FaciesZone)


#### Subsampling ####

# Function to subsample n observations from each facies level
subsample_by_FZ <- function(data, n_per_group = 20, seed = 42) {
  data %>%
    group_by(FaciesZone) %>%
    sample_n(size = min(n_per_group, n()), replace = FALSE) %>%
    ungroup()
}


# Run subsampling 
n_iterations <- 100
n_per_group <- 20
model_results <- list()
fitted_values_list <- list()
subsample_list <- list()

for (i in 1:n_iterations) {
  subsample_data <- subsample_by_FZ(data_combined, n_per_group = n_per_group)
  subsample_list[[i]] <- subsample_data
  
  model <- lm(PC1 ~ Length * FaciesZone, data = subsample_data)
  
  model_results[[i]] <- list(
    iteration = i,
    model = model,
    coefficients = coef(model),
    r_squared = summary(model)$r.squared,
    adj_r_squared = summary(model)$adj.r.squared,
    sigma = summary(model)$sigma
  )
  
  fitted_data <- subsample_data %>%
    select(ID, FaciesZone, Length, PC1) %>%
    mutate(
      iteration = i,
      fitted = fitted(model),
      residual = residuals(model)
    )
  fitted_values_list[[i]] <- fitted_data
}

# Combine all fitted values
all_fitted <- bind_rows(fitted_values_list)

# Extract coefficients across iterations
all_coefficients <- do.call(rbind, lapply(model_results, function(x) x$coefficients))
coef_names <- colnames(all_coefficients)

# Coefficient statistics
coef_summary <- data.frame(
  coefficient = coef_names,
  mean = colMeans(all_coefficients, na.rm = TRUE),
  sd = apply(all_coefficients, 2, sd, na.rm = TRUE),
  lower_95 = apply(all_coefficients, 2, quantile, 0.025, na.rm = TRUE),
  upper_95 = apply(all_coefficients, 2, quantile, 0.975, na.rm = TRUE),
  lower_50 = apply(all_coefficients, 2, quantile, 0.25, na.rm = TRUE),
  upper_50 = apply(all_coefficients, 2, quantile, 0.75, na.rm = TRUE)
)

# Raw coefficients 
print(coef_summary, digits = 4)

# Calculate slopes for each fz in each iteration
slope_data <- data.frame()
FZs <- levels(data_combined$FaciesZone)

for (i in 1:n_iterations) {
  coefs <- model_results[[i]]$coefficients
  
  # Base slope
  base_slope <- coefs["Length"]
  
  for (zone in FZs) {
    if (zone == FZs[1]) {
      # Reference category gets base slope
      slope <- base_slope
    } else {
      # Other zones get base + interaction
      interaction_name <- paste0("Length:FaciesZone", zone)
      if (interaction_name %in% names(coefs) && !is.na(coefs[interaction_name])) {
        slope <- base_slope + coefs[interaction_name]
      } else {
        slope <- base_slope
      }
    }
    slope_data <- rbind(slope_data, data.frame(
      iteration = i,
      FaciesZone = zone,
      slope = slope
    ))
  }
}

# Summarize slopes by FaciesZone
slope_summary <- slope_data %>%
  group_by(FaciesZone) %>%
  summarise(
    mean_slope = mean(slope, na.rm = TRUE),
    sd_slope = sd(slope, na.rm = TRUE),
    lower_95 = quantile(slope, 0.025, na.rm = TRUE),
    upper_95 = quantile(slope, 0.975, na.rm = TRUE),
    lower_50 = quantile(slope, 0.25, na.rm = TRUE),
    upper_50 = quantile(slope, 0.75, na.rm = TRUE),
    .groups = 'drop'
  )

# Export for publication
save(slope_summary, file="data/slope_summary.RData") # used in slope_summary.Rmd

##### Line plot for individual iterations #####

FZp1 <- ggplot(all_fitted, aes(x = Length, y = fitted)) +
  geom_smooth(aes(group = interaction(FaciesZone, iteration), color = FaciesZone), 
              method = "lm", se = FALSE, alpha = 0.1, size = 0.3) +
  geom_smooth(aes(color = FaciesZone), method = "lm", se = TRUE, size = 1.2) +
  scale_color_manual(values = colors_list$FaciesZone) +
  labs(
    x = "Length [µm]",
    y = "PC1",
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )



FZp2 <- ggplot(slope_data, aes(x = FaciesZone, y = slope, fill = FaciesZone)) +
  geom_boxplot(width = 0.2, alpha = 0.7, outlier.shape = NA) +
  geom_point(data = slope_summary, aes(y = mean_slope), 
             size = 3, color = "black", shape = 18) +
  geom_errorbar(data = slope_summary, 
                aes(y = mean_slope, ymin = lower_95, ymax = upper_95),
                width = 0.1, color = "black", size = 0.5) +
  scale_fill_manual(values = colors_list$FaciesZone) +
  labs(
    x = "Facies Zone",
    y = "Slope (Length effect on PC1)",
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
  )

# Diamond: mean, Error bars: 95% CI, Box: IQR, Line: full distribution
ggsave("supplementary_material/Fig.S5.jpg",FZp2, width=170, height=100, units="mm", dpi = 300)
##### Visualize slope distributions with confidence intervals #####

FZp3 <- ggplot(slope_data, aes(x = FaciesZone, y = slope, fill = FaciesZone)) +
  geom_violin(alpha = 0.4) +
  geom_boxplot(width = 0.2, alpha = 0.7, outlier.shape = NA) +
  geom_point(data = slope_summary, aes(y = mean_slope), 
             size = 3, color = "black", shape = 18) +
  geom_errorbar(data = slope_summary, 
                aes(y = mean_slope, ymin = lower_95, ymax = upper_95),
                width = 0.1, color = "black", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_fill_manual(values = colors_list$FaciesZone) +
  labs(
    x = "Facies Zone",
    y = "Slope (Length effect on PC1)",
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
  )

# Diamond: mean, Error bars: 95% CI, Box: IQR, Line: full distribution
ggsave("supplementary_material/Fig.S5.jpg",FZp3, width=170, height=100, units="mm", dpi = 300)

##### R² #####

# R-squared distribution
r_squared_values <- sapply(model_results, function(x) x$r_squared)
# Mean R²
round(mean(r_squared_values), 4)
# SD R²
round(sd(r_squared_values), 4)
# 95% CI R²
round(quantile(r_squared_values, 0.025), 4) 
round(quantile(r_squared_values, 0.975), 4)
