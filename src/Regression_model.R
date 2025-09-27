#### Dependencies ####
library(ggplot2)
library(dplyr)
library(visreg)

# Load and clean up data, set global variables
#source("src/import_data.R") 
load("data/processed_data.RData")
data_combined$FaciesZone <- as.factor(data_combined$FaciesZone)

length_model <- lm(PC1 ~ Length, data = data_combined)
summary(length_model)
visreg(length_model, gg = TRUE)
length_part_model <- lm(PC1 ~ Length * Part, data = data_combined)
summary(length_part_model)
visreg(length_part_model, gg = TRUE, "Length", by="Part")
length_country_model <- lm(PC1 ~ Length * Country, data = data_combined)
summary(length_country_model)
visreg(length_country_model, gg = TRUE)
length_section_model <- lm(PC1 ~ Length * Section, data = data_combined)
summary(length_section_model)
visreg(length_section_model, gg = TRUE, "Length", by="Section")
length_facies_model <- lm(PC1 ~ Length * FaciesZone + Section, data = data_combined)
summary(length_facies_model)
visreg(length_facies_model, "Length", by="FaciesZone", gg = TRUE, layout=c(3,1))

