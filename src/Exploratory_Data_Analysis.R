library(tidyr)
library(car)
library(ggplot2)
library(vegan)
source("src/import_data.R")

#### Length ####
##### Outlier detection #####
###### Section ###### 
ggplot(data = data_combined, aes(x = Section, y = Length)) +
  stat_boxplot(geom = 'errorbar', width = 0.3) + 
  geom_boxplot() +
  labs(x = "", y = "Length (µm)") +
  theme_bw()

######  Region / FZ ###### 
ggplot(data = data_combined, aes(x = Region, y = Length)) +
  stat_boxplot(geom = 'errorbar', width = 0.3) + 
  geom_boxplot() +
  labs(x = "", y = "Length (µm)") +
  theme_bw()

######  Chirality ###### 
ggplot(data = data_combined, aes(x = Chirality, y = Length)) +
  stat_boxplot(geom = 'errorbar', width = 0.3) + 
  geom_boxplot() +
  labs(x = "", y = "Length (µm)") +
  theme_bw()

##### Normality #####
######  section ###### 
data_combined %>% 
  ggplot(., aes(Length)) +
  geom_histogram() +
  labs(x = "Length (µm)", ylab = "") +
  theme_bw() +
  facet_grid(Section ~ .)


###### Region / FZ ###### 
data_combined %>% 
  ggplot(., aes(Length)) +
  geom_histogram() +
  labs(x = "Length (µm)", ylab = "") +
  theme_bw() +
  facet_grid(Region ~ .)

###### Chirality ###### 
data_combined %>% 
  ggplot(., aes(Length)) +
  geom_histogram() +
  labs(x = "Length (µm)", ylab = "") +
  theme_bw() +
  facet_grid(Chirality ~ .)



shapiro.test(data_combined$Length)

data_combined %>%
  group_by(Chirality) %>%
  summarise(
    W = shapiro.test(as.numeric(Length))$statistic,
    p_value = shapiro.test(as.numeric(Length))$p.value
  )

data_combined %>%
  group_by(Region) %>%
  summarise(
    W = shapiro.test(as.numeric(Length))$statistic,
    p_value = shapiro.test(as.numeric(Length))$p.value
  )

data_combined %>%
  group_by(Section) %>%
  summarise(
    W = shapiro.test(as.numeric(Length))$statistic,
    p_value = shapiro.test(as.numeric(Length))$p.value
  )




#### Informative PCs ####

#####  Chirality #####  
counts_Chirality <- data_combined %>%
  dplyr::group_by(Chirality) %>%
  dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
  dplyr::arrange(Chirality)

print(counts_Chirality)

set.seed(42)
bd_chirality <- betadisper(dist(PC_scores), group = data_combined$Chirality)

permutest(bd_chirality, permutations = 9999)

tapply(
  bd_chirality$distances,
  bd_chirality$group,
  mean
)

plot(bd_chirality, main = "Dispersion of elements by their chirality")

boxplot(
  bd_chirality$distances ~ bd_chirality$group,
  ylab = "Distance to centroid",
  xlab = "Section",
  main = "Within-chirality morphological disparity")

##### Region ##### 
counts_Region_chirality <- data_combined %>%
  dplyr::group_by(Region, Chirality) %>%
  dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
  dplyr::arrange(Region, Chirality)

print(counts_Region_chirality)

# Compare only Sinistral specimens: Western vs North-Eastern
bd_Sinistral_Region <- betadisper(
  dist(PC_scores[specimen_info_matched$Chirality == "Sinistral", ]),
  group = specimen_info_matched$Region[specimen_info_matched$Chirality == "Sinistral"]
)

permutest(bd_Sinistral_Region)
plot(bd_Sinistral_Region, main = "Dispersion of sinistral elements by region")
tapply(
  bd_Sinistral_Region$distances,
  bd_Sinistral_Region$group,
  mean
)
boxplot(
  bd_Sinistral_Region$distances ~ bd_Sinistral_Region$group,
  ylab = "Distance to centroid",
  xlab = "Section",
  main = "Within-region morphological disparity of sinistral elements")


# Compare only dextral specimens: Western vs North-Eastern

bd_Dextral_Region <- betadisper(
  dist(PC_scores[specimen_info_matched$Chirality == "Dextral", ]),
  group = specimen_info_matched$Region[specimen_info_matched$Chirality == "Dextral"]
)
set.seed(42)
permutest(bd_Dextral_Region)
plot(bd_Dextral_Region, main = "Dispersion of dextral elements by region")
tapply(
  bd_Dextral_Region$distances,
  bd_Dextral_Region$group,
  mean
)



##### Section ####

counts_Section <- data_combined %>%
dplyr::group_by(Section, Chirality) %>%
  dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
  dplyr::arrange(Section)

print(counts_Section)

# sinistral
bd_Sinistral_Section <- betadisper(
  dist(PC_scores[specimen_info_matched$Chirality == "Sinistral", ]),
  group = specimen_info_matched$Section[specimen_info_matched$Chirality == "Sinistral"]
)

permutest(bd_Sinistral_Section)
plot(bd_Sinistral_Section,
     main = "Dispersion of sinistral elements by section",
     label = FALSE)

legend("topleft",
       legend = levels(bd_Sinistral_Section$group),
       col = 1:length(levels(bd_Sinistral_Section$group)),
       pch = 19,
       bty = "n")

boxplot(
  bd_Sinistral_Section$distances ~ bd_Sinistral_Section$group,
  ylab = "Distance to centroid",
  xlab = "Section",
  main = "Within-section morphological disparity of sinistral elements")

tapply(
  bd_Sinistral_Section$distances,
  bd_Sinistral_Section$group,
  mean
)


# subset
sin_data <- subset(data_combined, Chirality == "Sinistral")
run_Section_dispersion <- function(df, n_sample = 10){
  
  replicate(500, {
    
    keep <- unlist(
      lapply(levels(df$Section), function(reg){
        
        idx <- which(df$Section == reg)
        
        if(length(idx) > n_sample){
          sample(idx, n_sample)
        } else {
          idx
        }
        
      })
    )
    
    bd <- betadisper(
      dist(df[keep, paste0("PC",1:7)]),
      group = df$Section[keep]
    )
    
    test <- permutest(
      bd,
      permutations = 999
    )
    
    test$tab$`Pr(>F)`[1]
    
  })
}



results_sin_disp_sect <- run_Section_dispersion(sin_data)

mean(results_sin_disp_sect < 0.05)


# dextral 

bd_Dextral_Section <- betadisper(
  dist(PC_scores[specimen_info_matched$Chirality == "Dextral", ]),
  group = specimen_info_matched$Section[specimen_info_matched$Chirality == "Dextral"]
)

permutest(bd_Dextral_Section)
plot(bd_Dextral_Section, main = "Dispersion of dextral elements by section",
     label = FALSE)

legend("topleft",
       legend = levels(bd_Dextral_Section$group),
       col = 1:length(levels(bd_Dextral_Section$group)),
       pch = 19,
       bty = "n")

boxplot(
  bd_Dextral_Section$distances ~ bd_Dextral_Section$group,
  ylab = "Distance to centroid",
  xlab = "Section",
  main = "Within-section morphological disparity of dextral elements")


tapply(
  bd_Dextral_Section$distances,
  bd_Dextral_Section$group,
  mean
)

#subset
dex_data <- subset(data_combined, Chirality == "Dextral")

run_Section_dispersion_dex <- function(df, n_sample = 8){
  
  replicate(500, {
    
    keep <- unlist(
      lapply(levels(df$Section), function(reg){
        
        idx <- which(df$Section == reg)
        
        if(length(idx) > n_sample){
          sample(idx, n_sample)
        } else {
          idx
        }
        
      })
    )
    
    bd <- betadisper(
      dist(df[keep, paste0("PC",1:7)]),
      group = df$Section[keep]
    )
    
    test <- permutest(
      bd,
      permutations = 999
    )
    
    test$tab$`Pr(>F)`[1]
    
  })
}



results_dex_disp_sect <- run_Section_dispersion_dex(dex_data)

mean(results_dex_disp_sect < 0.05)



