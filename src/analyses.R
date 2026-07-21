#### Dependencies ####
library(vegan)
library(geomorph)
library(dunn.test)
library(egg)
library(ggpubr)
library(dplyr)
library(pairwiseAdonis)
library(car)
library(tidyr)


#### Statistical Analysis and Visualization Script ####

# Load and clean up data, set global variables
source("src/import_data.R")


# Or load from saved file if import_data.R had run previously
# load("data/processed_data.RData")


#### Length ####

##### Chirality #####

stats::kruskal.test(Length ~ Chirality, data = data_combined)

aggregate(Length ~ Chirality, data = data_combined, mean)




##### Section #####

aggregate(Length ~ Section, data = data_combined, mean)

kruskal.test(Length ~ Section, data = data_combined)
dunn.test(data_combined$Length, g=data_combined$Section, method="bonferroni")




##### Region #####

aggregate(Length ~ Region, data = data_combined, mean)

kruskal.test(Length ~ Region, data = data_combined)




##### Length plots #####

source("src/plot_distance_boxplot.R")

FZ_length <- plot_distance_boxplot(data_combined, "FaciesZone", colors_list$FaciesZone,)

Region_length <- plot_distance_boxplot(data_combined, "Region", colors_list$Region)

Section_length <- plot_distance_boxplot(data_combined, "Section", colors_list$Section)

Chirality_length <- plot_distance_boxplot(data_combined, "Chirality", colors_list$Chirality)

ggarrange(Chirality_length, Region_length, Section_length,
          ncol=3, widths = c(1,1,2), 
          labels = c("A", "B", "C"))


data_combined %>%
  group_by(Chirality, Region, FaciesZone, Section) %>%
  summarise(
    Median = median(Length, na.rm = TRUE),
    Mean = mean(Length, na.rm = TRUE),
    Count = n(),
    .groups = "drop"
  )


source("src/plot_distance_violin.R")

FZ_length <- plot_distance_violin(data_combined, "FaciesZone", colors_list$FaciesZone, "Facies Zone"
)

Region_length <- plot_distance_violin(data_combined, "Region", colors_list$Region,"Region"
)

Section_length <- plot_distance_violin( data_combined, "Section", colors_list$Section, "Section"
)

Chirality_length <- plot_distance_violin(data_combined, "Chirality", colors_list$Chirality, "Chirality"
)


ggarrange(
  Chirality_length,
  Region_length,
  Section_length,
  ncol = 3,
  widths = c(1, 1, 1.8),
  labels = c("A", "B", "C"),
  align = "hv"
)

#### Statistical Analysis of first 7 PCs ####
##### PCA Analysis by Grouping Variables #####

pca1 <- ggplot(data_combined, aes(x = PC1, y = PC2,
                                  fill = Region,
                                  shape = Chirality)) +
  geom_point(size = 4, alpha = 0.7, color = "black") +
  scale_fill_manual(values = colors_list[["Region"]]) +
  scale_shape_manual(values = c("Sinistral" = 21,
                                "Dextral" = 24)) +
  labs(
    x = "PC1",
    y = "PC2",
    fill = "Region",
    shape = "Chirality"
  ) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21)),
    shape = guide_legend(override.aes = list(fill = "white"))
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )


pca2 <- ggplot(data_combined, aes(x = PC1, y = PC2,
                                  fill = Section,
                                  shape = Chirality)) +
  geom_point(size = 4, alpha = 0.7, color = "black") +
  scale_fill_manual(values = colors_list[["Section"]]) +
  scale_shape_manual(values = c("Sinistral" = 21,
                                "Dextral" = 24)) +
  labs(
    x = "PC1",
    y = "PC2",
    fill = "Section",
    shape = "Chirality"
  ) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21)),
    shape = guide_legend(override.aes = list(fill = "white"))
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )


pca3 <- ggplot(data_combined, aes(x = PC1, y = PC2,
                                  fill = Chirality,
                                  shape = Chirality)) +
  geom_point(size = 4, alpha = 0.7, color = "black") +
  scale_fill_manual(values = colors_list[["Chirality"]]) +
  scale_shape_manual(values = c("Sinistral" = 21,
                                "Dextral" = 24)) +
  labs(
    x = "PC1",
    y = "PC2",
    fill = "Chirality",
    shape = "Chirality"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )


ggarrange(pca3, pca1, pca2, 
          labels = c("A", "B", "C")) 



# sinistral by region
Region_Sinistral <- ggplot(subset(data_combined, Chirality == "Sinistral"),
                            aes(x = PC1, y = PC2,
                                fill = Region)) +
  geom_point(size = 4, alpha = 0.7, shape = 21, color = "black") +
  stat_ellipse(aes(color = Region), linewidth = 1) +
  scale_fill_manual(values = colors_list[["Region"]]) +
  scale_color_manual(values = colors_list[["Region"]]) +
  labs(
    x = "PC1",
    y = "PC2",
    fill = "Region",
    title = "Sinistral chirality by region"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

# dextral by region
Region_Dextral <- ggplot(subset(data_combined, Chirality == "Dextral"),
                          aes(x = PC1, y = PC2,
                              fill = Region)) +
  geom_point(size = 4, alpha = 0.7, shape = 21, color = "black") +
  stat_ellipse(aes(color = Region), linewidth = 1) +
  scale_fill_manual(values = colors_list[["Region"]]) +
  scale_color_manual(values = colors_list[["Region"]]) +
  labs(
    x = "PC1",
    y = "PC2",
    fill = "Region",
    title = "Dextral chirality by region"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")


ggarrange(Region_Sinistral, Region_Dextral, 
          labels = c("A", "B")) 




# sinistral by section
Section_Sinistral <- ggplot(subset(data_combined, Chirality == "Sinistral"),
       aes(x = PC1, y = PC2,
           fill = Section)) +
  geom_point(size = 4, alpha = 0.7, shape = 21, color = "black") +
  stat_ellipse(aes(color = Section), linewidth = 1) +
  scale_fill_manual(values = colors_list[["Section"]]) +
  scale_color_manual(values = colors_list[["Section"]]) +
  labs(
    x = "PC1",
    y = "PC2",
    fill = "Section",
    title = "Sinistral chirality"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

# dextral by section
Section_Dextral <- ggplot(subset(data_combined, Chirality == "Dextral"),
       aes(x = PC1, y = PC2,
           fill = Section)) +
  geom_point(size = 4, alpha = 0.7, shape = 21, color = "black") +
  stat_ellipse(aes(color = Section), linewidth = 1) +
  scale_fill_manual(values = colors_list[["Section"]]) +
  scale_color_manual(values = colors_list[["Section"]]) +
  labs(
    x = "PC1",
    y = "PC2",
    fill = "Section",
    title = "Dextral chirality"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")


ggarrange(Section_Sinistral, Section_Dextral, 
          labels = c("A", "B")) 




##### Extreme shape changes along first 7 PCs axes #####
par(mfrow = c(4, 4))

# PC1 shape changes
geomorph::plotRefToTarget(PCA$shapes$shapes.comp1$min, msho, method = "vector", 
                          main = "PC1 Minimum")
title(main = "PC1 Minimum")
geomorph::plotRefToTarget(PCA$shapes$shapes.comp1$max, msho, method = "vector",
                          main = "PC1 Maximum")
title(main = "PC1 Maximum")
# PC2 shape changes  
geomorph::plotRefToTarget(PCA$shapes$shapes.comp2$min, msho, method = "vector",
                          main = "PC2 Minimum")
title(main = "PC2 Minimum")
geomorph::plotRefToTarget(PCA$shapes$shapes.comp2$max, msho, method = "vector",
                          main = "PC2 Maximum")
title(main = "PC2 Maximum")
# PC3 shape changes  
geomorph::plotRefToTarget(PCA$shapes$shapes.comp3$min, msho, method = "vector",
                          main = "PC3 Minimum")
title(main = "PC3 Minimum")
geomorph::plotRefToTarget(PCA$shapes$shapes.comp3$max, msho, method = "vector",
                          main = "PC3 Maximum")
title(main = "PC3 Maximum")
# PC4 shape changes  
geomorph::plotRefToTarget(PCA$shapes$shapes.comp4$min, msho, method = "vector",
                          main = "PC4 Minimum")
title(main = "PC4 Minimum")
geomorph::plotRefToTarget(PCA$shapes$shapes.comp4$max, msho, method = "vector",
                          main = "PC4 Maximum")
title(main = "PC4 Maximum")
# PC5 shape changes  
geomorph::plotRefToTarget(PCA$shapes$shapes.comp5$min, msho, method = "vector",
                          main = "PC5 Minimum")
title(main = "PC5 Minimum")
geomorph::plotRefToTarget(PCA$shapes$shapes.comp5$max, msho, method = "vector",
                          main = "PC5 Maximum")
title(main = "PC5 Maximum")
# PC6 shape changes  
geomorph::plotRefToTarget(PCA$shapes$shapes.comp6$min, msho, method = "vector",
                          main = "PC6 Minimum")
title(main = "PC6 Minimum")
geomorph::plotRefToTarget(PCA$shapes$shapes.comp6$max, msho, method = "vector",
                          main = "PC6 Maximum")
title(main = "PC6 Maximum")


par(mfrow = c(1, 1))


##### Statistical test #####
###### Chirality ###### 

adonis2(
  PC_scores ~ Chirality,
  data = data_combined,
  permutations = 9999,
  method = "euclidean")


###### Region ###### 
# Compare only Sinistral specimens: Western vs North-Eastern

adonis_Sinistral <- adonis2(
  PC_scores[specimen_info_matched$Chirality == "Sinistral", ] ~ Region,
  data = data_combined[specimen_info_matched$Chirality == "Sinistral", ],
  method = "euclidean",
  permutations = 9999
)
adonis_Sinistral



# Compare only Dextral specimens: Western vs North-Eastern

adonis_Dextral <- adonis2(
  PC_scores[specimen_info_matched$Chirality == "Dextral", ] ~ Region,
  data = data_combined[specimen_info_matched$Chirality == "Dextral", ],
  method = "euclidean",
  permutations = 9999
)

adonis_Dextral




###### Section ###### 
# Sinistral

adonis_Sinistral <- adonis2(
  PC_scores[specimen_info_matched$Chirality == "Sinistral", ] ~ Section,
  data = data_combined[specimen_info_matched$Chirality == "Sinistral", ],
  method = "euclidean",
  permutations = 9999
)

pairwise_results_Section <- pairwise.adonis2(
  PC_scores[specimen_info_matched$Chirality == "Sinistral", ] ~ Section,
  data = data_combined[specimen_info_matched$Chirality == "Sinistral", ],
  method = "euclidean",
  permutations = 9999
)


# subset 

run_Section_test <- function(df, n_sample = 10){
  
  replicate(500, {
    
    keep <- unlist(
      lapply(levels(df$Section), function(sec){
        
        idx <- which(df$Section == sec)
        
        if(length(idx) > n_sample){
          sample(idx, n_sample)
        } else {
          idx
        }
        
      })
    )
    
    test <- adonis2(
      dist(df[keep, paste0("PC",1:7)]) ~ Section,
      data = df[keep, ],
      permutations = 999
    )
    
    test$`Pr(>F)`[1]
    
  })
}


results_sin_section <- run_Section_test(sin_data)

mean(results_sin_section < 0.05)



# Dextral

adonis_Sinistral <- adonis2(
  PC_scores[specimen_info_matched$Chirality == "Dextral", ] ~ Section,
  data = data_combined[specimen_info_matched$Chirality == "Dextral", ],
  method = "euclidean",
  permutations = 9999
)

pairwise_results_Section <- pairwise.adonis2(
  PC_scores[specimen_info_matched$Chirality == "Dextral", ] ~ Section,
  data = data_combined[specimen_info_matched$Chirality == "Dextral", ],
  method = "euclidean",
  permutations = 9999
)


#subset
dex_data <- subset(data_combined, Chirality == "Dextral")
run_Section_test_dex <- function(df, n_sample = 8){
  
  replicate(500, {
    
    keep <- unlist(
      lapply(levels(df$Section), function(sec){
        
        idx <- which(df$Section == sec)
        
        if(length(idx) > n_sample){
          sample(idx, n_sample)
        } else {
          idx
        }
        
      })
    )
    
    test <- adonis2(
      dist(df[keep, paste0("PC",1:7)]) ~ Section,
      data = df[keep, ],
      permutations = 999
    )
    
    test$`Pr(>F)`[1]
    
  })
}


results_dex_section <- run_Section_test_dex(dex_data)

mean(results_dex_section < 0.05)



##### Mean shapes #####
# Chirality
mean_shape_Sinistral <- geomorph::mshape(landmarks.gpa$coords[, , specimen_info_matched$Chirality == "Sinistral"])
mean_shape_Dextral <- geomorph::mshape(landmarks.gpa$coords[, , specimen_info_matched$Chirality == "Dextral"])

geomorph::plotRefToTarget(mean_shape_Sinistral, mean_shape_Dextral, method = "points",
                          main = "Mean Shape: Sinistral -> Dextral")
title(main = "Mean shape of sinistral and Dextral elements (grey dots represents mean shape of sinistral elements, and black dots mean shape of dextral elements")



plot(
  mean_shape_Sinistral[,1],
  mean_shape_Sinistral[,2],
  asp = 1,
  type = "n",
  xlab = "X",
  ylab = "Y",
  main = "Mean shapes of Sinistral and Dextral elements"
)

points(mean_shape_Sinistral[,1], mean_shape_Sinistral[,2],
       pch=21, bg="black")

points(mean_shape_Dextral[,1], mean_shape_Dextral[,2],
       pch=21, bg="grey")


legend(
  "topright",
  legend=c("Sinistral","Dextral"),
  pch=21,
  pt.bg=c("black","grey")
)





# Regions
mean_shape_W <- geomorph::mshape(landmarks.gpa$coords[, , specimen_info_matched$Region == "Western"])
mean_shape_N <- geomorph::mshape(landmarks.gpa$coords[, , specimen_info_matched$Region == "North-Eastern"]
)
geomorph::plotRefToTarget(mean_shape_W, mean_shape_N, method = "points",
                          main = "Mean Shape: Western -> North-Eastern region")
title(main = "Mean shape of Western and North-Eastern region elements (grey dots represents mean shape of elements from Weatern region, and black dots from Northeastern region")



mean_shape_W_Sinistral <- geomorph::mshape(
  landmarks.gpa$coords[, , specimen_info_matched$Region == "Western" &
                         specimen_info_matched$Chirality == "Sinistral"]
)

mean_shape_W_Dextral <- geomorph::mshape(
  landmarks.gpa$coords[, , specimen_info_matched$Region == "Western" &
                         specimen_info_matched$Chirality == "Dextral"]
)

mean_shape_N_Sinistral <- geomorph::mshape(
  landmarks.gpa$coords[, , specimen_info_matched$Region == "North-Eastern" &
                         specimen_info_matched$Chirality == "Sinistral"]
)

mean_shape_N_Dextral <- geomorph::mshape(
  landmarks.gpa$coords[, , specimen_info_matched$Region == "North-Eastern" &
                         specimen_info_matched$Chirality == "Dextral"]
)


geomorph::plotRefToTarget(
  mean_shape_W_Sinistral,
  mean_shape_N_Sinistral,
  method = "points",
  main = "Comparison of Sinistral elements between Western and North-Eastern region"
)
title(main = "Comparison of Sinistral elements between Western and North-Eastern region (grey dots represents mean shape of elements from Weatern region, and black dots from Northeastern region)")


geomorph::plotRefToTarget(
  mean_shape_W_Dextral,
  mean_shape_N_Dextral,
  method = "points",
  main = "Comparison of Sinistral elements between Western and North-Eastern region"
)
title(main = "Comparison of Dextral elements between Western and North-Eastern region (grey dots represents mean shape of elements from Weatern region, and black dots from Northeastern region)")





plot(
  mean_shape_W_Sinistral[,1],
  mean_shape_W_Sinistral[,2],
  asp = 1,
  type = "n",
  xlab = "X",
  ylab = "Y",
  main = "Mean shapes of sinistral elements across regions"
)
points(mean_shape_W_Sinistral[,1], mean_shape_W_Sinistral[,2],
       pch=21, bg="black")

points(mean_shape_N_Sinistral[,1], mean_shape_N_Sinistral[,2],
       pch=21, bg="grey")

legend(
  "topright",
  legend=c("Western region","Northeastern region"),
  pch=21,
  pt.bg=c("black","grey")
)



plot(
  mean_shape_W_Dextral[,1],
  mean_shape_W_Dextral[,2],
  asp = 1,
  type = "n",
  xlab = "X",
  ylab = "Y",
  main = "Mean shapes of dextral elements across regions"
)
points(mean_shape_W_Dextral[,1], mean_shape_W_Dextral[,2],
       pch=21, bg="black")

points(mean_shape_N_Dextral[,1], mean_shape_N_Dextral[,2],
       pch=21, bg="grey")

legend(
  "topright",
  legend=c("W region","NE region"),
  pch=21,
  pt.bg=c("black","grey")
)


# Section 
## Sinistral elements

mean_shape_Henarejos_L <- geomorph::mshape(
  landmarks.gpa$coords[, , specimen_info_matched$Section == "Henarejos" &
                         specimen_info_matched$Chirality == "Sinistral"]
)

mean_shape_Libros_L <- geomorph::mshape(
  landmarks.gpa$coords[, , specimen_info_matched$Section == "Libros" &
                         specimen_info_matched$Chirality == "Sinistral"]
)

mean_shape_Bugarra_L <- geomorph::mshape(
  landmarks.gpa$coords[, , specimen_info_matched$Section == "Bugarra" &
                         specimen_info_matched$Chirality == "Sinistral"]
)

mean_shape_Prikrnica_L <- geomorph::mshape(
  landmarks.gpa$coords[, , specimen_info_matched$Section == "Prikrnica" &
                         specimen_info_matched$Chirality == "Sinistral"]
)

mean_shape_Drežnica_L <- geomorph::mshape(
  landmarks.gpa$coords[, , specimen_info_matched$Section == "Drežnica" &
                         specimen_info_matched$Chirality == "Sinistral"]
)


plot(
  mean_shape_Henarejos_L[,1],
  mean_shape_Henarejos_L[,2],
  asp = 1,
  type = "n",
  xlab = "X",
  ylab = "Y",
  main = "Mean shapes of Sinistral elements by Section"
)

points(mean_shape_Henarejos_L[,1], mean_shape_Henarejos_L[,2],
       pch=21, bg="#fc8d59")

points(mean_shape_Libros_L[,1], mean_shape_Libros_L[,2],
       pch=21, bg="#fee090")

points(mean_shape_Bugarra_L[,1], mean_shape_Bugarra_L[,2],
       pch=21, bg="#ffffbf")

points(mean_shape_Prikrnica_L[,1], mean_shape_Prikrnica_L[,2],
       pch=21, bg="darkblue")

points(mean_shape_Drežnica_L[,1], mean_shape_Drežnica_L[,2],
       pch=21, bg="#4575b4")


legend(
  "topright",
  legend=c("Henarejos","Libros","Bugarra","Prikrnica","Drežnica"),
  pch=21,
  pt.bg=c("#fc8d59","#fee090","#ffffbf","darkblue","#4575b4")
)




## Dextral elements


mean_shape_Henarejos_R <- geomorph::mshape(
  landmarks.gpa$coords[, , specimen_info_matched$Section == "Henarejos" &
                         specimen_info_matched$Chirality == "Dextral"]
)

mean_shape_Libros_R <- geomorph::mshape(
  landmarks.gpa$coords[, , specimen_info_matched$Section == "Libros" &
                         specimen_info_matched$Chirality == "Dextral"]
)

mean_shape_Bugarra_R <- geomorph::mshape(
  landmarks.gpa$coords[, , specimen_info_matched$Section == "Bugarra" &
                         specimen_info_matched$Chirality == "Dextral"]
)

mean_shape_Prikrnica_R <- geomorph::mshape(
  landmarks.gpa$coords[, , specimen_info_matched$Section == "Prikrnica" &
                         specimen_info_matched$Chirality == "Dextral"]
)

mean_shape_Drežnica_R <- geomorph::mshape(
  landmarks.gpa$coords[, , specimen_info_matched$Section == "Drežnica" &
                         specimen_info_matched$Chirality == "Dextral"]
)


plot(
  mean_shape_Henarejos_R[,1],
  mean_shape_Henarejos_R[,2],
  asp = 1,
  type = "n",
  xlab = "X",
  ylab = "Y",
  main = "Mean shapes of Dextral elements by Section"
)

points(mean_shape_Henarejos_R[,1], mean_shape_Henarejos_R[,2],
       pch=21, bg="#fc8d59")

points(mean_shape_Libros_R[,1], mean_shape_Libros_R[,2],
       pch=21, bg="#fee090")

points(mean_shape_Bugarra_R[,1], mean_shape_Bugarra_R[,2],
       pch=21, bg="#ffffbf")

points(mean_shape_Prikrnica_R[,1], mean_shape_Prikrnica_R[,2],
       pch=21, bg="darkblue")

points(mean_shape_Drežnica_R[,1], mean_shape_Drežnica_R[,2],
       pch=21, bg="#4575b4")


legend(
  "topright",
  legend=c("Henarejos","Libros","Bugarra","Prikrnica","Drežnica"),
  pch=21,
  pt.bg=c("#fc8d59","#fee090","#ffffbf","darkblue","#4575b4")
)





