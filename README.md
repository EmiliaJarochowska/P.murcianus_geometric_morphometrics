# Morphological plasticity of *Pseudofurnishius murcianus*

Code and data supporting the manuscript "Morphological plasticity of the Middle Triassic conodont *Pseudofurnishius murcianus* in the Western Tethys"

## Authors

Katja Oselj
katja.oselj@geo-zs.si
(maintainer)

Emilia Jarochowska 
e.b.jarochowska@uu.nl

Luka Gale 
luka.gale@geo-zs.si
luka.gale@ntf.uni-lj.si

Tea Kolar-Jurkovšek 
tea.kolar-jurkovsek@geo-zs.si

Bogdan Jurkovšek

Gonçalo Silvério
gsilverio@uevora.pt

Carlos Martínez-Pérez 
carlos.martinez-perez@bristol.ac.uk
carlos.martinez-perez@uv.es

## Repository structure

```
P.murcianus_geometric_morphometrics/
├── data/
│   ├── All_sections.TPS          # Raw landmark coordinates (TPS format)
│   ├── Specimens_info.csv        # Specimen metadata (section, region, facies zone, etc.)
│   ├── Pr_samples.csv            # Prikrnica sample data
│   ├── Prikrnica_heights.csv     # Stratigraphic heights for the Prikrnica section
│   ├── curveslide.csv            # Semi-landmark sliding definitions
<<<<<<< HEAD
│   ├── processed_data.RData      # Pre-processed landmarks, GPA results, PCA 
├── src/
│   ├── import_data.R             # Read TPS file, run GPA and PCA, save processed_data.RData
│   ├── process_data.R            # Additional data wrangling utilities
│   ├── analyses.R                # Main statistical analyses (PERMANOVA, disparity, Kruskal-Wallis)
│   ├── Exploratory Data Analysis.R # Shapiro-Wilk normality tests, disparity test
│   ├── perform_pca_tests.R       # Kruskal-Wallis + Dunn post-hoc on PC scores
│   ├── calculate_distance.R      # Pairwise distance calculations
│   ├── mode_evolution_Pr.R       # paleoTS mode-of-evolution analysis
│   ├── plot_distance_boxplot.R   # Morphological distance boxplots
│   ├── plot_distance_violin.R   # Length
├── figs/                         # Output figures (Fig. 1–6)
├── supplementary_material/       # Supplementary figures (Fig. S1–S3) and table (S1)
├── Tab/                          # Manuscript tables (Word format)
├── RMarkdown report
=======
│   ├── Table_S4.csv              # Supplementary data table S4
│   ├── processed_data.RData      # Pre-processed landmarks, GPA results, PCA scores
│   ├── regression_part_summary.RData  # Regression model summaries
│   └── slope_summary.RData       # Slope estimates from subsampling
├── src/
│   ├── import_data.R             # Read TPS file, run GPA and PCA, save processed_data.RData
│   ├── process_data.R            # Additional data wrangling utilities
│   ├── analyses.R                # Main statistical analyses (MANOVA, disparity, Kruskal-Wallis)
│   ├── perform_normality_tests.R # Shapiro-Wilk normality tests
│   ├── perform_pca_tests.R       # Kruskal-Wallis + Dunn post-hoc on PC scores
│   ├── calculate_distance.R      # Pairwise distance calculations
│   ├── Variances.R               # Levene's tests for variance homogeneity
│   ├── Regression_model.R        # Linear regression (Length ~ PC1) with subsampling
│   ├── Regression_facies.R       # Regression stratified by facies zone
│   ├── mode_evolution_Pr.R       # paleoTS mode-of-evolution analysis
│   ├── plot_pca_scatter.R        # PCA scatter plots
│   ├── plot_pca_boxplot.R        # PC score boxplots by group
│   ├── plot_distance_boxplot.R   # Morphological distance boxplots
│   ├── plot_length_histogram.R   # Length frequency histograms
│   └── plot_correlation.R        # Correlation plots
├── figs/                         # Output figures (Fig. 1–9, mode_of_evolution.pdf)
├── supplementary_material/       # Supplementary figures (Fig. S1–S7) and tables (S1–S4)
├── Tab/                          # Manuscript tables (Word format)
>>>>>>> 443b2cd9fa0a39c50bd420c096615aa8bf3f36db
└── P.murcianus_geometric_morphometrics.Rproj
```
## Software requirements

R (≥ 4.1.0) with the following packages:

| Package | Purpose |
|---------|---------|
| `geomorph` | Landmark import, GPA, PCA, morphometric analyses |
| `paleoTS` | Paleontological time series / mode-of-evolution fitting |
| `vegan` | Multivariate community-ecology statistics |
| `car` | Levene's test for variance homogeneity |
| `dunn.test` | Post-hoc pairwise Dunn test |
| `ggplot2` | Plotting |
| `ggpubr` | Multi-panel figure assembly |
| `egg` | Additional ggplot2 extensions |
<<<<<<< HEAD
| `dplyr` | Data manipulation |
| `tidyr` | Data reshaping |
| `scales` | Axis formatting |
| `pairwiseAdonis` | pairwise Adonis |
| `knitr` | R Markdown rendering |


=======
| `visreg` | Regression visualisation |
| `dplyr` | Data manipulation |
| `tidyr` | Data reshaping |
| `purrr` | Functional programming / iteration |
| `scales` | Axis formatting |
| `knitr` | R Markdown rendering |



Install all packages at once:

```r
install.packages(c(
  "geomorph", "paleoTS", "vegan", "car", "dunn.test",
  "ggplot2", "ggpubr", "egg", "visreg",
  "dplyr", "tidyr", "purrr", "scales", "pairwiseAdonis", "knitr"
  "dplyr", "tidyr", "purrr", "scales", "knitr"
))
```

## Running instructions

Open `P.murcianus_geometric_morphometrics.Rproj` in RStudio (or set the project root as the working directory), then run the scripts in the following order:

1. **Import and pre-process data**
   ```r
   source("src/import_data.R")   # reads All_sections.TPS, runs GPA + PCA,
                                  # saves data/processed_data.RData
   source("src/process_data.R")
   ```

2. **Statistical analyses**
   ```r
   source("src/analyses.R")
   source("src/perform_pca_tests.R")
   source("src/calculate_distance.R")
   source("src/perform_normality_tests.R")
   source("src/perform_pca_tests.R")
   source("src/calculate_distance.R")
   source("src/Variances.R")
   source("src/Regression_model.R")
   source("src/Regression_facies.R")
   source("src/mode_evolution_Pr.R")
   ```

3. **Figures**
   ```r

   source("src/plot_distance_boxplot.R")
   source("src/plot_distance_violin.R.R")
   source("src/plot_pca_scatter.R")
   source("src/plot_pca_boxplot.R")
   source("src/plot_distance_boxplot.R")
   source("src/plot_length_histogram.R")
   source("src/plot_correlation.R")

   ```

Steps 2 and 3 load `data/processed_data.RData` (produced in step 1); pre-computed `.RData` files are included so steps 2–3 can be run independently.


## License 

Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

## Copyright

Copyright 2026 the Geological Survey of Slovenia, Ljubljana University, Utrecht University, Instituto de Ciências da Terra, University of Bristol and University of Valencia

## Funding

First author KO was funded by Slovenian Research and Innovation Agency (research core funding No. P1-0011), CMP by the Ministry of Science and Innovation of Spain (Research Project PID2020-117373GA-I00) and GS by Fundação para a Ciência e Tecnologia (PhD grant 2020.08450.BD, and ICT funds UIDB/04683/2020 and UIDP/04683/2020). EJ was funded by the European Union (ERC, MindTheGap, StG project no 101041077). Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Council. Neither the European Union nor the granting authority can be held responsible for them.