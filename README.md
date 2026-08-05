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
│   ├── Prikrnica_heights.csv     # Stratigraphic heights for the Prikrnica
│   ├── Prikrnica_mode_evolution.TPS     # Raw landmark coordinates (TPS format) with all adults specimens including upper part of Prikrnica
│   ├── Specimens_info_mode_evolution.csv     # Specimen metadata with all adults specimens including upper part of Prikrnica (section, region, facies zone, etc.)
section
│   ├── curveslide.csv            # Semi-landmark sliding definitions
│   ├── processed_data.RData      # Pre-processed landmarks, GPA results, PCA 
├── src/
│   ├── import_data.R             # Read TPS file, run GPA and PCA, save processed_data.RData
│   ├── process_data.R            # Additional data wrangling utilities
│   ├── analyses.R                # Main statistical analyses (MANOVA, disparity, Kruskal-Wallis)
│   ├── perform_pca_tests.R       # Kruskal-Wallis + Dunn post-hoc on PC scores
│   ├── calculate_distance.R      # Pairwise distance calculations
│   ├── Exploratory_Data_Analysis.R  # Exploratory Data Analysis (outlier detection, Normality, dispersion) 
│   ├── mode_evolution_Pr.R       # paleoTS mode-of-evolution analysis
│   ├── plot_distance_boxplot.R   # Morphological distance boxplots
│   ├── plot_distance_violin.R   # Morphological distance boxplots, overlayed violin plot
├── figs/                         # Output figures (Fig. 1–5)
├── supplementary_material/       # Supplementary figures (Fig. S1–S5) and table (S1)
├── Tab/                          # Manuscript tables (Word format, Tab. 1-4)
├── Ontogenetic_stages/           # All data required to generate the three models used to determine ontogenetic stages
├── P_murcianus.rmarkdown         # Report in RMarkdown
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
| `dplyr` | Data manipulation |
| `tidyr` | Data reshaping |
| `pairwiseAdonis` | pairwise Adonis |
| `visreg` | Regression visualisation |
| `purrr` | Functional programming / iteration |
| `scales` | Axis formatting |
| `knitr` | R Markdown rendering |
| `segmented` |  fits segmented regression models with one or more breakpoints  |
| `patchwork` |  combines multiple plots  |


Install all packages at once:

```r
install.packages(c(
  "geomorph", "paleoTS", "vegan", "car", "dunn.test",
  "ggplot2", "ggpubr", "egg", "visreg",
  "dplyr", "tidyr", "purrr", "scales", "pairwiseAdonis", "knitr", "segmented", "patchwork"
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
   source("src/perform_pca_tests.R")
   source("src/calculate_distance.R")
   source("src/mode_evolution_Pr.R")
   ```

3. **Figures**
   ```r

   source("src/plot_distance_boxplot.R")
   source("src/plot_distance_violin.R")

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