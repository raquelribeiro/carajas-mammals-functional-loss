# Assessing Habitat Loss Impacts on Mammal Functional Diversity and Ecosystem Services in Carajás

[![R Version](https://img.shields.io/badge/R-%E2%89%A5%204.3.2-blue.svg)](https://cran.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 1. Project Overview
This repository contains the computational pipeline and raw data to reproduce the study on the impacts of future climate scenarios and habitat loss on functional diversity and associated ecosystem services for mammals in Carajás, Pará State, Brazil.

By integrating Ecological Niche Modeling (MaxEnt via `ENMwizard`) with trait-based functional diversity metrics, this workflow quantifies how spatial habitat contractions under multiple Global Climate Models (CCSM4, MIROC6, MPI) and Shared Socioeconomic Pathways (SSP245 and SSP585) translate into the degradation of ecosystem services and functional values.

## 2. Repository Structure

    .
    ├── data/
    │   ├── spps_occ_env_thin_f.csv                 # Thinned occurrence coordinates for target species
    │   ├── traits.xlsx                             # Functional traits matrix
    │   └── shp_files.zip                           # Study area spatial boundaries (Shapefiles)
    ├── outputs/                                    # (Empty folder) Outputs are generated here
    │   ├── models/                                 # Spatial models (.rds) saved after execution
    │   └── tables/                                 # Analytical data frames and metrics(.csv) 
    ├── scripts/
    │   ├── 01_SDM_Maxent.R                         # Model calibration, tuning, and consensus projections
    │   └── 02_Functional_Diversity.R               # Trait scoring, functional values, and dbFD metrics
    └── README.md


## 3. Workflow Overview
The analytical pipeline executes in two sequential modules:
1. **Module 1 — Spatial Modeling (`01_SDM_Maxent.R`):**
   * Reads thinned occurrence coordinates (`spps_occ_env_thin_f.csv`) and spatial boundary shapefiles.
   * Buffers species occurrences and filters collinear environmental predictors (r > 0.75).
   * Calibrates and evaluates MaxEnt models across feature classes and regularization multipliers using `ENMeval`.
   * Projects suitability maps into future climate scenarios (2050 and 2070 under SSP245/SSP585) and builds multi-GCM consensus rasters.
2. **Module 2 — Functional Diversity & Ecosystem Services (`02_Functional_Diversity.R`):**
   * Processes species trait matrices (`traits.xlsx`, `functional.xlsx`, and `Eco_data_Mamiferos_Carajas_2019_12.csv`).
   * Calculates species-level baseline functional scores (`base_value`) based on diet, IUCN threat status, foraging stratum, and activity cycle.
   * Cross-references functional scores with current and projected habitat areas (in km²) to calculate relative functional loss.
   * Computes community-level functional diversity metrics (Functional Richness — FRic, Functional Dispersion — FDis, and Functional Redundancy) using `FD::dbFD()`.

## 4. Detailed Description of Each Script

### `01_SDM_Maxent.R`
* **Name:** `01_SDM_Maxent.R`
* **Objective:** Calibrate, evaluate, and project Maxent species distribution models across current and future climate scenarios.
* **Main Inputs:** `data/spps_occ_env_thin_f.csv`, shapefiles (`shp_files.zip`), baseline WorldClim rasters, and CMIP6 future GCM rasters.
* **Main Outputs:** Calibrated model lists and consensus projection rasters saved as `.rds` objects (`ENMeval.res.lst.rds`, `mxnt.mdls.preds.cf_cons.rds`, `mods_thrshld_lst.rds`).
* **Steps Executed:** Occurrence buffer creation (1.5 degrees), correlation cutoff filter (0.75), MaxEnt parameter tuning, block/jackknife spatial evaluation, and consensus ensemble generation.
* **Quality Control:** Spatial block evaluation for species with n > 15 occurrences; jackknife resampling for species with n <= 15.
* **Important Notes:** Environmental rasters must be downloaded separately due to file size constraints.

### `02_Functional_Diversity.R`
* **Name:** `02_Functional_Diversity.R`
* **Objective:** Compute intrinsic functional trait scores, habitat-weighted functional loss, and multidimensional community functional diversity metrics.
* **Main Inputs:** `data/traits.xlsx`, `data/functional.xlsx`, and `data/Eco_data_Mamiferos_Carajas_2019_12.csv`.
* **Main Outputs:** Cleaned analytical datasets exported to `outputs/tables/`:
  * `valores_funcionais_multitraco.csv` (Species-level functional values and loss)
  * `diversidade_funcional_multitraco.csv` (Community-level FRic, FDis, and redundancy)
* **Steps Executed:** String and name standardization (`janitor::clean_names`), trait score matrix construction, presence-absence community matrix reshaping, and multidimensional trait space ordination (`FD::dbFD()`).
* **Quality Control:** String sanitization (trimming whitespace and lowercase conversion) prevents mismatches during dataset joins.
* **Important Notes:** Relies on the spatial area outputs derived from Module 1.

## 5. Analytical Workflow
The methodological framework links geographical distribution dynamics directly with functional ecology:
1. **Data Preparation & Cleaning:** Occurrence records are thinned spatially. Trait attributes (body mass, dispersal distance, diet, IUCN threat, foraging stratum, activity pattern) are normalized and cleaned using `janitor`.
2. **Quality Control & Variable Selection:** Environmental predictors undergo Pearson correlation testing (r < 0.75) to strip redundant bioclimatic variables before modeling.
3. **Distribution Modeling:** MaxEnt algorithms evaluate species geographic responses. Spatial partitioning dynamically adapts between block validation (n > 15) and jackknife cross-validation (n <= 15). Multi-GCM ensemble consensus maps are built for 2050 and 2070.
4. **Calculations & Functional Metrics:** Continuous suitability rasters are thresholded into binary presence/absence. Range contractions (in km²) are multiplied by baseline trait scores (`base_value`) to compute relative functional loss per ecosystem service.
5. **Statistical & Ordination Analyses:** Community-level presence-absence matrices combined with functional trait spaces undergo distance-based Principal Coordinates Analysis (PCoA) via `FD::dbFD()` to derive Functional Richness (FRic), Functional Dispersion (FDis), and Functional Redundancy across regional management units (Pará State vs. Conservation Units).

## 6. Reproducibility
* **Execution Order:** 1. `scripts/01_SDM_Maxent.R`
  2. `scripts/02_Functional_Diversity.R`
* **Script Dependencies:** Script `02_Functional_Diversity.R` depends on the habitat area estimations derived from the spatial projections generated by `01_SDM_Maxent.R`.
* **Prerequisite Files:** Raw occurrence records (`spps_occ_env_thin_f.csv`), trait tables (`traits.xlsx`, `functional.xlsx`), and external climate rasters must exist before execution.
* **Auto-generated Files:** Model RDS files, binary threshold lists, and final CSV tables (`valores_funcionais_multitraco.csv`, `diversidade_funcional_multitraco.csv`) are generated automatically upon script execution.
* **Environmental Data Instructions:** Download baseline WorldClim v2.1 bioclimatic layers (2.5m resolution) and future CMIP6 projections (2050 and 2070 for CCSM4, MIROC6, and MPI-ESM-LR under SSP245/SSP585) from [WorldClim.org](https://www.worldclim.org/) and place them in your local working directory.

## 7. Software Requirements
* **Recommended R Version:** $\ge$ 4.3.2
* **Spatial & SDM Packages:** `ENMwizard`, `ENMeval` (v0.3.1), `dismo`, `raster`, `terra`, `rgdal`, `rgeos`, `spThin`, `rJava`.
* **Functional & Data Analytics Packages:** `FD`, `tidyverse`, `dplyr`, `tidyr`, `readxl`, `readr`, `janitor`.
* **External Dependencies:** Java Runtime Environment (JRE) configured for `rJava`/Maxent execution.

## 8. Outputs
The computational pipeline generates the following analytical artifacts:
* **Spatial Models & Objects:** Saved RDS files containing tuned Maxent model lists (`ENMeval.res.lst.rds`), multi-GCM consensus suitability models (`mxnt.mdls.preds.cf_cons.rds`), and binary thresholded projections (`mods_thrshld_lst.rds`).
* **CSV Data Tables:**
  * `valores_funcionais_multitraco.csv`: Species-level base values, habitat areas, and relative functional loss.
  * `diversidade_funcional_multitraco.csv`: Community-level Functional Richness (FRic), Functional Dispersion (FDis), and Functional Redundancy across scenarios and regions.

## 9. Quality Control
* **Predictor Screening:** Automated correlation filtering (r < 0.75) removes multicollinear environmental variables.
* **Sample-Size Validation Rules:** Algorithmic evaluation method dynamically switches between spatial block partitioning (n > 15) and jackknife cross-validation (n <= 15) to avoid over-optimistic evaluation metrics.
* **String Standardization:** Automated string trimming and lowercase conversion (`janitor::clean_names`) prevent silent failures during relational database joins.
* **Inspection vs. Final Scripts:** Script `01_SDM_Maxent.R` includes intermediate graphical inspection steps (plotting cropped predictors and thresholded maps) before outputting final RDS objects. Script `02_Functional_Diversity.R` produces the definitive analytical CSV outputs.

## 10. Reproducibility Notes
* **Computational Demands:** MaxEnt tuning (`ENMevaluate_b`) and consensus projections (`consensus_scn_b`) are highly memory-intensive and require multi-core processing (`numCores = 8`).
* **Package Deprecations:** The pipeline incorporates legacy spatial libraries (`rgdal`, `rgeos`, `raster`) alongside `terra`. Users on newer R environments should verify system-level JRE and `rgdal` compatibility.
* **Linearity Assumption:** Functional service loss calculation assumes a proportional relationship between spatial range contraction and functional contribution.

## 11. Citation
If you use this repository, code, or analytical framework, please cite both the scientific manuscript and this source code:

> **Ribeiro Da Silva, R., Alves Ferreira, G., Hemming, N., & Marinho-Filho, J. S.** (2025). *Projected range contractions threaten mammal-driven ecosystem functions in a conservation mosaic of southeastern Amazonia*. SSRN Preprint. DOI: [10.2139/ssrn.5432746](https://doi.org/10.2139/ssrn.5432746)
> 
> **Ribeiro Da Silva, R., Alves Ferreira, G., Hemming, N., & Marinho-Filho, J. S.** (2025). Repository: Data and computational pipeline for mammals functional diversity in Carajás [Source Code]. GitHub. https://github.com/raquelribeiro/carajas-mammals-functional-loss
