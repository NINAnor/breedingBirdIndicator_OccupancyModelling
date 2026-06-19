# Breeding Bird Occupancy Modelling for Ecosystem Indicators

## Overview

This repository contains code for fitting hierarchical occupancy models to Norwegian breeding bird data and generating ecosystem indicators based on predicted species occupancy.

The workflow integrates:

* Breeding bird point count data.
* Environmental covariates.
* Expert-defined species lists for ecosystems.
* Hierarchical occupancy models using **spOccupancy**.
* Spatial prediction of species occupancy.
* Estimation of ecosystem-level richness indicators.

Currently, the workflow supports:

* **Forest ecosystems**
* **Mountain ecosystems**

---

## Repository structure

```text
.
├── scripts/
│   ├── 00_masterScript.R
│   ├── 01_setPaths.R
│   └── 01_functions/
│       ├── 01_loadPackages.R
│       ├── dataProcessing/
│       └── modelFitting/
│   └── 02_pipeline/
│       ├── 01_dataPreparation/
│           ├── 01_getBreedingBirdData.R/
│           ├── 02_getPredictionCoords.R/
│           ├── 03_extractCovariates.R/
│       ├── 02_modelFit/
│           ├── 01_modelFit.R/
│           ├── 02_predictions.R/
│           ├── 03_estimating_indicators.R/
├── vignettes/
├── README.md
└── .Rproj
```

---

## Workflow

The analysis proceeds in four stages:

### 1. Data preparation

Breeding bird observations are loaded and filtered to the selected ecosystem and years.

Environmental covariates are extracted for:

* Survey locations
* Prediction locations

Prediction regions are generated for each ecosystem type.

---

### 2. Occupancy model fitting

Species occupancy models are fitted using **spOccupancy**.

Supported models include:

* Single-species occupancy models (`PGOcc`)
* Spatial occupancy models (`spPGOcc`)
* Multi-species occupancy models (optional)

Models account for:

* Detection probability
* Environmental covariates
* Route-level random effects
* Spatial autocorrelation (optional)

Check the vignette on the model fitting for the description of the model.

---

### 3. Prediction

Species occupancy probabilities are predicted across Norway.

Outputs include:

* Mean occupancy probability
* 95% credible intervals
* Species-specific raster maps
* GeoParquet prediction files

---

### 4. Indicator estimation

Predicted occupancies are aggregated to estimate ecosystem indicators.

Indicators are summarized by:

* Protected areas (reference condition)
* Biogeographic regions

Current indicators include:

* Scaled ecosystem richness

---

# Scripts

## Master script

```text
scripts/masterScript.R
```

Controls the complete workflow.

Responsible for:

* Loading user paths
* Defining model settings
* Running the pipeline
* Logging model runs

---

## Path configuration

```text
scripts/01_setPaths.R
```

Contains machine-specific paths to:

* Data folders
* Geospatial layers
* Output directories

Users should modify this file before running the workflow.

---

## Data preparation

### 01_getBreedingBirdData.R

Loads breeding bird observations from:

* SQL Server database (NINA users)
* Cached local files

Filters observations to the selected ecosystem and years.

---

### 02_getPredictionCoords.R

Creates prediction locations used during model prediction.

Forest ecosystems:

* Derived from land-cover classes.

Mountain ecosystems:

* Derived from mountain polygons.

Also loads:

* Norway boundaries
* County boundaries

---

### 03_extractCovariates.R

Creates a common raster template.

Builds a covariate stack and extracts environmental variables for:

* Survey points
* Prediction points

---

## Model fitting

### 01_modelFit.R

Fits occupancy models using **spOccupancy**.

Supports:

* Single-species models
* Spatial models
* Multi-species models

Outputs:

* RDS files containing fitted models

---

### 02_predictions.R

Generates occupancy predictions for all species.

Produces:

* Mean occupancy probability maps
* Lower and upper credible interval maps
* Species plots

Outputs are stored as GeoParquet files.

---

### 03_estimating_indicators.R

Loads:

* Biogeographic regions
* Protected areas

Calculates ecosystem indicators from predicted occupancies.

Outputs:

* Scaled richness indicators

---

# Input data

The repository itself does not contain all required data.

External data include:

### Breeding bird data

Survey observations used to fit occupancy models.

### Environmental covariates

Examples include:

* Elevation (DEM)
* Land cover
* Climate variables

### Expert species lists

Species defining each ecosystem.

### Geospatial layers

Including:

* Norway boundaries
* Protected areas
* Biogeographic regions
* Mountain polygons

Some datasets are available only internally at NINA.

External users should consult the project documentation or contact the corresponding authors.

---

# Outputs

Model outputs are organized as:

```text
results/
├── modelFit/
├── prediction_outputs/
├── indicators/
└── logs/
```

### modelFit

Contains fitted occupancy models (`.rds`).

### prediction_outputs

Contains:

* Occupancy prediction rasters
* GeoParquet files
* Figures

### indicators

Contains ecosystem indicators and summaries.

### logs

Contains logs from each model run.

---

# Reproducibility

The project uses:

* Relative paths via `here`
* Cached intermediate files
* GeoParquet outputs
* Logging of model runs

The workflow is designed to be reproducible across Windows, Linux, and macOS.

---

# Authors

Kwaku Adjei
Norwegian Institute for Nature Research (NINA)

## Correspondence

For questions regarding data access or model implementation, please contact:

**[kwaku.adjei@nina.no](mailto:kwaku.adjei@nina.no)**

---

# Citation

If using this repository or its outputs, please cite:

> Adjei et al. Breeding Bird Occupancy Modelling for Ecosystem Indicators. Norwegian Institute for Nature Research (NINA).
