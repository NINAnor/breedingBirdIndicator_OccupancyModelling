################################################################################
# Master script: breeding bird occupancy modelling
#
# This script controls the full modelling workflow:
#   1. Configure user-specific paths.
#   2. Check that required folders exist.
#   3. Load project functions.
#   4. Set analysis options.
#   5. Create output folders.
#   6. Run data preparation, model fitting, and optional prediction scripts.
#
# Before running:
#   - Edit scripts/01_setPaths.R for your local folder structure.
#   - Make sure required data folders are available.
################################################################################

################################################################################
# 1. User and path configuration
################################################################################

USER <- "kwaku" # options: all or kwaku or server
# If you are using the OneDrive, specify all. But if you are running on
# your own local disk, specify youur woth USER in the 01_setPaths.R script
if (USER %in% "server") {
  args <- commandArgs(trailingOnly = TRUE)
  print(args)
  tpr <- args[1]
  print(tpr)
}

source(here::here("scripts", "01_setPaths.R"))

dataFolder <- path_finder("dataFolder", user = USER)
geoSpatFolder <- path_finder("ninaGeospatialFiles", user = USER)
expertFolder <- path_finder("expertData", user = USER)
output_root <- normalizePath(
  path_finder("results", user = USER),
  mustWork = FALSE
)


# ==============================================================================
# 3. LOAD PACKAGES AND SET TERRA OPTIONS
# ==============================================================================

library(terra)
library(lme4)
library(emmeans)
library(ggplot2)
library(dplyr)

if (USER %in% "server") {
  terraOptions(
    memfrac = 0.5, # Fraction of available RAM
    memmax = 8, # Maximum memory allocation in GB
    threads = 4, # Number of processing threads
    parallel = TRUE,
    todisk = TRUE, # Prefer disk-based processing
    tempdir = "/data/scratch/kwaku",
    datatype = "FLT4S",
    progress = 3
  )
} else {
  terraOptions(
    memfrac = 0.5, # Fraction of available RAM
    memmax = 8, # Maximum memory allocation in GB
    threads = 4, # Number of processing threads
    parallel = TRUE,
    todisk = TRUE, # Prefer disk-based processing
    tempdir = "C:/terra_tmp",
    datatype = "FLT4S",
    progress = 3
  )
}

################################################################################
# 2. Check required folders
################################################################################

if (!dir.exists(dataFolder)) {
  stop(
    paste(
      "Project data directory not found:",
      dataFolder,
      "",
      "The data required to reproduce the analyses are not included in this repository.",
      "Please request access from the corresponding authors:",
      "",
      "    kwaku.adjei@nina.no",
      "",
      "After obtaining the files, place them in the project data directory",
      "and update the 'dataFolder' path if necessary.",
      sep = "\n"
    ),
    call. = FALSE
  )
}

# if (!dir.exists(geoSpatFolder)) {
#   stop(
#     paste(
#       "Geospatial files directory not found:",
#       geoSpatFolder,
#       "",
#       "NINA employees:",
#       "- Check that the R: drive is connected.",
#       "",
#       "External users:",
#       "- Download the raster and shape files described in the documentation.",
#       "- Store them locally.",
#       "- Update 'ninaGeospatialFiles' in scripts/01_setPaths.R.",
#       sep = "\n"
#     ),
#     call. = FALSE
#   )
# }

if (!dir.exists(expertFolder)) {
  stop(
    paste(
      "Expert information directory not found:",
      expertFolder,
      "",
      "The expert species lists and grouping files are required.",
      "",
      "Please check that the expert data folder exists and that",
      "'expertData' is correctly set in scripts/01_setPaths.R.",
      sep = "\n"
    ),
    call. = FALSE
  )
}


################################################################################
# 3. Locale and function loading
################################################################################

# Ensure Norwegian characters are read correctly in expert input files.
Sys.setlocale("LC_CTYPE", "Norwegian_Norway.utf8")

# Load all project functions.
function_files <- list.files(
  path_finder("functions", user = USER),
  pattern = "\\.R$",
  recursive = TRUE,
  full.names = TRUE
)

invisible(lapply(function_files, source))


################################################################################
# 4. Analysis settings
################################################################################

# Ecosystem to model.
ecosystemType <- "Mountain" # Options: "Forest", "Mountain"

# Version label used in output folder names.
# For testing, I use version 1 for route ID effect, and version 2 withour ruote ID
version <- 2
version_label <- paste0("version_0", version)
use_route_random_effect <- TRUE
if (version == 2) {
  use_route_random_effect <- FALSE
}


# Pipeline switches.
fitModel <- FALSE
runPredictions <- TRUE
runIndicators <- FALSE
runPostprocessing <- FALSE

# Model type switches.
fit_single_species <- FALSE # TRUE: PGOcc / spPGOcc per species
if (USER %in% "server") {
  fit_single_species <- tpr
}

fit_spatial_model <- FALSE # TRUE: spatial NNGP model

# Coordinate reference system used throughout the modelling workflow.
modelCrs <- 25833

# Years included in the model.
modelYear <- c(2020, 2024)

# Temporary directory.
TEMPDIR <- tempdir()


# Setting the resolution
pred_resolution <- 500 # in metres¨
aggregation_spacing <- 1000

# Which column of predictions should be used to estimate indicators
indicator_col <- "mean"

# How many groups should be used in grouping the species
maxSpeciesInGrps <- 24

modelChecks <- FALSE

################################################################################
# 5. Model covariates
################################################################################

covariateSelection <- readxl::read_excel(file.path(
  dataFolder,
  "Covariates",
  "covariatesCodingForModelling.xlsx"
))

if (ecosystemType %in% "Forest") {
  # Detection covariates affect the observation process.
  det_covs <- covariateSelection$covariate[
    covariateSelection$forest_det == TRUE
  ]

  # Occupancy covariates affect the ecological state process.
  occ_covs <- covariateSelection$covariate[
    covariateSelection$forest_occ == TRUE
  ]
} else if (ecosystemType %in% "Mountain") {
  # Detection covariates affect the observation process.
  det_covs <- covariateSelection$covariate[
    covariateSelection$mountain_det == TRUE
  ]

  # Occupancy covariates affect the ecological state process.
  occ_covs <- covariateSelection$covariate[
    covariateSelection$mountain_occ == TRUE
  ]
}


################################################################################
# 6. MCMC settings
################################################################################

n.batch <- 40
batch.length <- 500 #2000
n.burn <- 10000 #70000
n.thin <- 50 # use 10 for single species.
n.chains <- 3

n.neighbors <- 10
n.omp.threads <- 8


################################################################################
# 7. Output paths
################################################################################

species_dir <- if (fit_single_species) {
  "single_species_model"
} else {
  "multiple_species_model"
}

spatial_dir <- if (fit_spatial_model) {
  "spatialModel"
} else {
  "non_spatialModel"
}

path_to_results <- file.path(
  output_root,
  "modelFit",
  version_label,
  ecosystemType,
  species_dir,
  spatial_dir
)

path_to_predictions <- file.path(
  output_root,
  "prediction_outputs",
  version_label,
  ecosystemType,
  species_dir,
  spatial_dir
)

path_to_indicators <- file.path(
  output_root,
  "indicators_outputs",
  version_label,
  ecosystemType,
  species_dir,
  spatial_dir
)

# if runPostprocessing, specify whether you are interested in comparing all outputs from
# the various models within the same ecosystem or that you are just interested in postprocessing
# for the current configurations
compareMultipleOutputs <- FALSE

if (!compareMultipleOutputs) {
  path_to_postProcessing <- file.path(
    output_root,
    "postProcessing_outputs",
    version_label,
    ecosystemType,
    species_dir,
    spatial_dir
  )
} else {
  path_to_postProcessing <- file.path(
    output_root,
    "postProcessing_outputs",
    version_label,
    ecosystemType,
    "multipleOutputs"
  )
}

path_for_logs <- file.path(
  output_root,
  "logs",
  version_label,
  ecosystemType
)

dir.create(path_to_results, recursive = TRUE, showWarnings = FALSE)
dir.create(path_to_predictions, recursive = TRUE, showWarnings = FALSE)
dir.create(path_to_indicators, recursive = TRUE, showWarnings = FALSE)
dir.create(path_for_logs, recursive = TRUE, showWarnings = FALSE)
dir.create(path_to_postProcessing, recursive = TRUE, showWarnings = FALSE)

logfile <- file.path(
  path_for_logs,
  paste0("model_run_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
)


################################################################################
# 8. Run pipeline with logging
################################################################################

log_con <- file(logfile, open = "wt")

sink(log_con, split = TRUE)
sink(log_con, type = "message")

tryCatch(
  {
    cat("Started:", as.character(Sys.time()), "\n\n")

    cat("Selected ecosystem:\n")
    print(ecosystemType)

    cat("\nSelected years:\n")
    print(modelYear)

    cat("\nModel settings:\n")
    cat("Single-species model:", fit_single_species, "\n")
    cat("Spatial model:", fit_spatial_model, "\n")
    cat("Fit model:", fitModel, "\n")
    cat("Run predictions:", runPredictions, "\n")
    cat("Run indicators:", runIndicators, "\n\n")

    ##########################################################################
    # Data preparation
    ##########################################################################

    source(here::here(
      "scripts",
      "02_pipeline",
      "01_dataPreparation",
      "01_getBreedingBirdData.R"
    ))

    source(here::here(
      "scripts",
      "02_pipeline",
      "01_dataPreparation",
      "02_getPredictionCoords.R"
    ))

    source(here::here(
      "scripts",
      "02_pipeline",
      "01_dataPreparation",
      "03_extractCovariates.R"
    ))

    ##########################################################################
    # Model fitting
    ##########################################################################

    if (fitModel) {
      source(here::here(
        "scripts",
        "02_pipeline",
        "02_modelFit",
        "01_modelFit.R"
      ))
    }

    # Rest for 10 minutes
    Sys.sleep(1 * 60)

    ##########################################################################
    # Prediction
    ##########################################################################

    if (runPredictions) {
      source(here::here(
        "scripts",
        "02_pipeline",
        "02_modelFit",
        "02_predictions.R"
      ))
    }

    ##########################################################################
    # Indicator estimation
    ##########################################################################
    # Rest for 10 minutes
    Sys.sleep(1 * 60)

    if (runIndicators) {
      source(here::here(
        "scripts",
        "02_pipeline",
        "02_modelFit",
        "03_estimating_indicators.R"
      ))
    }

    if (runPostprocessing) {
      source(here::here(
        "scripts",
        "02_pipeline",
        "03_postProcessing",
        "01_process_one_path_for_results.R"
      ))
    }

    cat("\nFinished:", as.character(Sys.time()), "\n")
  },

  error = function(e) {
    cat("\nERROR:\n")
    cat(conditionMessage(e), "\n")
    stop(e)
  },

  finally = {
    sink(type = "message")
    sink()
    close(log_con)

    message("Log written to: ", logfile)
  }
)

files_to_remove <- list.files(
  TEMPDIR,
  full.names = TRUE,
  all.files = TRUE,
  no.. = TRUE
)

if (length(files_to_remove) > 0) {
  unlink(
    files_to_remove,
    recursive = TRUE,
    force = TRUE
  )
}
