# Master script

# Set directories

## Go here and edit the path for various data sources
# usethis::edit_r_environ()
## Alternatively we can also do this by creating a path for the various files.
## I will use the function, so that we can all use the same functions

## Load the project for the code to work seemlesly.

#First edit 01_setPaths.R to set the file paths for the scripts
# and folders needed to run this script. This load a path_finder script that returns
# a path of what you are looking for.
source("scripts/01_setPaths.R")

# Ensure we can read Norwegian text, especially from expert data input
Sys.setlocale("LC_CTYPE", "Norwegian_Norway.utf8")

# Load the functions
function_files <- list.files(
  file.path(path_finder("functions")),
  pattern = "\\.R$",
  recursive = TRUE,
  full.names = TRUE
)

invisible(
  lapply(function_files, source)
)


# Select the ecosystem of interest
ecosystemType <- "Forest" # "Mountain"
version <- 1
version_label <- paste0("version_0", version)


fitModel <- FALSE
runPredictions <- TRUE

# Path to results
path_to_results <- file.path(path_finder("results"),
                             "modelFit",
                             version_label)

if(!dir.exists(path_to_results)) dir.create(path_to_results,
                                            recursive = TRUE)


# create temporary directory
TEMPDIR <- tempdir()

# set model crs
modelCrs <- 25833

# Set the years to use data in the models
# Use the last seven full years

modelYear <- c(2021, 2025)

path_for_logs <- file.path(
  path_finder("results"),
  "logs",
  ecosystemType,
  version_label
)

dir.create(
  path_for_logs,
  recursive = TRUE,
  showWarnings = FALSE
)

logfile <- file.path(
  path_for_logs,
  paste0(
    "model_run_",
    format(Sys.time(), "%Y%m%d_%H%M%S"),
    ".log"
  )
)



# ============================================================
# 8. Run model with logfile
# ============================================================

log_con <- file(logfile, open = "wt")

sink(log_con, split = TRUE)
sink(log_con, type = "message")

tryCatch(

  {
cat("Started:", as.character(Sys.time()), "\n\n")
cat("Selected Ecosystem:\n")
print(ecosystemType)

# Specify the covariates for occupancy and detection
# This will be done automatically by Katherine in an excel sheet that
# will be called to provide the covariates here
det_covs <- c("DEM", "landCover")
occ_covs <- c("DEM", "landCover")

cat("Selected years:\n")
print(modelYear)

# Script for loading breeding bird data
source("scripts/02_pipeline/01_dataPreparation/01_getBreedingBirdData.R")

# Script for loading prediction for the ecosystem type
source("scripts/02_pipeline/01_dataPreparation/02_getPredictionCoords.R")

# Script for extracting covariates
source("scripts/02_pipeline/01_dataPreparation/03_extractCovariates.R")

# Run script to fit model for ecosystem type
if(fitModel){
source("scripts/02_pipeline/02_modelFit/01_modelFit.R")
}

# Run script to predict for ecosystem type
if(runPredictions){
  source("scripts/02_pipeline/02_modelFit/02_predictions.R")
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
}
)

