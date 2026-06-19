################################################################################
# Prepare environmental covariates for occupancy modelling
#
# This script:
#   1. Loads all environmental covariates required by the model.
#   2. Creates a common raster template to ensure consistent resolution and CRS.
#   3. Builds a covariate raster stack.
#   4. Extracts covariate values at bird survey locations.
#   5. Extracts covariate values at prediction locations.
#
# All covariates are aligned to the same grid before extraction to ensure
# consistency between model fitting and prediction.
################################################################################


################################################################################
# Define paths
################################################################################

# Folder containing all raster covariates used in the model.
# Keeping this path centralized improves portability and reproducibility.
cov_folder <- file.path(
  dataFolder,
  "Covariates"
)


################################################################################
# Check that model covariates have been defined
################################################################################

# Detection and occupancy covariates should have been specified in the master
# script. Stop execution if they are unavailable.

required_objects <- c("det_covs", "occ_covs")

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1),
    envir = .GlobalEnv
  )
]

if (length(missing_objects) > 0) {

  stop(
    paste0(
      "The following objects are missing from the global environment:\n",
      paste(missing_objects, collapse = ", "),
      "\n\nPlease run the master script before executing this section."
    ),
    call. = FALSE
  )

}


################################################################################
# Combine detection and occupancy covariates
################################################################################

# Create a unique list of all environmental variables needed for the model.

if (!exists("det_covs") || !exists("occ_covs")) {
  stop(
    "Objects 'det_covs' and 'occ_covs' must be defined. ",
    "Please run the master script first.",
    call. = FALSE
  )
}
all_covariates <- unique(c(det_covs, occ_covs))



################################################################################
# Create a common raster template
################################################################################
# Lets the load the protection regions # The file path is from NINA's R drive
# R:\GeoSpatialData\ProtectedSites\Norway_ProtectedAreas\Original\Natur_Vern_norge_med_svalbard_25833.gdb
# We load and crop to the boundary of Norway
# All covariates are resampled onto a common 1 km grid to ensure that
# resolution, extent, and coordinate reference system are identical.

baseRaster <- create_template_raster(
  study_region = norway_boundary,
  resolution   = 1000,       # 1 km resolution
  crs          = modelCrs,
  output_file  = file.path(
    dataFolder,
    "processedData",
    "prediction_grids",
    "base_raster.tif"
  ),
  overwrite = TRUE
)


################################################################################
# Build the covariate raster stack
################################################################################

# Load existing raster layers or recreate them if their geometry or variable
# names do not match the template raster.

cov_stack <- put_covariates_together(
  x              = ecosystemBreedingBird,
  modelCrs       = modelCrs,
  all_covariates = all_covariates,
  template       = baseRaster,
  output_file    = file.path(
    dataFolder,
    "processedData",
    "covariate_stack",
    "covariate_stack.tif"
  ),
  cov_folder     = cov_folder,
  data_folder    = dataFolder
)


################################################################################
# Extract covariates at bird survey locations
################################################################################

# Extract environmental values at each survey point.
# These data will be used for model fitting.

breeding_bird_data_for_modelling <-
  extract_model_covariates(
    x               = ecosystemBreedingBird %>%
      dplyr::select(-DEM),
    modelCrs        = modelCrs,
    all_covariates  = all_covariates,
    lon_col         = "Longitude",
    lat_col         = "Latitude",
    input_crs       = 4326,
    covariate_stack = cov_stack,
    return          = "data.frame"
  ) %>%

  # Split Route_point_id into route and point identifiers
  dplyr::mutate(
    route_id = sub("_.*", "", Route_point_id),
    point_id = sub(".*_", "", Route_point_id)
  )


################################################################################
# Extract covariates at prediction locations
################################################################################

# Environmental values are extracted at prediction points so that occupancy
# probabilities can be estimated across the landscape.

predictionPoints_for_modelling <-
  extract_model_covariates(
    x               = predictionPoints %>%
      dplyr::select(-any_of(all_covariates)),
    modelCrs        = modelCrs,
    all_covariates  = all_covariates,
    lon_col         = "Longitude",
    lat_col         = "Latitude",
    input_crs       = 4326,
    covariate_stack = cov_stack,
    return          = "data.frame"
  ) %>%

  # Ensure factor levels match those used during model fitting
  dplyr::mutate(
    landCover = factor(
      landCover,
      levels = levels(
        breeding_bird_data_for_modelling$landCover
      )
    )
  )
