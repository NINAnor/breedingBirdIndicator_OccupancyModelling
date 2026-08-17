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
  data_folder    = dataFolder,
  norway_boundary = norway_boundary,
  modelYear = modelYear,
  overwrite = TRUE
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

# --------------------------------------------------------
# Apply categorical transformations only when requested
# --------------------------------------------------------

if ("dominantTree" %in% all_covariates) {

  breeding_bird_data_for_modelling <-
    breeding_bird_data_for_modelling %>%
    dplyr::mutate(
      dominantTree = dplyr::if_else(
        is.na(dominantTree),
        4,
        dominantTree
      ),
      dominantTree = factor(
        dominantTree
      )
    )
}

if ("glad_landcover" %in% all_covariates) {

  breeding_bird_data_for_modelling <-
    breeding_bird_data_for_modelling %>%
    dplyr::mutate(
      glad_landcover = factor(
        glad_landCover
      )
    )
}

################################################################################
# Extract covariates at prediction locations
################################################################################

# Environmental values are extracted at prediction points so that occupancy
# probabilities can be estimated across the landscape.

# Environmental values are extracted at prediction points so that
# occupancy probabilities can be estimated across the landscape.

predictionPoints_for_modelling <-
  extract_model_covariates(
    x               = predictionPoints %>%
      dplyr::select(-dplyr::any_of(all_covariates)),
    modelCrs        = modelCrs,
    all_covariates  = all_covariates,
    lon_col         = "Longitude",
    lat_col         = "Latitude",
    input_crs       = 4326,
    covariate_stack = cov_stack,
    return          = "data.frame"
  )

# --------------------------------------------------------
# Apply dominantTree transformation only if requested
# --------------------------------------------------------

if ("dominantTree" %in% all_covariates) {

  predictionPoints_for_modelling <-
    predictionPoints_for_modelling %>%
    dplyr::mutate(
      dominantTree = dplyr::if_else(
        is.na(dominantTree),
        4,
        dominantTree
      ),
      dominantTree = factor(
        dominantTree,
        levels = levels(
          breeding_bird_data_for_modelling$dominantTree
        )
      )
    )
}

# --------------------------------------------------------
# Apply GLAD land-cover transformation only if requested
# --------------------------------------------------------

if ("glad_landcover" %in% all_covariates) {

  predictionPoints_for_modelling <-
    predictionPoints_for_modelling %>%
    dplyr::mutate(
      glad_landcover = factor(
        glad_landCover,
        levels = levels(
          breeding_bird_data_for_modelling$glad_landcover
        )
      )
    )
}
