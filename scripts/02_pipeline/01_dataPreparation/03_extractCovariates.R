########################################
# Extract covariates for the breeding birds data for the ecosystem
# and the prediction coordinates
#######################################

# Define covariate folder in a reproducible way.
# FAIR: avoids hard-coded absolute paths.
cov_folder <- file.path(
  path_finder("dataFolder"),
  "Covariates"
)

# Lets the load the protection regions
# The file path is from NINA's R drive
# R:\GeoSpatialData\ProtectedSites\Norway_ProtectedAreas\Original\Natur_Vern_norge_med_svalbard_25833.gdb
# We load and crop to the boundary of Norway


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
    "The following objects are missing from the global environment: ",
    paste(missing_objects, collapse = ", "),
    ". Please run the master script first.",
    call. = FALSE
  )
}

all_covariates <- unique(c(
  get("det_covs", envir = .GlobalEnv),
  get("occ_covs", envir = .GlobalEnv)
))

cov_folder <- file.path(
  path_finder("dataFolder"),
  "Covariates"
)


# Create template raster to ensure all covariates are on the same resolution
baseRaster <- create_template_raster(
  study_region = norway_boundary,
  resolution = 1000, # 1km for modelling. This needs to be changed maybe
  crs = modelCrs,
  output_file = file.path(
    path_finder("dataFolder"),
    "processedData",
    "prediction_grids",
    "base_raster.tif"
  ),
  overwrite = TRUE
)

# Load the stacked raster
# It automatically checks if the names and geometry matches the raster template,
# Else it will create a new raster.
cov_stack <- put_covariates_together(x = ecosystemBreedingBird,
                        modelCrs = modelCrs,
                        all_covariates = all_covariates,
                        template = baseRaster,
                        output_file = file.path(path_finder("dataFolder"),
                                                "processedData",
                                                "covariate_stack",
                                                "covariate_stack.tif"),
                        cov_folder = file.path(path_finder("dataFolder"), "Covariates"),
                        data_folder = path_finder("dataFolder"))


# Extract covariates for the Breeding bird data
breeding_bird_data_for_modelling <- extract_model_covariates(
    x = ecosystemBreedingBird %>% dplyr::select(-DEM),
    modelCrs,
    all_covariates,
    lon_col = "Longitude",
    lat_col = "Latitude",
    input_crs = 4326,
    covariate_stack = cov_stack,
    return = "data.frame"
) %>%
  dplyr::mutate(
    route_id  = sub("_.*", "", Route_point_id),
    point_id  = sub(".*_", "", Route_point_id)
  )


predictionPoints_for_modelling <- extract_model_covariates(
  x = predictionPoints %>% dplyr::select(-any_of(all_covariates)),
  modelCrs,
  all_covariates,
  lon_col = "Longitude",
  lat_col = "Latitude",
  input_crs = 4326,
  covariate_stack = cov_stack,
  return = "data.frame"
) %>%
  dplyr::mutate(landCover = factor(landCover, levels = levels(breeding_bird_data_for_modelling$landCover)))
