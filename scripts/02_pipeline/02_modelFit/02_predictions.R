################################################################################
# Generate species occupancy predictions
#
# This script:
#   1. Creates a fine-resolution raster template (100 m).
#   2. Loads all fitted species models from the model output directory.
#   3. Predicts occupancy probabilities at the prediction locations.
#   4. Interpolates predictions onto a raster grid.
#   5. Saves prediction rasters and optional diagnostic plots.
#
# Outputs:
#   - Mean occupancy probability rasters for each species.
#   - Lower and upper credible interval rasters (2.5% and 97.5%).
#   - Species-specific occupancy maps (optional).
################################################################################


################################################################################
# Create raster template for prediction maps
################################################################################

# Predictions are generated on a finer spatial resolution (100 m)
# than that used for model fitting (typically 1 km). This raster serves
# as a common grid onto which species occupancy probabilities are projected.

if(ecosystemType %in% "Mountain"){

  input_file <- file.path(
    geoSpatFolder,
    "LandCover",
    "Norway_MountainAreas",
    "Original",
    "Fjell_norge_publisering.gdb"
  )

  pred_boundary <- get_mountain_layer(
    input_file = input_file,
    layer = "fjell_grid50",
    norway_boundary = norway_boundary,
    modelCrs = modelCrs,
    output_file = file.path(
      dataFolder,
      "processedData",
      "prediction_regions",
      "mountain.parquet"
    )
  )

} else if (ecosystemType %in% "Forest"){

  legend <- readxl::read_excel(paste0(dataFolder,"/Covariates/Landcover_F_M/legend.xlsx"),
                               sheet = "NO_binnedClasses")
  classes <- legend$`lc_GLAD.tif value`[legend$`For Forest mask - predictions` == "T"]
  classes <- classes[!is.na(classes)]

  pred_boundary <- get_GLAD_landcover(
    input_file = file.path(
      dataFolder,
      "Covariates",
      "Landcover_F_M",
      "lc2020nor_mergedClasses.tif"
    ),
    model_crs = modelCrs,
    output_file = file.path(
      dataFolder,
      "processedData",
      "GLAD_lc",
      "LCglad.tif"
    )
  ) %>%
    resample_to_baseRaster(
      r = .,
      baseRaster = baseRaster,
      name = "glad_landCover"
    )
}

baseRasterForPreds <- create_template_raster(
  study_region = pred_boundary,
  resolution   = pred_resolution,     # 100-m prediction grid
  crs          = modelCrs,
  output_file  = file.path(
    dataFolder,
    "processedData",
    "prediction_grids",
    paste0("base_raster_for_predictions_",ecosystemType,"_100m.tif")
  ),
  overwrite = TRUE
)


################################################################################
# Generate occupancy predictions for all fitted species
################################################################################

# The function automatically:
#
#   - Reads all model objects (*.rds) in 'path_to_results'.
#   - Generates occupancy predictions for every species.
#   - Computes posterior summaries and credible intervals.
#   - Saves prediction rasters to 'path_to_predictions'.
#   - Optionally creates species-specific occupancy maps.
#
# Predictions are generated using the fixed effects only
# (ignore_RE = TRUE), which removes route-level random effects and
# provides smoother landscape-scale occupancy estimates.




prediction_results <- generate_occupancy_predictions_for_each_species_with_blocks(

  ##########################################################################
  # Input model directory
  ##########################################################################

  input_dir = path_to_results,


  ##########################################################################
  # Prediction data
  ##########################################################################

  prediction_points = predictionPoints_for_modelling,

  det_covariates = det_covs,
  occ_covariates = occ_covs,


  ##########################################################################
  # Spatial information
  ##########################################################################

  base_raster = baseRasterForPreds,
  model_crs   = modelCrs,
  ecosystem_type = ecosystemType,


  ##########################################################################
  # Output directory
  ##########################################################################

  output_dir = path_to_predictions,


  ##########################################################################
  # Output type
  ##########################################################################

  # Return predictions as raster layers
  return = "raster", #c("sf", "raster"), # options: sf or raster
  block_size = 10000,

  # Posterior credible intervals
  probs = c(0.025, 0.975),


  ##########################################################################
  # Random effects
  ##########################################################################

  # Route random effects were included during model fitting
  use_route_random_effect = TRUE,

  # Ignore route random effects during prediction
  # to obtain population-level occupancy estimates
  ignore_RE = TRUE,


  ##########################################################################
  # Plotting options
  ##########################################################################

  save_plots = TRUE,

  plot_width  = 8,
  plot_height = 6,
  plot_dpi    = 300
)



# prediction_results <- generate_occupancy_predictions_for_each_species(
#
#   ##########################################################################
#   # Input model directory
#   ##########################################################################
#
#   input_dir = path_to_results,
#
#
#   ##########################################################################
#   # Prediction data
#   ##########################################################################
#
#   prediction_points = predictionPoints_for_modelling[1:10000, ],
#
#   det_covariates = det_covs,
#   occ_covariates = occ_covs,
#
#
#   ##########################################################################
#   # Spatial information
#   ##########################################################################
#
#   base_raster = baseRasterForPreds,
#   model_crs   = modelCrs,
#   ecosystem_type = ecosystemType,
#
#
#   ##########################################################################
#   # Output directory
#   ##########################################################################
#
#   output_dir = path_to_predictions,
#
#
#   ##########################################################################
#   # Output type
#   ##########################################################################
#
#   # Return predictions as raster layers
#   return = "sf", # options: sf or raster
#
#   # Posterior credible intervals
#   probs = c(0.025, 0.975),
#
#
#   ##########################################################################
#   # Random effects
#   ##########################################################################
#
#   # Route random effects were included during model fitting
#   use_route_random_effect = TRUE,
#
#   # Ignore route random effects during prediction
#   # to obtain population-level occupancy estimates
#   ignore_RE = TRUE,
#
#
#   ##########################################################################
#   # Plotting options
#   ##########################################################################
#
#   save_plots = TRUE,
#
#   plot_width  = 8,
#   plot_height = 6,
#   plot_dpi    = 300
# )
