
#output is readRDS of all the files in the directory that end with .rds
# Write a function that generates the predictions for all the files in the directory and saves the output as well
# has an option to save the plot of each species occupancy probability for the
# ecosystem type
baseRasterForPreds <- create_template_raster(
  study_region = norway_boundary,
  resolution = 100,
  crs = modelCrs,
  output_file = file.path(
    path_finder("dataFolder"),
    "processedData",
    "prediction_grids",
    "base_raster_for_predictions_100m.tif"
  ),
  overwrite = TRUE
)


results <- generate_occupancy_predictions_for_each_species(
  input_dir = normalizePath(path_to_results, winslash = "/", mustWork = FALSE),
  prediction_points = predictionPoints_for_modelling,
  det_covariates = det_covs,
  occ_covariates = occ_covs,
  base_raster = baseRasterForPreds,
  model_crs = modelCrs,
  ecosystem_type = ecosystemType,
  output_dir = normalizePath(path_to_predictions, winslash = "/", mustWork = FALSE),
  return = "raster",
  probs = c(0.025, 0.975),
  use_route_random_effect = TRUE,
  ignore_RE = TRUE,
  save_plots = TRUE,
  plot_width = 8,
  plot_height = 6,
  plot_dpi = 300
)

