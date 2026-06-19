
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
  input_dir = file.path(
    path_finder("results"),
    "modelFit",
    version_label
  ),
  prediction_points = predictionPoints_for_modelling,
  det_covariates = det_covs,
  occ_covariates = occ_covs,
  base_raster = baseRasterForPreds,
  model_crs = modelCrs,
  ecosystem_type = ecosystemType,
  output_dir = file.path(
    path_finder("results"),
    "prediction_outputs"
  ),
  return = "sf",
  probs = c(0.025, 0.975),
  use_route_random_effect = FALSE,
  ignore_RE = TRUE,
  save_plots = TRUE,
  plot_width = 8,
  plot_height = 6,
  plot_dpi = 300
)

#
# model <- output$model_products$model
#
# det_covariates <- det_covs
# occ_covariates <- occ_covs
# use_route_random_effect <- FALSE
#
# make_rhs <- function(covariates, use_route_random_effect = FALSE) {
#   rhs <- paste(covariates, collapse = " + ")
#
#   if (use_route_random_effect) {
#     rhs <- paste(rhs, "+ (1 | route_id)")
#   }
#
#   rhs
# }
#
# # -------------------------------------------------------------------------
# # Formulas
# # -------------------------------------------------------------------------
#
# occ_rhs <- make_rhs(
#   covariates = occ_covariates,
#   use_route_random_effect = use_route_random_effect
# )
#
# det_rhs <- make_rhs(
#   covariates = det_covariates,
#   use_route_random_effect = FALSE
# )
#
# occ_formula <- as.formula(paste("~", occ_rhs))
# det_formula <- as.formula(paste("~", det_rhs))
#
# X.0 <- model.matrix(
#   occ_formula,
#   data = predictionPoints_for_modelling
# )
#
# X.0 <- X.0[, model$x.names, drop = FALSE]
#
# X.0 <- model.matrix(
#   reformulate(occ_terms),
#   data = predictionPoints_for_modelling
# )
#
# X.0 <- X.0[, model$x.names, drop = FALSE]
#
# coords.0 <- predictionPoints_for_modelling %>%
#   dplyr::select(longitude, latitude)
#
#
# pred <- predict(
#   model,
#   X.0 = X.0,
#   coords.0 = coords.0,
#   ignore.RE = TRUE,
#   verbose = TRUE
# )
#
#
# results <- summarize_occupancy_predictions(
#     pred,
#     coords = coords.0,
#     species_names = model$sp.names,
#     template = baseRaster,
#     return = c("sf"), #"raster", "sf"
#     crs = modelCrs,
#     probs = c(0.025, 0.975),
#     output_dir = file.path(path_finder("results"),
#                            "prediction_outputs",
#                            ecosystemType)
# )
#
# library(ggplot2)
#
# ggplot(
#   results %>%
#     filter(species == "Sylvia atricapilla")
# ) +
#   geom_sf(
#     aes(color = mean),
#     size = 0.15,
#     alpha = 0.8
#   ) +
#   scale_color_viridis_c(
#     name = "Occupancy\nprobability",
#     limits = c(0, 1)
#   ) +
#   coord_sf() +
#   theme_bw()
