################################################################################
# Generate prediction regions and prediction points
#
# This script creates the spatial regions used to estimate occupancy indicators
# and generates a set of prediction locations for each ecosystem type.
#
# Workflow:
#   1. Verify that required geospatial files are available.
#   2. Download or load Norway administrative boundaries.
#   3. Create prediction regions for the selected ecosystem.
#   4. Sample prediction points used during model prediction.
################################################################################


################################################################################
# Verify availability of geospatial files
################################################################################

# Several layers required for modelling (land cover, mountain regions, etc.)
# are distributed separately from this repository.
#
# NINA employees should ensure that the R: drive is mounted.
# External users should download the required files following the project
# documentation and set the corresponding path in set_paths.R.




################################################################################
# Load Norway administrative boundaries
################################################################################

# Download boundaries if needed and cache them locally.
# These boundaries are later used to clip and aggregate predictions.

norway_shape_files <- get_norway_boundaries(
  out_dir = file.path(
    dataFolder,
    "processedData",
    "norway_boundaries"
  ),
  temp_dir = TEMPDIR,
  source = "GADM",
  overwrite = FALSE
)

# Transform boundaries to the coordinate reference system used by the model.

norway_boundary <- norway_shape_files$national_boundary %>%
  sf::st_transform(modelCrs)

norway_county_boundary <- norway_shape_files$county_boundaries %>%
  sf::st_transform(modelCrs)



################################################################################
# Generate prediction regions
################################################################################

# Prediction points are generated differently for forest and mountain ecosystems.

if (ecosystemType == "Forest") {

  ##############################################################################
  # Forest ecosystem
  #
  # Forest prediction points are sampled from land-cover classes:
  #   311 = Broad-leaved forest
  #   312 = Coniferous forest
  #   313 = Mixed forest
  #
  # A regular sample of 10,000 points is generated and saved for reuse.
  ##############################################################################

  predictionPoints <- get_landcover_labels(
    ecosystemBreedingBird,
    modelCrs = modelCrs,
    cov_folder = file.path(
      dataFolder,
      "Covariates"
    ),
    output = "raster"
  ) %>%

    sample_prediction_points_from_raster(
      selected_codes = c(311, 312, 313),
      n_points = 10000,
      output_file = file.path(
        dataFolder,
        "processedData",
        "prediction_regions",
        "forest_region.parquet"
      ),
      metadata_file = NULL,
      overwrite = FALSE,
      seed = 123,
      sample_method = "regular"
    ) %>%

    dplyr::select(geometry)


} else if (ecosystemType == "Mountain") {

  ##############################################################################
  # Mountain ecosystem
  #
  # Prediction points are generated from the Norwegian mountain layer
  # (fjell_grid50) and a regular sample of 10,000 locations is drawn.
  ##############################################################################

  input_file <- file.path(
    geoSpatFolder,
    "LandCover",
    "Norway_MountainAreas",
    "Original",
    "Fjell_norge_publisering.gdb"
  )

  predictionPoints <- get_mountain_layer(
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
  ) %>%

    sf::st_sample(
      size = 10000,
      type = "regular"
    ) %>%

    sf::st_as_sf() %>%

    dplyr::rename(geometry = x)

}
