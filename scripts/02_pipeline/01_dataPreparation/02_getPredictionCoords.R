################################################################################
# Generate prediction regions and prediction points
################################################################################
#
# Purpose
# -------
# Create the spatial prediction regions and prediction locations used to
# estimate occupancy indicators for the selected ecosystem.
#
# Workflow
# --------
# 1. Load Norway administrative boundaries.
# 2. Create common raster templates for spatial processing.
# 3. Generate ecosystem-specific prediction locations:
#      - Forest: sample locations from selected GLAD land-cover classes.
#      - Mountain: sample locations from the Norwegian mountain layer and
#        remove locations classified as excluded land cover.
# 4. Save prediction locations for subsequent model prediction.
#
# Notes
# -----
# Several geospatial datasets required for modelling are distributed separately
# from this repository.
#
# NINA users:
#   Ensure that the R: drive / required NINA geospatial directories are
#   available before running this script.
#
# External users:
#   Download the required datasets according to the project documentation and
#   configure their locations in `set_paths.R`.
#
################################################################################

################################################################################
# 1. Load Norway administrative boundaries
################################################################################
#
# Administrative boundaries are downloaded when required and cached locally.
# The national boundary is used to define the study region and construct the
# prediction-grid templates.
#
################################################################################

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


# Transform the national and county boundaries to the model CRS.
norway_boundary <- norway_shape_files$national_boundary |>
  sf::st_transform(
    crs = modelCrs
  )

norway_county_boundary <- norway_shape_files$county_boundaries |>
  sf::st_transform(
    crs = modelCrs
  )


################################################################################
# 2. Create common raster templates
################################################################################
#
# Two raster templates are created:
#
#   baseRaster
#       100-m grid used for detailed spatial processing and covariate
#       preparation.
#
#   baseRaster_500m
#       500-m grid used for prediction locations and selected derived
#       covariates.
#
# All rasters used in the corresponding workflows should be aligned to the
# relevant template so that they have the same CRS, resolution, extent, and
# origin.
#
# The templates are cropped to the Norwegian national boundary.
#
################################################################################

# -------------------------------------------------------------------------
# 100-m raster template
# -------------------------------------------------------------------------

baseRaster <- create_template_raster(
  study_region = norway_boundary,
  resolution = 100,
  crs = modelCrs,
  output_file = file.path(
    dataFolder,
    "processedData",
    "prediction_grids",
    "base_raster.tif"
  ),
  overwrite = FALSE
)


# -------------------------------------------------------------------------
# 500-m raster template
# -------------------------------------------------------------------------

baseRaster_500m <- create_template_raster(
  study_region = norway_boundary,
  resolution = 500,
  crs = modelCrs,
  output_file = file.path(
    dataFolder,
    "processedData",
    "prediction_grids",
    "base_raster_500m.tif"
  ),
  overwrite = FALSE
)


################################################################################
# 3. Generate ecosystem-specific prediction locations
################################################################################
#
# Prediction locations are generated differently depending on the ecosystem
# being modelled.
#
#   Forest
#       Locations are sampled from selected GLAD land-cover classes.
#
#   Mountain
#       Locations are sampled from the Norwegian mountain-area layer and
#       locations classified as excluded land cover are removed.
#
################################################################################

if (ecosystemType == "Forest") {
  ##############################################################################
  # 3A. Forest ecosystem
  ##############################################################################
  #
  # Forest prediction locations are sampled from the GLAD land-cover raster.
  #
  # The eligible GLAD classes are defined in the project legend and selected
  # using the `For Forest mask - predictions` field.
  #
  # A regular sample of 100,000 locations is generated to provide broad spatial
  # coverage of the eligible forest area.
  #
  ##############################################################################

  message("Generating forest prediction locations...")

  # -------------------------------------------------------------------------
  # Read GLAD land-cover class legend
  # -------------------------------------------------------------------------

  legend <- readxl::read_excel(
    file.path(
      dataFolder,
      "Covariates",
      "Landcover_F_M",
      "legend.xlsx"
    ),
    sheet = "NO_binnedClasses"
  )

  # Identify GLAD classes eligible for forest predictions.
  classes <- legend$`lc_GLAD.tif value`[
    legend$`For Forest mask - predictions` == "T"
  ]

  classes <- classes[
    !is.na(classes)
  ]

  # -------------------------------------------------------------------------
  # Prepare GLAD land-cover raster and sample prediction locations
  # -------------------------------------------------------------------------

  predictionPoints <- get_GLAD_landcover(
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
  ) |>

    # Align GLAD land cover to the 100-m prediction grid.
    resample_to_baseRaster(
      r = _,
      baseRaster = baseRaster,
      name = "glad_landCover"
    ) |>

    # Sample prediction locations from eligible forest classes.
    sample_prediction_points_from_raster(
      raster = _,
      selected_codes = classes,
      n_points = 100000,
      output_file = file.path(
        dataFolder,
        "processedData",
        "prediction_regions",
        "forest_region.parquet"
      ),
      metadata_file = NULL,
      overwrite = TRUE,
      seed = 123,
      sample_method = "regular"
    ) |>

    # Only the spatial geometry is required for subsequent prediction.
    dplyr::select(
      geometry
    )

  message(
    "Forest prediction locations generated: ",
    nrow(predictionPoints)
  )
} else if (ecosystemType == "Mountain") {
  ##############################################################################
  # 3B. Mountain ecosystem
  ##############################################################################
  #
  # Mountain prediction locations are generated from the Norwegian mountain
  # layer (`fjell_grid50`).
  #
  # A regular sample of 150,000 locations is generated. These locations are
  # subsequently filtered using the GLAD land-cover raster to remove locations
  # classified as land-cover class 241.
  #
  ##############################################################################

  message("Generating mountain prediction locations...")

  # -------------------------------------------------------------------------
  # Define output file
  # -------------------------------------------------------------------------

  mountain_output <- file.path(
    dataFolder,
    "processedData",
    "prediction_regions",
    "mountain_regions.shp"
  )

  # -------------------------------------------------------------------------
  # Generate mountain prediction locations if they do not already exist
  # -------------------------------------------------------------------------

  if (!file.exists(mountain_output)) {
    # Path to the Norwegian mountain-area geodatabase.
    input_file <- file.path(
      geoSpatFolder,
      "LandCover",
      "Norway_MountainAreas",
      "Original",
      "Fjell_norge_publisering.gdb"
    )

    # -----------------------------------------------------------------------
    # Read and prepare mountain layer
    # -----------------------------------------------------------------------

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
    ) |>

      # Generate a regular spatial sample.
      sf::st_sample(
        size = 150000,
        type = "regular"
      ) |>

      # Convert the sampled geometries back to an sf object.
      sf::st_as_sf() |>

      # Rename the geometry column.
      dplyr::rename(
        geometry = x
      )

    # -----------------------------------------------------------------------
    # Load and prepare GLAD land cover
    # -----------------------------------------------------------------------

    landCoverRast <- get_GLAD_landcover(
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
    ) |>

      # Align GLAD land cover with the 100-m processing grid.
      resample_to_baseRaster(
        r = _,
        baseRaster = baseRaster,
        name = "glad_landCover"
      )

    # -----------------------------------------------------------------------
    # Extract land cover at mountain prediction locations
    # -----------------------------------------------------------------------

    # Extract GLAD land-cover values at prediction locations
    land_cover_values <- terra::extract(
      landCoverRast,
      terra::vect(predictionPoints)
    ) |>
      dplyr::select(
        -ID
      )

    # Combine extracted land-cover values with prediction points
    predictionPoints <- dplyr::bind_cols(
      predictionPoints,
      land_cover_values
    )

    # Remove locations belonging to excluded land-cover class
    predictionPoints <- predictionPoints |>
      dplyr::filter(
        glad_landCover != 241
      ) |>
      dplyr::select(
        geometry
      )

    # Remove intermediate object
    rm(land_cover_values)
    gc()

    # -----------------------------------------------------------------------
    # Save prediction locations
    # -----------------------------------------------------------------------

    sf::st_write(
      predictionPoints,
      mountain_output,
      delete_dsn = TRUE,
      quiet = TRUE
    )

    message(
      "Mountain prediction locations saved to: ",
      mountain_output
    )
  } else {
    # -----------------------------------------------------------------------
    # Load previously generated prediction locations
    # -----------------------------------------------------------------------

    message(
      "Reading existing mountain prediction locations..."
    )

    predictionPoints <- sf::st_read(
      mountain_output,
      quiet = TRUE
    )
  }
} else {
  ##############################################################################
  # Invalid ecosystem type
  ##############################################################################

  stop(
    "Unknown ecosystemType: '",
    ecosystemType,
    "'. Expected 'Forest' or 'Mountain'.",
    call. = FALSE
  )
}

################################################################################
# End of script
################################################################################
#
# Object returned by this workflow:
#
#   predictionPoints
#
# This is an `sf` POINT object containing the spatial locations used for
# subsequent model prediction.
#
################################################################################
