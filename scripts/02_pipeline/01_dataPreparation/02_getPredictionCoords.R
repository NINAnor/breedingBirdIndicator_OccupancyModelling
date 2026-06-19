# This script is used to generate the prediction coordinates as well as the
# regions needed to estimate the indicators

#Lets first download or load in the norway boundary
norway_shape_files <- get_norway_boundaries(
    out_dir = file.path(path_finder("dataFolder"),
                        "processedData",
                        "norway_boundaries"),
    temp_dir = TEMPDIR,
    source = "GADM",
    overwrite = FALSE
)

norway_boundary <- norway_shape_files$national_boundary %>%
  sf::st_transform(modelCrs)

norway_county_boundary <- norway_shape_files$county_boundaries %>%
  sf::st_transform(modelCrs)



#######################
## Get forest regions
########################

# Define covariate folder in a reproducible way.
# FAIR: avoids hard-coded absolute paths.
if(ecosystemType %in% "Forest"){
predictionPoints <- get_landcover_labels(ecosystemBreedingBird,
                                      modelCrs = modelCrs,
                                      cov_folder = file.path(
                                        path_finder("dataFolder"),
                                        "Covariates"
                                      ),
                                      output = "raster") %>%
  sample_prediction_points_from_raster(
    .,
    selected_codes = c(311,312, 313),
    n_points = 10000,
    output_file = file.path(path_finder("dataFolder"),
                            "processedData",
                            "prediction_regions",
                            "forest_region.parquet"),
    metadata_file = NULL,
    overwrite = FALSE,
    seed = 123,
    sample_method = c("regular")
  ) %>%
  dplyr::select(geometry)
} else if(ecosystemType %in% "Mountain"){

  input_file <- file.path(path_finder("ninaGeospatialFiles"),
                          "LandCover",
                          "Norway_MountainAreas",
                          "Original",
                          "Fjell_norge_publisering.gdb")

  predictionPoints <- get_mountain_layer(
    input_file = input_file,
    layer = "fjell_grid50",
    norway_boundary = norway_boundary,
    modelCrs = modelCrs,
    output_file = file.path(path_finder("dataFolder"),
                            "processedData",
                            "prediction_regions",
                            "mountain.parquet") )%>%
                  sf::st_sample(mountain,
                              size = 10000,
                              type = "regular") |>
      sf::st_as_sf() %>%
    dplyr::rename("geometry" = "x")

}

