#' Create Treeline Distance Raster
#'
#' Reads mountain and forest polygons from the Skoggrense geodatabase,
#' derives the distance of each 500-m prediction-grid cell to the nearest
#' mountain/treeline feature, assigns a signed distance according to whether
#' the location occurs within a mountain area, and rasterizes the result.
#'
#' Positive values indicate locations above/within the mountain area,
#' whereas negative values indicate locations below/outside the mountain
#' area. Distance is measured in metres.
#'
#' If the processed raster already exists and `overwrite = FALSE`, the
#' cached raster is read from disk instead of being recomputed.
#'
#' This function supports FAIR principles:
#' \itemize{
#'   \item Findable: output and metadata use descriptive filenames.
#'   \item Accessible: output is saved as GeoTIFF.
#'   \item Interoperable: GeoTIFF is readable by R, QGIS, GDAL, and Python.
#'   \item Reusable: metadata records the source geodatabase, layers,
#'         CRS, resolution, extent, processing method, and processing date.
#' }
#'
#' @param input_dir Character. Path to the Skoggrense file geodatabase,
#'   e.g. `"Skoggrense.gdb"`.
#' @param processing_raster Character. Path to the 500-m processing raster
#'   used as the template for the output.
#' @param output_file Character. Path to the output treeline-distance
#'   raster, e.g. `"treeline.tif"`.
#' @param target_crs Character or numeric. Target CRS for processing.
#'   Defaults to `"EPSG:25833"`.
#' @param metadata_file Character or `NULL`. Path to YAML metadata file.
#'   If `NULL`, metadata is written beside `output_file`.
#' @param mountain_layer Character. Name of the mountain polygon layer
#'   in the geodatabase. Defaults to `"mountain_areas"`.
#' @param forest_layer Character. Name of the forest polygon layer
#'   in the geodatabase. Defaults to `"forest_areas"`.
#' @param overwrite Logical. If `TRUE`, recompute the raster even if
#'   `output_file` already exists.
#'
#' @return A `terra::SpatRaster` containing signed treeline distance
#'   in metres.
#'
#' @export
get_Treeline <- function(
    input_dir,
    processing_raster,
    output_file,
    target_crs = "EPSG:25833",
    metadata_file = NULL,
    mountain_layer = "mountain_areas",
    forest_layer = "forest_areas",
    overwrite = FALSE
) {

  required_pkgs <- c(
    "terra",
    "sf",
    "yaml"
  )

  missing_pkgs <- required_pkgs[
    !vapply(
      required_pkgs,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  if (length(missing_pkgs) > 0) {

    stop(
      "Install missing packages first: ",
      paste(
        missing_pkgs,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  # --------------------------------------------------------
  # Normalise paths
  # --------------------------------------------------------

  input_dir <- normalizePath(
    input_dir,
    winslash = "/",
    mustWork = TRUE
  )

  processing_raster <- normalizePath(
    processing_raster,
    winslash = "/",
    mustWork = TRUE
  )

  output_file <- normalizePath(
    output_file,
    winslash = "/",
    mustWork = FALSE
  )

  # --------------------------------------------------------
  # Read cached output
  # --------------------------------------------------------

  if (
    file.exists(output_file) &&
    !overwrite
  ) {

    message(
      "Reading cached treeline distance raster: ",
      output_file
    )

    return(
      terra::rast(output_file)
    )
  }

  # --------------------------------------------------------
  # Create output directory
  # --------------------------------------------------------

  dir.create(
    dirname(output_file),
    recursive = TRUE,
    showWarnings = FALSE
  )

  # --------------------------------------------------------
  # Target CRS
  # --------------------------------------------------------

  if (is.numeric(target_crs)) {

    target_crs <- paste0(
      "EPSG:",
      target_crs
    )
  }

  target_crs <- as.character(
    target_crs
  )

  # --------------------------------------------------------
  # Read processing raster
  # --------------------------------------------------------

  message(
    "Reading 500-m processing raster..."
  )

  processing_rast <- terra::rast(
    processing_raster
  )

  # --------------------------------------------------------
  # Read mountain and forest layers
  # --------------------------------------------------------

  message(
    "Reading mountain areas..."
  )

  mountain <- sf::st_read(
    dsn = input_dir,
    layer = mountain_layer,
    quiet = TRUE
  )

  message(
    "Reading forest areas..."
  )

  forest <- sf::st_read(
    dsn = input_dir,
    layer = forest_layer,
    quiet = TRUE
  )

  # --------------------------------------------------------
  # Project spatial data
  # --------------------------------------------------------

  message(
    "Projecting treeline layers to ",
    target_crs,
    "..."
  )

  mountain <- sf::st_transform(
    mountain,
    crs = target_crs
  )

  forest <- sf::st_transform(
    forest,
    crs = target_crs
  )

  # --------------------------------------------------------
  # Convert mountain polygons to lines
  # --------------------------------------------------------

  message(
    "Converting mountain areas to boundary lines..."
  )

  mountain_line <- sf::st_cast(
    mountain,
    "MULTILINESTRING"
  )

  # --------------------------------------------------------
  # Create prediction points
  # --------------------------------------------------------

  message(
    "Creating prediction points from the 500-m raster..."
  )

  pts_sf <- terra::as.data.frame(
    processing_rast,
    xy = TRUE,
    na.rm = TRUE
  ) |>
    sf::st_as_sf(
      coords = c("x", "y"),
      crs = terra::crs(processing_rast)
    ) |>
    sf::st_transform(
      crs = target_crs
    )

  # --------------------------------------------------------
  # Find nearest mountain line
  # --------------------------------------------------------

  message(
    "Finding nearest mountain line for ",
    nrow(pts_sf),
    " prediction locations..."
  )

  nearest_mountain <- sf::st_nearest_feature(
    pts_sf,
    mountain_line
  )

  # --------------------------------------------------------
  # Calculate distance to mountain line
  # --------------------------------------------------------

  message(
    "Calculating distance to nearest mountain line..."
  )

  pts_sf$dist_mountainline <- as.numeric(
    sf::st_distance(
      pts_sf,
      mountain_line[
        nearest_mountain,
        ,
        drop = FALSE
      ],
      by_element = TRUE
    )
  )

  # --------------------------------------------------------
  # Determine whether points occur within mountain areas
  # --------------------------------------------------------

  message(
    "Determining treeline side..."
  )

  pts_sf$TL_above <- lengths(
    sf::st_intersects(
      pts_sf,
      mountain
    )
  ) > 0

  # --------------------------------------------------------
  # Calculate signed treeline distance
  # --------------------------------------------------------

  pts_sf$TL_dist <- ifelse(
    pts_sf$TL_above,
    pts_sf$dist_mountainline,
    -pts_sf$dist_mountainline
  )

  # --------------------------------------------------------
  # Rasterize signed treeline distance
  # --------------------------------------------------------

  message(
    "Rasterizing treeline distance..."
  )

  TL_dist_rast <- terra::rasterize(
    terra::vect(pts_sf),
    processing_rast,
    field = "TL_dist",
    fun = "mean",
    na.rm = TRUE
  )

  names(TL_dist_rast) <- "TL_dist"

  # --------------------------------------------------------
  # Write output raster
  # --------------------------------------------------------

  message(
    "Writing treeline distance raster: ",
    output_file
  )

  terra::writeRaster(
    TL_dist_rast,
    output_file,
    datatype = "FLT4S",
    overwrite = TRUE
  )

  # --------------------------------------------------------
  # Metadata path
  # --------------------------------------------------------

  if (is.null(metadata_file)) {

    metadata_file <- file.path(
      dirname(output_file),
      paste0(
        tools::file_path_sans_ext(
          basename(output_file)
        ),
        "_metadata.yml"
      )
    )
  }

  metadata_file <- normalizePath(
    metadata_file,
    winslash = "/",
    mustWork = FALSE
  )

  dir.create(
    dirname(metadata_file),
    recursive = TRUE,
    showWarnings = FALSE
  )

  # --------------------------------------------------------
  # Metadata
  # --------------------------------------------------------

  metadata <- list(

    title = "Signed treeline distance",

    description = paste(
      "Signed distance to the nearest mountain/treeline boundary.",
      "Positive values indicate locations within mountain areas,",
      "whereas negative values indicate locations outside mountain",
      "areas. Distance is expressed in metres."
    ),

    created = as.character(
      Sys.time()
    ),

    source = list(

      dataset = "Skoggrense",

      source_directory = input_dir,

      mountain_layer = mountain_layer,

      forest_layer = forest_layer,

      variables = c(
        "mountain areas",
        "forest areas"
      )
    ),

    output = list(

      file = output_file,

      file_format = "GeoTIFF",

      variable = "TL_dist",

      units = "metres",

      data_type = "continuous",

      signed_distance = TRUE
    ),

    processing_method = list(

      prediction_grid = processing_raster,

      nearest_feature = paste(
        "Nearest mountain area boundary"
      ),

      distance_function = "sf::st_distance()",

      nearest_feature_function = "sf::st_nearest_feature()",

      rasterization_function = "terra::rasterize()",

      rasterization_summary = "mean",

      positive_values = paste(
        "Locations intersecting mountain_areas"
      ),

      negative_values = paste(
        "Locations outside mountain_areas"
      ),

      target_crs = target_crs
    ),

    raster_properties = list(

      output_crs = terra::crs(
        TL_dist_rast,
        describe = TRUE
      ),

      extent = as.character(
        terra::ext(TL_dist_rast)
      ),

      resolution = as.numeric(
        terra::res(TL_dist_rast)
      ),

      dimensions = c(
        nrow = terra::nrow(TL_dist_rast),
        ncol = terra::ncol(TL_dist_rast),
        nlyr = terra::nlyr(TL_dist_rast)
      ),

      names = names(
        TL_dist_rast
      )
    ),

    FAIR = list(

      Findable = paste(
        "Descriptive output and metadata filenames are used."
      ),

      Accessible = paste(
        "Processed raster is stored as GeoTIFF."
      ),

      Interoperable = paste(
        "GeoTIFF is compatible with R terra, QGIS, GDAL,",
        "and Python rasterio."
      ),

      Reusable = paste(
        "Metadata records source layers, processing method,",
        "CRS, resolution, extent, dimensions, units, and",
        "processing date."
      )
    )
  )

  # --------------------------------------------------------
  # Write metadata
  # --------------------------------------------------------

  tryCatch(

    yaml::write_yaml(
      metadata,
      metadata_file
    ),

    error = function(e) {

      warning(
        "Treeline raster was saved, but metadata could not ",
        "be written: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  message(
    "Metadata saved to: ",
    metadata_file
  )

  # --------------------------------------------------------
  # Clean up large intermediate objects
  # --------------------------------------------------------

  rm(
    pts_sf,
    nearest_mountain,
    mountain_line
  )

  gc()

  # --------------------------------------------------------
  # Return raster
  # --------------------------------------------------------

  TL_dist_rast
}
