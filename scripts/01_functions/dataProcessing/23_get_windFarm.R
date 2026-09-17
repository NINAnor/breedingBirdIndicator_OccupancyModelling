#' Create Wind Farm Distance Raster
#'
#' Reads wind farm features from the Norwegian wind power geodatabase,
#' calculates the distance of each 500-m prediction-grid cell to the
#' nearest wind farm boundary, assigns a signed distance according to
#' whether the location occurs within a wind farm, and rasterizes the
#' result.
#'
#' Negative values indicate locations within wind farm areas, whereas
#' positive values indicate locations outside wind farm areas.
#' Distance is expressed in metres.
#'
#' If the processed raster already exists and `overwrite = FALSE`, the
#' cached raster is read from disk instead of being recomputed.
#'
#' This function supports FAIR principles:
#' \itemize{
#'   \item Findable: output and metadata use descriptive filenames.
#'   \item Accessible: output is saved as GeoTIFF.
#'   \item Interoperable: GeoTIFF is readable by R, QGIS, GDAL, and Python.
#'   \item Reusable: metadata records the source geodatabase, source layer,
#'         CRS, resolution, extent, processing method, and processing date.
#' }
#'
#' @param input_dir Character. Path to the directory containing the
#'   `Vindkraft.gdb` file geodatabase.
#' @param processing_raster Character. Path to the 500-m processing raster
#'   used as the template for the output.
#' @param output_file Character. Path to the output wind farm distance
#'   raster, e.g. `"windFarm.tif"`.
#' @param target_crs Character or numeric. Target CRS for processing.
#'   Defaults to `"EPSG:25833"`.
#' @param metadata_file Character or `NULL`. Path to YAML metadata file.
#'   If `NULL`, metadata is written beside `output_file`.
#' @param gdb_name Character. Name of the wind power file geodatabase.
#'   Defaults to `"Vindkraft.gdb"`.
#' @param windfarm_layer Character. Name of the wind farm layer in the
#'   geodatabase. Defaults to `"VindkraftanleggF"`.
#' @param overwrite Logical. If `TRUE`, recompute the raster even if
#'   `output_file` already exists.
#'
#' @return A `terra::SpatRaster` containing signed wind farm distance
#'   in metres.
#'
#' @export
get_WindFarm <- function(
    input_dir,
    processing_raster,
    output_file,
    target_crs = "EPSG:25833",
    metadata_file = NULL,
    gdb_name = "Vindkraft.gdb",
    windfarm_layer = "VindkraftanleggF",
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
  # Wind power geodatabase
  # --------------------------------------------------------

  gdb <- file.path(
    input_dir,
    gdb_name
  )

  if (!dir.exists(gdb)) {

    stop(
      "Wind power geodatabase not found: ",
      gdb,
      call. = FALSE
    )
  }

  # --------------------------------------------------------
  # Read cached output
  # --------------------------------------------------------

  if (
    file.exists(output_file) &&
    !overwrite
  ) {

    message(
      "Reading cached wind farm distance raster: ",
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
  # Read wind farm layer
  # --------------------------------------------------------

  message(
    "Reading wind farm layer..."
  )

  windFarm <- sf::st_read(
    dsn = gdb,
    layer = windfarm_layer,
    quiet = TRUE
  )

  # --------------------------------------------------------
  # Project wind farm features
  # --------------------------------------------------------

  message(
    "Projecting wind farm features to ",
    target_crs,
    "..."
  )

  windFarm <- sf::st_transform(
    windFarm,
    crs = target_crs
  )

  # --------------------------------------------------------
  # Convert wind farm polygons/features to lines
  # --------------------------------------------------------

  message(
    "Converting wind farm features to boundary lines..."
  )

  windFarm_line <- sf::st_cast(
    windFarm,
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
  # Find nearest wind farm feature
  # --------------------------------------------------------

  message(
    "Finding nearest wind farm feature for ",
    nrow(pts_sf),
    " prediction locations..."
  )

  nearest_windFarm <- sf::st_nearest_feature(
    pts_sf,
    windFarm_line
  )

  # --------------------------------------------------------
  # Calculate distance to nearest wind farm
  # --------------------------------------------------------

  message(
    "Calculating distance to nearest wind farm..."
  )

  pts_sf$dist_windFarmline <- as.numeric(
    sf::st_distance(
      pts_sf,
      windFarm_line[
        nearest_windFarm,
        ,
        drop = FALSE
      ],
      by_element = TRUE
    )
  )

  # --------------------------------------------------------
  # Determine whether points occur within wind farms
  # --------------------------------------------------------

  message(
    "Determining whether locations occur within wind farms..."
  )

  pts_sf$WF_in <- lengths(
    sf::st_intersects(
      pts_sf,
      windFarm
    )
  ) > 0

  # --------------------------------------------------------
  # Calculate signed wind farm distance
  # --------------------------------------------------------

  # Negative = inside wind farm
  # Positive = outside wind farm

  pts_sf$WF_dist <- ifelse(
    pts_sf$WF_in,
    -pts_sf$dist_windFarmline,
    pts_sf$dist_windFarmline
  )

  # --------------------------------------------------------
  # Rasterize signed distance
  # --------------------------------------------------------

  message(
    "Rasterizing wind farm distance..."
  )

  WF_dist_rast <- terra::rasterize(
    terra::vect(pts_sf),
    processing_rast,
    field = "WF_dist",
    fun = "mean",
    na.rm = TRUE
  )

  names(WF_dist_rast) <- "WF_dist"

  # --------------------------------------------------------
  # Write output raster
  # --------------------------------------------------------

  message(
    "Writing wind farm distance raster: ",
    output_file
  )

  terra::writeRaster(
    WF_dist_rast,
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

    title = "Signed wind farm distance",

    description = paste(
      "Signed distance to the nearest wind farm feature.",
      "Negative values indicate locations within wind farm",
      "areas, whereas positive values indicate locations",
      "outside wind farm areas. Distance is expressed in metres."
    ),

    created = as.character(
      Sys.time()
    ),

    source = list(

      dataset = "Norway Windpower",

      provider = "Norwegian governmental geospatial data",

      source_geodatabase = gdb,

      source_layer = windfarm_layer,

      source_version = "versjon2022"
    ),

    output = list(

      file = output_file,

      file_format = "GeoTIFF",

      variable = "WF_dist",

      units = "metres",

      data_type = "continuous",

      signed_distance = TRUE,

      negative_values = "Inside wind farm",

      positive_values = "Outside wind farm"
    ),

    processing_method = list(

      prediction_grid = processing_raster,

      nearest_feature = "Nearest wind farm boundary",

      nearest_feature_function =
        "sf::st_nearest_feature()",

      distance_function =
        "sf::st_distance()",

      rasterization_function =
        "terra::rasterize()",

      rasterization_summary = "mean",

      target_crs = target_crs
    ),

    raster_properties = list(

      output_crs = terra::crs(
        WF_dist_rast,
        describe = TRUE
      ),

      extent = as.character(
        terra::ext(WF_dist_rast)
      ),

      resolution = as.numeric(
        terra::res(WF_dist_rast)
      ),

      dimensions = c(
        nrow = terra::nrow(WF_dist_rast),
        ncol = terra::ncol(WF_dist_rast),
        nlyr = terra::nlyr(WF_dist_rast)
      ),

      names = names(
        WF_dist_rast
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
        "Metadata records source geodatabase, source layer,",
        "processing method, CRS, resolution, extent,",
        "dimensions, units, and processing date."
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
        "Wind farm raster was saved, but metadata could ",
        "not be written: ",
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
    nearest_windFarm,
    windFarm_line,
    windFarm
  )

  gc()

  # --------------------------------------------------------
  # Return raster
  # --------------------------------------------------------

  WF_dist_rast
}
