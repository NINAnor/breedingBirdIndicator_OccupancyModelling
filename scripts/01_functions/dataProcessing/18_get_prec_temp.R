#' Calculate Mean Breeding-Season Temperature and Precipitation
#'
#' Reads daily temperature and precipitation layers from annual NetCDF files,
#' extracts the requested breeding-season days, calculates the mean across
#' the breeding season, and saves the resulting rasters as GeoTIFF files.
#'
#' If the output rasters already exist and `overwrite = FALSE`, the cached
#' rasters are read from disk instead of being recomputed.
#'
#' This function supports FAIR principles:
#' \itemize{
#'   \item Findable: output and metadata use descriptive filenames.
#'   \item Accessible: outputs are saved as GeoTIFF.
#'   \item Interoperable: GeoTIFF is readable by R, QGIS, GDAL, and Python.
#'   \item Reusable: metadata records source files, year, selected DOYs,
#'         CRS, resolution, extent, processing date, and methods.
#' }
#'
#' @param nc_dir Character. Directory containing annual NetCDF files.
#' @param output_dir Character. Directory where processed rasters are saved.
#' @param year Integer. Year to process.
#' @param nc_pattern Character. Filename pattern containing one `%s` or
#'   equivalent year placeholder, e.g. `"seNorge_%s.nc"`.
#' @param tg_sub Character or integer. NetCDF subdataset/layer specification
#'   for temperature.
#' @param rr_sub Character or integer. NetCDF subdataset/layer specification
#'   for precipitation.
#' @param breeding_season Integer vector. Day-of-year values defining the
#'   breeding season.
#' @param raster_crs Character. CRS assigned to the input rasters,
#'   e.g. `"EPSG:25833"`.
#' @param metadata_file Character or `NULL`. Path to YAML metadata file.
#'   If `NULL`, metadata is written inside `output_dir`.
#' @param overwrite Logical. If `TRUE`, recompute outputs even if they exist.
#'
#' @return A named list containing two `terra::SpatRaster` objects:
#'   `temperature` and `precipitation`.
#'
#' @export
get_prec_temp <- function(
    nc_dir,
    output_dir,
    year,
    nc_pattern,
    tg_sub,
    rr_sub,
    breeding_season,
    raster_crs = "EPSG:25833",
    metadata_file = NULL,
    overwrite = FALSE
) {

  required_pkgs <- c("terra", "yaml")

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
      paste(missing_pkgs, collapse = ", "),
      call. = FALSE
    )
  }

  # --------------------------------------------------------
  # Validate inputs
  # --------------------------------------------------------

  if (
    length(year) != 1L ||
    is.na(year)
  ) {
    stop(
      "`year` must be a single valid year.",
      call. = FALSE
    )
  }

  year <- as.integer(year)

  if (length(breeding_season) == 0L) {
    stop(
      "`breeding_season` must contain at least one day-of-year.",
      call. = FALSE
    )
  }

  breeding_season <- as.integer(
    breeding_season
  )

  nc_dir <- normalizePath(
    nc_dir,
    winslash = "/",
    mustWork = TRUE
  )

  output_dir <- normalizePath(
    output_dir,
    winslash = "/",
    mustWork = FALSE
  )

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  # --------------------------------------------------------
  # Output files
  # --------------------------------------------------------

  temp_output <- file.path(
    output_dir,
    paste0(
      "mean_temperature_breeding_season_",
      year,
      ".tif"
    )
  )

  prec_output <- file.path(
    output_dir,
    paste0(
      "mean_precipitation_breeding_season_",
      year,
      ".tif"
    )
  )

  # --------------------------------------------------------
  # Read cached outputs
  # --------------------------------------------------------

  if (
    file.exists(temp_output) &&
    file.exists(prec_output) &&
    !overwrite
  ) {

    message(
      "Reading cached temperature and precipitation rasters for ",
      year
    )

    return(
      list(
        temperature = terra::rast(
          temp_output
        ),
        precipitation = terra::rast(
          prec_output
        )
      )
    )
  }

  # --------------------------------------------------------
  # Input NetCDF file
  # --------------------------------------------------------

  nc_file <- file.path(
    nc_dir,
    sprintf(
      nc_pattern,
      year
    )
  )

  if (!file.exists(nc_file)) {

    stop(
      "NetCDF file not found: ",
      nc_file,
      call. = FALSE
    )
  }

  message(
    "Processing: ",
    year
  )

  message(
    "Input file: ",
    nc_file
  )

  # --------------------------------------------------------
  # Read temperature
  # --------------------------------------------------------

  message(
    "Reading temperature layers..."
  )

  tg_rast <- terra::rast(
    nc_file,
    sub = tg_sub
  )

  # --------------------------------------------------------
  # Read precipitation
  # --------------------------------------------------------

  message(
    "Reading precipitation layers..."
  )

  rr_rast <- terra::rast(
    nc_file,
    sub = rr_sub
  )

  # --------------------------------------------------------
  # Assign CRS
  # --------------------------------------------------------

  terra::crs(tg_rast) <- raster_crs
  terra::crs(rr_rast) <- raster_crs

  # --------------------------------------------------------
  # Number of available days
  # --------------------------------------------------------

  n_days_tg <- terra::nlyr(tg_rast)
  n_days_rr <- terra::nlyr(rr_rast)

  if (n_days_tg != n_days_rr) {
    stop(
      "Temperature and precipitation have different numbers ",
      "of layers: temperature = ",
      n_days_tg,
      ", precipitation = ",
      n_days_rr,
      ".",
      call. = FALSE
    )
  }

  n_days <- n_days_tg

  valid_doys <- breeding_season[
    breeding_season >= 1 &
      breeding_season <= n_days
  ]

  if (length(valid_doys) == 0L) {
    stop(
      "No valid breeding-season layers available for ",
      year,
      ".",
      call. = FALSE
    )
  }

  if (length(valid_doys) < length(breeding_season)) {
    warning(
      length(breeding_season) -
        length(valid_doys),
      " requested breeding-season DOYs are outside ",
      "the available layers for ",
      year,
      ".",
      call. = FALSE
    )
  }

  # Select breeding-season layers
  tg_season <- tg_rast[[valid_doys]]
  rr_season <- rr_rast[[valid_doys]]

  # Calculate breeding-season means
  tg_mean <- terra::mean(
    tg_season,
    na.rm = TRUE
  )

  rr_mean <- terra::mean(
    rr_season,
    na.rm = TRUE
  )

  names(tg_mean) <- paste0(
    "mean_temperature_",
    year
  )

  names(rr_mean) <- paste0(
    "mean_precipitation_",
    year
  )

  # --------------------------------------------------------
  # Write output rasters
  # --------------------------------------------------------

  message(
    "Writing temperature raster..."
  )

  terra::writeRaster(
    tg_mean,
    temp_output,
    overwrite = TRUE
  )

  message(
    "Writing precipitation raster..."
  )

  terra::writeRaster(
    rr_mean,
    prec_output,
    overwrite = TRUE
  )

  # --------------------------------------------------------
  # Metadata
  # --------------------------------------------------------

  if (is.null(metadata_file)) {

    metadata_file <- file.path(
      output_dir,
      paste0(
        "prec_temp_",
        year,
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

  metadata <- list(

    title = paste(
      "Mean breeding-season temperature and precipitation",
      year
    ),

    description = paste(
      "Annual mean temperature and precipitation calculated",
      "across the specified breeding-season day-of-year layers."
    ),

    created = as.character(
      Sys.time()
    ),

    input = list(

      directory = nc_dir,

      file = nc_file,

      year = year,

      temperature_subdataset = tg_sub,

      precipitation_subdataset = rr_sub,

      available_days = n_days,

      requested_breeding_season_DOY =
        breeding_season,

      valid_breeding_season_DOY =
        valid_doys
    ),

    output = list(

      temperature = list(
        file = temp_output,
        variable = "mean temperature",
        unit = "as provided by source dataset",
        file_format = "GeoTIFF"
      ),

      precipitation = list(
        file = prec_output,
        variable = "mean precipitation",
        unit = "as provided by source dataset",
        file_format = "GeoTIFF"
      )
    ),

    processing_method = list(

      "function" = "terra::mean()",

      aggregation = paste(
        "Arithmetic mean across the valid breeding-season",
        "daily raster layers."
      ),

      na_rm = TRUE,

      CRS = raster_crs
    ),

    raster_properties = list(

      temperature = list(
        crs = terra::crs(
          tg_mean,
          describe = TRUE
        ),
        extent = as.character(
          terra::ext(tg_mean)
        ),
        resolution = as.numeric(
          terra::res(tg_mean)
        ),
        dimensions = c(
          nrow = terra::nrow(tg_mean),
          ncol = terra::ncol(tg_mean),
          nlyr = terra::nlyr(tg_mean)
        )
      ),

      precipitation = list(
        crs = terra::crs(
          rr_mean,
          describe = TRUE
        ),
        extent = as.character(
          terra::ext(rr_mean)
        ),
        resolution = as.numeric(
          terra::res(rr_mean)
        ),
        dimensions = c(
          nrow = terra::nrow(rr_mean),
          ncol = terra::ncol(rr_mean),
          nlyr = terra::nlyr(rr_mean)
        )
      )
    ),

    FAIR = list(

      Findable = paste(
        "Descriptive output and metadata filenames are used."
      ),

      Accessible = paste(
        "Temperature and precipitation outputs are stored as",
        "GeoTIFF files."
      ),

      Interoperable = paste(
        "GeoTIFF is compatible with R terra, QGIS, GDAL,",
        "and Python rasterio."
      ),

      Reusable = paste(
        "Metadata records the input file, year, selected",
        "breeding-season DOYs, CRS, processing method,",
        "raster properties, and creation date."
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
        "Temperature and precipitation rasters were saved, ",
        "but metadata could not be written: ",
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
  # Return
  # --------------------------------------------------------

  list(
    temperature = tg_mean,
    precipitation = rr_mean
  )
}
