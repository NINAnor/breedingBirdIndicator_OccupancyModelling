#' Calculate Annual Snowmelt Date from seNorge Snow Depth
#'
#' Calculates the last day with snow before the first sustained snow-free
#' period from daily seNorge snow-depth data. A sustained snow-free period is
#' defined as `run_threshold` consecutive days with snow depth less than or
#' equal to `snow_threshold`.
#'
#' The function reads daily seNorge snow-depth NetCDF files for a specified
#' year, calculates the snowmelt day-of-year (DOY) for each raster cell, and
#' saves the resulting raster as a GeoTIFF. If the output already exists and
#' `overwrite = FALSE`, the cached raster is returned instead of being
#' recomputed.
#'
#' The seNorge snow model provides spatially distributed estimates of snow
#' conditions across Norway.
#'
#' This function supports FAIR principles:
#' \itemize{
#'   \item Findable: output and metadata use descriptive filenames.
#'   \item Accessible: output is saved as GeoTIFF.
#'   \item Interoperable: GeoTIFF is readable by R, QGIS, GDAL, and Python.
#'   \item Reusable: metadata records the input directory, source files,
#'         year, CRS, resolution, extent, processing date, threshold,
#'         snow-free run length, and processing method.
#' }
#'
#' @param snow_dir Character. Directory containing the daily seNorge
#'   snow-depth files.
#' @param output_file Character. Path to the output snowmelt raster,
#'   e.g. `"snowmelt_2021.tif"`.
#' @param year Integer. Year for which snowmelt should be calculated.
#' @param pattern Character or `NULL`. Pattern used to identify files for
#'   the specified year. If `NULL`, the year is used as the search pattern.
#' @param snow_threshold Numeric. Snow-depth threshold defining a snow-free
#'   day. Default is `0`.
#' @param run_threshold Integer. Number of consecutive snow-free days required
#'   to define sustained snowmelt. Default is `14`.
#' @param target_crs Character or numeric. Target CRS. If `NULL`, the CRS of
#'   the input raster is retained. Default is `"EPSG:25833"`.
#' @param metadata_file Character or `NULL`. Path to YAML metadata file.
#'   If `NULL`, metadata is written beside `output_file`.
#' @param overwrite Logical. If `TRUE`, recompute the raster even if
#'   `output_file` already exists.
#'
#' @return A `terra::SpatRaster` containing the last snow day as
#'   day-of-year.
#'
#' @references
#' Saloranta, T. (2016). Operational snow mapping with simplified data
#' assimilation using the seNorge snow model. Journal of Hydrology, 538,
#' 314--325. https://doi.org/10.1016/j.jhydrol.2016.03.06
#'
#' @export
get_snowmelt <- function(
    snow_dir,
    output_file,
    year,
    pattern = NULL,
    snow_threshold = 0,
    run_threshold = 14,
    target_crs = "EPSG:25833",
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
    is.na(year) ||
    !is.numeric(year)
  ) {
    stop(
      "`year` must be a single numeric value.",
      call. = FALSE
    )
  }

  year <- as.integer(year)

  if (
    length(snow_threshold) != 1L ||
    is.na(snow_threshold)
  ) {
    stop(
      "`snow_threshold` must be a single numeric value.",
      call. = FALSE
    )
  }

  if (
    length(run_threshold) != 1L ||
    is.na(run_threshold) ||
    run_threshold < 1
  ) {
    stop(
      "`run_threshold` must be a positive integer.",
      call. = FALSE
    )
  }

  run_threshold <- as.integer(run_threshold)

  # --------------------------------------------------------
  # Validate directories and paths
  # --------------------------------------------------------

  snow_dir <- normalizePath(
    snow_dir,
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

  if (file.exists(output_file) && !overwrite) {

    message(
      "Reading cached snowmelt raster: ",
      output_file
    )

    return(
      terra::rast(output_file)
    )
  }

  dir.create(
    dirname(output_file),
    recursive = TRUE,
    showWarnings = FALSE
  )

  # --------------------------------------------------------
  # Find daily snow-depth files
  # --------------------------------------------------------

  if (is.null(pattern)) {
    pattern <- as.character(year)
  }

  files <- list.files(
    snow_dir,
    pattern = pattern,
    full.names = TRUE,
    recursive = TRUE
  )

  if (length(files) == 0L) {

    stop(
      "No snow-depth files found for year ",
      year,
      " in: ",
      snow_dir,
      call. = FALSE
    )
  }

  files <- normalizePath(
    files,
    winslash = "/",
    mustWork = TRUE
  )

  files <- sort(files)

  message(
    "Found ",
    length(files),
    " snow-depth files for ",
    year,
    "."
  )

  # --------------------------------------------------------
  # Read daily snow-depth data
  # --------------------------------------------------------

  message("Reading daily snow-depth data...")

  r <- terra::rast(files)

  message(
    terra::nlyr(r),
    " daily raster layers loaded."
  )

  # --------------------------------------------------------
  # CRS
  # --------------------------------------------------------

  input_crs <- terra::crs(
    r,
    describe = TRUE
  )

  # Set target CRS if supplied
  if (!is.null(target_crs)) {

    if (is.numeric(target_crs)) {
      target_crs <- paste0(
        "EPSG:",
        target_crs
      )
    }

    target_crs <- as.character(target_crs)

    if (
      !identical(
        terra::crs(r),
        target_crs
      )
    ) {

      message(
        "Reprojecting snow-depth data to ",
        target_crs,
        "..."
      )

      r <- terra::project(
        r,
        target_crs,
        method = "bilinear"
      )
    }
  }

  # --------------------------------------------------------
  # Day of year
  # --------------------------------------------------------

  start_date <- as.Date(
    paste0(year, "-01-01")
  )

  doy <- as.integer(
    format(
      seq.Date(
        start_date,
        by = "day",
        length.out = terra::nlyr(r)
      ),
      "%j"
    )
  )

  # --------------------------------------------------------
  # Calculate last snow day
  # --------------------------------------------------------

  message(
    "Calculating last snow day before a ",
    run_threshold,
    "-day snow-free period..."
  )

  last_snow <- terra::app(

    r,

    fun = function(x) {

      # Identify snow-free days
      snow_free <- x <= snow_threshold

      # Remove missing values only if necessary
      if (all(is.na(snow_free))) {
        return(NA_real_)
      }

      # Run-length encoding
      rr <- rle(snow_free)

      # Start and end positions of each run
      ends <- cumsum(rr$lengths)

      starts <- ends - rr$lengths + 1

      # Find snow-free runs of sufficient duration
      keep <- which(
        rr$values &
          rr$lengths >= run_threshold
      )

      # No sustained snow-free period
      if (length(keep) == 0L) {
        return(NA_real_)
      }

      # First sustained snow-free period
      first_run <- keep[1]

      # First day of sustained snow-free period
      melt_day <- starts[first_run]

      # If already snow-free on the first day,
      # there is no preceding snow-covered period
      if (melt_day == 1L) {
        return(NA_real_)
      }

      # Last day with snow
      doy[melt_day - 1L]
    }
  )

  names(last_snow) <- paste0(
    "LastSnowDOY_",
    year
  )

  # --------------------------------------------------------
  # Write output raster
  # --------------------------------------------------------

  message(
    "Writing snowmelt raster: ",
    output_file
  )

  terra::writeRaster(
    last_snow,
    output_file,
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

    title = paste(
      "Annual snowmelt date derived from seNorge snow depth",
      year
    ),

    description = paste(
      "Last day with snow before the first sustained",
      run_threshold,
      "-day snow-free period."
    ),

    created = as.character(
      Sys.time()
    ),

    input = list(

      directory = snow_dir,

      year = year,

      n_input_files = length(files),

      input_files = basename(files),

      dataset = "seNorge snow model",

      variable = "snow depth",

      source = "Norwegian Meteorological Institute (MET Norway)",

      source_url = paste0(
        "https://thredds.met.no/thredds/",
        "catalog/senorge/seNorge_snow/catalog.html"
      ),

      documentation_url = paste0(
        "https://github.com/metno/seNorge_docs/",
        "wiki/senorgesnow_variables"
      )
    ),

    output = list(

      file = output_file,

      file_format = "GeoTIFF",

      variable = "last snow day",

      unit = "day of year",

      year = year
    ),

    processing_method = list(

      "function" = "terra::app()",

      snow_threshold = snow_threshold,

      snow_free_definition = paste(
        "Snow depth <=",
        snow_threshold
      ),

      run_threshold_days = run_threshold,

      method = paste(
        "First run of",
        run_threshold,
        "or more consecutive snow-free days."
      ),

      output_definition = paste(
        "Day of year corresponding to the final day with snow",
        "immediately before the first sustained snow-free period."
      ),

      missing_value_definition = paste(
        "NA when no sustained snow-free period is detected",
        "or when the year begins with a sustained snow-free period."
      )
    ),

    raster_properties = list(

      input_crs = input_crs,

      output_crs = terra::crs(
        last_snow,
        describe = TRUE
      ),

      extent = as.character(
        terra::ext(last_snow)
      ),

      resolution = as.numeric(
        terra::res(last_snow)
      ),

      dimensions = c(
        nrow = terra::nrow(last_snow),
        ncol = terra::ncol(last_snow),
        nlyr = terra::nlyr(last_snow)
      ),

      names = names(last_snow)
    ),

    FAIR = list(

      Findable = paste(
        "Descriptive output and metadata filenames are used."
      ),

      Accessible = paste(
        "Snowmelt raster is stored as GeoTIFF."
      ),

      Interoperable = paste(
        "GeoTIFF is compatible with R terra, QGIS, GDAL,",
        "and Python rasterio."
      ),

      Reusable = paste(
        "Metadata records the source files, year, snow-depth",
        "threshold, snow-free run length, processing method,",
        "CRS, extent, resolution, dimensions, and creation date."
      )
    ),

    citation = paste(
      "Saloranta, T. (2016). Operational snow mapping with",
      "simplified data assimilation using the seNorge snow model.",
      "Journal of Hydrology, 538, 314-325.",
      "https://doi.org/10.1016/j.jhydrol.2016.03.06"
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
        "Snowmelt raster was saved, but metadata could not be written: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  message(
    "Snowmelt metadata saved to: ",
    metadata_file
  )

  # --------------------------------------------------------
  # Return raster
  # --------------------------------------------------------

  last_snow
}
