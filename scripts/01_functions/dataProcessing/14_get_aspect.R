#' Calculate Aspect from a DEM and Cache Output
#'
#' Calculates terrain aspect in degrees from a digital elevation model (DEM)
#' using `terra::terrain()`. If the aspect raster already exists and
#' `overwrite = FALSE`, the cached raster is read from disk instead of being
#' recomputed.
#'
#' This function supports FAIR principles:
#' \itemize{
#'   \item Findable: output and metadata use descriptive filenames.
#'   \item Accessible: output is saved as GeoTIFF.
#'   \item Interoperable: GeoTIFF is readable by R, QGIS, GDAL, and Python.
#'   \item Reusable: metadata records the source DEM, CRS, resolution,
#'         extent, processing date, method, and raster properties.
#' }
#'
#' @param dem A `terra::SpatRaster` containing the digital elevation model.
#' @param output_file Character. Path to the output aspect raster,
#'   e.g. `"aspect.tif"`.
#' @param metadata_file Character or `NULL`. Path to YAML metadata file.
#'   If `NULL`, metadata is written beside `output_file`.
#' @param overwrite Logical. If `TRUE`, recompute the aspect raster even if
#'   `output_file` already exists.
#'
#' @return A `terra::SpatRaster` containing aspect in degrees.
#'
#' @export
get_aspect <- function(
    dem,
    output_file,
    metadata_file = NULL,
    overwrite = FALSE
) {

  required_pkgs <- c("terra", "yaml")

  missing_pkgs <- required_pkgs[
    !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_pkgs) > 0) {
    stop(
      "Install missing packages first: ",
      paste(missing_pkgs, collapse = ", "),
      call. = FALSE
    )
  }

  # --------------------------------------------------------
  # Validate DEM
  # --------------------------------------------------------

  if (!inherits(dem, "SpatRaster")) {
    stop(
      "`dem` must be a terra::SpatRaster returned by get_DEM().",
      call. = FALSE
    )
  }

  if (terra::nlyr(dem) != 1L) {
    warning(
      "`dem` contains ", terra::nlyr(dem),
      " layers. Using the first layer only.",
      call. = FALSE
    )

    dem <- dem[[1]]
  }

  # --------------------------------------------------------
  # Validate output path
  # --------------------------------------------------------

  output_file <- normalizePath(
    output_file,
    winslash = "/",
    mustWork = FALSE
  )

  if (file.exists(output_file) && !overwrite) {
    message("Reading cached aspect raster: ", output_file)
    return(terra::rast(output_file))
  }

  dir.create(
    dirname(output_file),
    recursive = TRUE,
    showWarnings = FALSE
  )

  # --------------------------------------------------------
  # Calculate aspect
  # --------------------------------------------------------

  message("Calculating aspect from DEM...")

  aspect <- terra::terrain(
    dem,
    v = "aspect",
    unit = "degrees"
  )

  names(aspect) <- "aspect_degrees"

  # --------------------------------------------------------
  # Save raster
  # --------------------------------------------------------

  terra::writeRaster(
    aspect,
    output_file,
    overwrite = TRUE
  )

  message("Aspect raster saved to: ", output_file)

  # --------------------------------------------------------
  # Metadata
  # --------------------------------------------------------

  if (is.null(metadata_file)) {
    metadata_file <- file.path(
      dirname(output_file),
      paste0(
        tools::file_path_sans_ext(basename(output_file)),
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

  # Determine source DEM information where available
  dem_source <- tryCatch(
    terra::sources(dem),
    error = function(e) character(0)
  )

  metadata <- list(
    title = "Terrain aspect derived from DEM",

    description = paste(
      "Aspect calculated from a digital elevation model using",
      "terra::terrain() and expressed in degrees."
    ),

    created = as.character(Sys.time()),

    input = list(
      type = "Digital Elevation Model",
      source_file = if (length(dem_source) > 0) {
        dem_source[1]
      } else {
        NA_character_
      },
      raster_name = names(dem)
    ),

    output = list(
      file = output_file,
      file_format = "GeoTIFF",
      variable = "aspect",
      unit = "degrees"
    ),

    processing_method = list(
      "function" = "terra::terrain()",
      terrain_variable = "aspect",
      unit = "degrees"
    ),

    raster_properties = list(
      crs = terra::crs(aspect),
      extent = as.character(terra::ext(aspect)),
      resolution = as.numeric(terra::res(aspect)),
      dimensions = c(
        nrow = terra::nrow(aspect),
        ncol = terra::ncol(aspect),
        nlyr = terra::nlyr(aspect)
      ),
      names = names(aspect)
    ),

    FAIR = list(
      Findable = "Descriptive output and metadata filenames are used.",
      Accessible = "Aspect raster is stored as GeoTIFF.",
      Interoperable = paste(
        "GeoTIFF is compatible with R terra, QGIS, GDAL,",
        "and Python rasterio."
      ),
      Reusable = paste(
        "Metadata records the input DEM, processing method,",
        "CRS, extent, resolution, dimensions, output format,",
        "and creation date."
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
        "Aspect raster was saved, but metadata could not be written: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  message("Aspect metadata saved to: ", metadata_file)

  aspect
}
