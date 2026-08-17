#' Create Latitude Raster
#'
#' Creates a raster containing geographic latitude in decimal degrees for
#' each raster cell centre. The raster is created using a supplied raster
#' template so that its extent, resolution, dimensions, and CRS match the
#' model raster.
#'
#' Latitude is calculated by transforming raster cell-centre coordinates
#' from the model CRS to geographic coordinates (EPSG:4326). Values are
#' rounded to four decimal places.
#'
#' If the output raster already exists and `overwrite = FALSE`, the cached
#' raster is loaded and returned without recomputation.
#'
#' This function supports FAIR principles:
#' \itemize{
#'   \item Findable: output and metadata use descriptive filenames.
#'   \item Accessible: latitude is saved as GeoTIFF.
#'   \item Interoperable: GeoTIFF is readable by R, QGIS, GDAL, and Python.
#'   \item Reusable: metadata records CRS, resolution, extent, processing
#'         method, precision, and creation date.
#' }
#'
#' @param boundary An `sf` polygon containing the spatial extent for the
#'   latitude raster.
#' @param model_crs Character or numeric. CRS used by the model,
#'   e.g. `"EPSG:25833"` or `25833`.
#' @param template Character or `terra::SpatRaster`. Raster template defining
#'   the extent, resolution, dimensions, and CRS of the output.
#' @param output_file Character. Path to the output GeoTIFF.
#' @param metadata_file Character or `NULL`. Path to YAML metadata file.
#'   If `NULL`, metadata is written beside `output_file`.
#' @param overwrite Logical. If `TRUE`, recreate the raster even if it exists.
#'
#' @return A `terra::SpatRaster` containing latitude in decimal degrees.
#'
#' @export
get_latitude <- function(
    boundary,
    model_crs,
    template,
    output_file,
    metadata_file = NULL,
    overwrite = FALSE
) {

  # --------------------------------------------------------
  # Required packages
  # --------------------------------------------------------

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
  # Validate model CRS
  # --------------------------------------------------------

  if (is.numeric(model_crs)) {

    model_crs <- paste0(
      "EPSG:",
      model_crs
    )
  }

  model_crs <- as.character(
    model_crs
  )

  tryCatch(

    terra::crs(model_crs),

    error = function(e) {

      stop(
        "Invalid model CRS: ",
        model_crs,
        "\n",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  # --------------------------------------------------------
  # Validate boundary
  # --------------------------------------------------------

  if (!inherits(boundary, "sf")) {

    stop(
      "`boundary` must be an sf object.",
      call. = FALSE
    )
  }

  boundary_crs <- sf::st_crs(
    boundary
  )

  if (is.na(boundary_crs)) {

    stop(
      "`boundary` does not have a valid CRS.",
      call. = FALSE
    )
  }

  # --------------------------------------------------------
  # Transform boundary to model CRS if necessary
  # --------------------------------------------------------

  if (
    boundary_crs$wkt !=
    sf::st_crs(model_crs)$wkt
  ) {

    message(
      "Transforming boundary to model CRS: ",
      model_crs
    )

    boundary <- sf::st_transform(
      boundary,
      model_crs
    )
  }

  # --------------------------------------------------------
  # Normalise output paths
  # --------------------------------------------------------

  output_file <- normalizePath(
    output_file,
    winslash = "/",
    mustWork = FALSE
  )

  dir.create(
    dirname(output_file),
    recursive = TRUE,
    showWarnings = FALSE
  )

  # --------------------------------------------------------
  # Read cached raster first
  # --------------------------------------------------------

  if (
    file.exists(output_file) &&
    !overwrite
  ) {

    message(
      "Reading cached latitude raster: ",
      output_file
    )

    cached_latitude <- terra::rast(
      output_file
    )

    # ------------------------------------------------------
    # Check cached raster CRS
    # ------------------------------------------------------

    if (
      terra::same.crs(
        cached_latitude,
        terra::crs(model_crs)
      )
    ) {

      message(
        "Cached latitude raster already uses model CRS."
      )

      return(
        cached_latitude
      )
    }

    # ------------------------------------------------------
    # Cached raster has unexpected CRS
    # ------------------------------------------------------

    message(
      "Cached latitude raster has a different CRS."
    )

    message(
      "Cached CRS: ",
      terra::crs(cached_latitude)
    )

    message(
      "Expected CRS: ",
      model_crs
    )

    message(
      "Recreating latitude raster from the template."
    )
  }

  # --------------------------------------------------------
  # Read raster template
  # --------------------------------------------------------

  if (inherits(template, "SpatRaster")) {

    template_raster <- template

  } else if (is.character(template)) {

    template <- normalizePath(
      template,
      winslash = "/",
      mustWork = TRUE
    )

    template_raster <- terra::rast(
      template
    )

  } else {

    stop(
      "`template` must be either a terra::SpatRaster or ",
      "a path to a raster file.",
      call. = FALSE
    )
  }

  # --------------------------------------------------------
  # Check template CRS
  # --------------------------------------------------------

  if (
    !terra::same.crs(
      template_raster,
      terra::crs(model_crs)
    )
  ) {

    stop(
      "The supplied raster template does not use model_crs.\n",
      "Template CRS: ",
      terra::crs(template_raster),
      "\n",
      "Model CRS: ",
      model_crs,
      call. = FALSE
    )
  }

  # --------------------------------------------------------
  # Create output raster from template
  # --------------------------------------------------------

  latitude_raster <- template_raster[[1]]

  names(latitude_raster) <- "latitude"

  # --------------------------------------------------------
  # Get raster cell-centre coordinates
  # --------------------------------------------------------

  message(
    "Calculating latitude for ",
    terra::ncell(latitude_raster),
    " raster cells..."
  )

  cell_ids <- seq_len(
    terra::ncell(latitude_raster)
  )

  xy <- terra::xyFromCell(
    latitude_raster,
    cell_ids
  )

  # --------------------------------------------------------
  # Transform cell centres to geographic coordinates
  # --------------------------------------------------------

  xy_sf <- sf::st_as_sf(

    data.frame(
      x = xy[, 1],
      y = xy[, 2]
    ),

    coords = c(
      "x",
      "y"
    ),

    crs = model_crs
  )

  xy_lonlat <- sf::st_transform(
    xy_sf,
    "EPSG:4326"
  )

  # --------------------------------------------------------
  # Extract latitude and round to four decimal places
  # --------------------------------------------------------

  latitude_values <- round(
    sf::st_coordinates(
      xy_lonlat
    )[, 2],
    4
  )

  # --------------------------------------------------------
  # Populate raster
  # --------------------------------------------------------

  latitude_raster <- terra::setValues(
    latitude_raster,
    latitude_values
  )

  names(latitude_raster) <- "latitude"

  # --------------------------------------------------------
  # Mask to Norway boundary
  # --------------------------------------------------------

  message(
    "Masking latitude raster to boundary."
  )

  latitude_raster <- terra::mask(
    latitude_raster,
    terra::vect(boundary)
  )

  # --------------------------------------------------------
  # Write output raster
  # --------------------------------------------------------

  message(
    "Writing latitude raster: ",
    output_file
  )

  terra::writeRaster(
    latitude_raster,
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

    title = "Latitude raster",

    description = paste(
      "Raster containing geographic latitude at each raster",
      "cell centre, calculated from the model CRS."
    ),

    created = as.character(
      Sys.time()
    ),

    input = list(

      boundary_crs =
        sf::st_crs(boundary)$input,

      model_crs =
        model_crs,

      geographic_crs =
        "EPSG:4326",

      template =
        if (inherits(template, "SpatRaster")) {
          "terra::SpatRaster"
        } else {
          template
        }
    ),

    output = list(

      file =
        output_file,

      file_format =
        "GeoTIFF",

      variable =
        "latitude",

      units =
        "decimal degrees",

      decimal_places =
        4
    ),

    processing_method = list(

      "coordinate_transformation" =
        "sf::st_transform()",

      "cell_coordinates" =
        "terra::xyFromCell()",

      "geographic_crs" =
        "EPSG:4326",

      "rounding" =
        "Latitude values rounded to 4 decimal places.",

      "masking" =
        "Raster masked to the supplied boundary."
    ),

    raster_properties = list(

      crs =
        terra::crs(
          latitude_raster,
          describe = TRUE
        ),

      extent =
        as.character(
          terra::ext(
            latitude_raster
          )
        ),

      resolution =
        as.numeric(
          terra::res(
            latitude_raster
          )
        ),

      dimensions = c(

        nrow =
          terra::nrow(
            latitude_raster
          ),

        ncol =
          terra::ncol(
            latitude_raster
          ),

        nlyr =
          terra::nlyr(
            latitude_raster
          )
      ),

      names =
        names(
          latitude_raster
        )
    ),

    FAIR = list(

      Findable =
        "Descriptive output and metadata filenames are used.",

      Accessible =
        "Latitude raster is stored as GeoTIFF.",

      Interoperable =
        "GeoTIFF is compatible with R terra, QGIS, GDAL, and Python.",

      Reusable =
        paste(
          "Metadata records the model CRS, geographic CRS,",
          "template, resolution, extent, processing method,",
          "precision, and creation date."
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
        "Latitude raster was saved, but metadata could not ",
        "be written: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  message(
    "Latitude metadata saved to: ",
    metadata_file
  )

  # --------------------------------------------------------
  # Return raster
  # --------------------------------------------------------

  latitude_raster
}
