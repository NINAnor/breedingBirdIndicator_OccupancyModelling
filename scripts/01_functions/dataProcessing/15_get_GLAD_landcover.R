#' Get GLAD Land-Cover Raster
#'
#' Reads the GLAD 2020 land-cover raster and ensures that it is in the
#' specified model CRS. If a processed raster already exists and
#' `overwrite = FALSE`, the existing raster is loaded first and its CRS
#' is checked against `model_crs`. Reprojection is performed only when
#' the CRS differs.
#'
#' Nearest-neighbour interpolation is used because the GLAD raster
#' contains categorical land-cover classes.
#'
#' @param input_file Character. Path to the original GLAD raster.
#' @param model_crs Character or numeric. Target model CRS,
#'   e.g. `"EPSG:25833"` or `25833`.
#' @param output_file Character. Path to processed output raster.
#' @param metadata_file Character or `NULL`. Path to YAML metadata file.
#' @param overwrite Logical. If `TRUE`, recompute the output raster.
#'
#' @return A `terra::SpatRaster`.
#'
#' @export
get_GLAD_landcover <- function(
    input_file = file.path(
      dataFolder,
      "Covariates",
      "Landcover_F_M",
      "lc2020nor_mergedClasses.tif"
    ),
    model_crs,
    output_file,
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
  # Validate CRS
  # --------------------------------------------------------

  if (is.numeric(model_crs)) {
    model_crs <- paste0(
      "EPSG:",
      model_crs
    )
  }

  model_crs <- as.character(model_crs)

  # Validate target CRS without projecting a raster
  tryCatch(
    terra::crs(model_crs),
    error = function(e) {
      stop(
        "Invalid target CRS: ",
        model_crs,
        "\n",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  # --------------------------------------------------------
  # Normalise paths
  # --------------------------------------------------------

  input_file <- normalizePath(
    input_file,
    winslash = "/",
    mustWork = TRUE
  )

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
  # Check existing processed raster FIRST
  # --------------------------------------------------------

  if (file.exists(output_file) && !overwrite) {

    message(
      "Reading existing GLAD raster: ",
      output_file
    )

    GLADRaster <- terra::rast(
      output_file
    )

    existing_crs <- terra::crs(
      GLADRaster
    )

    # ------------------------------------------------------
    # Existing raster already has model CRS
    # ------------------------------------------------------

    if (
      terra::same.crs(
        GLADRaster,
        terra::crs(model_crs)
      )
    ) {

      message(
        "Existing GLAD raster already uses model CRS: ",
        model_crs
      )

      return(GLADRaster)
    }

    # ------------------------------------------------------
    # Existing raster has different CRS
    # ------------------------------------------------------

    message(
      "Existing GLAD raster has a different CRS."
    )

    message(
      "Current CRS: ",
      existing_crs
    )

    message(
      "Projecting existing GLAD raster to: ",
      model_crs
    )

    GLADRaster_projected <- terra::project(
      GLADRaster,
      model_crs,
      method = "near"
    )

    names(GLADRaster_projected) <- names(
      GLADRaster
    )

    terra::writeRaster(
      GLADRaster_projected,
      output_file,
      overwrite = TRUE
    )

    return(
      GLADRaster_projected
    )
  }

  # --------------------------------------------------------
  # No cached output: read original raster
  # --------------------------------------------------------

  message(
    "Processed GLAD raster not found."
  )

  message(
    "Reading original GLAD raster: ",
    input_file
  )

  GLADRaster <- terra::rast(
    input_file
  )

  # --------------------------------------------------------
  # Check original CRS
  # --------------------------------------------------------

  input_crs <- terra::crs(
    GLADRaster,
    describe = TRUE
  )

  message(
    "Input CRS: ",
    terra::crs(GLADRaster)
  )

  # --------------------------------------------------------
  # Project only if necessary
  # --------------------------------------------------------

  if (
    terra::same.crs(
      GLADRaster,
      terra::crs(model_crs)
    )
  ) {

    message(
      "Input GLAD raster already uses model CRS."
    )

    GLADRaster_projected <- GLADRaster

  } else {

    message(
      "Projecting GLAD raster to: ",
      model_crs
    )

    GLADRaster_projected <- terra::project(
      GLADRaster,
      model_crs,
      method = "near"
    )
  }

  # --------------------------------------------------------
  # Set descriptive layer names
  # --------------------------------------------------------

  if (
    terra::nlyr(GLADRaster_projected) == 1L
  ) {

    names(GLADRaster_projected) <-
      "GLAD_landcover"
  }

  # --------------------------------------------------------
  # Save processed raster
  # --------------------------------------------------------

  message(
    "Writing processed GLAD raster: ",
    output_file
  )

  terra::writeRaster(
    GLADRaster_projected,
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
  # Class summary
  # --------------------------------------------------------

  class_summary <- tryCatch(
    terra::freq(
      GLADRaster_projected
    ),
    error = function(e) {
      NULL
    }
  )

  # --------------------------------------------------------
  # Metadata
  # --------------------------------------------------------

  metadata <- list(

    title = "GLAD 2020 land-cover raster for Norway",

    description = paste(
      "GLAD 2020 land-cover raster harmonised to the model CRS.",
      "Reprojection is performed only when the input or cached",
      "raster CRS differs from the requested model CRS."
    ),

    created = as.character(
      Sys.time()
    ),

    input = list(

      file = input_file,

      dataset = paste(
        "Global 2000-2020 Land Cover and Land Use Change",
        "dataset (GLAD)"
      ),

      year = 2020,

      source = paste(
        "University of Maryland / Global Land Analysis",
        "and Discovery Laboratory"
      )
    ),

    output = list(

      file = output_file,

      file_format = "GeoTIFF",

      variable = "land_cover",

      target_crs = model_crs
    ),

    processing_method = list(

      "function" = "terra::project()",

      reprojection = paste(
        "Raster was projected to the model CRS only when",
        "the source CRS differed from the model CRS."
      ),

      resampling = "nearest neighbour",

      resampling_method = "near",

      rationale = paste(
        "Nearest-neighbour interpolation preserves categorical",
        "land-cover class values and prevents creation of",
        "non-existent intermediate classes."
      )
    ),

    raster_properties = list(

      input_crs = input_crs,

      output_crs = terra::crs(
        GLADRaster_projected
      ),

      extent = as.character(
        terra::ext(
          GLADRaster_projected
        )
      ),

      resolution = as.numeric(
        terra::res(
          GLADRaster_projected
        )
      ),

      dimensions = c(
        nrow = terra::nrow(
          GLADRaster_projected
        ),
        ncol = terra::ncol(
          GLADRaster_projected
        ),
        nlyr = terra::nlyr(
          GLADRaster_projected
        )
      ),

      names = names(
        GLADRaster_projected
      )
    ),

    class_summary = class_summary,

    FAIR = list(

      Findable =
        "Descriptive output and metadata filenames are used.",

      Accessible =
        "Processed land-cover raster is stored as GeoTIFF.",

      Interoperable =
        "GeoTIFF is compatible with R terra, QGIS, GDAL, and Python rasterio.",

      Reusable =
        "Metadata records the source raster, dataset, year, CRS,",
      "resolution, extent, dimensions, resampling method, and",
      "processing date."
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
        "GLAD raster was saved, but metadata could not be written: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  message(
    "GLAD metadata saved to: ",
    metadata_file
  )

  # --------------------------------------------------------
  # Return
  # --------------------------------------------------------

  GLADRaster_projected
}
