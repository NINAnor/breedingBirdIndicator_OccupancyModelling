#' Reproject GLAD Land-Cover Raster and Cache Output
#'
#' Reads the GLAD 2020 land-cover raster for Norway, reprojects it to the
#' specified model CRS using nearest-neighbour interpolation, and saves the
#' result as a GeoTIFF. If the processed raster already exists and
#' `overwrite = FALSE`, the cached raster is read from disk instead of being
#' recomputed.
#'
#' The GLAD raster contains land-cover classes and associated vegetation
#' information derived from the Global 2000--2020 Land Cover and Land Use
#' Change dataset.
#'
#' Nearest-neighbour interpolation is used because the input raster contains
#' categorical land-cover classes and therefore should not be interpolated
#' using bilinear or cubic methods.
#'
#' This function supports FAIR principles:
#' \itemize{
#'   \item Findable: output and metadata use descriptive filenames.
#'   \item Accessible: output is saved as GeoTIFF.
#'   \item Interoperable: GeoTIFF is readable by R, QGIS, GDAL, and Python.
#'   \item Reusable: metadata records the input raster, source, CRS,
#'         resolution, extent, processing date, method, and raster properties.
#' }
#'
#' @param input_file Character. Path to the original GLAD raster.
#'   If `NULL`, the default GLAD raster within `dataFolder` is used.
#' @param model_crs Character or numeric. Target model CRS, e.g.
#'   `"EPSG:25833"` or `25833`.
#' @param output_file Character. Path to the processed output raster,
#'   e.g. `"LCglad.tif"`.
#' @param metadata_file Character or `NULL`. Path to YAML metadata file.
#'   If `NULL`, metadata is written beside `output_file`.
#' @param overwrite Logical. If `TRUE`, recompute the raster even if
#'   `output_file` already exists.
#'
#' @return A `terra::SpatRaster`.
#'
#' @references
#' Potapov, P., Hansen, M.C., Pickens, A., Hernandez-Serna, A.,
#' Tyukavina, A., Turubanova, S., Zalles, V., Li, X., Khan, A.,
#' Stolle, F., Harris, N., Song, X.-P., Baggett, A., Kommareddy, I.,
#' & Kommareddy, A. (2022). The global 2000--2020 land cover and land
#' use change dataset derived from the Landsat archive: first results.
#' Frontiers in Remote Sensing, 3, 856903.
#' https://doi.org/10.3389/frsen.2022.856903
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
  # Validate input
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

  # --------------------------------------------------------
  # Validate target CRS
  # --------------------------------------------------------

  if (is.numeric(model_crs)) {
    model_crs <- paste0("EPSG:", model_crs)
  }

  model_crs <- as.character(model_crs)

  tryCatch(
    {
      terra::crs(
        terra::rast(input_file)
      )

      terra::crs(
        terra::rast(input_file),
        describe = TRUE
      )

      # Validate target CRS independently
      terra::crs(
        terra::project(
          terra::rast(input_file)[[1]],
          model_crs,
          method = "near"
        ),
        describe = TRUE
      )
    },
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
  # Read cached output
  # --------------------------------------------------------

  if (file.exists(output_file) && !overwrite) {

    message(
      "Reading cached GLAD land-cover raster: ",
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
  # Read GLAD raster
  # --------------------------------------------------------

  message(
    "Reading GLAD land-cover raster: ",
    input_file
  )

  GLADRaster <- terra::rast(input_file)

  # --------------------------------------------------------
  # Check input CRS
  # --------------------------------------------------------

  if (terra::is.lonlat(GLADRaster)) {
    message(
      "Input GLAD raster is in geographic coordinates."
    )
  }

  input_crs <- terra::crs(
    GLADRaster,
    describe = TRUE
  )

  # --------------------------------------------------------
  # Reproject
  # --------------------------------------------------------

  message(
    "Reprojecting GLAD raster to ",
    model_crs,
    "..."
  )

  GLADRaster_projected <- terra::project(
    GLADRaster,
    model_crs,
    method = "near"
  )

  # --------------------------------------------------------
  # Set descriptive layer names
  # --------------------------------------------------------

  if (terra::nlyr(GLADRaster_projected) == 1L) {
    names(GLADRaster_projected) <- "GLAD_landcover"
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
    terra::freq(GLADRaster_projected),
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
      "GLAD 2020 land-cover raster reprojected to the model CRS",
      "using nearest-neighbour interpolation."
    ),

    created = as.character(Sys.time()),

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
      ),

      source_url = paste0(
        "https://storage.googleapis.com/",
        "earthenginepartners-hansen/",
        "GLCLU2000-2020/v2/download.html"
      ),

      STAC_catalogue = paste0(
        "https://stac.openlandmap.org/",
        "lc_glad.glcluc/collection.json"
      ),

      description = paste(
        "Norwegian GLAD land-cover raster prepared from",
        "the original GLAD dataset. Tree-height classes",
        "have been merged where appropriate."
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
        "Input GLAD raster reprojected to the model CRS."
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
        terra::ext(GLADRaster_projected)
      ),

      resolution = as.numeric(
        terra::res(GLADRaster_projected)
      ),

      dimensions = c(
        nrow = terra::nrow(GLADRaster_projected),
        ncol = terra::ncol(GLADRaster_projected),
        nlyr = terra::nlyr(GLADRaster_projected)
      ),

      names = names(
        GLADRaster_projected
      )
    ),

    class_summary = class_summary,

    FAIR = list(

      Findable = paste(
        "Descriptive output and metadata filenames are used."
      ),

      Accessible = paste(
        "Processed land-cover raster is stored as GeoTIFF."
      ),

      Interoperable = paste(
        "GeoTIFF is compatible with R terra, QGIS, GDAL,",
        "and Python rasterio."
      ),

      Reusable = paste(
        "Metadata records the source raster, dataset, year,",
        "CRS, resolution, extent, dimensions, resampling",
        "method, and processing date."
      )
    ),

    citation = paste(
      "Potapov et al. (2022). The global 2000-2020 land cover",
      "and land use change dataset derived from the Landsat",
      "archive: first results. Frontiers in Remote Sensing,",
      "3, 856903."
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
  # Return raster
  # --------------------------------------------------------

  GLADRaster_projected
}
