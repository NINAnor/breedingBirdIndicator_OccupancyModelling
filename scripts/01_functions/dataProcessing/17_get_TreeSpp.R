#' Read and Project SR16 Dominant Tree Species Raster
#'
#' Reads the SR16 dominant tree species raster prepared for Norway and
#' projects it to the specified target CRS. Because dominant tree species
#' is a categorical variable, nearest-neighbour resampling is used.
#'
#' The raster values 1, 2 and 3 are spruce, pine and deciduous trees. The
#' value 0 does not indicate a tree type and is set to `NA`.
#'
#' If the processed raster already exists and `overwrite = FALSE`, the
#' cached raster is read from disk instead of being recomputed (use
#' `overwrite = TRUE` to apply the recoding of 0 to a raster created earlier).
#'
#' This function supports FAIR principles:
#' \itemize{
#'   \item Findable: output and metadata use descriptive filenames.
#'   \item Accessible: output is saved as GeoTIFF.
#'   \item Interoperable: GeoTIFF is readable by R, QGIS, GDAL, and Python.
#'   \item Reusable: metadata records the input file, CRS, resolution,
#'         extent, processing date, resampling method, and source.
#' }
#'
#' @param input_file Character. Path to the SR16 dominant tree species
#'   raster, e.g. `"dominantTreeSpecies.tif"`.
#' @param output_file Character. Path to the processed output raster,
#'   e.g. `"TreeSpp.tif"`.
#' @param target_crs Character or numeric. Target CRS for the output raster.
#'   For example, `"EPSG:25833"` or `25833`.
#' @param metadata_file Character or `NULL`. Path to YAML metadata file.
#'   If `NULL`, metadata is written beside `output_file`.
#' @param overwrite Logical. If `TRUE`, recompute the raster even if
#'   `output_file` already exists.
#'
#' @return A `terra::SpatRaster`.
#'
#' @references
#' NIBIO (2026). The Norwegian Forest Resource Map (Skogressurskart)
#' SR16 and SR16-beta Raster Dataset.
#'
#' Hauglin, M., Breidenbach, J., McLean, J. P., Astrup, R., & Gobakken, T.
#' (2021). Large scale mapping of forest attributes using heterogeneous
#' sets of airborne laser scanning and National Forest Inventory data.
#' Forest Ecosystems, 8(1), 65.
#'
#' @export
get_TreeSpp <- function(
  input_file,
  output_file,
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

  # --------------------------------------------------------
  # Read cached output
  # --------------------------------------------------------

  if (file.exists(output_file) && !overwrite) {
    message(
      "Reading cached SR16 dominant tree species raster: ",
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

  target_crs <- as.character(target_crs)

  # --------------------------------------------------------
  # Read source raster
  # --------------------------------------------------------

  message(
    "Reading SR16 dominant tree species raster..."
  )

  tree_spp <- terra::rast(
    input_file
  )

  input_crs <- terra::crs(
    tree_spp,
    describe = TRUE
  )

  input_resolution <- terra::res(
    tree_spp
  )

  input_extent <- terra::ext(
    tree_spp
  )

  # --------------------------------------------------------
  # Project raster
  # --------------------------------------------------------

  message(
    "Projecting dominant tree species raster to ",
    target_crs,
    "..."
  )

  TreeSppRaster <- terra::project(
    tree_spp,
    target_crs,
    method = "near"
  )

  # --------------------------------------------------------
  # Class 0 does not indicate a tree type: set it to NA
  #
  # Values 1, 2 and 3 are spruce, pine and deciduous trees.
  # --------------------------------------------------------

  TreeSppRaster <- terra::classify(
    TreeSppRaster,
    cbind(0, NA)
  )

  # --------------------------------------------------------
  # Name output layer
  # --------------------------------------------------------

  names(TreeSppRaster) <- "dominantTreeSpecies"

  # --------------------------------------------------------
  # Write output
  # --------------------------------------------------------

  message(
    "Writing processed raster: ",
    output_file
  )

  terra::writeRaster(
    TreeSppRaster,
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
    title = "SR16 dominant tree species",

    description = paste(
      "Dominant tree species raster derived from the Norwegian",
      "Forest Resource Map (SR16). The raster is projected to",
      "the specified model CRS using nearest-neighbour resampling."
    ),

    created = as.character(
      Sys.time()
    ),

    source = list(
      dataset = paste(
        "Norwegian Forest Resource Map",
        "(Skogressurskart) SR16"
      ),

      provider = "Norwegian Institute of Bioeconomy Research (NIBIO)",

      variable = "Dominant tree species",

      source_file = input_file,

      baseline_year = 2021
    ),

    output = list(
      file = output_file,

      file_format = "GeoTIFF",

      variable = "dominant tree species",

      categorical = TRUE
    ),

    processing_method = list(
      "function" = "terra::project()",

      resampling = "nearest neighbour",

      rationale = paste(
        "Nearest-neighbour resampling preserves categorical",
        "tree-species class values."
      ),

      target_crs = target_crs
    ),

    raster_properties = list(
      input_crs = input_crs,

      output_crs = terra::crs(
        TreeSppRaster,
        describe = TRUE
      ),

      input_extent = as.character(
        input_extent
      ),

      output_extent = as.character(
        terra::ext(TreeSppRaster)
      ),

      input_resolution = as.numeric(
        input_resolution
      ),

      output_resolution = as.numeric(
        terra::res(TreeSppRaster)
      ),

      dimensions = c(
        nrow = terra::nrow(TreeSppRaster),
        ncol = terra::ncol(TreeSppRaster),
        nlyr = terra::nlyr(TreeSppRaster)
      ),

      names = names(TreeSppRaster)
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
        "Metadata records source, baseline year, input and",
        "output CRS, resolution, extent, dimensions,",
        "resampling method, and processing date."
      )
    ),

    citation = list(
      NIBIO = paste(
        "NIBIO (2026). The Norwegian Forest Resource Map",
        "(Skogressurskart) SR16 and SR16-beta Raster Dataset",
        "(16x16 meter). Norwegian Institute of Bioeconomy."
      ),

      Hauglin_et_al = paste(
        "Hauglin, M., Breidenbach, J., McLean, J. P., Astrup, R.,",
        "& Gobakken, T. (2021). Large scale mapping of forest",
        "attributes using heterogeneous sets of airborne laser",
        "scanning and National Forest Inventory data.",
        "Forest Ecosystems, 8(1), 65."
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
        "Tree species raster was saved, but metadata could not",
        " be written: ",
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
  # Return raster
  # --------------------------------------------------------

  TreeSppRaster
}
