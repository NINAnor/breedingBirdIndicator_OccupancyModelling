#' Merge TIF Raster Tiles and Cache Output
#'
#' Merges all `.tif` files in a directory into one raster using `terra::mosaic()`.
#' If the merged raster already exists and `overwrite = FALSE`, the saved raster
#' is read from disk instead of being recomputed.
#'
#' This function supports FAIR principles:
#' \itemize{
#'   \item Findable: output and metadata use descriptive filenames.
#'   \item Accessible: output is saved as GeoTIFF.
#'   \item Interoperable: GeoTIFF is readable by R, QGIS, GDAL, and Python.
#'   \item Reusable: metadata records input directory, source files, CRS,
#'         resolution, extent, processing date, and method.
#' }
#'
#' @param input_dir Character. Directory containing `.tif` files.
#' @param output_file Character. Path to merged output raster, e.g. `"merged_dem.tif"`.
#' @param metadata_file Character or `NULL`. Path to YAML metadata file.
#'   If `NULL`, metadata is written beside `output_file`.
#' @param overwrite Logical. If `TRUE`, recompute even if `output_file` exists.
#'
#' @return A `terra::SpatRaster`.
#'
#' @export
get_DEM <- function(
    input_dir = file.path(path_finder("ninaGeospatialFiles"),
                          "Elevation",
                          "Norway_DEM_50m_Kartverket",
                          "Original",
                          "DTM50_UTM33_20210416"),
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

  input_dir <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
  output_file <- normalizePath(output_file, winslash = "/", mustWork = FALSE)

  if (file.exists(output_file) && !overwrite) {
    message("Reading cached merged raster: ", output_file)
    return(terra::rast(output_file))
  }

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

  tif_files <- list.files(
    input_dir,
    pattern = "\\.tif$",
    full.names = TRUE,
    recursive = TRUE
  )

  if (length(tif_files) == 0) {
    stop("No .tif files found in input_dir: ", input_dir, call. = FALSE)
  }

  tif_files <- normalizePath(tif_files, winslash = "/", mustWork = TRUE)

  message("Found ", length(tif_files), " tif files.")
  message("Merging rasters...")

  r_collection <- terra::sprc(
    lapply(tif_files, terra::rast)
  )

  merged_raster <- terra::mosaic(r_collection)

  terra::writeRaster(
    merged_raster,
    output_file,
    overwrite = TRUE
  )

  if (is.null(metadata_file)) {
    metadata_file <- file.path(
      dirname(output_file),
      paste0(
        tools::file_path_sans_ext(basename(output_file)),
        "_metadata.yml"
      )
    )
  }

  metadata_file <- normalizePath(metadata_file, winslash = "/", mustWork = FALSE)
  dir.create(dirname(metadata_file), recursive = TRUE, showWarnings = FALSE)

  metadata <- list(
    title = "Merged raster tiles",
    description = "Multiple GeoTIFF raster tiles merged into one GeoTIFF raster.",
    created = as.character(Sys.time()),
    input_directory = input_dir,
    output_file = output_file,
    file_format = "GeoTIFF",
    n_input_files = length(tif_files),
    input_files = basename(tif_files),
    processing_method = "terra::sprc(lapply(tif_files, terra::rast)) followed by terra::mosaic()",
    raster_properties = list(
      crs = terra::crs(merged_raster),
      extent = as.character(terra::ext(merged_raster)),
      resolution = as.numeric(terra::res(merged_raster)),
      dimensions = c(
        nrow = terra::nrow(merged_raster),
        ncol = terra::ncol(merged_raster),
        nlyr = terra::nlyr(merged_raster)
      ),
      names = names(merged_raster)
    ),
    FAIR = list(
      Findable = "Descriptive output and metadata filenames are used.",
      Accessible = "Merged raster is stored as GeoTIFF.",
      Interoperable = "GeoTIFF is compatible with R terra, QGIS, GDAL, and Python rasterio.",
      Reusable = "Metadata records inputs, method, CRS, extent, resolution, dimensions, and creation date."
    )
  )

  tryCatch(
    yaml::write_yaml(metadata, metadata_file),
    error = function(e) {
      warning(
        "Merged raster was saved, but metadata could not be written: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  merged_raster
}
