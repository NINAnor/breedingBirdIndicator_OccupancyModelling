#' Process and Cache Naturverngrense Layer
#'
#' Reads the `Naturverngrense` layer from a File Geodatabase, transforms it to
#' the model CRS, clips it to the Norway boundary, and caches the processed
#' result. If `output_file` already exists and `overwrite = FALSE`, the cached
#' file is read instead.
#'
#' @param input_file Character. Path to the `.gdb` directory.
#' @param norway_boundary An `sf` object containing the Norway boundary.
#' @param modelCrs CRS object, EPSG code, or proj string used for transformation.
#' @param output_file Character. Path to output GeoParquet file.
#' @param layer Character. Layer name in the geodatabase.
#' @param overwrite Logical. If `TRUE`, reprocess and overwrite cached file.
#'
#' @return An `sf` object containing the processed Naturverngrense layer.
#'
#' @export
get_protected_areas <- function(
    input_file,
    norway_boundary,
    modelCrs,
    output_file,
    layer = "Naturverngrense",
    overwrite = FALSE
) {
  required_pkgs <- c("sf", "dplyr", "sfarrow")

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

  output_file <- normalizePath(output_file,
                               winslash = "/",
                               mustWork = FALSE)

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

  if (file.exists(output_file) && !overwrite) {
    message("Reading cached Naturverngrense file.")

    naturverngrense <- sfarrow::st_read_parquet(output_file)

  } else {
    message("Processing Naturverngrense layer.")

    norway_boundary <- norway_boundary |>
      sf::st_transform(modelCrs) |>
      sf::st_make_valid() |>
      sf::st_union()

    naturverngrense <- sf::st_read(
      dsn = input_file,
      layer = layer,
      quiet = TRUE
    ) |>
      sf::st_transform(modelCrs) |>
      sf::st_make_valid() |>
      sf::st_intersection(norway_boundary)

    sfarrow::st_write_parquet(
      naturverngrense,
      output_file
    )
  }

  naturverngrense
}
