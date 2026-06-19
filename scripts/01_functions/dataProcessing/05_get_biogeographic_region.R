#' Process and Cache Biogeographic Regions
#'
#' Reads a biogeographic region layer, clips it to the Norway boundary,
#' transforms it to the model coordinate reference system, and caches the
#' result as a GeoParquet file. If the cached file exists and
#' `overwrite = FALSE`, the cached file is returned.
#'
#' @param input_file Character. Path to the input shapefile or vector file.
#' @param norway_boundary An `sf` object containing the Norway boundary.
#' @param modelCrs Coordinate reference system to transform the output to.
#' @param output_file Character. Path to the cached GeoParquet file.
#' @param overwrite Logical. If `TRUE`, overwrite existing output.
#'
#' @return An `sf` object with one column:
#' \describe{
#'   \item{bioGeographic_region}{Biogeographic region identifier}
#' }
#'
#' @details
#' This function follows FAIR principles:
#' \itemize{
#'   \item Findable: descriptive file names.
#'   \item Accessible: GeoParquet output.
#'   \item Interoperable: returns an `sf` object.
#'   \item Reusable: processing steps are deterministic.
#' }
#'
#' @examples
#' \dontrun{
#' region <- process_biogeographic_region(
#'   input_file = "biogeo.shp",
#'   norway_boundary = norway_boundary,
#'   modelCrs = modelCrs,
#'   output_file = "biogeographic_region.parquet"
#' )
#' }
#'
#' @export
get_biogeographic_region <- function(
    input_file,
    norway_boundary,
    modelCrs,
    output_file,
    overwrite = FALSE
) {

  required_pkgs <- c(
    "sf",
    "dplyr",
    "sfarrow"
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
      "Missing packages: ",
      paste(missing_pkgs, collapse = ", "),
      call. = FALSE
    )
  }

  output_file <- normalizePath(output_file,
                               winslash = "/",
                               mustWork = FALSE)

  dir.create(
    dirname(output_file),
    recursive = TRUE,
    showWarnings = FALSE
  )

  if (file.exists(output_file) && !overwrite) {

    message("Reading cached biogeographic regions.")

    region <- sfarrow::st_read_parquet(
      output_file
    )

  } else {

    message("Processing biogeographic regions.")

    # Ensure common CRS
    norway_boundary <- sf::st_transform(
      norway_boundary,
      modelCrs
    ) |>
      sf::st_union() |>
      sf::st_make_valid()

    region <- sf::st_read(
      input_file,
      quiet = TRUE
    ) |>
      dplyr::select(
        biogeo
      ) |>
      sf::st_transform(
        modelCrs
      ) |>
      sf::st_make_valid() |>
      sf::st_intersection(
        norway_boundary
      ) |>
      dplyr::rename(
        bioGeographic_region = biogeo
      ) |>
      dplyr::select(
        bioGeographic_region
      )

    sfarrow::st_write_parquet(
      region,
      output_file
    )

  }

  return(region)
}
