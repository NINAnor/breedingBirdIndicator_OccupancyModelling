#' Process and Cache Mountain Layer Clipped to Norway
#'
#' Reads a vector layer, transforms it to the model CRS, makes geometries valid,
#' clips it to the Norway boundary, and caches the processed result as
#' GeoParquet. If the cached file already exists and `overwrite = FALSE`, the
#' file is loaded directly.
#'
#' This function supports FAIR principles:
#' \itemize{
#'   \item Findable: descriptive output and metadata file names.
#'   \item Accessible: cached output is stored as GeoParquet.
#'   \item Interoperable: output is returned as an `sf` object with CRS.
#'   \item Reusable: metadata records source, layer, CRS, processing steps, and date.
#' }
#'
#' @param input_file Character. Path to input data source, for example a `.gdb`.
#' @param layer Character. Layer name to read from `input_file`.
#' @param norway_boundary An `sf` object containing Norway's boundary.
#' @param modelCrs CRS object, EPSG code, or proj string used for transformation.
#' @param output_file Character. Path to output GeoParquet file.
#' @param metadata_file Character or `NULL`. Path to YAML metadata file.
#'   If `NULL`, metadata is written beside `output_file`.
#' @param overwrite Logical. If `TRUE`, reprocess even if output exists.
#'
#' @return An `sf` object containing the processed and clipped layer.
#'
#' @export
get_mountain_layer <- function(
    input_file,
    layer,
    norway_boundary,
    modelCrs,
    output_file,
    metadata_file = NULL,
    overwrite = FALSE
) {
  required_pkgs <- c("sf", "sfarrow", "yaml")

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

  output_file <- normalizePath(
    output_file,
    winslash = "/",
    mustWork = FALSE
  )

  output_exists <- file.exists(output_file) || dir.exists(output_file)

  if (output_exists && !overwrite) {
    message("Reading cached processed layer: ", output_file)

    return(
      sfarrow::st_read_parquet(output_file)
    )
  }

  if (!inherits(norway_boundary, "sf")) {
    stop("`norway_boundary` must be an sf object.", call. = FALSE)
  }

  dir.create(
    dirname(output_file),
    recursive = TRUE,
    showWarnings = FALSE
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

  metadata_file <- normalizePath(
    metadata_file,
    winslash = "/",
    mustWork = FALSE
  )

  message("Processing layer: ", layer)

  norway_boundary_processed <- norway_boundary |>
    sf::st_transform(modelCrs) |>
    sf::st_make_valid() |>
    sf::st_union()

  mountain <- sf::st_read(
    dsn = input_file,
    layer = layer,
    quiet = TRUE
  ) |>
    sf::st_transform(modelCrs) |>
    sf::st_make_valid() |>
    sf::st_intersection(norway_boundary_processed)

  sfarrow::st_write_parquet(
    mountain,
    output_file
  )

  metadata <- list(
    title = paste("Processed", layer, "layer clipped to Norway"),
    description = "Vector layer transformed to model CRS, validated, clipped to Norway boundary, and cached as GeoParquet.",
    created = as.character(Sys.time()),
    input_file = input_file,
    input_layer = layer,
    output_file = output_file,
    file_format = "GeoParquet",
    processing_steps = list(
      "Transform Norway boundary to model CRS",
      "Make Norway boundary geometry valid",
      "Union Norway boundary",
      "Read input layer",
      "Transform input layer to model CRS",
      "Make input layer geometry valid",
      "Intersect input layer with Norway boundary",
      "Write processed layer to GeoParquet"
    ),
    coordinate_reference_system = list(
      output = sf::st_crs(mountain)$input
    ),
    n_features = nrow(mountain),
    FAIR = list(
      Findable = "Output and metadata use descriptive names.",
      Accessible = "Output is stored in open GeoParquet format.",
      Interoperable = "Output is an sf-compatible GeoParquet layer with CRS information.",
      Reusable = "Metadata records input source, layer name, CRS, processing steps, creation date, and feature count."
    )
  )

  dir.create(
    dirname(metadata_file),
    recursive = TRUE,
    showWarnings = FALSE
  )

  tryCatch(
    yaml::write_yaml(
      metadata,
      metadata_file
    ),
    error = function(e) {
      warning(
        "Processed layer was saved, but metadata could not be written: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  mountain
}
