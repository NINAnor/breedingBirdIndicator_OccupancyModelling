#' Sample Prediction Points from Raster Classes
#'
#' Samples prediction points from a `terra::SpatRaster`, optionally restricted
#' to selected raster cell values such as land-cover codes. Results are cached
#' as GeoParquet and accompanied by YAML metadata to support FAIR data use.
#'
#' @param raster A `terra::SpatRaster`.
#' @param selected_codes Optional vector of raster values to sample from.
#'   If `NULL`, all non-NA raster cells are eligible.
#' @param n_points Integer. Number of prediction points to sample.
#' @param output_file Character. Path to output GeoParquet file.
#' @param metadata_file Character or `NULL`. Path to YAML metadata file.
#'   If `NULL`, metadata is written beside `output_file`.
#' @param overwrite Logical. If `TRUE`, recompute even if cached output exists.
#' @param seed Integer or `NULL`. Random seed for reproducible sampling.
#' @param sample_method Character. `"random"` or `"regular"`.
#'
#' @return An `sf` point object with `site_id`, `cell`, `x`, `y`, and raster value.
#'
#' @export
sample_prediction_points_from_raster <- function(
    raster,
    selected_codes = NULL,
    n_points = 10000,
    output_file,
    metadata_file = NULL,
    overwrite = FALSE,
    seed = 123,
    sample_method = c("random", "regular")
) {
  sample_method <- match.arg(sample_method)

  required_pkgs <- c("terra", "sf", "dplyr", "tibble", "sfarrow", "yaml")

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
    message("Reading cached prediction points: ", output_file)
    return(sfarrow::st_read_parquet(output_file))
  }

  if (!inherits(raster, "SpatRaster")) {
    stop("`raster` must be a terra::SpatRaster.", call. = FALSE)
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

  values_vec <- terra::values(raster, mat = FALSE)

  eligible_cells <- which(!is.na(values_vec))

  if (!is.null(selected_codes)) {
    eligible_cells <- which(
      !is.na(values_vec) & values_vec %in% selected_codes
    )
  }

  if (length(eligible_cells) == 0) {
    stop("No eligible raster cells found.", call. = FALSE)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (length(eligible_cells) <= n_points) {
    selected_cells <- eligible_cells
    warning(
      "Requested more points than available cells. Returning all eligible cells."
    )
  } else if (sample_method == "random") {
    selected_cells <- sample(eligible_cells, n_points)
  } else {
    selected_cells <- eligible_cells[
      round(seq(1, length(eligible_cells), length.out = n_points))
    ]
  }

  xy <- terra::xyFromCell(raster, selected_cells)
  raster_values <- values_vec[selected_cells]

  pts <- tibble::tibble(
    site_id = seq_along(selected_cells),
    cell = selected_cells,
    x = xy[, 1],
    y = xy[, 2],
    raster_value = raster_values
  ) |>
    sf::st_as_sf(
      coords = c("x", "y"),
      crs = terra::crs(raster)
    ) |>
    dplyr::mutate(
      x = sf::st_coordinates(geometry)[, 1],
      y = sf::st_coordinates(geometry)[, 2]
    ) |>
    dplyr::select(
      site_id,
      cell,
      x,
      y,
      raster_value,
      geometry
    )

  sfarrow::st_write_parquet(
    pts,
    output_file
  )

  metadata <- list(
    title = "Prediction points sampled from raster",
    description = "Spatial prediction points sampled from eligible raster cells.",
    created = as.character(Sys.time()),
    output_file = output_file,
    file_format = "GeoParquet",
    raster_layer = names(raster),
    raster_crs = terra::crs(raster),
    raster_extent = as.character(terra::ext(raster)),
    raster_resolution = terra::res(raster),
    selected_codes = if (is.null(selected_codes)) {
      "all non-NA cells"
    } else {
      selected_codes
    },
    n_requested = n_points,
    n_returned = nrow(pts),
    sample_method = sample_method,
    seed = seed,
    FAIR = list(
      Findable = "Descriptive file and metadata names are used.",
      Accessible = "Prediction points are stored in open GeoParquet format.",
      Interoperable = "Output is readable as sf, GeoParquet, and Arrow-compatible data.",
      Reusable = "Metadata records raster properties, selected codes, CRS, sampling method, and date."
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
        "Prediction points were saved, but metadata could not be written: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  pts
}
