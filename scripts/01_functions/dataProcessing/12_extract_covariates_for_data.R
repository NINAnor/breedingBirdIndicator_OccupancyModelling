#' Extract Model Covariates to Prediction or Survey Points
#'
#' Extracts raster covariates to point locations supplied either as an `sf`
#' object or as a data frame with longitude and latitude columns.
#'
#' @param x An `sf` point object or data frame containing coordinates.
#' @param modelCrs CRS used by the model and covariate rasters.
#' @param ecosystemBreedingBird Object required by `get_landcover_labels()`.
#' @param lon_col Character. Longitude/x-coordinate column if `x` is a data frame.
#' @param lat_col Character. Latitude/y-coordinate column if `x` is a data frame.
#' @param input_crs CRS of `lon_col` and `lat_col` if `x` is a data frame.
#'   Default is `4326`.
#' @param return Character. Either `"sf"` or `"data.frame"`.
#'
#' @return An `sf` object or data frame containing original variables,
#' coordinates, and extracted covariates.
#'
#' @export
extract_model_covariates <- function(
    x,
    modelCrs,
    all_covariates,
    lon_col = "Longitude",
    lat_col = "Latitude",
    input_crs = 4326,
    covariate_stack = NULL,
    return = c("sf", "data.frame")
) {
  return <- match.arg(return)

  required_pkgs <- c("sf", "terra", "dplyr")

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




  if (inherits(x, "sf")) {
    points_sf <- x |>
      sf::st_transform(modelCrs)
  } else {
    if (!all(c(lon_col, lat_col) %in% names(x))) {
      stop(
        "The data frame must contain columns: ",
        lon_col,
        " and ",
        lat_col,
        call. = FALSE
      )
    }

    missing_coords <- is.na(x[[lon_col]]) | is.na(x[[lat_col]])

    if (any(missing_coords)) {
      warning(
        sum(missing_coords),
        " rows have missing coordinates and will be removed.",
        call. = FALSE
      )

      x <- x[!missing_coords, , drop = FALSE]
    }

    points_sf <- sf::st_as_sf(
      x,
      coords = c(lon_col, lat_col),
      crs = input_crs,
      remove = FALSE
    ) |>
      sf::st_transform(modelCrs)
  }

  extracted <- terra::extract(
    covariate_stack,
    terra::vect(points_sf)
  )

  extracted <- extracted |>
    dplyr::select(-ID)

  coords <- sf::st_coordinates(points_sf)

  out <- points_sf |>
    dplyr::mutate(
      long = coords[, 1],
      lat = coords[, 2]
    ) |>
    dplyr::bind_cols(extracted)

  if (return == "data.frame") {
    out <- sf::st_drop_geometry(out)
  }

  out
}
