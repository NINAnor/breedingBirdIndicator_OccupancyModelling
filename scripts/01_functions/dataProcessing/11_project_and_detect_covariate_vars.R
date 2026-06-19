#' Detect Whether a Raster is Categorical
#'
#' @param r A terra::SpatRaster.
#' @param max_unique Maximum number of unique values allowed for categorical data.
#' @param sample_size Number of raster cells to sample.
#'
#' @return TRUE if raster is likely categorical, otherwise FALSE.
is_categorical_raster <- function(
    r,
    max_unique = 100,
    sample_size = 10000
) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required.", call. = FALSE)
  }

  if (!inherits(r, "SpatRaster")) {
    stop("`r` must be a terra::SpatRaster.", call. = FALSE)
  }

  if (terra::is.factor(r)) {
    return(TRUE)
  }

  vals <- terra::spatSample(
    r,
    size = sample_size,
    method = "regular",
    na.rm = TRUE,
    values = TRUE,
    as.df = TRUE
  )

  if (nrow(vals) == 0) {
    return(FALSE)
  }

  vals <- vals[[1]]
  vals <- vals[!is.na(vals)]

  if (length(vals) == 0) {
    return(FALSE)
  }

  integer_like <- all(vals == round(vals), na.rm = TRUE)
  #few_unique <- length(unique(vals)) <= max_unique

  integer_like #&& few_unique
}

#' Resample a Raster to a Base Template Raster
#'
#' Aligns a raster to a base raster by matching CRS, extent, resolution,
#' origin, and dimensions. The resampling method is selected automatically:
#' nearest-neighbour for categorical rasters and bilinear interpolation for
#' continuous rasters.
#'
#' This supports FAIR processing by ensuring all covariates are spatially
#' interoperable and reproducibly aligned to a documented raster template.
#'
#' @param r A `terra::SpatRaster` or object coercible to `SpatRaster`.
#' @param baseRaster A `terra::SpatRaster` used as the spatial template.
#' @param name Character or `NULL`. Optional output layer name.
#' @param max_unique Integer. Passed to `is_categorical_raster()`.
#' @param sample_size Integer. Passed to `is_categorical_raster()`.
#'
#' @return A `terra::SpatRaster` aligned to `baseRaster`.
#'
#' @export
resample_to_baseRaster <- function(
    r,
    baseRaster,
    name = NULL,
    max_unique = 100,
    sample_size = 10000,
    scale_continuous = TRUE
) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required.", call. = FALSE)
  }

  if (!inherits(r, "SpatRaster")) {
    r <- terra::rast(r)
  }

  if (!inherits(baseRaster, "SpatRaster")) {
    stop("`baseRaster` must be a terra::SpatRaster.", call. = FALSE)
  }

  is_categorical <- is_categorical_raster(
    r,
    max_unique = max_unique,
    sample_size = sample_size
  )

  method <- if (is_categorical) "near" else "bilinear"

  if (!terra::same.crs(r, baseRaster)) {
    r <- terra::project(r, baseRaster, method = method)
  }

  if (!terra::compareGeom(r, baseRaster, stopOnError = FALSE)) {
    r <- terra::resample(r, baseRaster, method = method)
  }


  if (scale_continuous && !is_categorical) {
    message("Scaling ", ifelse(is.null(name), names(r)[1], name), " covariate")

    r_mean <- terra::global(r, "mean", na.rm = TRUE)[1, 1]
    r_sd   <- terra::global(r, "sd", na.rm = TRUE)[1, 1]

    if (is.na(r_sd) || r_sd == 0) {
      warning("Raster could not be scaled because SD is zero or NA.", call. = FALSE)
    } else {
      r <- (r - r_mean) / r_sd
    }
  }

  if (!is.null(name)) {
    names(r) <- name
  }

  r
}

