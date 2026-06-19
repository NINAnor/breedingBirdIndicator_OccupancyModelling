#' Build and Align Model Covariate Raster Stack
#'
#' Loads requested model covariates, aligns them to a common template raster,
#' and returns a single `terra::SpatRaster` stack.
#'
#' This function supports FAIR principles by:
#' \itemize{
#'   \item \strong{Findable}: covariates are loaded from explicit file paths.
#'   \item \strong{Accessible}: outputs are standard `terra::SpatRaster` objects.
#'   \item \strong{Interoperable}: all rasters are aligned to one CRS, extent,
#'         resolution, and origin.
#'   \item \strong{Reusable}: the function records deterministic processing steps
#'         through explicit arguments and consistent layer names.
#' }
#'
#' @param x Object passed to `get_landcover_labels()`, usually the ecosystem or
#'   study object used to derive land-cover covariates.
#' @param modelCrs CRS used by the model and covariate rasters.
#' @param all_covariates Character vector of requested covariate names, e.g.
#'   `c("DEM", "LABEL1")`.
#' @param template A `terra::SpatRaster` used as the base raster for CRS,
#'   extent, resolution, and origin.
#' @param cov_folder Character. Folder containing land-cover covariate files.
#' @param data_folder Character. Base data folder. If `NULL`, uses
#'   `path_finder("dataFolder")`.
#'
#' @return A `terra::SpatRaster` containing aligned covariate layers.
#'
#' @export
put_covariates_together <- function(
    x,
    modelCrs,
    all_covariates,
    template,
    output_file,
    cov_folder = file.path(path_finder("dataFolder"), "Covariates"),
    data_folder = path_finder("dataFolder"),
    overwrite = FALSE,
    scale_vars = TRUE
) {
  required_pkgs <- c("sf", "terra", "dplyr", "tidyterra")

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

  if (!inherits(template, "SpatRaster")) {
    stop("`template` must be a terra::SpatRaster.", call. = FALSE)
  }

  if (!is.character(all_covariates)) {
    stop("`all_covariates` must be a character vector.", call. = FALSE)
  }

  output_file <- normalizePath(output_file, winslash = "/", mustWork = FALSE)

  if (file.exists(output_file) && !overwrite) {
    message("Found cached covariate stack: ", output_file)

    cached_stack <- terra::rast(output_file)

    same_as_template <- terra::compareGeom(
      cached_stack,
      template,
      stopOnError = FALSE
    )

    matched_cov_names <- all(all_covariates %in% names(cached_stack))

    if (same_as_template && matched_cov_names) {
      message("Cached covariate stack matches template and requested covariates.")
      return(cached_stack[[all_covariates]])
    }

    if (!same_as_template) {
      message("Cached covariate stack does not match template.")
    }

    if (!matched_cov_names) {
      missing_covariates <- setdiff(all_covariates, names(cached_stack))
      message(
        "Cached covariate stack is missing: ",
        paste(missing_covariates, collapse = ", ")
      )
    }

    message("Rebuilding covariate stack.")
  }

  covariate_stack <- NULL

  if ("DEM" %in% all_covariates) {
    dem <- get_DEM(
      output_file = file.path(
        data_folder,
        "processedData",
        "DEM",
        "DEM.tif"
      )
    )

    dem <- resample_to_baseRaster(
      r = dem,
      baseRaster = template,
      name = "DEM",
      scale_continuous = scale_vars
    )

    covariate_stack <- dem
  }

  if ("landCover" %in% all_covariates) {
    landCoverVars <- get_landcover_labels(
      x,
      modelCrs = modelCrs,
      cov_folder = cov_folder,
      output = "raster"
    ) |>
      tidyterra::select(LABEL1)

    landCoverVars <- resample_to_baseRaster(
      r = landCoverVars,
      baseRaster = template,
      name = "landCover"
    )

    covariate_stack <- if (is.null(covariate_stack)) {
      landCoverVars
    } else {
      c(covariate_stack, landCoverVars)
    }
  }

  if (is.null(covariate_stack)) {
    stop(
      "No supported covariates found in `all_covariates`. Supported values currently include: DEM, LABEL1.",
      call. = FALSE
    )
  }

  #covariate_stack <- terra::rast(covariate_stack)
  names(covariate_stack) <- make.unique(names(covariate_stack))

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

  terra::writeRaster(
    covariate_stack,
    output_file,
    overwrite = TRUE
  )

  covariate_stack
}
