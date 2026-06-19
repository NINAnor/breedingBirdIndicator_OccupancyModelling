#' Create a Template Raster for a Study Region
#'
#' Creates, saves, or loads a template raster covering a study area with a
#' specified resolution and CRS.
#'
#' @param study_region An sf polygon object.
#' @param resolution Numeric. Cell size in map units. Default 1000 m.
#' @param crs CRS of output raster. Defaults to the CRS of `study_region`.
#' @param output_file Optional file path for saving/loading the template raster.
#' @param overwrite Logical. If `TRUE`, recreate and overwrite `output_file`.
#'   Default `FALSE`.
#'
#' @return A `terra::SpatRaster`.
#'
#' @export
create_template_raster <- function(
    study_region,
    resolution = 1000,
    crs = NULL,
    output_file = NULL,
    overwrite = FALSE
) {
  if (!inherits(study_region, "sf")) {
    stop("`study_region` must be an sf object.", call. = FALSE)
  }

  if (is.null(crs)) {
    crs <- sf::st_crs(study_region)
  } else {
    crs <- sf::st_crs(crs)
  }

  if (is.na(crs)) {
    stop("CRS is missing or invalid.", call. = FALSE)
  }

  if (!is.numeric(resolution) || length(resolution) != 1 || resolution <= 0) {
    stop("`resolution` must be a single positive numeric value.", call. = FALSE)
  }

  if (!is.null(output_file) && file.exists(output_file) && !overwrite) {
    template <- terra::rast(output_file)

    raster_res <- terra::res(template)

    if (!all(abs(raster_res - resolution) < .Machine$double.eps^0.5)) {
      stop(
        "`output_file` already exists, but its resolution is ",
        paste(raster_res, collapse = " × "),
        "; expected ",
        resolution,
        ". Use `overwrite = TRUE` to recreate it.",
        call. = FALSE
      )
    }

    return(template)
  }

  study_region <- study_region |>
    sf::st_transform(crs) |>
    sf::st_make_valid()

  template <- terra::rast(
    ext = terra::ext(terra::vect(study_region)),
    resolution = resolution,
    crs = crs$wkt
  )

  terra::values(template) <- 1

  template <- terra::mask(
    template,
    terra::vect(study_region)
  )

  names(template) <- "template"

  if (!is.null(output_file)) {
    dir.create(
      dirname(output_file),
      recursive = TRUE,
      showWarnings = FALSE
    )

    terra::writeRaster(
      template,
      output_file,
      overwrite = overwrite
    )
  }

  template
}
