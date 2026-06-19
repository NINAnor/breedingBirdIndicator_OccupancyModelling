#' Summarize posterior occupancy predictions by species
#'
#' @description
#' Converts posterior occupancy prediction samples into spatial summaries by
#' species and location. Outputs can be returned as an `sf` object or as a named
#' list of `terra::SpatRaster` objects. When `output_dir` is supplied, results
#' are written to species-specific subdirectories to support FAIR data handling.
#'
#' @details
#' FAIR-oriented practices supported by this function include:
#'
#' **Findable**: outputs are organized by species name in predictable
#' subdirectories and use descriptive filenames.
#'
#' **Accessible**: outputs are written in common geospatial formats:
#' GeoParquet for vector data and GeoTIFF for raster data.
#'
#' **Interoperable**: spatial outputs include coordinate reference system
#' information through `sf` or `terra`.
#'
#' **Reusable**: outputs include standardized summary fields:
#' `mean`, `median`, `sd`, `lower`, and `upper`, with quantile bounds controlled
#' by `probs`.
#'
#' @param pred A prediction object containing `psi.0.samples`, a 3-dimensional
#'   numeric array with dimensions samples × species × locations.
#' @param coords A two-column matrix or data frame of prediction coordinates.
#'   The first column is interpreted as x-coordinate and the second as
#'   y-coordinate. Must contain one row per prediction location.
#' @param species_names Optional character vector of species names. Must have
#'   one entry per species in `pred$psi.0.samples`. If `NULL`, names are generated
#'   as `"species_1"`, `"species_2"`, and so on.
#' @param template Optional `terra::SpatRaster` template used when
#'   `return = "raster"`. Defines the raster geometry, resolution, extent, and CRS.
#' @param return Character. Output type, either `"sf"` or `"raster"`.
#' @param crs Coordinate reference system used when `return = "sf"`. Required
#'   for `sf` output.
#' @param probs Numeric vector of length 2 giving lower and upper posterior
#'   quantiles. Defaults to `c(0.025, 0.975)`.
#' @param output_dir Optional directory where results are written. If supplied,
#'   a subdirectory is created for each species.
#'
#' @return
#' If `return = "sf"`, an `sf` object containing one row per species-location
#' combination.
#'
#' If `return = "raster"`, a named list of `terra::SpatRaster` objects, with one
#' raster stack per species.
#'
#' @section Output files:
#' When `output_dir` is provided, the function writes:
#'
#' - `output_dir/<species>/<species>_occupancy_predictions.parquet` for `sf`
#'   output.
#' - `output_dir/<species>/<species>_occupancy_predictions.tif` for raster output.
#'
#' Species names are sanitized for directory and file names by replacing
#' non-alphanumeric characters with underscores.
#'
#' @section Recommended metadata:
#' For stronger FAIR compliance, store an accompanying metadata file alongside
#' each output documenting:
#'
#' - model name and version
#' - prediction date
#' - coordinate reference system
#' - meaning of each summary field
#' - quantile probabilities used
#' - software versions for `sf`, `terra`, and R
#' - source data citation or DOI, where available
#'
#' @examples
#' \dontrun{
#' occ_sf <- summarize_occupancy_predictions(
#'   pred = pred,
#'   coords = coords,
#'   species_names = species_names,
#'   return = "sf",
#'   crs = 4326,
#'   output_dir = "outputs/occupancy_predictions"
#' )
#'
#' occ_rasters <- summarize_occupancy_predictions(
#'   pred = pred,
#'   coords = coords,
#'   species_names = species_names,
#'   template = template,
#'   return = "raster",
#'   output_dir = "outputs/occupancy_predictions"
#' )
#'
#' occ_rasters[["Species name"]]
#' }
#'
#' @export
summarize_occupancy_predictions <- function(
    pred,
    coords,
    species_names = NULL,
    template = NULL,
    return = c("sf", "raster"),
    crs = NULL,
    probs = c(0.025, 0.975),
    output_dir = NULL
) {
  return <- match.arg(return)

  psi <- pred$psi.0.samples

  if (length(dim(psi)) != 3) {
    stop("`pred$psi.0.samples` must be samples × species × locations.", call. = FALSE)
  }

  n_species <- dim(psi)[2]
  n_locations <- dim(psi)[3]

  if (is.null(species_names)) {
    species_names <- paste0("species_", seq_len(n_species))
  }

  if (length(species_names) != n_species) {
    stop("`species_names` must have one name per species.", call. = FALSE)
  }

  coords <- as.data.frame(coords)

  if (nrow(coords) != n_locations) {
    stop("`coords` must have one row per prediction location.", call. = FALSE)
  }

  if (!is.numeric(probs) || length(probs) != 2 || any(probs < 0 | probs > 1)) {
    stop("`probs` must be a numeric vector of length 2 between 0 and 1.", call. = FALSE)
  }

  safe_species_names <- gsub("[^A-Za-z0-9_-]+", "_", species_names)
  safe_species_names <- gsub("^_+|_+$", "", safe_species_names)

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  summaries <- lapply(seq_len(n_species), function(i) {
    m <- psi[, i, , drop = FALSE]
    m <- matrix(m, nrow = dim(psi)[1], ncol = n_locations)

    q <- apply(m, 2, stats::quantile, probs = probs, na.rm = TRUE)

    data.frame(
      species = species_names[i],
      location_id = seq_len(n_locations),
      x = coords[[1]],
      y = coords[[2]],
      mean = colMeans(m, na.rm = TRUE),
      median = apply(m, 2, stats::median, na.rm = TRUE),
      sd = apply(m, 2, stats::sd, na.rm = TRUE),
      lower = q[1, ],
      upper = q[2, ],
      stringsAsFactors = FALSE
    )
  })

  out_df <- do.call(rbind, summaries)

  if (return == "sf") {
    if (is.null(crs)) {
      stop("Provide `crs` when returning sf.", call. = FALSE)
    }

    out_sf <- sf::st_as_sf(
      out_df,
      coords = c("x", "y"),
      crs = crs,
      remove = FALSE
    )

    if (!is.null(output_dir)) {
      for (i in seq_along(species_names)) {
        sp <- species_names[i]
        sp_safe <- safe_species_names[i]

        sp_dir <- file.path(output_dir, sp_safe)
        dir.create(sp_dir, recursive = TRUE, showWarnings = FALSE)

        sp_sf <- out_sf[out_sf$species == sp, ]

        sfarrow::st_write_parquet(
          obj = sp_sf,
          dsn = file.path(sp_dir, paste0(sp_safe, "_occupancy_predictions.parquet")),
          compression = "snappy"
        )
      }
    }

    return(out_sf)
  }

  if (return == "raster") {
    if (is.null(template)) {
      stop("Provide `template` when returning raster.", call. = FALSE)
    }

    rasters <- lapply(seq_along(species_names), function(i) {
      sp <- species_names[i]
      sp_df <- out_df[out_df$species == sp, ]

      sp_vect <- terra::vect(
        sp_df,
        geom = c("x", "y"),
        crs = terra::crs(template)
      )

      r <- terra::rast(
        lapply(
          c("mean", "median", "sd", "lower", "upper"),
          function(field) {
            terra::rasterize(sp_vect, template, field = field)
          }
        )
      )

      names(r) <- c("mean", "median", "sd", "lower", "upper")
      r
    })

    names(rasters) <- species_names

    if (!is.null(output_dir)) {
      for (i in seq_along(species_names)) {
        sp <- species_names[i]
        sp_safe <- safe_species_names[i]

        sp_dir <- file.path(output_dir, sp_safe)
        dir.create(sp_dir, recursive = TRUE, showWarnings = FALSE)

        terra::writeRaster(
          rasters[[sp]],
          filename = file.path(sp_dir, paste0(sp_safe, "_occupancy_predictions.tif")),
          overwrite = TRUE
        )
      }
    }

    return(rasters)
  }
}
