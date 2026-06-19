#' Summarize posterior occupancy predictions by species
#'
#' @description
#' Converts posterior occupancy prediction samples into spatial summaries by
#' species and location. Handles output from both single-species models
#' (`PGOcc`, `spPGOcc`) and multi-species models (`msPGOcc`, `spMsPGOcc`).
#' Outputs can be returned as an `sf` object or as a named list of
#' `terra::SpatRaster` objects. When `output_dir` is supplied, results are
#' written to species-specific subdirectories to support FAIR data handling.
#'
#' @details
#'
#' ## Single-species vs multi-species input
#'
#' The shape of `pred$psi.0.samples` determines which path is taken:
#'
#' - **Multi-species** (`msPGOcc`, `spMsPGOcc`): 3-D array with dimensions
#'   samples × species × locations. `species_names` must have one entry per
#'   species (second dimension).
#' - **Single-species** (`PGOcc`, `spPGOcc`): 2-D matrix with dimensions
#'   samples × locations. Exactly one species name must be provided via
#'   `species_names`.
#'
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
#' `mean`, `median`, `sd`, `lower`, and `upper`, with quantile bounds
#' controlled by `probs`.
#'
#' @param pred A prediction object. Must contain `psi.0.samples`, which is
#'   either:
#'   - a 3-D array (samples × species × locations) for multi-species models, or
#'   - a 2-D matrix (samples × locations) for single-species models.
#' @param coords A two-column matrix or data frame of prediction coordinates.
#'   The first column is interpreted as x and the second as y. Must contain
#'   one row per prediction location.
#' @param species_names Character vector of species names.
#'   - For multi-species input: must have one entry per species (second
#'     dimension of `psi.0.samples`). If `NULL`, names are generated as
#'     `"species_1"`, `"species_2"`, and so on.
#'   - For single-species input: must be a length-1 character vector. There is
#'     no dimension to infer from, so this is required. An error is raised if
#'     it is missing or has length other than 1.
#' @param template Optional `terra::SpatRaster` template used when
#'   `return = "raster"`. Defines raster geometry, resolution, extent, and CRS.
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
#' combination, with columns `species`, `location_id`, `x`, `y`, `mean`,
#' `median`, `sd`, `lower`, and `upper`.
#'
#' If `return = "raster"`, a named list of `terra::SpatRaster` objects, with
#' one raster stack per species containing layers `mean`, `median`, `sd`,
#' `lower`, and `upper`.
#'
#' @section Output files:
#' When `output_dir` is provided, the function writes:
#'
#' - `output_dir/<species>/<species>_occupancy_predictions.parquet` for `sf`
#'   output.
#' - `output_dir/<species>/<species>_occupancy_predictions.tif` for raster
#'   output.
#'
#' Species names are sanitized for use in directory and file names by replacing
#' non-alphanumeric characters with underscores.
#'
#' @section Recommended metadata:
#' For stronger FAIR compliance, store an accompanying metadata file alongside
#' each output documenting:
#'
#' - model name, type (single/multi-species, spatial/non-spatial), and version
#' - prediction date
#' - coordinate reference system
#' - meaning of each summary field
#' - quantile probabilities used
#' - software versions for `sf`, `terra`, and R
#' - source data citation or DOI, where available
#'
#' @examples
#' \dontrun{
#' # Multi-species spatial model — psi.0.samples is samples x species x locations
#' occ_sf <- summarize_occupancy_predictions(
#'   pred         = pred,
#'   coords       = coords,
#'   species_names = species_names,
#'   return       = "sf",
#'   crs          = 4326,
#'   output_dir   = "outputs/occupancy_predictions"
#' )
#'
#' # Single-species model — psi.0.samples is samples x locations
#' occ_sf_sp <- summarize_occupancy_predictions(
#'   pred          = pred_single,
#'   coords        = coords,
#'   species_names = "Turdus merula",
#'   return        = "sf",
#'   crs           = 4326,
#'   output_dir    = "outputs/occupancy_predictions"
#' )
#'
#' # Raster output
#' occ_rasters <- summarize_occupancy_predictions(
#'   pred          = pred,
#'   coords        = coords,
#'   species_names = species_names,
#'   template      = template,
#'   return        = "raster",
#'   output_dir    = "outputs/occupancy_predictions"
#' )
#'
#' occ_rasters[["Turdus merula"]]
#' }
#'
#' @export
summarize_occupancy_predictions <- function(
    pred,
    coords,
    species_names = NULL,
    template      = NULL,
    return        = c("sf", "raster"),
    crs           = NULL,
    probs         = c(0.025, 0.975),
    output_dir    = NULL
) {
  return <- match.arg(return)

  psi <- pred$psi.0.samples

  if (is.null(psi)) {
    stop("`pred$psi.0.samples` is NULL. Check that prediction was run successfully.",
         call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Detect model type from array dimensions and coerce to a common 3-D shape
  # ---------------------------------------------------------------------------
  #
  # Multi-species (msPGOcc / spMsPGOcc): psi.0.samples is samples x species x locations
  # Single-species (PGOcc / spPGOcc)   : psi.0.samples is samples x locations
  #
  # After this block, `psi` is always 3-D: samples x species x locations,
  # and `n_species`, `n_samples`, `n_locations` are set.

  nd <- length(dim(psi))

  if (nd == 3) {
    # --- Multi-species ---
    n_samples   <- dim(psi)[1]
    n_species   <- dim(psi)[2]
    n_locations <- dim(psi)[3]

    if (is.null(species_names)) {
      species_names <- paste0("species_", seq_len(n_species))
    }

    if (length(species_names) != n_species) {
      stop(
        "`species_names` must have one name per species in `psi.0.samples` ",
        "(", n_species, " species found).",
        call. = FALSE
      )
    }

  } else if (nd == 2) {
    # --- Single-species ---
    n_samples   <- dim(psi)[1]
    n_locations <- dim(psi)[2]
    n_species   <- 1L

    if (is.null(species_names)) {
      stop(
        "`pred$psi.0.samples` is 2-D (single-species output). ",
        "Provide `species_names` as a length-1 character vector, e.g. ",
        "`species_names = \"Turdus merula\"`.",
        call. = FALSE
      )
    }

    if (length(species_names) != 1L) {
      stop(
        "`pred$psi.0.samples` is 2-D (single-species output) but ",
        "`species_names` has length ", length(species_names), ". ",
        "Supply exactly one species name.",
        call. = FALSE
      )
    }

    # Promote to 3-D so downstream code is uniform: 1 x samples x locations -> reorder
    psi <- array(psi, dim = c(n_samples, 1L, n_locations))

  } else {
    stop(
      "`pred$psi.0.samples` must be a 2-D matrix (single-species) or a ",
      "3-D array (multi-species), but has ", nd, " dimensions.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Validate coordinates
  # ---------------------------------------------------------------------------

  coords <- as.data.frame(coords)

  if (nrow(coords) != n_locations) {
    stop(
      "`coords` has ", nrow(coords), " rows but `psi.0.samples` has ",
      n_locations, " locations. They must match.",
      call. = FALSE
    )
  }

  if (!is.numeric(probs) || length(probs) != 2 || any(probs < 0 | probs > 1)) {
    stop("`probs` must be a numeric vector of length 2 with values in [0, 1].",
         call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Sanitize species names for file paths
  # ---------------------------------------------------------------------------

  safe_species_names <- gsub("[^A-Za-z0-9_-]+", "_", species_names)
  safe_species_names <- gsub("^_+|_+$", "", safe_species_names)

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Compute per-species summaries
  # psi is now always 3-D: samples x species x locations
  # ---------------------------------------------------------------------------

  summaries <- lapply(seq_len(n_species), function(i) {

    # Extract samples x locations matrix for species i
    m <- matrix(psi[, i, ], nrow = n_samples, ncol = n_locations)

    q <- apply(m, 2, stats::quantile, probs = probs, na.rm = TRUE)

    data.frame(
      species     = species_names[i],
      location_id = seq_len(n_locations),
      x           = coords[[1]],
      y           = coords[[2]],
      mean        = colMeans(m, na.rm = TRUE),
      median      = apply(m, 2, stats::median, na.rm = TRUE),
      sd          = apply(m, 2, stats::sd,     na.rm = TRUE),
      lower       = q[1, ],
      upper       = q[2, ],
      stringsAsFactors = FALSE
    )
  })

  out_df <- do.call(rbind, summaries)

  # ---------------------------------------------------------------------------
  # Return sf
  # ---------------------------------------------------------------------------
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  print(output_dir)
  if (return == "sf") {

    if (is.null(crs)) {
      stop("Provide `crs` when `return = \"sf\"`.", call. = FALSE)
    }

    out_sf <- sf::st_as_sf(
      out_df,
      coords = c("x", "y"),
      crs    = crs,
      remove = FALSE
    )

    if (!is.null(output_dir)) {
      for (i in seq_along(species_names)) {
        sp      <- species_names[i]
        sp_safe <- safe_species_names[i]
        sp_dir  <- file.path(output_dir,
                             sp_safe)
        dir.create(sp_dir, recursive = TRUE, showWarnings = FALSE)

        sp_sf <- out_sf[out_sf$species == sp, ]

        sfarrow::st_write_parquet(
          sp_sf,
          file.path(
            sp_dir,
            paste0(sp_safe, "_occupancy_predictions.parquet")
          ),
          compression = "snappy"
        )
      }
    }

    return(out_sf)
  }

  # ---------------------------------------------------------------------------
  # Return raster
  # ---------------------------------------------------------------------------

  if (return == "raster") {

    if (is.null(template)) {
      stop("Provide `template` when `return = \"raster\"`.", call. = FALSE)
    }

    rasters <- lapply(seq_along(species_names), function(i) {
      sp    <- species_names[i]
      sp_df <- out_df[out_df$species == sp, ]

      sp_vect <- terra::vect(
        sp_df,
        geom = c("x", "y"),
        crs  = terra::crs(template)
      )

      r <- terra::rast(
        lapply(c("mean", "median", "sd", "lower", "upper"), function(field) {
          terra::rasterize(sp_vect, template, field = field)
        })
      )

      names(r) <- c("mean", "median", "sd", "lower", "upper")
      r
    })

    names(rasters) <- species_names

    if (!is.null(output_dir)) {
      for (i in seq_along(species_names)) {
        sp      <- species_names[i]
        sp_safe <- safe_species_names[i]
        sp_dir  <- file.path(output_dir, sp_safe)
        dir.create(sp_dir, recursive = TRUE, showWarnings = FALSE)

        terra::writeRaster(
          rasters[[sp]],
          filename  = file.path(sp_dir,
                                paste0(sp_safe, "_occupancy_predictions.tif")),
          overwrite = TRUE
        )
      }
    }

    return(rasters)
  }
}
