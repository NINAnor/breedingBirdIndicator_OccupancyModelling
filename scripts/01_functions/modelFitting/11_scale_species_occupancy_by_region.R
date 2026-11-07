#' Scale species occupancy by a region-specific reference value
#'
#' For each species, the reference value of a region is the
#' `reference_quantile` (from `reference_data`) of that species' occupancy
#' among all prediction locations in the region. Occupancy at every location
#' is then scaled with [ecTools::ec_normalise()], using that reference value
#' as `x100`.
#'
#' @param results Named list with one prediction object per species. Each
#'   element is either an `sf` object or a `SpatRaster` containing
#'   `value_col`. All species must share the same prediction locations, in the
#'   same order. Names may use underscores (`Anthus_pratensis`).
#' @param reference_data Data frame with columns `Species` and
#'   `reference_quantile` (e.g. from [get_reference_quantiles()]). Species
#'   without a reference are skipped with a warning.
#' @param categorisation_shape `sf` polygons defining the regions.
#' @param region_col Character. Column of `categorisation_shape` with the
#'   region names.
#' @param value_col Character. Occupancy column to scale.
#' @param truncation Logical. Passed to [ecTools::ec_normalise()]; if `TRUE`
#'   (default) scaled values are limited to the range 0 to 1.
#' @param ... Further arguments for [ecTools::ec_normalise()], e.g. `x0`,
#'   `fun`, `convex_exponent`. `variable` and `x100` are set by this function.
#'
#' @return A list with
#'   \describe{
#'     \item{scaled}{`sf` with one scaled column per species (named by species
#'       with underscores), a `bio_region` column and the geometry.}
#'     \item{mean_raw}{Numeric vector: unscaled mean of `value_col` across
#'       the species, one value per location.}
#'     \item{reference_values}{Data frame with `bio_region`, `species`,
#'       `reference_quantile`, `reference_value` and `n_locations`.}
#'   }
#'
#' @details Locations outside all regions get `NA` for the region and the
#'   scaled values. If a location falls in several regions the first is used.
#'   Regions whose reference value is missing or not positive give `NA`
#'   scaled values (with a warning).
#' @export
scale_species_occupancy_by_region <- function(
  results,
  reference_data,
  categorisation_shape,
  region_col,
  value_col = "occupancy_mean",
  truncation = TRUE,
  ...
) {
  # ---- Validate ---------------------------------------------------------------
  if (!is.list(results) || is.null(names(results)) || !length(results)) {
    stop(
      "`results` must be a named list of species predictions.",
      call. = FALSE
    )
  }
  if (!all(c("Species", "reference_quantile") %in% names(reference_data))) {
    stop(
      "`reference_data` needs columns `Species` and `reference_quantile`.",
      call. = FALSE
    )
  }
  if (!requireNamespace("ecTools", quietly = TRUE)) {
    stop("Package `ecTools` is required.", call. = FALSE)
  }
  if (!inherits(categorisation_shape, "sf")) {
    stop("`categorisation_shape` must be an sf object.", call. = FALSE)
  }
  if (!region_col %in% names(categorisation_shape)) {
    stop("`region_col` not found in `categorisation_shape`.", call. = FALSE)
  }
  q <- reference_data$reference_quantile
  if (anyNA(q) || any(q <= 0 | q >= 1)) {
    stop("Reference quantiles must be between 0 and 1.", call. = FALSE)
  }

  # ---- Match species to reference quantiles ------------------------------------
  species_label <- gsub("_", " ", names(results), fixed = TRUE)
  ref_q <- reference_data$reference_quantile[
    match(species_label, reference_data$Species)
  ]

  if (anyNA(ref_q)) {
    warning(
      "No reference quantile for (skipped): ",
      paste(species_label[is.na(ref_q)], collapse = ", "),
      call. = FALSE
    )
  }
  keep <- !is.na(ref_q)
  if (!any(keep)) {
    stop("No species in `results` have a reference quantile.", call. = FALSE)
  }
  results <- results[keep]
  species_label <- species_label[keep]
  ref_q <- ref_q[keep]

  # ---- Values and geometry -----------------------------------------------------
  get_values <- function(x, nm) {
    if (inherits(x, "SpatRaster")) {
      if (!value_col %in% names(x)) {
        stop("Layer `", value_col, "` missing for ", nm, call. = FALSE)
      }
      terra::values(x[[value_col]], mat = FALSE)
    } else if (inherits(x, "sf")) {
      if (!value_col %in% names(x)) {
        stop("Column `", value_col, "` missing for ", nm, call. = FALSE)
      }
      x[[value_col]]
    } else {
      stop("Unsupported prediction class for ", nm, call. = FALSE)
    }
  }

  values <- Map(get_values, results, names(results))

  n_loc <- unique(lengths(values))
  if (length(n_loc) != 1L) {
    stop(
      "Species have different numbers of prediction locations.",
      call. = FALSE
    )
  }

  first <- results[[1]]
  geometry <- if (inherits(first, "SpatRaster")) {
    sf::st_geometry(sf::st_as_sf(terra::as.points(first[[1]], na.rm = FALSE)))
  } else {
    sf::st_geometry(first)
  }

  # ---- Assign each location to a region (once for all species) -----------------
  regions <- categorisation_shape[region_col]
  regions <- sf::st_transform(regions, sf::st_crs(geometry))

  within_idx <- sf::st_within(geometry, regions)
  region_row <- vapply(
    within_idx,
    function(i) if (length(i)) i[1L] else NA_integer_,
    integer(1)
  )
  region <- as.character(regions[[region_col]])[region_row]

  # ---- Reference value per region and species ----------------------------------
  scaled <- vector("list", length(values))
  ref_tables <- vector("list", length(values))
  species_names <- gsub(" ", "_", species_label, fixed = TRUE)

  for (k in seq_along(values)) {
    v <- values[[k]]

    ref_by_region <- tapply(
      v,
      region,
      stats::quantile,
      probs = ref_q[k],
      na.rm = TRUE,
      names = FALSE,
      type = 7
    )
    ref_by_region <- stats::setNames(
      as.numeric(ref_by_region),
      names(ref_by_region)
    )
    ref_by_region[is.finite(ref_by_region) & ref_by_region <= 0] <- NA_real_

    bad <- names(ref_by_region)[!is.finite(ref_by_region)]
    if (length(bad)) {
      warning(
        species_label[k],
        ": no usable reference value in region(s): ",
        paste(bad, collapse = ", "),
        call. = FALSE
      )
    }

    # ec_normalise() takes a single x100, so scale region by region
    scaled_k <- rep(NA_real_, length(v))
    for (r in names(ref_by_region)) {
      x100 <- ref_by_region[[r]]
      if (!is.finite(x100)) {
        next
      }
      idx <- which(region == r)
      scaled_k[idx] <- ecTools::ec_normalise(
        variable = v[idx],
        x100 = x100,
        truncation = truncation,
        ...
      )
    }
    scaled[[k]] <- scaled_k

    ref_tables[[k]] <- data.frame(
      bio_region = names(ref_by_region),
      species = species_label[k],
      reference_quantile = ref_q[k],
      reference_value = unname(ref_by_region),
      n_locations = as.integer(table(region)[names(ref_by_region)]),
      stringsAsFactors = FALSE
    )
  }

  # ---- Assemble outputs --------------------------------------------------------
  scaled_df <- as.data.frame(stats::setNames(scaled, species_names))
  scaled_df$bio_region <- region

  scaled_sf <- sf::st_sf(scaled_df, geometry = geometry)

  # Unscaled mean occupancy across the species used (for mapping)
  mean_raw <- rowMeans(do.call(cbind, values), na.rm = TRUE)
  mean_raw[is.nan(mean_raw)] <- NA_real_

  list(
    scaled = scaled_sf,
    reference_values = do.call(rbind, ref_tables),
    mean_raw = mean_raw
  )
}
