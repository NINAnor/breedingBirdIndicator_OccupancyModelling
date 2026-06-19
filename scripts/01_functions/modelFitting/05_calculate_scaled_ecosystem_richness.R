#' Calculate Scaled Ecosystem Richness from Species Predictions
#'
#' Computes scaled ecosystem richness by aggregating species prediction outputs
#' (raster or vector) and normalizing expected richness against region-specific
#' reference conditions derived from spatial intersections. The function supports
#' both raster (`terra`) and vector (`sf`/GeoParquet) inputs, performs spatial joins,
#' and produces FAIR-compliant outputs and metadata.
#'
#' This function supports FAIR principles:
#' \itemize{
#'   \item \strong{Findable}: Outputs use structured filenames and species-specific files.
#'   \item \strong{Accessible}: Supports GeoTIFF, GeoParquet, CSV, and RDS formats.
#'   \item \strong{Interoperable}: Preserves CRS, uses standard spatial formats, and
#'         harmonized species naming.
#'   \item \strong{Reusable}: Records parameters, processing steps, and metadata.
#' }
#'
#' @param input_dir Character. Directory containing species prediction outputs
#'   (subdirectories with `.tif` or `.parquet` files).
#' @param reference_shape An \code{sf} object defining reference locations
#'   used to estimate baseline richness.
#' @param categorisation_shape An \code{sf} object containing region categories
#'   (e.g., biogeographic regions).
#' @param ecosystem_type Character. Ecosystem type used for filtering reference
#'   data (e.g., `"Forest"`).
#' @param value_col Character. Name of the column or raster layer containing
#'   prediction values. Default is `"mean"`.
#' @param point_spacing Numeric. Sampling spacing (in meters) when extracting
#'   raster predictions. Default is `1000`.
#' @param region_col Character. Column name defining regions used for scaling.
#'   Default is `"bioGeographic_region"`.
#' @param forest_only_strict Logical. If `TRUE`, restricts forest analyses to
#'   protected forest areas. Default is `TRUE`.
#' @param na.rm Logical. If `TRUE`, removes missing values in aggregation steps.
#' @param output_dir Character or `NULL`. Directory to save outputs.
#'   Required if `save_outputs = TRUE`.
#' @param save_outputs Logical. If `TRUE`, saves outputs and metadata to disk.
#' @param overwrite Logical. If `TRUE`, overwrite existing raster outputs.
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{results}: List of species prediction objects.
#'   \item \code{reference_shape_with_categories}: Reference features with
#'         assigned region categories.
#'   \item \code{reference_predictions}: Extracted species predictions at reference locations.
#'   \item \code{reference_predictions_expected_richness}: Regional summary of
#'         expected richness for reference conditions.
#'   \item \code{prediction_values_scaled_richness}: Per-point scaled richness values.
#'   \item \code{scaled_richness_summary}: Regional summary statistics.
#'   \item \code{lower_bound}: Minimum expected richness used in scaling.
#'   \item \code{regional_reference}: Maximum expected richness per region.
#' }
#'
#' @details
#' The function performs the following workflow:
#' \enumerate{
#'   \item Validates inputs and ensures consistent coordinate reference systems.
#'   \item Intersects \code{reference_shape} with \code{categorisation_shape}.
#'   \item Optionally filters reference data for forest ecosystems.
#'   \item Assigns each reference unit to a dominant region category.
#'   \item Loads species prediction outputs:
#'     \itemize{
#'       \item GeoTIFF files (\pkg{terra})
#'       \item GeoParquet files (\pkg{sfarrow})
#'     }
#'   \item Extracts species values at reference locations.
#'   \item Computes expected species richness per reference point.
#'   \item Aggregates richness by region to determine reference benchmarks.
#'   \item Converts all species predictions to point-based `sf` objects.
#'   \item Joins predictions and assigns regional categories spatially.
#'   \item Computes expected richness across species.
#'   \item Scales richness values between a lower bound and regional reference:
#'     \deqn{
#'       scaled = (expected - lower) / (reference - lower)
#'     }
#'   \item Constrains scaled values between 0 and 1.
#'   \item Produces regional summary statistics.
#'   \item Optionally writes outputs (GeoParquet, GeoTIFF, CSV, RDS) and JSON metadata.
#' }
#'
#' Scaling interpretation:
#' \itemize{
#'   \item Values close to 0 indicate low richness relative to baseline conditions.
#'   \item Values close to 1 indicate richness approaching regional reference levels.
#' }
#'
#' Output files (if \code{save_outputs = TRUE}) include:
#' \itemize{
#'   \item Reference and prediction datasets (GeoParquet / GeoTIFF)
#'   \item Regional summaries (CSV)
#'   \item Species-specific outputs
#'   \item Metadata file (`metadata.json`)
#'   \item Serialized results (`all_results.rds`)
#' }
#'
#' @examples
#' \dontrun{
#' result <- calculate_scaled_ecosystem_richness(
#'   input_dir = "data/species_predictions",
#'   reference_shape = ref_sf,
#'   categorisation_shape = region_sf,
#'   ecosystem_type = "Forest",
#'   output_dir = "outputs/richness",
#'   save_outputs = TRUE
#' )
#'
#' result$scaled_richness_summary
#' }
#'
#' @seealso
#' [terra::rast()], [sf::st_intersection()], [sfarrow::st_read_parquet()],
#' [dplyr::summarise()]
#'
#' @export


calculate_scaled_ecosystem_richness <- function(
    input_dir,
    reference_shape,
    categorisation_shape,
    ecosystem_type,
    value_col = "mean",
    point_spacing = 1000,
    region_col = "bioGeographic_region",
    forest_only_strict = TRUE,
    na.rm = TRUE,
    output_dir = NULL,
    save_outputs = FALSE,
    overwrite = TRUE
) {

  stopifnot(inherits(reference_shape, "sf"))
  stopifnot(inherits(categorisation_shape, "sf"))

  if (!dir.exists(input_dir)) {
    stop("`input_dir` does not exist: ", input_dir, call. = FALSE)
  }

  if (save_outputs && is.null(output_dir)) {
    stop("`output_dir` must be supplied when `save_outputs = TRUE`.", call. = FALSE)
  }

  if (save_outputs) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  if (sf::st_crs(reference_shape) != sf::st_crs(categorisation_shape)) {
    categorisation_shape <- sf::st_transform(
      categorisation_shape,
      sf::st_crs(reference_shape)
    )
  }

  intersections <- sf::st_intersection(reference_shape, categorisation_shape)
  intersections$overlap_area <- sf::st_area(intersections)

  intersections_for_ecosystem <- intersections

  if (ecosystem_type %in% c("Forest", "Forests", "forest", "forests")) {
    if (forest_only_strict) {
      intersections_for_ecosystem <- intersections_for_ecosystem |>
        dplyr::filter(skogvern %in% "ja") |>
        dplyr::filter(verneform %in% c("nasjonalpark", "naturreservat"))
    } else {
      intersections_for_ecosystem <- intersections_for_ecosystem |>
        dplyr::filter(skogvern %in% "ja")
    }
  }

  if (nrow(intersections_for_ecosystem) == 0) {
    stop("No reference features remain after ecosystem filtering.", call. = FALSE)
  }

  if (!"identifikasjon_lokalid" %in% names(intersections_for_ecosystem)) {
    intersections_for_ecosystem$.reference_id <- seq_len(nrow(intersections_for_ecosystem))
    reference_id_col <- ".reference_id"
  } else {
    reference_id_col <- "identifikasjon_lokalid"
  }

  reference_shape_with_categories <- intersections_for_ecosystem |>
    dplyr::group_by(.data[[reference_id_col]]) |>
    dplyr::slice_max(overlap_area, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(.data[[region_col]])) |>
    dplyr::select(dplyr::all_of(region_col))

  species_dirs <- list.dirs(input_dir, recursive = FALSE, full.names = TRUE)

  results <- setNames(
    lapply(species_dirs, function(sp_dir) {
      files <- list.files(sp_dir, full.names = TRUE)

      tif_files <- files[grepl("\\.tif$", files)]
      parquet_files <- files[grepl("\\.parquet$", files)]

      if (length(tif_files) > 0) {
        terra::rast(tif_files[1])
      } else if (length(parquet_files) > 0) {
        sfarrow::st_read_parquet(parquet_files[1])
      } else {
        NULL
      }
    }),
    basename(species_dirs)
  )

  results <- Filter(Negate(is.null), results)

  if (length(results) == 0) {
    stop("No valid .tif or .parquet species outputs were found.", call. = FALSE)
  }

  species_cols <- safe_name(names(results))

  reference_predictions <- extract_species_values_at_reference(
    allResults = results,
    reference_shape_with_regions = reference_shape_with_categories,
    value_col = value_col,
    summarise_fun = mean,
    na.rm = na.rm,
    keep_fair_metadata = FALSE,
    point_spacing = point_spacing,
    region_col = region_col
  )

  reference_predictions_expected_richness <- reference_predictions |>
    dplyr::mutate(
      expected_richness = rowSums(
        sf::st_drop_geometry(
          dplyr::select(reference_predictions, dplyr::all_of(species_cols))
        ),
        na.rm = TRUE
      )
    ) |>
    sf::st_drop_geometry() |>
    dplyr::filter(!is.na(.data[[region_col]])) |>
    dplyr::group_by(.data[[region_col]]) |>
    dplyr::summarise(
      n_reference_points = dplyr::n(),
      min_expected_richness = min(expected_richness, na.rm = TRUE),
      mean_expected_richness = mean(expected_richness, na.rm = TRUE),
      max_expected_richness = max(expected_richness, na.rm = TRUE),
      .groups = "drop"
    )

  regional_reference <- reference_predictions_expected_richness |>
    dplyr::select(
      dplyr::all_of(region_col),
      reference_richness = max_expected_richness
    )

  species_points <- lapply(names(results), function(sp) {
    x <- results[[sp]]
    sp_safe <- safe_name(sp)

    if (inherits(x, "sf")) {
      if (!value_col %in% names(x)) {
        warning("Skipping ", sp, ": column not found: ", value_col)
        return(NULL)
      }

      x <- x |>
        dplyr::select(dplyr::all_of(value_col), geometry)

      names(x)[names(x) == value_col] <- sp_safe
      return(x)
    }

    if (inherits(x, "SpatRaster")) {
      if (!value_col %in% names(x)) {
        warning("Skipping ", sp, ": raster layer not found: ", value_col)
        return(NULL)
      }

      pts <- terra::as.points(x[[value_col]], na.rm = FALSE)
      pts_sf <- sf::st_as_sf(pts)

      names(pts_sf)[names(pts_sf) == value_col] <- sp_safe
      return(pts_sf)
    }

    warning("Skipping ", sp, ": unsupported object type.")
    NULL
  })

  species_points <- Filter(Negate(is.null), species_points)

  if (length(species_points) == 0) {
    stop("No valid species prediction point objects could be created.", call. = FALSE)
  }

  prediction_values <- Reduce(function(a, b) {
    if (sf::st_crs(a) != sf::st_crs(b)) {
      b <- sf::st_transform(b, sf::st_crs(a))
    }

    sf::st_join(a, b, join = sf::st_equals, left = TRUE)
  }, species_points)

  bio_regions_for_join <- categorisation_shape |>
    dplyr::select(dplyr::all_of(region_col))

  if (sf::st_crs(bio_regions_for_join) != sf::st_crs(prediction_values)) {
    bio_regions_for_join <- sf::st_transform(
      bio_regions_for_join,
      sf::st_crs(prediction_values)
    )
  }

  idx <- sf::st_intersects(prediction_values, bio_regions_for_join)

  prediction_values[[region_col]] <- vapply(
    idx,
    function(i) {
      if (length(i) == 0) NA_character_
      else bio_regions_for_join[[region_col]][i[1]]
    },
    character(1)
  )

  prediction_values <- prediction_values |>
    dplyr::filter(!is.na(.data[[region_col]]))

  species_cols_present <- intersect(species_cols, names(prediction_values))

  prediction_values_scaled_richness <- prediction_values |>
    dplyr::mutate(
      expected_richness = rowSums(
        sf::st_drop_geometry(
          dplyr::select(prediction_values, dplyr::all_of(species_cols_present))
        ),
        na.rm = TRUE
      )
    )

  lower_bound <- prediction_values_scaled_richness |>
    sf::st_drop_geometry() |>
    dplyr::summarise(
      lower_bound = min(expected_richness, na.rm = TRUE)
    ) |>
    dplyr::pull(lower_bound)

  prediction_values_scaled_richness <- prediction_values_scaled_richness |>
    dplyr::left_join(regional_reference, by = region_col) |>
    dplyr::mutate(
      lower_bound = lower_bound,
      scaled_richness =
        (expected_richness - lower_bound) /
        (reference_richness - lower_bound),
      scaled_richness = dplyr::case_when(
        is.na(scaled_richness) ~ NA_real_,
        is.infinite(scaled_richness) ~ NA_real_,
        TRUE ~ pmax(0, pmin(1, scaled_richness))
      )
    )

  scaled_richness_summary <- prediction_values_scaled_richness |>
    sf::st_drop_geometry() |>
    dplyr::filter(!is.na(.data[[region_col]])) |>
    dplyr::group_by(.data[[region_col]]) |>
    dplyr::summarise(
      n_points = dplyr::n(),
      lower_bound = dplyr::first(lower_bound),
      reference_richness = dplyr::first(reference_richness),
      min_expected_richness = min(expected_richness, na.rm = TRUE),
      mean_expected_richness = mean(expected_richness, na.rm = TRUE),
      max_expected_richness = max(expected_richness, na.rm = TRUE),
      min_scaled_richness = min(scaled_richness, na.rm = TRUE),
      mean_scaled_richness = mean(scaled_richness, na.rm = TRUE),
      max_scaled_richness = max(scaled_richness, na.rm = TRUE),
      .groups = "drop"
    )

  output <- list(
    results = results,
    reference_shape_with_categories = reference_shape_with_categories,
    reference_predictions = reference_predictions,
    reference_predictions_expected_richness = reference_predictions_expected_richness,
    prediction_values_scaled_richness = prediction_values_scaled_richness,
    scaled_richness_summary = scaled_richness_summary,
    lower_bound = lower_bound,
    regional_reference = regional_reference
  )

  if (save_outputs) {

    save_object <- function(x, filename_base, dir = output_dir) {
      if (inherits(x, "sf")) {
        sfarrow::st_write_parquet(
          obj = x,
          dsn = file.path(dir, paste0(filename_base, ".parquet")),
          compression = "snappy"
        )
      } else if (inherits(x, "SpatRaster")) {
        terra::writeRaster(
          x,
          filename = file.path(dir, paste0(filename_base, ".tif")),
          overwrite = overwrite
        )
      } else if (inherits(x, "data.frame")) {
        utils::write.csv(
          x,
          file = file.path(dir, paste0(filename_base, ".csv")),
          row.names = FALSE
        )
      } else {
        saveRDS(
          x,
          file = file.path(dir, paste0(filename_base, ".rds"))
        )
      }
    }

    save_object(reference_shape_with_categories, "reference_shape_with_categories")
    save_object(reference_predictions, "reference_predictions")
    save_object(reference_predictions_expected_richness, "reference_predictions_expected_richness")
    save_object(prediction_values_scaled_richness, "prediction_values_scaled_richness")
    save_object(scaled_richness_summary, "scaled_richness_summary")
    save_object(regional_reference, "regional_reference")
    save_object(lower_bound, "lower_bound")

    species_output_dir <- file.path(output_dir, "species_results")
    dir.create(species_output_dir, recursive = TRUE, showWarnings = FALSE)

    for (sp in names(results)) {
      save_object(
        results[[sp]],
        safe_name(sp),
        dir = species_output_dir
      )
    }

    metadata <- list(
      title = "Scaled ecosystem richness indicators",
      created = as.character(Sys.time()),
      ecosystem_type = ecosystem_type,
      input_dir = normalizePath(input_dir, mustWork = FALSE),
      output_dir = normalizePath(output_dir, mustWork = FALSE),
      value_column = value_col,
      point_spacing = point_spacing,
      region_column = region_col,
      forest_only_strict = forest_only_strict,
      na_rm = na.rm,
      lower_bound = lower_bound,
      species = names(results),
      species_columns = safe_name(names(results)),
      file_format_rules = list(
        sf = "GeoParquet",
        SpatRaster = "GeoTIFF",
        data_frame = "CSV",
        other = "RDS"
      ),
      FAIR_notes = list(
        findable = "Outputs use descriptive filenames and species-specific files.",
        accessible = "Vector outputs are GeoParquet, raster outputs are GeoTIFF, tables are CSV.",
        interoperable = "Spatial outputs preserve CRS information.",
        reusable = "Metadata records parameters, sources, species names, software versions, and format rules."
      ),
      software_versions = list(
        R = R.version.string,
        sf = as.character(utils::packageVersion("sf")),
        terra = as.character(utils::packageVersion("terra")),
        dplyr = as.character(utils::packageVersion("dplyr")),
        sfarrow = as.character(utils::packageVersion("sfarrow")),
        jsonlite = as.character(utils::packageVersion("jsonlite"))
      )
    )

    jsonlite::write_json(
      metadata,
      path = file.path(output_dir, "metadata.json"),
      pretty = TRUE,
      auto_unbox = TRUE
    )

    saveRDS(
      output,
      file = file.path(output_dir, "all_results.rds")
    )
  }

  output
}
